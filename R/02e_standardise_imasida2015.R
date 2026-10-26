#!/usr/bin/env Rscript
# =============================================================================
# Standardise Mozambique IMASIDA 2015 (a DHS-format AIDS Indicator Survey, not
# a PHIA) to the round-0 schema, and add Mozambique to the combined two-round
# survey-undiagnosed file (Z_hat = logit P(HIV+ and unaware) per province).
# Reads:  data/raw/phia/mozambique/imasida2015/*.DTA (not in the release) and
#         data/processed/phia/{phia_mozambique, survey_undiagnosed_round0_round1}.rds
# Writes: data/processed/phia/phia_mozambique_2015.rds; updates the combined
#         survey_undiagnosed_round0_round1.{rds,csv} in data/processed/phia/
# Run:    Rscript R/02e_standardise_imasida2015.R  (from the repository root)
# =============================================================================

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(survey)
})
options(survey.lonely.psu = "adjust")

PROJ <- "."
IM   <- file.path(PROJ, "data/raw/phia/mozambique/imasida2015")
OUT  <- file.path(PROJ, "data/processed/phia")

zap <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

region_map <- c("1"="Niassa","2"="Cabo Delgado","3"="Nampula","4"="Zambézia",
                "5"="Tete","6"="Manica","7"="Sofala","8"="Inhambane","9"="Gaza",
                "10"="Maputo","11"="Maputo")  # 10 & 11 both -> Maputo

# ---- individual recodes (women + men): keys, region, testing history --------
read_ind <- function(path, pfx, sex_label) {
  v <- function(s) paste0(pfx, s)
  d <- read_dta(path, col_select = any_of(c(
    v("001"), v("002"), v("003"), v("005"), v("012"),
    v("022"), v("023"), v("024"), v("025"), v("781"), v("828"))))
  tibble(
    cluster = zap(d[[v("001")]]),
    hh      = zap(d[[v("002")]]),
    line    = zap(d[[v("003")]]),
    age     = zap(d[[v("012")]]),
    region  = zap(d[[v("024")]]),
    urban   = zap(d[[v("025")]]),
    stratum = if (v("023") %in% names(d)) zap(d[[v("023")]])
              else if (v("022") %in% names(d)) zap(d[[v("022")]]) else zap(d[[v("024")]]),
    ever_tested  = zap(d[[v("781")]]),
    got_result   = if (v("828") %in% names(d)) zap(d[[v("828")]]) else NA_real_,
    sex = sex_label)
}

message("=== IMASIDA 2015: reading individual recodes ===")
ind <- bind_rows(
  read_ind(file.path(IM,"MZIR71DT/MZIR71FL.DTA"), "v",  "Female"),
  read_ind(file.path(IM,"MZMR71DT/MZMR71FL.DTA"), "mv", "Male")
)

# ---- AR biomarker recode ----------------------------------------------------
# hiv03: 0 = negative, 1 to 3 = positive, 7 and 9 = indeterminate or inconclusive (NA); hiv05: HIV weight, stored x 1e6
message("=== IMASIDA 2015: reading AR biomarker ===")
ar <- read_dta(file.path(IM,"MZAR72DT/MZAR72FL.DTA"),
               col_select = any_of(c("hivclust","hivnumb","hivline","hiv03","hiv05")))
ar <- tibble(cluster=zap(ar$hivclust), hh=zap(ar$hivnumb), line=zap(ar$hivline),
             hiv03=zap(ar$hiv03), hivweight=zap(ar$hiv05)/1e6)

# ---- merge biomarker onto individuals --------------------------------------
# No awareness question in a DHS-AIS; diagnosed proxy = ever tested (v781 == 1) and received the last result (v828 == 1)
phia_moz0 <- ind %>%
  inner_join(ar, by = c("cluster","hh","line")) %>%
  mutate(
    hivstatus = case_when(hiv03 %in% c(1,2,3) ~ 1L, hiv03 == 0 ~ 0L, TRUE ~ NA_integer_),
    diagnosed = ever_tested == 1 & got_result == 1,
    diagnosed = ifelse(is.na(diagnosed), FALSE, diagnosed),
    aware_bin = case_when(hivstatus != 1 ~ NA, diagnosed ~ TRUE, TRUE ~ FALSE),
    province  = recode(as.character(region), !!!region_map, .default = NA_character_),
    age_group = cut(age, breaks=c(15,20,25,50,Inf),
                    labels=c("15-19","20-24","25-49","50+"), right=FALSE)
  ) %>%
  transmute(country="Mozambique", wave="round0", year=2015L,
            personid = paste(cluster, hh, line, sep="-"),
            province, district = NA_character_, age, sex,
            hivstatus, aware_bin,
            awareselfreported = NA_real_, arvstatus = NA_real_,
            btwt0 = hivweight, indweight = hivweight,
            varstrat = stratum, varunit = cluster) %>%
  filter(age >= 15, !is.na(hivstatus), !is.na(btwt0), btwt0 > 0)

wp <- with(filter(phia_moz0, age<=49), weighted.mean(hivstatus, btwt0))
cat(sprintf("IMASIDA 2015: rows=%d  HIV+=%d  weighted prev 15-49=%.1f%%  provinces=%d (NA geo=%d)\n",
            nrow(phia_moz0), sum(phia_moz0$hivstatus==1), 100*wp,
            dplyr::n_distinct(phia_moz0$province, na.rm=TRUE), sum(is.na(phia_moz0$province))))
cat("Diagnosed proxy among HIV+:",
    round(100*weighted.mean(phia_moz0$aware_bin[phia_moz0$hivstatus==1], phia_moz0$btwt0[phia_moz0$hivstatus==1], na.rm=TRUE),1), "%\n")

saveRDS(phia_moz0, file.path(OUT, "phia_mozambique_2015.rds"))

# =============================================================================
# Z_hat_p = logit P(HIV+ & unaware) per province (same method as 02d)
# =============================================================================
compute_undiagnosed <- function(phia, area_col="province") {
  d <- phia %>%
    filter(!is.na(.data[[area_col]]), !is.na(btwt0), btwt0>0) %>%
    mutate(hiv_unaware = case_when(
      hivstatus==0 ~ 0L,
      hivstatus==1 & aware_bin==FALSE ~ 1L,
      hivstatus==1 & aware_bin==TRUE  ~ 0L, TRUE ~ NA_integer_)) %>%
    filter(!is.na(hiv_unaware))
  d$area <- d[[area_col]]
  des <- svydesign(ids=~varunit, strata=~varstrat, weights=~btwt0, data=d, nest=TRUE)
  est <- svyby(~hiv_unaware, ~area, des, svymean, na.rm=TRUE)
  ns <- d %>% group_by(area) %>%
    summarise(n=n(), n_hiv_pos=sum(hivstatus==1), n_unaware=sum(hiv_unaware==1), .groups="drop")
  tibble(area=as.character(est$area), p_undiag=as.numeric(est$hiv_unaware), se_p=as.numeric(est$se)) %>%
    mutate(p_undiag=pmin(pmax(p_undiag,1e-6),1-1e-6), Z_hat=qlogis(p_undiag),
           se_Z=se_p/(p_undiag*(1-p_undiag))) %>%
    left_join(ns, by="area") %>% arrange(desc(n_hiv_pos))
}

moz_r1 <- readRDS(file.path(OUT,"phia_mozambique.rds"))   # INSIDA 2021 (true PHIA)
moz_undiag <- bind_rows(
  compute_undiagnosed(phia_moz0) %>% mutate(country="Mozambique", wave="round0", year=2015L),
  compute_undiagnosed(moz_r1)   %>% mutate(country="Mozambique", wave="round1", year=2021L)
) %>% select(country, wave, year, area, n, n_hiv_pos, n_unaware, p_undiag, se_p, Z_hat, se_Z)

# ---- append Mozambique to the combined four-country file --------------------
comb_path <- file.path(OUT, "survey_undiagnosed_round0_round1.rds")
comb <- readRDS(comb_path) %>% filter(country != "Mozambique")
comb <- bind_rows(comb, moz_undiag) %>% arrange(country, wave, desc(n_hiv_pos))
saveRDS(comb, comb_path)
readr::write_csv(comb, file.path(OUT, "survey_undiagnosed_round0_round1.csv"))

cat("\n=== Mozambique Z_hat_p by wave/province ===\n")
print(as.data.frame(moz_undiag %>% mutate(across(c(p_undiag,se_p,Z_hat,se_Z), ~round(.,3)))), row.names=FALSE)
cat("\nCombined file now covers countries:", paste(sort(unique(comb$country)), collapse=", "), "\n")
message("=== Done (Mozambique IMASIDA 2015) ===")

#!/usr/bin/env Rscript
# =============================================================================
# Standardise the round-0 PHIA surveys (MPHIA 2015-16, LePHIA 2016-17, ZIMPHIA
# 2015-16) to the round-1 schema, then compute Z_hat = logit P(HIV+ and unaware)
# per area with design-based SE for both rounds (Mozambique IMASIDA 2015: 02e).
# Reads:  data/raw/phia/{malawi,lesotho,zimbabwe}/ round-0 CSVs (not in the
#         release); data/processed/phia/phia_{malawi,lesotho,zimbabwe}.rds
# Writes: data/processed/phia/phia_{malawi_2015,lesotho_2016,zimbabwe_2015}.rds
#         and survey_undiagnosed_round0_round1.{rds,csv} (same folder)
# Run:    Rscript R/02d_standardise_phia_round0.R  (from the repository root)
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(survey)
})
options(survey.lonely.psu = "adjust")

PROJ <- "."
PHIA <- file.path(PROJ, "data", "raw", "phia")
OUT  <- file.path(PROJ, "data", "processed", "phia")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))

# ---- geography maps ----------------------------------------------------------
les_district_map <- c("1"="Berea","2"="Butha Buthe","3"="Leribe","4"="Mafeteng",
                      "5"="Maseru","6"="Mohale's Hoek","7"="Mokhotlong",
                      "8"="Qacha's Nek","9"="Quthing","10"="Thaba Tseka")
zim_province_map <- c("0"="Bulawayo","1"="Manicaland","2"="Mashonaland Central",
                      "3"="Mashonaland East","4"="Mashonaland West",
                      "5"="Matabeleland North","6"="Matabeleland South",
                      "7"="Midlands","8"="Masvingo","9"="Harare")
mwi_region_map   <- c("1"="Northern","2"="Central","3"="Southern")

# round-1 Malawi 5-zone -> 3-region crosswalk (so both waves compare at region level)
mwi_zone_to_region <- c(
  "Northern Zone"     = "Northern",
  "Central East Zone" = "Central", "Central West Zone" = "Central",
  "South East Zone"   = "Southern", "South West Zone"   = "Southern")

# =============================================================================
# Recode one round-0 PHIA survey to the standardised schema
# =============================================================================
recode_phia_round0 <- function(country, year, bio_file, geo_file, geo_key, geo_var,
                               geo_map, geo_is_district = FALSE) {
  message("\n=== ", country, " ", year, " ===")
  bio <- read_csv(bio_file, show_col_types = FALSE, progress = FALSE,
                  col_select = any_of(c("personid","householdid","centroidid",
                    "gender","age","hivstatusfinal","aware","awareselfreported",
                    "art","arvstatus","btwt0","varstrat","varunit")))

  # geography join (region from hh by householdid, or district/province from ind by personid)
  geo <- read_csv(geo_file, show_col_types = FALSE, progress = FALSE,
                  col_select = any_of(c(geo_key, geo_var)))
  geo <- geo %>% distinct(across(all_of(geo_key)), .keep_all = TRUE)
  bio <- left_join(bio, geo, by = geo_key)

  area_lab <- recode(as.character(num(bio[[geo_var]])), !!!geo_map, .default = NA_character_)

  # hivstatusfinal: 1 = positive, 2 = negative; aware: 1 = aware, 2 = unaware (HIV+ only)
  out <- tibble(
    country  = country,
    wave     = "round0",
    year     = year,
    personid = bio$personid,
    province = area_lab,
    district = if (geo_is_district) area_lab else NA_character_,
    age      = num(bio$age),
    sex      = case_when(num(bio$gender)==1 ~ "Male", num(bio$gender)==2 ~ "Female", TRUE ~ NA_character_),
    hivstatus  = case_when(num(bio$hivstatusfinal)==1 ~ 1L, num(bio$hivstatusfinal)==2 ~ 0L, TRUE ~ NA_integer_),
    aware_bin  = case_when(num(bio$aware)==1 ~ TRUE, num(bio$aware)==2 ~ FALSE, TRUE ~ NA),
    awareselfreported = if ("awareselfreported" %in% names(bio)) num(bio$awareselfreported) else NA_real_,
    arvstatus  = if ("arvstatus" %in% names(bio)) num(bio$arvstatus) else NA_real_,
    btwt0    = num(bio$btwt0),
    indweight = num(bio$btwt0),
    varstrat = num(bio$varstrat),
    varunit  = num(bio$varunit)
  ) %>%
    mutate(age_group = cut(age, breaks=c(15,20,25,50,Inf),
                           labels=c("15-19","20-24","25-49","50+"), right=FALSE)) %>%
    filter(age >= 15, !is.na(hivstatus))

  # sanity
  wp <- with(filter(out, !is.na(btwt0), btwt0>0, age<=64),
             weighted.mean(hivstatus, btwt0))
  cat(sprintf("  rows=%d  HIV+=%d  weighted prev 15-64=%.1f%%  areas=%d (NA geo=%d)\n",
              nrow(out), sum(out$hivstatus==1), 100*wp,
              dplyr::n_distinct(out$province, na.rm=TRUE), sum(is.na(out$province))))
  out
}

# =============================================================================
# Survey undiagnosed-prevalence term:  Z_hat_p = logit P(HIV+ & unaware) per area
# design-based mean + SE, delta-method SE on the logit. Shared by both waves.
# =============================================================================
compute_undiagnosed <- function(phia, area_col = "province") {
  d <- phia %>%
    filter(!is.na(.data[[area_col]]), !is.na(btwt0), btwt0 > 0) %>%
    mutate(hiv_unaware = case_when(
      hivstatus == 0 ~ 0L,
      hivstatus == 1 & aware_bin == FALSE ~ 1L,
      hivstatus == 1 & aware_bin == TRUE  ~ 0L,
      TRUE ~ NA_integer_)) %>%
    filter(!is.na(hiv_unaware))
  d$area <- d[[area_col]]

  des <- svydesign(ids=~varunit, strata=~varstrat, weights=~btwt0, data=d, nest=TRUE)
  est <- svyby(~hiv_unaware, ~area, des, svymean, na.rm = TRUE)

  ns <- d %>% group_by(area) %>%
    summarise(n=n(), n_hiv_pos=sum(hivstatus==1), n_unaware=sum(hiv_unaware==1), .groups="drop")

  tibble(area = as.character(est$area),
         p_undiag = as.numeric(est$hiv_unaware),
         se_p     = as.numeric(est$se)) %>%
    mutate(
      p_undiag = pmin(pmax(p_undiag, 1e-6), 1-1e-6),
      Z_hat    = qlogis(p_undiag),
      se_Z     = se_p / (p_undiag * (1 - p_undiag))   # delta method
    ) %>%
    left_join(ns, by = "area") %>%
    arrange(desc(n_hiv_pos))
}

# =============================================================================
# RUN: round-0 recodes
# =============================================================================
phia_mwi0 <- recode_phia_round0("Malawi", 2015L,
  file.path(PHIA,"malawi/mphia2015adultbio.csv"),
  file.path(PHIA,"malawi/mphia2015hh.csv"), "householdid", "region", mwi_region_map, FALSE)

phia_lso0 <- recode_phia_round0("Lesotho", 2016L,
  file.path(PHIA,"lesotho/lephia2016adultbio.csv"),
  file.path(PHIA,"lesotho/lephia2016adultind.csv"), "personid", "district", les_district_map, TRUE)

phia_zwe0 <- recode_phia_round0("Zimbabwe", 2015L,
  file.path(PHIA,"zimbabwe/zimphia2015adultbio.csv"),
  file.path(PHIA,"zimbabwe/zimphia2015adultind.csv"), "personid", "province", zim_province_map, FALSE)

saveRDS(phia_mwi0, file.path(OUT, "phia_malawi_2015.rds"))
saveRDS(phia_lso0, file.path(OUT, "phia_lesotho_2016.rds"))
saveRDS(phia_zwe0, file.path(OUT, "phia_zimbabwe_2015.rds"))

# =============================================================================
# RUN: survey undiagnosed-prevalence, BOTH waves, identical method
# Malawi harmonised to region (3) in both waves.
# =============================================================================
message("\n=== Survey undiagnosed prevalence Z_hat_p per area (both waves) ===")

# round-1 processed
r1 <- list(
  Malawi   = readRDS(file.path(OUT,"phia_malawi.rds")) %>%
               mutate(province = recode(province, !!!mwi_zone_to_region, .default = province)),
  Lesotho  = readRDS(file.path(OUT,"phia_lesotho.rds")),
  Zimbabwe = readRDS(file.path(OUT,"phia_zimbabwe.rds"))
)
# round-0 (Malawi already at region)
r0 <- list(Malawi = phia_mwi0, Lesotho = phia_lso0, Zimbabwe = phia_zwe0)
r1_year <- c(Malawi=2020L, Lesotho=2020L, Zimbabwe=2020L)
r0_year <- c(Malawi=2015L, Lesotho=2016L, Zimbabwe=2015L)

undiag <- bind_rows(lapply(names(r0), function(c) {
  bind_rows(
    compute_undiagnosed(r0[[c]]) %>% mutate(country=c, wave="round0", year=r0_year[[c]]),
    compute_undiagnosed(r1[[c]]) %>% mutate(country=c, wave="round1", year=r1_year[[c]])
  )
})) %>% select(country, wave, year, area, n, n_hiv_pos, n_unaware, p_undiag, se_p, Z_hat, se_Z)

saveRDS(undiag, file.path(OUT, "survey_undiagnosed_round0_round1.rds"))
write_csv(undiag, file.path(OUT, "survey_undiagnosed_round0_round1.csv"))

cat("\n=== Z_hat_p by country/wave/area ===\n")
print(as.data.frame(undiag %>% mutate(across(c(p_undiag,se_p,Z_hat,se_Z), ~round(.,3)))), row.names=FALSE)

cat("\nSaved:\n  phia_{malawi_2015,lesotho_2016,zimbabwe_2015}.rds\n",
    "  survey_undiagnosed_round0_round1.{rds,csv}\n", sep="")
message("=== Done (MWI/LSO/ZWE) ===")

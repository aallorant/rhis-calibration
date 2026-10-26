#!/usr/bin/env Rscript
# =============================================================================
# 04d_temporal_lambda.R -- two-round stability of the calibration offset lambda_area =
# z_area - Z_hat_area (routine logit positivity minus survey logit undiagnosed prevalence) at
# round 0 (2015-16 surveys) and round 1 (2020-21 surveys), for areas with routine data at both
# waves: level shift mean(Delta lambda), area drift sd(Delta lambda), rank stability
# cor(lambda_r0, lambda_r1), and a composition probe against the change in area HIV prevalence.
# Routine source: PEPFAR MER HTS_TST and HTS_TST_POS (same source both waves); survey source: 02d/02e.
# Reads:  data/raw/covariates/pepfar-program-data/hts_positivity_4country.csv, data/processed/phia/*.rds
# Writes: output/temporal_lambda_{by_area_wave,paired,summary}.csv (+ .rds).  Usage: Rscript R/04d_temporal_lambda.R
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

PROJ <- "."
OUT  <- file.path(PROJ, "data/processed/phia")
HTS  <- file.path(PROJ, "data/raw/covariates/pepfar-program-data/hts_positivity_4country.csv")
RES  <- file.path(PROJ, "output"); dir.create(RES, showWarnings = FALSE, recursive = TRUE)

# ---- routine fiscal-year windows per country and wave (2-FY windows around survey fieldwork) ----
WIN <- list(
  Malawi     = list(round0 = c(2015,2016), round1 = c(2019,2020), area = "snu1_raw"),
  Lesotho    = list(round0 = c(2016,2017), round1 = c(2019,2020), area = "district"),
  Zimbabwe   = list(round0 = c(2015,2016), round1 = c(2020,2021), area = "snu1_raw"),
  Mozambique = list(round0 = c(2015,2016), round1 = c(2020,2021), area = "snu1_raw")
)
# amfAR area label -> survey-area label
AREA_FIX <- c(
  "Central Region"="Central","Northern Region"="Northern","Southern Region"="Southern",
  "Cidade De Maputo"="Maputo","Zambezia"="Zambézia")
fix_area <- function(x) ifelse(x %in% names(AREA_FIX), AREA_FIX[x], x)

# ---- 1. routine positivity z per area per wave ----
hts <- read_csv(HTS, show_col_types = FALSE) %>%
  filter(!is.na(snu1_raw), !grepl("above PSNU", snu1_raw, ignore.case = TRUE))

routine_z <- function(country) {
  cfg <- WIN[[country]]; akey <- cfg$area
  one_wave <- function(years, wave) {
    hts %>%
      filter(country == !!country, year %in% years, indicator %in% c("HTS_TST","HTS_TST_POS")) %>%
      mutate(area = fix_area(.data[[akey]])) %>%
      group_by(area, indicator) %>% summarise(v = sum(value, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = indicator, values_from = v) %>%
      filter(!is.na(HTS_TST), !is.na(HTS_TST_POS), HTS_TST_POS > 0, HTS_TST > HTS_TST_POS) %>%
      transmute(country, wave, area,
                TST = HTS_TST, POS = HTS_TST_POS,
                positivity = POS / TST,
                z = log(POS / (TST - POS)),
                se_z = sqrt(1/POS + 1/(TST - POS)))   # nominal binomial (lower bound)
  }
  bind_rows(one_wave(cfg$round0, "round0"), one_wave(cfg$round1, "round1"))
}
routine <- bind_rows(lapply(names(WIN), routine_z))

# ---- 2. survey undiagnosed Z_hat (from 02d/02e) and area HIV prevalence (composition) ----
undiag <- readRDS(file.path(OUT, "survey_undiagnosed_round0_round1.rds")) %>%
  rename(area = area)

# area HIV prevalence per wave, weighted (composition probe)
mwi_zone_to_region <- c("Northern Zone"="Northern","Central East Zone"="Central",
  "Central West Zone"="Central","South East Zone"="Southern","South West Zone"="Southern")
prev_one <- function(path, country, wave, year, collapse_mwi = FALSE) {
  d <- readRDS(path)
  if (collapse_mwi) d <- d %>% mutate(province = recode(province, !!!mwi_zone_to_region, .default = province))
  d %>% filter(!is.na(province), !is.na(btwt0), btwt0 > 0) %>%
    group_by(area = province) %>%
    summarise(prev = weighted.mean(hivstatus, btwt0, na.rm = TRUE), .groups = "drop") %>%
    mutate(country = country, wave = wave)
}
prev <- bind_rows(
  prev_one(file.path(OUT,"phia_malawi_2015.rds"),   "Malawi","round0",2015),
  prev_one(file.path(OUT,"phia_malawi.rds"),        "Malawi","round1",2020, collapse_mwi=TRUE),
  prev_one(file.path(OUT,"phia_lesotho_2016.rds"),  "Lesotho","round0",2016),
  prev_one(file.path(OUT,"phia_lesotho.rds"),       "Lesotho","round1",2020),
  prev_one(file.path(OUT,"phia_zimbabwe_2015.rds"), "Zimbabwe","round0",2015),
  prev_one(file.path(OUT,"phia_zimbabwe.rds"),      "Zimbabwe","round1",2020),
  prev_one(file.path(OUT,"phia_mozambique_2015.rds"),"Mozambique","round0",2015),
  prev_one(file.path(OUT,"phia_mozambique.rds"),    "Mozambique","round1",2021)
)

# ---- 3. lambda = z - Z_hat per area per wave ----
lambda <- undiag %>%
  inner_join(routine, by = c("country","wave","area")) %>%
  left_join(prev, by = c("country","wave","area")) %>%
  mutate(lambda = z - Z_hat,
         se_lambda = sqrt(se_z^2 + se_Z^2)) %>%
  select(country, wave, year, area, n_hiv_pos, TST, POS, positivity,
         z, se_z, p_undiag, Z_hat, se_Z, prev, lambda, se_lambda)

write_csv(lambda, file.path(RES, "temporal_lambda_by_area_wave.csv"))
saveRDS(lambda, file.path(RES, "temporal_lambda_by_area_wave.rds"))

# ---- 4. pair waves: Delta lambda, stability and composition probe ----
wide <- lambda %>%
  select(country, area, wave, lambda, se_lambda, z, Z_hat, prev) %>%
  pivot_wider(names_from = wave, values_from = c(lambda, se_lambda, z, Z_hat, prev)) %>%
  filter(!is.na(lambda_round0), !is.na(lambda_round1)) %>%
  mutate(
    d_lambda = lambda_round1 - lambda_round0,
    se_d     = sqrt(se_lambda_round0^2 + se_lambda_round1^2),
    z_score  = d_lambda / se_d,
    d_z      = z_round1 - z_round0,            # routine positivity change
    d_Zhat   = Z_hat_round1 - Z_hat_round0,    # survey undiagnosed change
    d_prev   = prev_round1 - prev_round0)      # composition change

write_csv(wide, file.path(RES, "temporal_lambda_paired.csv"))
saveRDS(wide, file.path(RES, "temporal_lambda_paired.rds"))

summ <- wide %>% group_by(country) %>% summarise(
  n_areas      = n(),
  mean_lambda0 = mean(lambda_round0),
  mean_lambda1 = mean(lambda_round1),
  mean_dlambda = mean(d_lambda),              # global level shift
  sd_dlambda   = sd(d_lambda),                # area-specific drift (key)
  cor_lambda   = if (n() > 2) cor(lambda_round0, lambda_round1) else NA_real_,
  frac_sig     = mean(abs(z_score) > 1.96),   # share of areas with sig. drift
  cor_dl_dprev = if (n() > 2) cor(d_lambda, d_prev, use="complete.obs") else NA_real_,
  .groups = "drop")

cat("\n================ lambda by area and wave ================\n")
print(as.data.frame(lambda %>% mutate(across(where(is.numeric), ~round(.,3)))), row.names = FALSE)
cat("\n================ paired Delta lambda (round1 - round0) ================\n")
print(as.data.frame(wide %>% transmute(country, area,
        lambda_r0=round(lambda_round0,2), lambda_r1=round(lambda_round1,2),
        d_lambda=round(d_lambda,2), z=round(z_score,1),
        d_z=round(d_z,2), d_Zhat=round(d_Zhat,2), d_prev=round(d_prev,3))), row.names = FALSE)
cat("\n================ stability summary by country ================\n")
print(as.data.frame(summ %>% mutate(across(where(is.numeric), ~round(.,3)))), row.names = FALSE)

cat("\nInterpretation:\n",
"  mean_dlambda = uniform (re-anchorable) shift in the offset between rounds\n",
"  sd_dlambda   = area-specific drift -- the part a single-survey calibration cannot carry over\n",
"  cor_lambda   = rank stability of the offset across rounds\n",
"  cor_dl_dprev = composition confounding probe (Delta lambda vs Delta prevalence)\n", sep="")
write_csv(summ, file.path(RES, "temporal_lambda_summary.csv"))
message("\n=== Done: temporal lambda ===")

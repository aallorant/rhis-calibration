#!/usr/bin/env Rscript
# R/03a_moz_covariates.R
# Build the Mozambique district covariate matrix x_d (159 canonical districts):
# PEPFAR priority tier and self-test volume, raw and per RHIS test in 2021.
# Reads:  data/covariates/pepfar_covariates_harmonised.csv (from 07), the demographic
#         estimation results in output/mozambique/fits/ (from 03) and the canonical
#         district list in output/crosscountry/fits/two_stage_results.rds
# Writes: output/mozambique/fits/moz_covariate_matrix.rds
# Run:    Rscript R/03a_moz_covariates.R  (from the repository root)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

PROJ_DIR <- normalizePath(".")
HARM_PATH <- file.path(PROJ_DIR, "data/covariates/pepfar_covariates_harmonised.csv")
DEMO_PATH <- file.path(PROJ_DIR, "output/mozambique/fits/estimation_demographic.rds")
CANON_PATH <- file.path(PROJ_DIR, "output/crosscountry/fits/two_stage_results.rds")
OUT_PATH  <- file.path(PROJ_DIR, "output/mozambique/fits/moz_covariate_matrix.rds")

# ---- 1. canonical district list -----------------------------------------
ts <- readRDS(CANON_PATH)
canonical_moz <- ts$Mozambique$district_pred |>
  as_tibble() |>
  distinct(province, district) |>
  arrange(province, district) |>
  rename(canonical_province = province, canonical_district = district)
stopifnot(nrow(canonical_moz) == 159L)

# ---- 2. harmonised pepfar covariates ------------------------------------
harm <- read_csv(HARM_PATH, comment = "#", show_col_types = FALSE) |>
  filter(country == "Mozambique") |>
  select(province, district, pepfar_priority_3lvl,
         hts_self_sum_2022_2024, pepfar_covered)
stopifnot(nrow(harm) == 159L)

# ---- 3. rhis test volume per district (2021 demographic data) -----------
est_dg <- readRDS(DEMO_PATH) |> as_tibble()

# Sum n_d across age x sex within each district (using demographic dataset's
# area_name + province as the join key into canonical).
rhis_totals <- est_dg |>
  filter(has_rhis) |>
  group_by(province, area_name) |>
  summarise(rhis_n_d_2021 = sum(n_d, na.rm = TRUE), .groups = "drop")

# Diagnose the join before applying it: demographic area_id count vs canonical.
n_demo <- nrow(rhis_totals)
n_canon <- nrow(canonical_moz)
cat(sprintf("demographic data districts: %d\ncanonical districts: %d\n",
            n_demo, n_canon))

# Quick exact-name check
demo_keys <- paste(rhis_totals$province, rhis_totals$area_name, sep = " | ")
canon_keys <- paste(canonical_moz$canonical_province,
                    canonical_moz$canonical_district, sep = " | ")
unmatched_canon <- setdiff(canon_keys, demo_keys)
unmatched_demo  <- setdiff(demo_keys, canon_keys)
cat(sprintf("canonical districts without a demographic match: %d\n",
            length(unmatched_canon)))
cat(sprintf("demographic districts without a canonical match: %d\n",
            length(unmatched_demo)))
if (length(unmatched_canon) > 0) {
  cat("Canonical without demo match:\n")
  for (k in unmatched_canon) cat("  -", k, "\n")
}
if (length(unmatched_demo) > 0) {
  cat("Demo without canonical match:\n")
  for (k in unmatched_demo) cat("  -", k, "\n")
}

# ---- 4. build the wide-format matrix ------------------------------------
covmat <- canonical_moz |>
  left_join(
    harm |> rename(canonical_province = province,
                   canonical_district = district),
    by = c("canonical_province", "canonical_district")
  ) |>
  left_join(
    rhis_totals |> rename(canonical_province = province,
                          canonical_district = area_name),
    by = c("canonical_province", "canonical_district")
  ) |>
  mutate(
    pepfar_priority_scaleup  = as.integer(pepfar_priority_3lvl == "scaleup"),
    pepfar_priority_notset   = as.integer(pepfar_priority_3lvl == "not_set"),
    pepfar_selftest_log_raw  = log(hts_self_sum_2022_2024 + 1),
    pepfar_selftest_log_norm = log((hts_self_sum_2022_2024 + 1) /
                                   (coalesce(rhis_n_d_2021, 0L) + 1)),
    pepfar_selftest_binary   = as.integer(hts_self_sum_2022_2024 > 0)
  )

# ---- 5. NA check + report ------------------------------------------------
covar_cols <- c("pepfar_priority_scaleup", "pepfar_priority_notset",
                "pepfar_selftest_log_raw", "pepfar_selftest_log_norm",
                "pepfar_selftest_binary")
na_check <- covmat |> select(all_of(covar_cols)) |>
  summarise(across(everything(), ~ sum(is.na(.))))
cat("\nNA counts per covariate column:\n"); print(na_check)
if (any(unlist(na_check) > 0)) {
  cat("\nDistricts with NA covariates (raw rhis_n_d_2021 missing):\n")
  print(covmat |> filter(is.na(rhis_n_d_2021)) |>
        select(canonical_province, canonical_district))
}

# ---- 6. summary report ---------------------------------------------------
cat("\n=== Step 1 report ===\n")
cat(sprintf("Districts in covariate matrix: %d\n", nrow(covmat)))

cat("\nContinuous covariates (mean, sd, range, n zeros):\n")
cont_stats <- function(x, name) {
  data.frame(
    covariate = name,
    mean = mean(x, na.rm = TRUE),
    sd   = sd(x, na.rm = TRUE),
    min  = min(x, na.rm = TRUE),
    max  = max(x, na.rm = TRUE),
    n_zero_raw_hts = sum(x == 0 & name == "pepfar_selftest_log_raw", na.rm = TRUE)
  )
}
print(rbind(
  cont_stats(covmat$pepfar_selftest_log_raw, "pepfar_selftest_log_raw"),
  cont_stats(covmat$pepfar_selftest_log_norm, "pepfar_selftest_log_norm")
))
cat(sprintf("districts with hts_self_sum_2022_2024 == 0: %d / %d\n",
            sum(covmat$hts_self_sum_2022_2024 == 0),
            nrow(covmat)))

cat("\nPriority level counts:\n")
print(covmat |> count(pepfar_priority_3lvl))

cat("\nRHIS test-volume denominator (rhis_n_d_2021) summary:\n")
print(summary(covmat$rhis_n_d_2021))

if (!dir.exists(dirname(OUT_PATH))) dir.create(dirname(OUT_PATH), recursive = TRUE)
saveRDS(list(table = covmat, covar_cols = covar_cols), OUT_PATH)
cat(sprintf("\nWrote: %s\n", OUT_PATH))

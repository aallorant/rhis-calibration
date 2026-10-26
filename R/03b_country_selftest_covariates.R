#!/usr/bin/env Rscript
# R/03b_country_selftest_covariates.R
# District covariate matrices for Lesotho, Malawi and Zimbabwe: PEPFAR priority
# tier plus self-test kits from national programme workbooks, raw and per RHIS test.
# Reads:  data/raw/rhis/{malawi,lesotho,zimbabwe}/*.xls(x) (not in the release),
#         data/covariates/pepfar_covariates_harmonised.csv (from 07), canonical
#         districts and RHIS test volumes from output/*/fits/ (see each section)
# Writes: output/{lesotho,malawi,zimbabwe}/fits/country_covariate_matrix.rds
# Run:    Rscript R/03b_country_selftest_covariates.R  (from the repository root)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(tidyr)
  library(stringr)
})

PROJ <- normalizePath(".")
HARM <- file.path(PROJ, "data/covariates/pepfar_covariates_harmonised.csv")
CANON <- file.path(PROJ, "output/crosscountry/fits/two_stage_results.rds")

# ---- canonical district keys --------------------------------------------
ts <- readRDS(CANON)
canon_country <- function(name) {
  ts[[name]]$district_pred |>
    as_tibble() |>
    distinct(province, district) |>
    arrange(province, district)
}
canon_lso <- canon_country("Lesotho")
canon_mwi <- canon_country("Malawi")
canon_zwe <- canon_country("Zimbabwe")

# ---- harmonised pepfar priorities ---------------------------------------
harm <- read_csv(HARM, comment = "#", show_col_types = FALSE) |>
  select(country, province, district,
         pepfar_priority_3lvl, hts_self_sum_2022_2024, pepfar_covered)

# ---- helper: build covariate matrix --------------------------------------
build_country_covmat <- function(country, canon, selftest_df, rhis_volume_df,
                                 join_key = c("province", "district")) {
  pep <- harm |> filter(country == !!country) |>
    select(-country)

  if (country == "Lesotho") {
    # Lesotho joins on district only (province is NA in harmonised; canonical
    # has province == district by override). Use district-only join.
    out <- canon |>
      left_join(pep |> select(-province), by = "district")
  } else {
    out <- canon |>
      left_join(pep, by = c("province", "district"))
  }

  out <- out |>
    left_join(selftest_df,
              by = setNames(c("province", "district"), join_key)) |>
    left_join(rhis_volume_df,
              by = c("province", "district"))

  # If selftest is missing for a canonical entry, treat as 0 (no kits recorded).
  out <- out |>
    mutate(
      selftest_local      = tidyr::replace_na(selftest_local, 0),
      rhis_n_d            = tidyr::replace_na(rhis_n_d, 0L),
      pepfar_priority_scaleup = as.integer(pepfar_priority_3lvl == "scaleup"),
      pepfar_priority_notset  = as.integer(pepfar_priority_3lvl == "not_set"),
      pepfar_priority_scaleup = tidyr::replace_na(pepfar_priority_scaleup, 0L),
      pepfar_priority_notset  = tidyr::replace_na(pepfar_priority_notset, 0L),
      selftest_log_raw    = log(selftest_local + 1),
      selftest_log_norm   = log((selftest_local + 1) / (rhis_n_d + 1)),
      selftest_binary     = as.integer(selftest_local > 0),
      pepfar_covered      = tidyr::replace_na(pepfar_covered, TRUE)
    )
  out
}

report_join <- function(country, canon, selftest_df,
                        canon_keys = c("province","district"),
                        st_keys = c("province","district")) {
  ck <- paste(canon[[canon_keys[1]]], canon[[canon_keys[2]]], sep = " | ")
  sk <- paste(selftest_df[[st_keys[1]]], selftest_df[[st_keys[2]]], sep = " | ")
  matched <- sum(ck %in% sk)
  pct <- 100 * matched / length(ck)
  cat(sprintf("%s: %d / %d canonical matched (%.1f%%)\n",
              country, matched, length(ck), pct))
  missing_canon <- setdiff(ck, sk)
  extra_st <- setdiff(sk, ck)
  if (length(missing_canon)) {
    cat("  canonical entries without selftest data:\n")
    for (x in missing_canon) cat("    -", x, "\n")
  }
  if (length(extra_st)) {
    cat("  selftest entries without canonical match:\n")
    for (x in extra_st) cat("    -", x, "\n")
  }
  list(matched = matched, total = length(ck), pct = pct,
       missing_canon = missing_canon, extra_st = extra_st)
}

# =========================================================================
# MWI
# =========================================================================
cat("\n===== MWI =====\n")
# MWI: kits = Self + Sex partner + Other, summed over all quarters and sub-groups
mwi_xlsx <- file.path(PROJ, "data/raw/rhis/malawi",
                      "Malawi HTS HST ART program data 2019-2022.xlsx")
mwi_raw <- read_excel(mwi_xlsx, sheet = "HST_distribution", col_types = "text")
mwi_raw$Self          <- suppressWarnings(as.numeric(mwi_raw$Self))
mwi_raw$`Sex partner` <- suppressWarnings(as.numeric(mwi_raw$`Sex partner`))
mwi_raw$Other         <- suppressWarnings(as.numeric(mwi_raw$Other))
mwi_raw$total_kits    <- rowSums(mwi_raw[, c("Self","Sex partner","Other")], na.rm = TRUE)

# Map Zone label to canonical naming (already match: 'Central East Zone' etc.)
mwi_selftest <- mwi_raw |>
  group_by(Zone, district) |>
  summarise(selftest_local = sum(total_kits, na.rm = TRUE), .groups = "drop") |>
  rename(province = Zone)

# RHIS test-volume denominator from district_predictions_demo (sex-pooled n_tests)
mwi_demo <- readRDS(file.path(PROJ, "output/malawi/fits/district_predictions_demo.rds")) |>
  as_tibble()
mwi_rhis_vol <- mwi_demo |>
  group_by(province, district) |>
  summarise(rhis_n_d = sum(n_tests, na.rm = TRUE), .groups = "drop")

mwi_join <- report_join("MWI", canon_mwi, mwi_selftest)
mwi_covmat <- build_country_covmat("Malawi", canon_mwi, mwi_selftest, mwi_rhis_vol)

# =========================================================================
# LSO
# =========================================================================
cat("\n===== LSO =====\n")
# LSO: kits = Total Distribution summed per district
lso_xls <- file.path(PROJ, "data/raw/rhis/lesotho", "HTS_TST, POS & HTS_SELF-v2.xls")
lso_raw <- read_excel(lso_xls, sheet = "FY20-FY23 Q1 HTS_SELF",
                      col_types = "text", .name_repair = "minimal")
# Drop the trailing 'Total' row (District is NA) and 'Unknown Age' + 'Total' Age rows
lso_raw <- lso_raw |>
  filter(!is.na(District), !Age %in% c("Total", "Unknown Age"))
lso_raw$total_dist <- suppressWarnings(as.numeric(lso_raw$`Total Distribution`))

lso_selftest <- lso_raw |>
  group_by(District) |>
  summarise(selftest_local = sum(total_dist, na.rm = TRUE), .groups = "drop") |>
  rename(district = District) |>
  mutate(province = district)   # LSO canonical has province == district

# RHIS test volume denominator
lso_demo <- readRDS(file.path(PROJ, "output/lesotho/fits/district_predictions_demo.rds")) |>
  as_tibble()
lso_rhis_vol <- lso_demo |>
  group_by(province, district) |>
  summarise(rhis_n_d = sum(n_tests, na.rm = TRUE), .groups = "drop")

lso_join <- report_join("LSO", canon_lso, lso_selftest,
                        canon_keys = c("province","district"),
                        st_keys = c("province","district"))
lso_covmat <- build_country_covmat("Lesotho", canon_lso, lso_selftest, lso_rhis_vol)

# =========================================================================
# ZWE
# =========================================================================
cat("\n===== ZWE =====\n")
zwe_xlsx <- file.path(PROJ, "data/raw/rhis/zimbabwe",
                      "HIVST Distribution Data 2019-2022.mn.xlsx")

# Load each sheet; the 2022 sheet carries a two-letter province-code prefix
# ('bu ', 'ha ', ...) on Province and District, stripped here.
read_zwe_sheet <- function(sheet_name) {
  df <- read_excel(zwe_xlsx, sheet = sheet_name,
                   col_types = "text", .name_repair = "minimal")
  if (sheet_name == "2022") {
    df$Province <- sub("^[a-z]{2} ", "", df$Province)
    df$District <- sub("^[a-z]{2} ", "", df$District)
  }
  df
}
zwe_2019_2020 <- read_zwe_sheet("20219 and 2020")
zwe_2021      <- read_zwe_sheet("2021")
zwe_2022      <- read_zwe_sheet("2022")
zwe_raw <- bind_rows(zwe_2019_2020, zwe_2021, zwe_2022)

# Roll up hospital and city labels into their canonical province (as in 07).
zwe_raw$Province <- dplyr::recode(zwe_raw$Province,
  "Bulawayo Central Hospitals" = "Bulawayo",
  "Harare Central Hospitals"   = "Harare",
  "Chitungwiza City"           = "Harare"
)
# The canonical district list keeps the older name of this hospital.
zwe_raw$District <- dplyr::recode(zwe_raw$District,
  "Sally Mugabe Central Hospital" = "Harare Central Hospital"
)

# ZWE: kits = the 16 age x sex columns summed over all quarters 2019-2022
age_sex_cols <- grep(",", names(zwe_raw), value = TRUE)
stopifnot(length(age_sex_cols) == 16L)
zwe_raw[, age_sex_cols] <- lapply(zwe_raw[, age_sex_cols], function(x)
  suppressWarnings(as.numeric(x)))
zwe_raw$total_kits <- rowSums(zwe_raw[, age_sex_cols], na.rm = TRUE)

zwe_selftest <- zwe_raw |>
  group_by(Province, District) |>
  summarise(selftest_local = sum(total_kits, na.rm = TRUE), .groups = "drop") |>
  rename(province = Province, district = District)

# RHIS test-volume denominator
zwe_demo <- readRDS(file.path(PROJ, "output/zimbabwe/fits/district_predictions_demo.rds")) |>
  as_tibble()
zwe_rhis_vol <- zwe_demo |>
  group_by(province, district) |>
  summarise(rhis_n_d = sum(n_tests, na.rm = TRUE), .groups = "drop")

zwe_join <- report_join("ZWE", canon_zwe, zwe_selftest)
zwe_covmat <- build_country_covmat("Zimbabwe", canon_zwe, zwe_selftest, zwe_rhis_vol)

# =========================================================================
# 10% mismatch gate
# =========================================================================
fails <- c()
for (j in list(mwi_join, lso_join, zwe_join)) {
  if (j$pct < 90) fails <- c(fails, sprintf("%.1f%%", j$pct))
}
if (length(fails)) {
  cat("\n*** ONE OR MORE COUNTRIES MISSED 90% MATCH THRESHOLD ***\n")
  cat("Match rates below 90%: ", paste(fails, collapse = ", "), "\n")
  cat("Stopping per spec.\n")
  quit(status = 1)
}

# =========================================================================
# Persist outputs
# =========================================================================
write_country <- function(country, covmat, key = c("province","district"),
                          fit_dir) {
  out <- list(
    table = covmat,
    covar_cols = c("pepfar_priority_scaleup", "pepfar_priority_notset",
                   "selftest_log_raw", "selftest_log_norm", "selftest_binary",
                   "pepfar_covered")
  )
  saveRDS(out, file.path(PROJ, "output", tolower(country), "fits",
                         "country_covariate_matrix.rds"))
}
write_country("Lesotho",  lso_covmat)
write_country("Malawi",   mwi_covmat)
write_country("Zimbabwe", zwe_covmat)

cat("\n===== Summary =====\n")
for (sym in c("lso_covmat","mwi_covmat","zwe_covmat")) {
  cm <- get(sym)
  ctry <- sub("_covmat", "", sym)
  cat(sprintf("\n%s (%d districts)\n", toupper(ctry), nrow(cm)))
  cat(sprintf("  selftest_local: median=%.0f  max=%.0f  zeros=%d\n",
              median(cm$selftest_local), max(cm$selftest_local),
              sum(cm$selftest_local == 0)))
  cat(sprintf("  rhis_n_d:       median=%.0f  max=%.0f  zeros=%d\n",
              median(cm$rhis_n_d), max(cm$rhis_n_d), sum(cm$rhis_n_d == 0)))
  cat("  priority:\n"); print(table(cm$pepfar_priority_3lvl, useNA="ifany"))
  if ("pepfar_covered" %in% names(cm)) {
    cat("  pepfar_covered:\n"); print(table(cm$pepfar_covered, useNA="ifany"))
  }
}

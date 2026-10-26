# 03_calibration_model_moz.R: Mozambique district calibration of routine HIV-testing
# positivity against the 2021 survey. Survey estimates (district + province), RHIS
# log-odds, district crosswalk, estimation dataset, descriptive figures, weighted
# regression, INLA fits (fixed / iid / ICAR / BYM2), age x sex extension, stacked model.
# Reads:  data/all_moz_data.rda (or ./all_moz_data.rda), output/df_long_ajuda.rds,
#         output/phia_insida2021_with_awareness.rds
# Writes: output/*.rds, output/district_crosswalk.csv, output/figures/*.png and *.pdf
# Run:    Rscript R/03_calibration_model_moz.R   (from the repository root)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(sf)
  library(stringr)
  library(lubridate)
  library(patchwork)
})

# ---- Configuration ----
out_dir <- "output"
fig_dir <- file.path(out_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

age_breaks <- c(15, 20, 25, 50, Inf)
age_labels <- c("15-19", "20-24", "25-49", "50+")

# assumed survey design effect
deff_default <- 2

# Toggle INLA fits
run_inla <- TRUE

# Manual overrides for district name mismatches (normalized names)
manual_district_map <- tibble::tibble(
  district_norm = c(
    "xai xai extinto",
    "nacala porto",
    "chimbonila",
    "doa",
    "kamaxakene",
    "kampfumo",
    "murrupula",
    "nacaroa",
    "namarroi",
    "nlhamankulu",
    "sussundenga",
    "vilankulos",
    "zumbo"
  ),
  area_name = c(
    "Xai-Xai",
    "Nacala",
    "Chimbunila",
    "Dôa",
    "KaMaxakeni",
    "KaMpfumu",
    "Morrupula",
    "Nacarôa",
    "Namarói",
    "Nhlamankulu",
    "Sussendenga",
    "Vilankulo",
    "Zumbu"
  )
)

# ---- Helpers ----
normalize_name <- function(x) {
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[`']", "", x)
  x <- gsub("[^a-z0-9 ]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x <- gsub("^cidade\\s+(de|da|do)\\s+", "", x)
  x <- gsub("^distrito\\s+municipal\\s+(de|da|do)\\s+", "", x)
  x <- gsub("^distrito\\s+(de|da|do)\\s+", "", x)
  x <- gsub("^municipio\\s+(de|da|do)\\s+", "", x)
  x <- trimws(x)
  x
}

safe_qlogis <- function(p, n_eff) {
  # Guard against 0/1 using a small correction based on n_eff
  n_eff <- ifelse(is.na(n_eff) | n_eff <= 0, NA_real_, n_eff)
  eps <- ifelse(is.na(n_eff), 1e-6, pmin(0.5 / n_eff, 1e-3))
  p_adj <- pmin(pmax(p, eps), 1 - eps)
  qlogis(p_adj)
}

safe_log_odds_counts <- function(pos, neg) {
  # Continuity correction for zeros
  pos_adj <- pos + 0.5
  neg_adj <- neg + 0.5
  log(pos_adj / neg_adj)
}

# ---- 1) Load data ----

message("=== Loading data ===")

rda_candidates <- c("all_moz_data.rda", file.path("data", "all_moz_data.rda"))
rda_path <- rda_candidates[file.exists(rda_candidates)][1]
if (is.na(rda_path)) {
  stop("Could not find all_moz_data.rda in working directory or data/. Please place the file and re-run.")
}

load(rda_path) # expects area_merged

if (!exists("area_merged")) stop("area_merged not found in all_moz_data.rda")

# Convert area_merged to sf (if geometry exists)
if ("geometry" %in% names(area_merged)) {
  area_merged_sf <- sf::st_as_sf(area_merged, sf_column_name = "geometry", crs = 4326)
} else {
  area_merged_sf <- NULL
}

area_merged_df <- if (!is.null(area_merged_sf)) sf::st_drop_geometry(area_merged_sf) else area_merged

# Load RHIS data
rhis_path <- file.path(out_dir, "df_long_ajuda.rds")
if (!file.exists(rhis_path)) stop("Missing RHIS file: output/df_long_ajuda.rds")

df_long <- readRDS(rhis_path)

# ---- 2) Prepare PHIA data + awareness ----

message("=== Preparing PHIA data ===")

phia_path <- file.path(out_dir, "phia_insida2021_with_awareness.rds")
if (!file.exists(phia_path)) stop("Missing PHIA file: output/phia_insida2021_with_awareness.rds")

phia_moz <- readRDS(phia_path)

cat("PHIA columns:\n")
cat(paste(names(phia_moz), collapse = "\n"), "\n")

stopifnot("aware_bin" %in% names(phia_moz))
stopifnot("area_id" %in% names(phia_moz))
stopifnot("indweight" %in% names(phia_moz))

phia_moz <- phia_moz %>% mutate(aware = aware_bin)

# Harmonize demographics
phia_moz <- phia_moz %>%
  mutate(
    age_group = cut(age, breaks = age_breaks, labels = age_labels, right = FALSE),
    sex_label = ifelse(sex == 1, "Male", "Female")
  )

# ---- 3) Compute Z_hat_d and psi_d (district + province) ----

message("=== Computing survey estimates (district + province) ===")

# Add province and district names from area_merged
provinces <- area_merged_df %>%
  filter(area_level == 1) %>%
  select(province_id = area_id, province = area_name)

districts <- area_merged_df %>%
  filter(area_level == 2) %>%
  select(area_id, area_name, parent_area_id)

districts <- districts %>%
  left_join(provinces, by = c("parent_area_id" = "province_id"))

needs_join <- !("area_name" %in% names(phia_moz) && "province" %in% names(phia_moz))
if (needs_join) {
  phia_moz <- phia_moz %>%
    left_join(districts %>% select(area_id, area_name, province), by = "area_id")
}

# Standardize area/province columns if duplicated by joins
if (!"area_name" %in% names(phia_moz)) {
  if (all(c("area_name.x", "area_name.y") %in% names(phia_moz))) {
    phia_moz$area_name <- dplyr::coalesce(phia_moz$area_name.x, phia_moz$area_name.y)
  } else if ("area_name.x" %in% names(phia_moz)) {
    phia_moz$area_name <- phia_moz$area_name.x
  } else if ("area_name.y" %in% names(phia_moz)) {
    phia_moz$area_name <- phia_moz$area_name.y
  }
}

if (!"province" %in% names(phia_moz)) {
  if ("province_name" %in% names(phia_moz)) {
    phia_moz$province <- phia_moz$province_name
  } else if (all(c("province.x", "province.y") %in% names(phia_moz))) {
    phia_moz$province <- dplyr::coalesce(phia_moz$province.x, phia_moz$province.y)
  } else if ("province.x" %in% names(phia_moz)) {
    phia_moz$province <- phia_moz$province.x
  } else if ("province.y" %in% names(phia_moz)) {
    phia_moz$province <- phia_moz$province.y
  }
}

# Exclude diagnosed positives from target population
phia_undx <- phia_moz %>%
  filter(!(hivstatus == 1 & aware == TRUE))

# District level
phia_district <- phia_undx %>%
  group_by(area_id, area_name, province_name) %>%
  summarise(
    n = n(),
    w_total = sum(indweight, na.rm = TRUE),
    w_pos = sum(indweight * hivstatus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(province = province_name) %>%
  mutate(
    p_hat_raw = ifelse(w_total > 0, w_pos / w_total, NA_real_),
    n_eff = n / deff_default,
    p_hat_var = ifelse(!is.na(n_eff),
                       pmin(pmax(p_hat_raw, 0.5 / n_eff), 1 - 0.5 / n_eff),
                       NA_real_),
    Z_hat = safe_qlogis(p_hat_raw, n_eff),
    psi = ifelse(!is.na(n_eff) & !is.na(p_hat_var),
                 1 / (n_eff * p_hat_var * (1 - p_hat_var)), NA_real_),
    p_hat_raw_is_boundary = (w_pos == 0) | (abs(w_pos - w_total) < 0.001)
  )

# Province level
phia_province <- phia_undx %>%
  group_by(province_name) %>%
  summarise(
    n = n(),
    w_total = sum(indweight, na.rm = TRUE),
    w_pos = sum(indweight * hivstatus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(province = province_name) %>%
  mutate(
    p_hat_raw = ifelse(w_total > 0, w_pos / w_total, NA_real_),
    n_eff = n / deff_default,
    p_hat_var = ifelse(!is.na(n_eff),
                       pmin(pmax(p_hat_raw, 0.5 / n_eff), 1 - 0.5 / n_eff),
                       NA_real_),
    Z_hat = safe_qlogis(p_hat_raw, n_eff),
    psi = ifelse(!is.na(n_eff) & !is.na(p_hat_var),
                 1 / (n_eff * p_hat_var * (1 - p_hat_var)), NA_real_)
  )

saveRDS(phia_district, file.path(out_dir, "phia_district_estimates.rds"))
saveRDS(phia_province, file.path(out_dir, "phia_province_estimates.rds"))

# ---- 4) Compute RHIS z_d (district + optional finer resolutions) ----

message("=== Computing RHIS log-odds (2021) ===")

# Ensure year is available
if (!"year" %in% names(df_long)) {
  if ("period" %in% names(df_long)) {
    df_long <- df_long %>% mutate(year = lubridate::year(period))
  } else {
    stop("df_long_ajuda.rds must contain 'year' or 'period'.")
  }
}

rhis_2021 <- df_long %>%
  filter(
    indicator %in% c("hts_tests", "hts_positive"),
    age_group %in% age_labels,
    sex %in% c("Male", "Female"),
    year == 2021
  )

# Aggregate to district
rhis_district <- rhis_2021 %>%
  group_by(district, indicator) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = indicator, values_from = value, values_fill = 0) %>%
  rename(n_d = hts_tests, ny_d = hts_positive) %>%
  mutate(
    n_neg_d = n_d - ny_d,
    z_d = ifelse(n_d > 0, safe_log_odds_counts(ny_d, n_neg_d), NA_real_),
    v_d = ifelse(ny_d > 0 & n_neg_d > 0, 1 / ny_d + 1 / n_neg_d, NA_real_)
  )

# District x age x sex
rhis_district_age_sex <- rhis_2021 %>%
  group_by(district, age_group, sex, indicator) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = indicator, values_from = value, values_fill = 0) %>%
  rename(n_d = hts_tests, ny_d = hts_positive) %>%
  mutate(
    n_neg_d = n_d - ny_d,
    z_d = ifelse(n_d > 0, safe_log_odds_counts(ny_d, n_neg_d), NA_real_),
    v_d = ifelse(ny_d > 0 & n_neg_d > 0, 1 / ny_d + 1 / n_neg_d, NA_real_)
  )

# District x month
rhis_district_month <- rhis_2021 %>%
  group_by(district, period, indicator) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = indicator, values_from = value, values_fill = 0) %>%
  rename(n_d = hts_tests, ny_d = hts_positive) %>%
  mutate(
    n_neg_d = n_d - ny_d,
    z_d = ifelse(n_d > 0, safe_log_odds_counts(ny_d, n_neg_d), NA_real_),
    v_d = ifelse(ny_d > 0 & n_neg_d > 0, 1 / ny_d + 1 / n_neg_d, NA_real_)
  )

saveRDS(rhis_district, file.path(out_dir, "rhis_district_2021.rds"))
saveRDS(rhis_district_age_sex, file.path(out_dir, "rhis_district_age_sex_2021.rds"))
saveRDS(rhis_district_month, file.path(out_dir, "rhis_district_month_2021.rds"))

# ---- 5) Build district crosswalk (RHIS district -> area_id) ----

message("=== Building district crosswalk ===")

rhis_districts <- rhis_district %>%
  distinct(district) %>%
  mutate(district_norm = normalize_name(district))

area_districts <- districts %>%
  mutate(area_name_norm = normalize_name(area_name))

# Manual overrides first (join on normalized area_name for robustness)
manual_lookup <- manual_district_map %>%
  mutate(area_name_norm = normalize_name(area_name)) %>%
  left_join(area_districts %>% select(area_id, area_name, province, area_name_norm),
            by = "area_name_norm")
if ("area_name.y" %in% names(manual_lookup) || "area_name.x" %in% names(manual_lookup)) {
  manual_lookup <- manual_lookup %>%
    mutate(area_name = dplyr::coalesce(.data$area_name.y, .data$area_name.x))
}
if (!"province" %in% names(manual_lookup) && "province.y" %in% names(manual_lookup)) {
  manual_lookup <- manual_lookup %>%
    mutate(province = .data$province.y)
}

crosswalk <- rhis_districts %>%
  left_join(manual_lookup %>% select(district_norm, area_id, area_name, province),
            by = "district_norm") %>%
  left_join(area_districts, by = c("district_norm" = "area_name_norm"),
            suffix = c("_manual", ""))

crosswalk <- crosswalk %>%
  mutate(
    area_id = coalesce(area_id_manual, area_id),
    area_name = coalesce(area_name_manual, area_name),
    province = coalesce(province_manual, province)
  ) %>%
  select(district, district_norm, area_id, area_name, province)

# Fuzzy match for unmatched (do not auto-accept)
unmatched <- crosswalk %>% filter(is.na(area_id))
review_path <- file.path(out_dir, "district_crosswalk_NEEDS_REVIEW.csv")

if (nrow(unmatched) > 0) {
  message("  Unmatched districts (exact): ", nrow(unmatched), " -> computing fuzzy suggestions")
  adist_mat <- adist(unmatched$district_norm, area_districts$area_name_norm)

  top_k <- 3
  suggestions <- lapply(seq_len(nrow(adist_mat)), function(i) {
    ord <- order(adist_mat[i, ])
    data.frame(
      rhis_district = unmatched$district[i],
      rhis_district_normalized = unmatched$district_norm[i],
      suggested_area_id = area_districts$area_id[ord[1]],
      suggested_area_name = area_districts$area_name[ord[1]],
      edit_distance = adist_mat[i, ord[1]],
      suggestion_2 = area_districts$area_name[ord[2]],
      suggestion_2_dist = adist_mat[i, ord[2]],
      suggestion_3 = area_districts$area_name[ord[3]],
      suggestion_3_dist = adist_mat[i, ord[3]],
      accepted = FALSE,
      stringsAsFactors = FALSE
    )
  })

  review_df <- bind_rows(suggestions)
  write.csv(review_df, review_path, row.names = FALSE)
  message("  Review required: output/district_crosswalk_NEEDS_REVIEW.csv")
}

crosswalk <- crosswalk %>%
  select(district, district_norm, area_id, area_name, province)

write.csv(crosswalk, file.path(out_dir, "district_crosswalk.csv"), row.names = FALSE)

# Crosswalk quality report
total_tests <- sum(rhis_district$n_d, na.rm = TRUE)
matched <- crosswalk %>% filter(!is.na(area_id))
unmatched <- crosswalk %>% filter(is.na(area_id))
matched_tests <- rhis_district %>%
  filter(district %in% matched$district) %>%
  summarise(n = sum(n_d, na.rm = TRUE)) %>%
  pull(n)

cat("\n=== District crosswalk summary ===\n")
cat("Total RHIS districts:", nrow(crosswalk), "\n")
cat("Matched:", nrow(matched),
    sprintf("(%.1f%% of districts, %.1f%% of tests)\n",
            100 * nrow(matched) / nrow(crosswalk),
            100 * matched_tests / total_tests))
cat("Unmatched:", nrow(unmatched),
    sprintf("(%.1f%% of districts, %.1f%% of tests)\n",
            100 * nrow(unmatched) / nrow(crosswalk),
            100 * (total_tests - matched_tests) / total_tests))

if (nrow(unmatched) > 0) {
  cat("\nUnmatched districts:\n")
  for (i in seq_len(nrow(unmatched))) {
    cat("  -", unmatched$district[i], "\n")
  }
  cat("\n*** Review output/district_crosswalk_NEEDS_REVIEW.csv ***\n")
  cat("*** Add corrections to manual_district_map at top of script ***\n\n")
}

# ---- 6) Assemble estimation dataset (full scaffold) ----

message("=== Assembling estimation dataset ===")

rhis_district_xw <- rhis_district %>%
  left_join(crosswalk, by = "district") %>%
  mutate(log_tests = ifelse(n_d > 0, log(n_d), NA_real_))

# Report duplicated area_id mappings (expected when RHIS names collapse)
dup_area_ids <- crosswalk %>%
  filter(!is.na(area_id)) %>%
  filter(duplicated(area_id) | duplicated(area_id, fromLast = TRUE))
if (nrow(dup_area_ids) > 0) {
  cat("\n=== Duplicated area_id mappings in crosswalk ===\n")
  print(dup_area_ids)
}

# Aggregate RHIS to area_id to avoid duplicate joins
rhis_by_area <- rhis_district_xw %>%
  filter(!is.na(area_id)) %>%
  group_by(area_id) %>%
  summarise(
    n_d = sum(n_d, na.rm = TRUE),
    ny_d = sum(ny_d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_neg_d = n_d - ny_d,
    z_d = ifelse(n_d > 0, safe_log_odds_counts(ny_d, n_neg_d), NA_real_),
    v_d = ifelse(ny_d > 0 & n_neg_d > 0, 1 / ny_d + 1 / n_neg_d, NA_real_),
    log_tests = ifelse(n_d > 0, log(n_d), NA_real_)
  )

# Build full district scaffold for INLA alignment
if (!is.null(area_merged_sf)) {
  all_districts_sf <- area_merged_sf %>%
    filter(area_level == 2) %>%
    st_make_valid()
} else {
  stop("area_merged has no geometry; required for INLA scaffold")
}

scaffold <- all_districts_sf %>%
  st_drop_geometry() %>%
  select(area_id, area_name, parent_area_id) %>%
  mutate(district_idx = row_number())

provinces_df <- st_drop_geometry(area_merged_sf) %>%
  filter(area_level == 1) %>%
  select(province_id = area_id, province = area_name)

scaffold <- scaffold %>%
  left_join(provinces_df, by = c("parent_area_id" = "province_id"))

est_dat <- scaffold %>%
  left_join(
    rhis_by_area %>% select(area_id, n_d, ny_d, n_neg_d, z_d, v_d, log_tests),
    by = "area_id"
  ) %>%
  left_join(
    phia_district %>% select(area_id, Z_hat = Z_hat, psi = psi, n_phia = n, p_hat_raw_is_boundary),
    by = "area_id"
  ) %>%
  rename(Z_hat_d = Z_hat, psi_d = psi) %>%
  mutate(
    has_rhis = !is.na(z_d),
    # Only include informative PHIA districts (exclude boundary p_hat_raw)
    in_Ds = !is.na(Z_hat_d) & !is.na(psi_d) & is.finite(psi_d) & !p_hat_raw_is_boundary,
    has_both = has_rhis & in_Ds,
    y = ifelse(has_both, z_d - Z_hat_d, NA_real_),
    prec = ifelse(!is.na(y) & !is.na(v_d) & !is.na(psi_d),
                  1 / (v_d + psi_d), NA_real_),
    prec_inla = ifelse(!is.na(y), 1 / (v_d + psi_d), 1)
  )

cat("\n=== Estimation dataset ===\n")
cat("Total districts (scaffold):", nrow(est_dat), "\n")
cat("  With RHIS data:", sum(est_dat$has_rhis), "\n")
cat("  With PHIA data (informative):", sum(est_dat$in_Ds), "\n")
cat("  With both (usable for calibration):", sum(est_dat$has_both), "\n")
cat("  RHIS only (predict using model):", sum(est_dat$has_rhis & !est_dat$in_Ds), "\n")
cat("  Neither (predict from spatial prior only):", sum(!est_dat$has_rhis & !est_dat$in_Ds), "\n")
cat("  Excluded from calibration (PHIA boundary p_hat):",
    sum(!is.na(est_dat$Z_hat_d) & !est_dat$in_Ds), "\n")

saveRDS(est_dat, file.path(out_dir, "estimation_dataset.rds"))

# Province-level fallback for excluded districts
excluded_phia <- est_dat %>%
  filter(!is.na(Z_hat_d) & !in_Ds) %>%
  left_join(
    phia_province %>% select(province, Z_hat_prov = Z_hat, psi_prov = psi),
    by = "province"
  )

if (nrow(excluded_phia) > 0) {
  cat("\n=== Province-level fallback for excluded districts ===\n")
  print(excluded_phia %>% select(area_name, area_id, province, Z_hat_d, Z_hat_prov, psi_prov))
  saveRDS(excluded_phia, file.path(out_dir, "phia_excluded_with_province_fallback.rds"))
}

# ---- 7) Descriptive analysis ----

message("=== Descriptive analysis (sampled districts) ===")

sampled <- est_dat %>% filter(in_Ds)

sampled <- sampled %>%
  mutate(log_lambda_hat = z_d - Z_hat_d)

saveRDS(sampled, file.path(out_dir, "sampled_districts.rds"))

# Plotting dataset (drop non-finite values)
sampled_plot <- sampled %>%
  filter(is.finite(log_lambda_hat), is.finite(Z_hat_d), is.finite(z_d))

# 1) Histogram of log_lambda_hat
p_hist <- ggplot(sampled_plot, aes(x = log_lambda_hat)) +
  geom_histogram(bins = 30, fill = "#4C78A8", color = "white") +
  labs(x = "log-RTP (z_d - Z_hat_d)", y = "Count",
       title = "Distribution of observed log-RTP") +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "log_rtp_hist.png"), p_hist, width = 7, height = 4)
ggsave(file.path(fig_dir, "log_rtp_hist.pdf"), p_hist, width = 7, height = 4)

# 2) Scatter Z_hat vs z_d
p_scatter <- ggplot(sampled_plot, aes(x = Z_hat_d, y = z_d)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(x = "Survey log-odds (Z_hat_d)", y = "RHIS log-odds (z_d)",
       title = "Survey vs RHIS log-odds") +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "z_hat_vs_z_d.png"), p_scatter, width = 6, height = 5)
ggsave(file.path(fig_dir, "z_hat_vs_z_d.pdf"), p_scatter, width = 6, height = 5)

# 3) log_lambda_hat vs Z_hat_d
p_scatter2 <- ggplot(sampled_plot, aes(x = Z_hat_d, y = log_lambda_hat)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(x = "Survey log-odds (Z_hat_d)", y = "log-RTP (z_d - Z_hat_d)",
       title = "Selection bias vs prevalence") +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "log_rtp_vs_z_hat.png"), p_scatter2, width = 6, height = 5)
ggsave(file.path(fig_dir, "log_rtp_vs_z_hat.pdf"), p_scatter2, width = 6, height = 5)

# 4) Map of log_lambda_hat (if geometry available)
if (!is.null(area_merged_sf)) {
  districts_sf <- area_merged_sf %>% filter(area_level == 2)
  map_df <- districts_sf %>%
    left_join(sampled %>% select(area_id, log_lambda_hat), by = "area_id")

  p_map <- ggplot(map_df) +
    geom_sf(aes(fill = log_lambda_hat), color = "white", size = 0.1) +
    scale_fill_viridis_c(option = "C", na.value = "gray90") +
    labs(fill = "log-RTP", title = "Observed log-RTP by district") +
    theme_void(base_size = 12)

  ggsave(file.path(fig_dir, "log_rtp_map.png"), p_map, width = 7, height = 6)
  ggsave(file.path(fig_dir, "log_rtp_map.pdf"), p_map, width = 7, height = 6)
}

# ---- 8) Weighted regression ----

message("=== Weighted regression (no measurement error model) ===")

sampled <- sampled %>%
  mutate(weight = 1 / (v_d + psi_d))

# Simple model: log_lambda_hat ~ province + log_tests
# If province is missing, drop those rows
lm_data <- sampled %>% filter(!is.na(province), !is.na(weight))

if (nrow(lm_data) > 0) {
  fit_lm <- lm(log_lambda_hat ~ province + log_tests, data = lm_data, weights = weight)
  saveRDS(fit_lm, file.path(out_dir, "wls_fit.rds"))
}

# ---- 9) INLA scaffolding (optional; run if INLA installed) ----

if (run_inla && requireNamespace("INLA", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(INLA)
    library(spdep)
  })

  if (!is.null(area_merged_sf)) {
    all_districts_sf <- area_merged_sf %>%
      filter(area_level == 2) %>%
      st_make_valid()

    D <- nrow(all_districts_sf)
    cat("Total districts in Mozambique:", D, "\n")

    inla_ok <- TRUE
    adj_graph <- NULL
    nb <- NULL
    tryCatch({
      nb <- spdep::poly2nb(all_districts_sf, queen = TRUE)
      island_idx <- which(spdep::card(nb) == 0)
      if (length(island_idx) > 0) {
        cat("WARNING:", length(island_idx), "island districts with no neighbours:\n")
        cat(paste(" ", all_districts_sf$area_name[island_idx], collapse = "\n"), "\n")
        for (i in island_idx) {
          other_idx <- setdiff(seq_len(D), i)
          dists <- sf::st_distance(
            sf::st_centroid(all_districts_sf[i, ]),
            sf::st_centroid(all_districts_sf[other_idx, ])
          )
          nearest_other <- other_idx[which.min(dists)]
          nb[[i]] <- as.integer(nearest_other)
          nb[[nearest_other]] <- sort(unique(c(nb[[nearest_other]], as.integer(i))))
        }
        cat("Islands connected to nearest neighbour.\n")
      }

      cat("Neighbours per district: min =", min(spdep::card(nb)),
          ", median =", median(spdep::card(nb)),
          ", max =", max(spdep::card(nb)), "\n")

      inla_graph_path <- file.path(out_dir, "district_adj.graph")
      spdep::nb2INLA(inla_graph_path, nb)
      adj_graph <- INLA::inla.read.graph(inla_graph_path)
      stopifnot(adj_graph$n == D)
    }, error = function(e) {
      message("INLA adjacency graph failed: ", e$message)
      inla_ok <<- FALSE
    })

    # Scaffold alignment checks
    stopifnot(nrow(est_dat) == D)
    stopifnot(all(est_dat$district_idx == 1:D))
    stopifnot(!any(duplicated(est_dat$area_id)))

    est_dat <- est_dat %>%
      mutate(
        prec_inla = case_when(
          is.na(y) ~ 1,
          is.na(prec) | prec <= 0 | is.infinite(prec) ~ NA_real_,
          TRUE ~ prec
        ),
        y = ifelse(is.na(prec_inla), NA_real_, y),
        prec_inla = ifelse(is.na(y), 1, prec_inla)
      )

    dropped <- est_dat %>% filter(has_both & is.na(y))
    if (nrow(dropped) > 0) {
      cat("Districts dropped due to invalid precision:\n")
      print(dropped %>% select(area_name, ny_d, n_neg_d, v_d, n_phia, psi_d))
    }
    est_dat$province_factor <- as.factor(ifelse(is.na(est_dat$province), "Unknown", est_dat$province))
    est_dat$district_idx2 <- est_dat$district_idx

    if (any(is.na(est_dat$prec_inla))) {
      cat("ERROR: prec_inla still has NA values. Rows:\n")
      print(est_dat %>% filter(is.na(prec_inla)) %>% select(area_name, y, prec, prec_inla))
      stop("prec_inla contains NA values")
    }

    if (inla_ok) {
      tryCatch({
        formula_1a <- y ~ 1 + province_factor
        fit_1a <- INLA::inla(
          formula_1a,
          family = "gaussian",
          data = est_dat,
          scale = est_dat$prec_inla,
          control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
          control.predictor = list(compute = TRUE),
          control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
        )

        # Diagnostic: QQ plot after calibration (fit_1a)
        corrected_z <- sampled$z_d - fit_1a$summary.fitted.values$mean[sampled$district_idx]
        ks_res <- ks.test(corrected_z, sampled$Z_hat_d)

        qq_probs <- seq(0.01, 0.99, by = 0.01)
        qq_df <- data.frame(
          q = qq_probs,
          x = as.numeric(quantile(sampled$Z_hat_d, qq_probs, na.rm = TRUE)),
          y = as.numeric(quantile(corrected_z, qq_probs, na.rm = TRUE))
        )

        p_scatter_raw <- ggplot(sampled_plot, aes(x = Z_hat_d, y = z_d)) +
          geom_point(alpha = 0.7) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
          labs(x = "Survey log-odds (Z_hat_d)", y = "RHIS log-odds (z_d)",
               title = "A. Raw RHIS vs survey") +
          theme_minimal(base_size = 12)

        p_qq <- ggplot(qq_df, aes(x = x, y = y)) +
          geom_point(alpha = 0.7) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
          labs(x = "Survey log-odds quantiles", y = "Corrected RHIS quantiles",
               title = "B. Post-calibration QQ plot") +
          theme_minimal(base_size = 12)

        fig_diag <- p_scatter_raw + p_qq +
          plot_annotation(
            title = "Diagnostic: RHIS vs PHIA",
            caption = paste0(
              "Left: raw RHIS vs survey log-odds. Right: QQ plot after calibration using fit_1a. ",
              "KS test p = ", signif(ks_res$p.value, 3), "."
            )
          )

        ggsave(file.path(fig_dir, "fig_diagnostics.png"), fig_diag, width = 10, height = 4.5)
        ggsave(file.path(fig_dir, "fig_diagnostics.pdf"), fig_diag, width = 10, height = 4.5)

        formula_1b <- y ~ 1 + province_factor + f(district_idx, model = "iid")
        fit_1b <- INLA::inla(
          formula_1b,
          family = "gaussian",
          data = est_dat,
          scale = est_dat$prec_inla,
          control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
          control.predictor = list(compute = TRUE),
          control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
        )

        formula_1c <- y ~ 1 + province_factor + f(district_idx, model = "besag", graph = adj_graph)
        fit_1c <- INLA::inla(
          formula_1c,
          family = "gaussian",
          data = est_dat,
          scale = est_dat$prec_inla,
          control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
          control.predictor = list(compute = TRUE),
          control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
        )

        fit_1d <- NULL
        # BYM2 with PC priors
        formula_1d <- y ~ 1 + province_factor +
          f(district_idx, model = "bym2", graph = adj_graph,
            hyper = list(
              prec = list(prior = "pc.prec", param = c(1, 0.01)),
              phi = list(prior = "pc", param = c(0.5, 0.5))
            ))
        fit_1d <- tryCatch({
          INLA::inla(
            formula_1d,
            family = "gaussian",
            data = est_dat,
            scale = est_dat$prec_inla,
            control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
            control.predictor = list(compute = TRUE),
            control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
          )
        }, error = function(e) e)

        # fallback: adaptive integration strategy
        if (inherits(fit_1d, "error")) {
          fit_1d <- tryCatch({
            INLA::inla(
              formula_1d,
              family = "gaussian",
              data = est_dat,
              scale = est_dat$prec_inla,
              control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
              control.inla = list(strategy = "adaptive", int.strategy = "eb"),
              control.predictor = list(compute = TRUE),
              control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
            )
          }, error = function(e) e)
        }

        # fallback: Besag + iid if BYM2 fails
        if (inherits(fit_1d, "error")) {
          formula_1d_alt <- y ~ 1 + province_factor +
            f(district_idx, model = "besag", graph = adj_graph) +
            f(district_idx2, model = "iid")
          fit_1d <- tryCatch({
            INLA::inla(
              formula_1d_alt,
              family = "gaussian",
              data = est_dat,
              scale = est_dat$prec_inla,
              control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
              control.predictor = list(compute = TRUE),
              control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
            )
          }, error = function(e) e)
        }

        model_comparison <- data.frame(
          model = c("Fixed only", "iid RE", "ICAR", "BYM2/Alt"),
          dic = c(fit_1a$dic$dic, fit_1b$dic$dic, fit_1c$dic$dic,
                  if (!inherits(fit_1d, "error")) fit_1d$dic$dic else NA_real_),
          waic = c(fit_1a$waic$waic, fit_1b$waic$waic, fit_1c$waic$waic,
                   if (!inherits(fit_1d, "error")) fit_1d$waic$waic else NA_real_)
        )
        print(model_comparison)

        saveRDS(list(fit_1a = fit_1a, fit_1b = fit_1b, fit_1c = fit_1c, fit_1d = fit_1d),
                file.path(out_dir, "inla_fits.rds"))

        # Predictions for all districts
        best_fit <- if (!inherits(fit_1d, "error")) fit_1d else fit_1a
        est_dat$log_rtp_fitted <- best_fit$summary.fitted.values$mean
        est_dat$log_rtp_sd <- best_fit$summary.fitted.values$sd
        est_dat <- est_dat %>%
          mutate(
            p_hat = ifelse(has_rhis, plogis(z_d - log_rtp_fitted), NA_real_),
            p_hat_lower = ifelse(has_rhis,
                                 plogis(z_d - (log_rtp_fitted + 1.96 * log_rtp_sd)), NA_real_),
            p_hat_upper = ifelse(has_rhis,
                                 plogis(z_d - (log_rtp_fitted - 1.96 * log_rtp_sd)), NA_real_)
          )

        saveRDS(est_dat, file.path(out_dir, "estimation_results.rds"))

        # Maps
        map_dat <- all_districts_sf %>%
          left_join(est_dat %>% select(area_id, log_rtp_fitted, log_rtp_sd,
                                       p_hat, has_rhis, in_Ds, has_both),
                    by = "area_id")

        p_map_phat <- ggplot(map_dat) +
          geom_sf(aes(fill = p_hat), color = "white", size = 0.1) +
          scale_fill_viridis_c(option = "C", na.value = "gray90", labels = scales::percent) +
          labs(fill = "Undiagnosed\npositivity", title = "Calibrated undiagnosed positivity by district") +
          theme_void(base_size = 12)

        p_map_uncertainty <- ggplot(map_dat) +
          geom_sf(aes(fill = log_rtp_sd), color = "white", size = 0.1) +
          scale_fill_viridis_c(option = "B", na.value = "gray90", direction = -1) +
          labs(fill = "Posterior SD\n(log-RTP)",
               title = "Uncertainty in selection bias estimate\n(higher = priority for audit)") +
          theme_void(base_size = 12)

        map_dat$data_status <- dplyr::case_when(
          map_dat$has_both ~ "Both RHIS + PHIA",
          map_dat$has_rhis ~ "RHIS only",
          map_dat$in_Ds ~ "PHIA only",
          TRUE ~ "Neither"
        )

        p_map_data <- ggplot(map_dat) +
          geom_sf(aes(fill = data_status), color = "white", size = 0.1) +
          scale_fill_manual(values = c("Both RHIS + PHIA" = "#2166ac",
                                        "RHIS only" = "#67a9cf",
                                        "PHIA only" = "#ef8a62",
                                        "Neither" = "gray90")) +
          labs(fill = "Data\navailability", title = "District data coverage") +
          theme_void(base_size = 12)

        ggsave(file.path(fig_dir, "calibrated_positivity_map.png"), p_map_phat, width = 7, height = 6)
        ggsave(file.path(fig_dir, "calibrated_positivity_map.pdf"), p_map_phat, width = 7, height = 6)
        ggsave(file.path(fig_dir, "uncertainty_map.png"), p_map_uncertainty, width = 7, height = 6)
        ggsave(file.path(fig_dir, "uncertainty_map.pdf"), p_map_uncertainty, width = 7, height = 6)
        ggsave(file.path(fig_dir, "data_availability_map.png"), p_map_data, width = 7, height = 6)
        ggsave(file.path(fig_dir, "data_availability_map.pdf"), p_map_data, width = 7, height = 6)
      }, error = function(e) {
        message("INLA failed: ", e$message)
      })
    } else {
      message("Skipping INLA: adjacency graph failed.")
    }
  }
}

# ---- 10) EXTENSION 1: Demographic group-specific RTP (age x sex x district) ----

message("=== Extension 1: Demographic group-specific RTP ===")

rhis_dg_path <- file.path(out_dir, "rhis_district_age_sex_2021.rds")
if (file.exists(rhis_dg_path)) {
  rhis_dg <- readRDS(rhis_dg_path)

  rhis_dg <- rhis_dg %>%
    left_join(crosswalk %>% select(district, area_id), by = "district") %>%
    filter(!is.na(area_id)) %>%
    group_by(area_id, age_group, sex) %>%
    summarise(
      n_d = sum(n_d, na.rm = TRUE),
      ny_d = sum(ny_d, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      n_neg_d = n_d - ny_d,
      z_d = safe_log_odds_counts(ny_d, n_neg_d),
      v_d = ifelse(ny_d > 0 & n_neg_d > 0, 1 / ny_d + 1 / n_neg_d, NA_real_)
    )

  phia_moz_dg <- phia_moz %>%
    mutate(
      age_group = cut(age, breaks = age_breaks, labels = age_labels, right = FALSE),
      sex_label = ifelse(sex == 1, "Male", "Female")
    )

  phia_undx_dg <- phia_moz_dg %>%
    filter(!(hivstatus == 1 & aware_bin == TRUE))

  phia_dg <- phia_undx_dg %>%
    group_by(area_id, age_group, sex_label) %>%
    summarise(
      n_phia = n(),
      w_total = sum(indweight, na.rm = TRUE),
      w_pos = sum(indweight * hivstatus, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(sex = sex_label) %>%
    mutate(
      p_hat_raw = ifelse(w_total > 0, w_pos / w_total, NA_real_),
      n_eff = n_phia / deff_default,
      p_hat_var = pmin(pmax(p_hat_raw, 0.5 / n_eff), 1 - 0.5 / n_eff),
      Z_hat = safe_qlogis(p_hat_raw, n_eff),
      psi = ifelse(!is.na(n_eff) & !is.na(p_hat_var),
                   1 / (n_eff * p_hat_var * (1 - p_hat_var)), NA_real_),
      p_hat_raw_is_boundary = (w_pos == 0) | (abs(w_pos - w_total) < 0.001)
    )

  cat("RHIS cells:", nrow(rhis_dg), "\n")
  cat("PHIA cells:", nrow(phia_dg), "\n")
  cat("PHIA informative:", sum(!phia_dg$p_hat_raw_is_boundary, na.rm = TRUE), "\n")
  cat("PHIA cells with n < 10:", sum(phia_dg$n_phia < 10, na.rm = TRUE), "\n")
  print(summary(phia_dg$n_phia))

  if (!is.null(all_districts_sf)) {
    all_districts <- st_drop_geometry(all_districts_sf) %>%
      select(area_id, area_name, parent_area_id) %>%
      left_join(provinces_df, by = c("parent_area_id" = "province_id"))

    scaffold_dg <- tidyr::expand_grid(
      area_id = all_districts$area_id,
      age_group = age_labels,
      sex = c("Male", "Female")
    ) %>%
      left_join(all_districts, by = "area_id") %>%
      mutate(
        cell_idx = row_number(),
        district_idx = match(area_id, all_districts$area_id),
        group = paste(age_group, sex, sep = "_"),
        group_idx = as.integer(as.factor(group))
      )

    est_dg <- scaffold_dg %>%
      left_join(rhis_dg, by = c("area_id", "age_group", "sex")) %>%
      left_join(phia_dg %>% select(area_id, age_group, sex,
                                   Z_hat, psi, n_phia, p_hat_raw_is_boundary),
                by = c("area_id", "age_group", "sex")) %>%
      rename(Z_hat_d = Z_hat, psi_d = psi) %>%
      mutate(
        has_rhis = !is.na(z_d),
        in_Ds = !is.na(Z_hat_d) & !is.na(psi_d) & is.finite(psi_d) & !p_hat_raw_is_boundary,
        has_both = has_rhis & in_Ds,
        y = ifelse(has_both, z_d - Z_hat_d, NA_real_),
        prec = ifelse(!is.na(y), 1 / (v_d + psi_d), NA_real_),
        prec_inla = ifelse(is.na(y), 1, pmin(prec, 1e6)),
        province_factor = as.factor(ifelse(is.na(province), "Unknown", province)),
        group_factor = as.factor(group)
      )

    cat("\n=== Demographic estimation dataset ===\n")
    cat("Total cells (scaffold):", nrow(est_dg), "\n")
    cat("  With RHIS:", sum(est_dg$has_rhis), "\n")
    cat("  With PHIA (informative):", sum(est_dg$in_Ds), "\n")
    cat("  With both:", sum(est_dg$has_both), "\n")
    cat("  Groups:", paste(unique(est_dg$group), collapse = ", "), "\n")

    if (requireNamespace("INLA", quietly = TRUE)) {
      fit_dg <- INLA::inla(
        y ~ -1 + group_factor + province_factor,
        family = "gaussian",
        data = est_dg,
        scale = est_dg$prec_inla,
        control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
        control.predictor = list(compute = TRUE),
        control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
      )

      group_effects <- fit_dg$summary.fixed[grep("group_factor", rownames(fit_dg$summary.fixed)), ]
      group_effects$group <- gsub("group_factor", "", rownames(group_effects))
      group_effects$OR <- exp(group_effects$mean)
      group_effects$ci_lo <- exp(group_effects$`0.025quant`)
      group_effects$ci_hi <- exp(group_effects$`0.975quant`)

      print(group_effects %>% select(group, mean, sd, OR))

      saveRDS(fit_dg, file.path(out_dir, "inla_fit_demographic.rds"))
      saveRDS(est_dg, file.path(out_dir, "estimation_demographic.rds"))

      p_group <- ggplot(group_effects, aes(x = reorder(group, OR), y = OR)) +
        geom_point(size = 2) +
        geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
        coord_flip() +
        labs(x = NULL, y = "Implied OR (relative test propensity)",
             title = "Selection bias by demographic group") +
        theme_minimal(base_size = 11)

      ggsave(file.path(fig_dir, "or_by_demographic_group.png"), p_group, width = 7, height = 5, dpi = 300)
      ggsave(file.path(fig_dir, "or_by_demographic_group.pdf"), p_group, width = 7, height = 5)

      if (exists("fit_1a")) {
        cat("Pooled DIC:", fit_1a$dic$dic, "\n")
        cat("Group-specific DIC:", fit_dg$dic$dic, "\n")
      }
    } else {
      message("INLA not available; skipping demographic model fit.")
    }
  }
} else {
  message("Missing rhis_district_age_sex_2021.rds; skipping Extension 1.")
}

# ---- 11) EXTENSION 2: Full measurement error model (stacked INLA) ----

message("=== Extension 2: Full measurement error model ===")

if (requireNamespace("INLA", quietly = TRUE)) {
  # RHIS rows
  rhis_rows <- est_dat %>%
    filter(has_rhis) %>%
    transmute(
      response = z_d,
      precision = 1 / v_d,
      district_idx = district_idx,
      province_factor = province_factor,
      source = "rhis",
      idx_Z = district_idx,
      idx_lambda = district_idx
    )

  # Survey rows
  survey_rows <- est_dat %>%
    filter(in_Ds) %>%
    transmute(
      response = Z_hat_d,
      precision = 1 / psi_d,
      district_idx = district_idx,
      province_factor = province_factor,
      source = "survey",
      idx_Z = district_idx,
      idx_lambda = NA_integer_
    )

  stacked <- bind_rows(rhis_rows, survey_rows)
  stacked$is_rhis <- as.numeric(stacked$source == "rhis")
  stacked$province_lambda <- ifelse(stacked$source == "rhis",
                                    as.character(stacked$province_factor),
                                    NA_character_)
  stacked$province_lambda <- as.factor(stacked$province_lambda)

  cat("Stacked dataset:", nrow(stacked), "rows\n")
  cat("  RHIS rows:", sum(stacked$source == "rhis"), "\n")
  cat("  Survey rows:", sum(stacked$source == "survey"), "\n")

  formula_full_a <- response ~ -1 +
    is_rhis +
    f(idx_Z, model = "iid") +
    province_lambda

  fit_full_a <- INLA::inla(
    formula_full_a,
    family = "gaussian",
    data = stacked,
    scale = stacked$precision,
    control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
  )

  saveRDS(fit_full_a, file.path(out_dir, "inla_fit_full_a.rds"))

  if (exists("adj_graph")) {
    formula_full_b <- response ~ -1 +
      is_rhis +
      f(idx_Z, model = "besag", graph = adj_graph) +
      province_lambda

    fit_full_b <- INLA::inla(
      formula_full_b,
      family = "gaussian",
      data = stacked,
      scale = stacked$precision,
      control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
      control.predictor = list(compute = TRUE),
      control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
    )
    saveRDS(fit_full_b, file.path(out_dir, "inla_fit_full_b.rds"))
  }
} else {
  message("INLA not available; skipping Extension 2.")
}

message("=== Done ===")

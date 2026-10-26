# =============================================================================
# 02_standardise_phia_moz.R
# Build the Mozambique INSIDA 2021 (PHIA) person-level dataset: HIV status and
# awareness from the biomarker file, district by spatial join of the cluster
# centroids to district boundaries, interview covariates merged by personid.
# Reads:  data/insida2021{adultbio,centroids,adultind}.csv and
#         data/all_moz_data.rda (district boundaries); not in the release
# Writes: output/phia_insida2021_with_awareness.rds, output/phia_awareness_summary_by_province.csv
# Run:    Rscript R/02_standardise_phia_moz.R  (from the repository root)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Step 1: Load biomarker data and report aware coding ----
message("=== Step 1: Load biomarker data ===")

bio <- read.csv("data/insida2021adultbio.csv", stringsAsFactors = FALSE)

cat("=== aware variable ===\n")
cat("Values:", paste(sort(unique(bio$aware)), collapse = ", "), "\n")
print(table(bio$aware, bio$hivstatusfinal, useNA = "always"))

bio_pos <- bio[bio$hivstatusfinal == 1, ]

if ("arv" %in% names(bio)) {
  cat("\n=== aware vs arv (among HIV+) ===\n")
  print(table(bio_pos$aware, bio_pos$arv, useNA = "always"))
}

if ("awareselfreported" %in% names(bio)) {
  cat("\n=== aware vs awareselfreported (among HIV+) ===\n")
  print(table(bio_pos$aware, bio_pos$awareselfreported, useNA = "always"))
}

cat("\n=== Summary ===\n")
cat("Total in bio:", nrow(bio), "\n")
cat("HIV+ (hivstatusfinal == 1):", sum(bio$hivstatusfinal == 1, na.rm = TRUE), "\n")
cat("HIV- (hivstatusfinal == 2):", sum(bio$hivstatusfinal == 2, na.rm = TRUE), "\n")

# ---- Step 2: Recode to clean binary variables ----
message("=== Step 2: Recode and clean ===")

bio_clean <- bio %>%
  filter(age >= 15) %>%
  mutate(
    hivstatus = case_when(
      hivstatusfinal == 1 ~ 1L,
      hivstatusfinal == 2 ~ 0L,
      TRUE ~ NA_integer_
    ),
    aware_bin = case_when(
      hivstatus != 1 ~ NA,
      aware == 1 ~ TRUE,
      aware == 2 ~ FALSE,
      TRUE ~ NA
    ),
    sex = gender,
    age_group = cut(age, breaks = c(15, 20, 25, 50, Inf),
                    labels = c("15-19", "20-24", "25-49", "50+"), right = FALSE),
    sex_label = ifelse(gender == 1, "Male", "Female"),
    indweight = btwt0
  ) %>%
  filter(!is.na(hivstatus))

cat("Clean bio dataset:", nrow(bio_clean), "rows\n")
cat("HIV+:", sum(bio_clean$hivstatus == 1), "\n")
cat("  Aware:", sum(bio_clean$aware_bin == TRUE, na.rm = TRUE), "\n")
cat("  Not aware:", sum(bio_clean$aware_bin == FALSE, na.rm = TRUE), "\n")
cat("  Awareness missing:", sum(is.na(bio_clean$aware_bin) & bio_clean$hivstatus == 1), "\n")

# ---- Step 3: Add geographic coordinates from centroids ----
message("=== Step 3: Join centroids ===")

centroids <- read.csv("data/insida2021centroids.csv", stringsAsFactors = FALSE)

bio_clean <- bio_clean %>%
  left_join(centroids, by = c("centroidid" = "CentroidID"))

cat("Individuals with coordinates:", sum(!is.na(bio_clean$longitude)), "\n")
cat("Individuals WITHOUT coordinates:", sum(is.na(bio_clean$longitude)), "\n")

# ---- Step 4: Spatial join to districts ----
message("=== Step 4: Spatial join to districts ===")

load("data/all_moz_data.rda")

if (!exists("area_merged")) stop("area_merged not found in all_moz_data.rda")

# Prepare district boundaries
if (!"geometry" %in% names(area_merged)) stop("area_merged has no geometry")

districts_sf <- st_as_sf(area_merged, sf_column_name = "geometry", crs = 4326) %>%
  filter(area_level == 2) %>%
  st_make_valid()

# Convert survey points to sf
bio_sf <- bio_clean %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Point-in-polygon join to districts
bio_matched <- st_join(
  bio_sf,
  districts_sf %>% select(area_id, area_name, parent_area_id),
  join = st_within
)

n_unmatched <- sum(is.na(bio_matched$area_id))
cat("Points not matching any district:", n_unmatched, "\n")

if (n_unmatched > 0) {
  unmatched_idx <- which(is.na(bio_matched$area_id))
  nearest <- st_nearest_feature(bio_matched[unmatched_idx, ], districts_sf)
  bio_matched$area_id[unmatched_idx] <- districts_sf$area_id[nearest]
  bio_matched$area_name[unmatched_idx] <- districts_sf$area_name[nearest]
  bio_matched$parent_area_id[unmatched_idx] <- districts_sf$parent_area_id[nearest]
}

bio_final <- st_drop_geometry(bio_matched)

# ---- Step 5: Add province names ----
message("=== Step 5: Add province names ===")

provinces_df <- st_drop_geometry(area_merged) %>%
  filter(area_level == 1) %>%
  select(province_id = area_id, province_name = area_name)

bio_final <- bio_final %>%
  left_join(provinces_df, by = c("parent_area_id" = "province_id"))

# ---- Step 6: Add back individuals without coordinates ----
message("=== Step 6: Add back missing coords ===")

no_coords <- bio_clean %>%
  filter(is.na(longitude) | is.na(latitude)) %>%
  mutate(
    area_id = NA_character_,
    area_name = NA_character_,
    parent_area_id = NA_character_,
    province_name = NA_character_
  )

bio_final <- bind_rows(bio_final, no_coords)

cat("Individuals with no geographic assignment:", sum(is.na(bio_final$area_id)), "\n")

# ---- Step 7: Merge adultd covariates ----
message("=== Step 7: Merge adultd covariates ===")

adultd <- read.csv("data/insida2021adultind.csv", stringsAsFactors = FALSE)

adultd_select <- adultd %>%
  select(
    personid,
    urban,
    hivtstever, hivtesty, hivtstrslt,
    arvscurrent, arvstakenev,
    any_of(c("education", "wealthquintile", "wealthscorecont"))
  )

bio_final <- bio_final %>%
  left_join(adultd_select, by = "personid")

# ---- Step 8: Report and save ----
message("=== Step 8: Report and save ===")

cat("\n=== Final dataset ===\n")
cat("Total individuals:", nrow(bio_final), "\n")
cat("Districts covered:", dplyr::n_distinct(bio_final$area_id, na.rm = TRUE), "\n")
cat("Provinces covered:", dplyr::n_distinct(bio_final$province_name, na.rm = TRUE), "\n")

summary_tbl <- bio_final %>%
  filter(!is.na(area_id)) %>%
  group_by(province_name) %>%
  summarise(
    n_total = n(),
    n_hiv_pos = sum(hivstatus == 1, na.rm = TRUE),
    n_aware = sum(aware_bin == TRUE, na.rm = TRUE),
    n_not_aware = sum(aware_bin == FALSE, na.rm = TRUE),
    pct_aware = round(100 * n_aware / (n_aware + n_not_aware), 1),
    n_districts = n_distinct(area_id),
    .groups = "drop"
  )

print(summary_tbl)

bio_pos_final <- bio_final %>% filter(hivstatus == 1)

if ("arvscurrent" %in% names(bio_pos_final)) {
  cat("ART coverage among HIV+ (self-reported arvscurrent): ",
      weighted.mean(bio_pos_final$arvscurrent == 1, bio_pos_final$indweight, na.rm = TRUE),
      "\n", sep = "")
}

if ("awareselfreported" %in% names(bio_pos_final)) {
  cat("Non-disclosure rate (aware but self-report unaware): ",
      mean(bio_pos_final$awareselfreported[bio_pos_final$aware_bin == TRUE] == 2, na.rm = TRUE),
      "\n", sep = "")
}

saveRDS(bio_final, file.path(out_dir, "phia_insida2021_with_awareness.rds"))
write.csv(summary_tbl, file.path(out_dir, "phia_awareness_summary_by_province.csv"), row.names = FALSE)

message("=== Done ===")

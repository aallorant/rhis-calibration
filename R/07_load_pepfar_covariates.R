#!/usr/bin/env Rscript
# R/07_load_pepfar_covariates.R
# Harmonise the PEPFAR MER programme-data extracts of the four countries into
# one district-level covariate file: modal priority tier over 2022-2024
# (scaleup / non_scaleup / not_set) and HTS_SELF volume, with manual crosswalks
# to the canonical district keys; assumptions are repeated in the csv header.
# Reads:  data/raw/covariates/pepfar-program-data/{moz,lso,mwi,zwe}.csv (not in the release)
# Writes: data/covariates/pepfar_covariates_harmonised.csv
# Run:    Rscript R/07_load_pepfar_covariates.R  (from the repository root)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

PEPFAR_DIR <- "data/raw/covariates/pepfar-program-data"
OUT_PATH   <- "data/covariates/pepfar_covariates_harmonised.csv"
WINDOW     <- 2022:2024

# Tie-break order (highest first).
TIER_RANK <- c(
  "Scale-Up Aggressive"  = 1,
  "Scale-Up Saturation"  = 2,
  "Attained"             = 3,
  "Centrally Supported"  = 4,
  "Sustained"            = 5,
  "Not Set"              = 6
)

COUNTRY_FILES <- c(
  Mozambique = "moz.csv",
  Lesotho    = "lso.csv",
  Malawi     = "mwi.csv",
  Zimbabwe   = "zwe.csv"
)

load_one <- function(country, file) {
  raw <- read_csv(file.path(PEPFAR_DIR, file), show_col_types = FALSE)

  # Filter PSNU meta-rows from both district and province columns.
  has_psnu <- function(x) stringr::str_detect(coalesce(x, ""),
                                              regex("PSNU", ignore_case = TRUE))
  raw <- raw |>
    filter(!has_psnu(district), !has_psnu(province)) |>
    filter(!is.na(district), nzchar(trimws(district)))

  raw <- raw |>
    mutate(
      year    = as.integer(year),
      quarter = suppressWarnings(as.integer(quarter))
    )

  # Country-specific province override (lossless; raw cols preserved).
  if (country == "Malawi") {
    raw <- raw |> mutate(province = snu2_raw)
  } else if (country == "Lesotho") {
    raw <- raw |> mutate(province = district)
  }

  raw |>
    select(country, province, district,
           year, quarter, indicator, priority, value)
}

modal_priority <- function(priorities) {
  # Modal tier within a district across 2022-2024; ties broken by TIER_RANK.
  tab <- table(priorities)
  if (length(tab) == 0) return(NA_character_)
  top_count <- max(tab)
  candidates <- names(tab)[tab == top_count]
  ranks <- TIER_RANK[candidates]
  ranks[is.na(ranks)] <- 99L  # unknown tiers sort last
  candidates[which.min(ranks)]
}

collapse_to_3lvl <- function(modal) {
  dplyr::case_when(
    is.na(modal)                                          ~ NA_character_,
    stringr::str_detect(modal, "Scale-Up")                ~ "scaleup",
    modal %in% c("Attained", "Centrally Supported",
                 "Sustained")                             ~ "non_scaleup",
    modal == "Not Set"                                    ~ "not_set",
    TRUE                                                  ~ NA_character_
  )
}

# ---- crosswalks (pepfar -> canonical), set by hand ------------------------

MWI_PROV_RENAME <- c(
  "Central Eastern Zone" = "Central East Zone",
  "Central Western Zone" = "Central West Zone",
  "South Eastern Zone"   = "South East Zone",
  "South Western Zone"   = "South West Zone"
)

MOZ_PROV_RENAME <- c(
  "Cidade De Maputo" = "Maputo",
  "Zambezia"         = "Zambézia"
)

# (post-province-rename province, pepfar district) -> canonical district
MOZ_DIST_RENAME <- tibble::tribble(
  ~province,    ~old,                  ~new,
  "Maputo",     "Kamavota",            "KaMavota",
  "Maputo",     "Kamaxakeni",          "KaMaxakeni",
  "Maputo",     "Kampfumu",            "KaMpfumu",
  "Maputo",     "Kamubukwana",         "KaMubukwana",
  "Maputo",     "Kanyaka",             "KaNyaka",
  "Maputo",     "Katembe",             "KaTembe",
  "Maputo",     "Nlhamankulu",         "Nhlamankulu",
  "Gaza",       "Chonguene",           "Chongoene",
  "Gaza",       "Chókwè",              "Chókwé",
  "Inhambane",  "Homoine",             "Homoíne",
  "Manica",     "Barue",               "Báruè",
  "Manica",     "Sussundenga",         "Sussendenga",
  "Nampula",    "Erati",               "Eráti",
  "Nampula",    "Ilha De Moçambique",  "Ilha de Moçambique",
  "Nampula",    "Liupo",               "Liúpo",
  "Nampula",    "Mecuburi",            "Mecubúri",
  "Nampula",    "Murrupula",           "Morrupula",
  "Nampula",    "Nacala-A-Velha",      "Nacala-a-Velha",
  "Nampula",    "Nacaroa",             "Nacarôa",
  "Nampula",    "Ribaue",              "Ribáuè",
  "Niassa",     "Chimbonila",          "Chimbunila",
  "Niassa",     "Ngauma",              "Ngaúma",
  "Sofala",     "Buzi",                "Búzi",
  "Sofala",     "Maringue",            "Maríngué",
  "Tete",       "Chiuta",              "Chiúta",
  "Tete",       "Doa",                 "Dôa",
  "Tete",       "Magoe",               "Mágoè",
  "Tete",       "Maravia",             "Marávia",
  "Zambézia",   "Alto Molocue",        "Alto Molócuè",
  "Zambézia",   "Gile",                "Gilé",
  "Zambézia",   "Gurue",               "Guruè",
  "Zambézia",   "Maganja Da Costa",    "Maganja da Costa",
  "Zambézia",   "Namarroi",            "Namarói"
)

MOZ_DROP <- tibble::tribble(
  ~province,      ~district,
  "Cabo Delgado", "Mocimboa Da Praia",
  "Cabo Delgado", "Quissanga"
)

ZWE_DIST_RENAME <- tibble::tribble(
  ~province,             ~old,        ~new,
  "Mashonaland Central", "Mt Darwin", "Mount Darwin",
  "Mashonaland East",    "Murehwa",   "Murewa"
)

ZWE_WARDS <- tibble::tribble(
  ~province,  ~ward,                          ~rollup_source,
  "Bulawayo", "Emakhandeni",                  "Bulawayo",
  "Bulawayo", "Nkulumane",                    "Bulawayo",
  "Bulawayo", "Northern Suburbs",             "Bulawayo",
  "Harare",   "Chitungwiza Urban",            "Harare",
  "Harare",   "Eastern",                      "Harare",
  "Harare",   "Harare Central South Eastern", "Harare",
  "Harare",   "North Western",                "Harare",
  "Harare",   "Northern",                     "Harare",
  "Harare",   "South Western",                "Harare",
  "Harare",   "Southern",                     "Harare",
  "Harare",   "West South West",              "Harare",
  "Harare",   "Western",                      "Harare"
)

ZWE_NONCOVERED <- tibble::tribble(
  ~province,             ~district,
  "Manicaland",          "Chimanimani",
  "Manicaland",          "Mutare City",
  "Manicaland",          "Nyanga",
  "Mashonaland Central", "Bindura",
  "Mashonaland Central", "Centenary",
  "Mashonaland Central", "Mbire",
  "Mashonaland Central", "Rushinga",
  "Mashonaland East",    "Chikomba",
  "Mashonaland East",    "Hwedza",
  "Mashonaland East",    "Mudzi",
  "Mashonaland East",    "Mutoko",
  "Mashonaland East",    "Uzumba Maramba Pfungwe",
  "Mashonaland West",    "Kariba",
  "Mashonaland West",    "Mhondoro",
  "Mashonaland West",    "Sanyati",
  "Midlands",            "Chirumhanzu",
  "Midlands",            "Gokwe North",
  "Midlands",            "Zvishavane",
  "Masvingo",            "Bikita",
  "Matabeleland North",  "Binga"
)

apply_crosswalks <- function(harm) {
  # MWI: rename zones, then duplicate Mzimba into Mzimba North / Mzimba South.
  mwi <- harm |> dplyr::filter(country == "Malawi") |>
    dplyr::mutate(province = dplyr::recode(province, !!!MWI_PROV_RENAME))
  mzimba <- mwi |> dplyr::filter(district == "Mzimba")
  if (nrow(mzimba) != 1L) {
    stop("MWI: expected exactly one Mzimba row, got ", nrow(mzimba))
  }
  mwi_split <- dplyr::bind_rows(
    mzimba |> dplyr::mutate(district = "Mzimba North"),
    mzimba |> dplyr::mutate(district = "Mzimba South")
  )
  mwi <- mwi |> dplyr::filter(district != "Mzimba") |>
    dplyr::bind_rows(mwi_split)

  # MOZ: rename provinces, then rename districts, then drop two Cabo Delgado.
  moz <- harm |> dplyr::filter(country == "Mozambique") |>
    dplyr::mutate(province = dplyr::recode(province, !!!MOZ_PROV_RENAME))
  moz <- moz |>
    dplyr::left_join(MOZ_DIST_RENAME |>
                       dplyr::rename(district = old, district_new = new),
                     by = c("province", "district")) |>
    dplyr::mutate(district = dplyr::coalesce(district_new, district)) |>
    dplyr::select(-district_new)
  moz <- moz |>
    dplyr::anti_join(MOZ_DROP, by = c("province", "district"))

  # ZWE: rename two districts, then add ward roll-ups + non-covered entries.
  zwe <- harm |> dplyr::filter(country == "Zimbabwe")
  zwe <- zwe |>
    dplyr::left_join(ZWE_DIST_RENAME |>
                       dplyr::rename(district = old, district_new = new),
                     by = c("province", "district")) |>
    dplyr::mutate(district = dplyr::coalesce(district_new, district)) |>
    dplyr::select(-district_new)

  # Roll up: each ward gets the values from its city's pepfar row.
  city_rows <- zwe |>
    dplyr::filter((province == "Bulawayo" & district == "Bulawayo") |
                  (province == "Harare"   & district == "Harare"))
  if (nrow(city_rows) != 2L) {
    stop("ZWE: expected Bulawayo|Bulawayo and Harare|Harare city rows, got ",
         nrow(city_rows))
  }
  ward_rows <- ZWE_WARDS |>
    dplyr::left_join(
      city_rows |>
        dplyr::select(rollup_source = district,
                      pepfar_priority_modal,
                      pepfar_priority_3lvl,
                      hts_self_sum_2022_2024,
                      n_pepfar_rows_2022_2024),
      by = "rollup_source"
    ) |>
    dplyr::transmute(
      country = "Zimbabwe",
      province,
      district = ward,
      pepfar_priority_modal,
      pepfar_priority_3lvl,
      hts_self_sum_2022_2024,
      n_pepfar_rows_2022_2024
    )
  # Drop the city-level pepfar rows now that they've been propagated to wards.
  zwe <- zwe |>
    dplyr::filter(!((province == "Bulawayo" & district == "Bulawayo") |
                    (province == "Harare"   & district == "Harare")))
  # Append wards.
  zwe <- dplyr::bind_rows(zwe, ward_rows)
  # Append non-covered admin_2_districts.
  noncov_rows <- ZWE_NONCOVERED |>
    dplyr::mutate(
      country = "Zimbabwe",
      pepfar_priority_modal   = NA_character_,
      pepfar_priority_3lvl    = "non_scaleup",
      hts_self_sum_2022_2024  = 0,
      n_pepfar_rows_2022_2024 = 0L
    )
  zwe <- dplyr::bind_rows(zwe, noncov_rows)

  # LSO: untouched.
  lso <- harm |> dplyr::filter(country == "Lesotho")

  out <- dplyr::bind_rows(lso, mwi, moz, zwe) |>
    dplyr::arrange(country, province, district)

  # pepfar_covered: TRUE everywhere except the 20 ZWE non-covered rows.
  noncov_keys <- paste(ZWE_NONCOVERED$province, ZWE_NONCOVERED$district, sep = " | ")
  out <- out |>
    dplyr::mutate(
      pepfar_covered = !(country == "Zimbabwe" &
                         paste(province, district, sep = " | ") %in% noncov_keys)
    )

  out
}

main <- function() {
  long <- purrr::imap_dfr(COUNTRY_FILES, ~ load_one(country = .y, file = .x))

  windowed <- long |> filter(year %in% WINDOW)

  # Priority: modal per (country, province, district), regardless of indicator.
  prio <- windowed |>
    filter(!is.na(priority)) |>
    distinct(country, province, district, year, quarter, priority) |>
    group_by(country, province, district) |>
    summarise(
      pepfar_priority_modal = modal_priority(priority),
      .groups = "drop"
    ) |>
    mutate(pepfar_priority_3lvl = collapse_to_3lvl(pepfar_priority_modal))

  # HTS_SELF volume: sum per (country, province, district) in window.
  hts <- windowed |>
    filter(indicator == "HTS_SELF") |>
    group_by(country, province, district) |>
    summarise(
      hts_self_sum_2022_2024 = sum(as.numeric(value), na.rm = TRUE),
      n_pepfar_rows_2022_2024 = dplyr::n(),
      .groups = "drop"
    )

  # All (country, province, district) keys observed in 2022-2024.
  keys <- windowed |>
    distinct(country, province, district)

  harmonised <- keys |>
    left_join(prio, by = c("country", "province", "district")) |>
    left_join(hts,  by = c("country", "province", "district")) |>
    mutate(
      hts_self_sum_2022_2024  = tidyr::replace_na(hts_self_sum_2022_2024, 0),
      n_pepfar_rows_2022_2024 = tidyr::replace_na(n_pepfar_rows_2022_2024, 0L)
    ) |>
    arrange(country, province, district)

  # Apply the manual crosswalks defined above.
  harmonised <- apply_crosswalks(harmonised)
  harmonised <- harmonised |>
    mutate(pepfar_priority_3lvl = factor(pepfar_priority_3lvl,
                                         levels = c("non_scaleup", "scaleup", "not_set")))

  # File header (csv comments).
  header <- c(
    sprintf("# pepfar covariates harmonised - generated %s", Sys.Date()),
    "# source: data/raw/covariates/pepfar-program-data/{moz,lso,mwi,zwe}.csv",
    "# window: 2022-2024 (results level only)",
    "# assumptions:",
    "#   - districts with zero pepfar rows in 2022-2024 get hts_self_sum=0",
    "#   - modal priority tier across 2022-2024 collapses time variation",
    "#   - ties broken: aggressive > saturation > attained > centrally supported > sustained > not set",
    "#   - 'not set' retained as its own categorical level",
    "# join keys downstream:",
    "#   - moz, mwi, zwe: (province, district) composite",
    "#   - lso: district alone (province = district by override; functionally district-only)",
    "# overrides relative to scrape:",
    "#   - malawi province = snu2_raw (pepfar zone) to match canonical district keys",
    "#   - lesotho province = district (cosmetic, to match canonical)",
    "# crosswalks applied to align with canonical district keys:",
    "#   - mwi: zone names renamed (eastern -> east, western -> west); mzimba",
    "#     duplicated into mzimba north + mzimba south, both inheriting the",
    "#     full pepfar mzimba value (geographic split unknown).",
    "#   - moz: province renames (cidade de maputo -> maputo, zambezia -> zambézia);",
    "#     33 district renames for accents/case/orthography; 2 districts dropped",
    "#     (cabo delgado | mocimboa da praia, cabo delgado | quissanga) as",
    "#     genuinely missing from canonical (both hts_self = 0).",
    "#   - zwe: mt darwin -> mount darwin; murehwa -> murewa.",
    "#     bulawayo|bulawayo and harare|harare are city-level pepfar rollups;",
    "#     each is replaced by its constituent canonical wards (3 in bulawayo,",
    "#     9 in harare), with each ward inheriting the city's priority and",
    "#     hts_self_sum. 7 facility-type canonical entries dropped (mpilo,",
    "#     ingutsheni, ubh, etc.). 20 canonical admin_2_districts that have no",
    "#     pepfar coverage at psnu level (binga, bikita, etc.) are added with",
    "#     pepfar_covered=false, priority=non_scaleup, hts_self_sum=0.",
    "# pepfar_covered column: TRUE except for the 20 ZWE non-covered rows.",
    ""
  )

  if (!dir.exists(dirname(OUT_PATH))) dir.create(dirname(OUT_PATH), recursive = TRUE)
  writeLines(header, OUT_PATH)
  write.table(harmonised, OUT_PATH, sep = ",", row.names = FALSE,
              quote = TRUE, append = TRUE, col.names = TRUE)

  cat("Wrote ", nrow(harmonised), " rows to ", OUT_PATH, "\n", sep = "")
  cat("Per country:\n")
  print(harmonised |> count(country))
  invisible(harmonised)
}

if (!interactive()) {
  out <- main()
}

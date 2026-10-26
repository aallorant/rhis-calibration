# ============================================================
# Standardise the round-1 PHIA microdata of Lesotho (LePHIA 2020) and
# Zimbabwe (ZIMPHIA 2020) to one person-level schema (province, district,
# age, sex, hivstatus, aware_bin, weights, strata, PSU), with diagnostics.
# Reads:  data/phia/{lesotho,zimbabwe}/*adultbio*.csv and *adultind*.csv
#         (not in the release); data/rhis/standardised/rhis_{country}.rds
# Writes: data/phia/standardised/phia_{lesotho,zimbabwe}.rds
# Run:    Rscript R/02_standardise_phia.R (runs standardise_zimbabwe_phia();
#         call standardise_lesotho_phia() first to write phia_lesotho.rds)
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(survey)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(path) > 0) {
    return(dirname(normalizePath(path)))
  }
  normalizePath(getwd())
}

find_project_dir <- function() {
  script_dir <- get_script_dir()
  candidates <- unique(c(
    script_dir,
    dirname(script_dir),
    getwd(),
    file.path(getwd(), "rhis-calibration"),
    file.path(script_dir, "rhis-calibration"),
    file.path(dirname(script_dir), "rhis-calibration")
  ))

  marker_rel <- file.path("data", "phia", "lesotho", "lephia2020adultbio.csv")

  for (cand in candidates) {
    cand_norm <- tryCatch(normalizePath(cand, mustWork = FALSE), error = function(e) cand)
    if (file.exists(file.path(cand_norm, marker_rel))) {
      return(cand_norm)
    }
  }

  for (cand in candidates) {
    cand_norm <- tryCatch(normalizePath(cand, mustWork = FALSE), error = function(e) cand)
    if (dir.exists(file.path(cand_norm, "data", "phia"))) {
      return(cand_norm)
    }
  }

  stop("Could not locate the project root containing data/phia.")
}

project_dir <- find_project_dir()
phia_dir <- file.path(project_dir, "data", "phia")
options(survey.lonely.psu = "adjust")

summarise_weight_vars <- function(df, label) {
  wt_vars <- grep("wt|weight", names(df), ignore.case = TRUE, value = TRUE)
  if (length(wt_vars) == 0) {
    cat("  No matches for", label, "\n")
    return(invisible(NULL))
  }

  rep_pattern <- "^(btwt|intwt)\\d{3}$"
  rep_vars <- wt_vars[grepl(rep_pattern, wt_vars, ignore.case = TRUE)]
  main_vars <- setdiff(wt_vars, rep_vars)

  cat("  Matches for", label, ":\n")
  for (v in main_vars) {
    cat("   -", v, "\n")
  }
  if (length(rep_vars) > 0) {
    cat(sprintf("   - %d replicate weights %s-%s\n",
                length(rep_vars), min(rep_vars), max(rep_vars)))
  }
  cat("\n")

  for (v in main_vars) {
    cat("  >>>", v, "\n")
    x <- df[[v]]
    if (is.numeric(x)) {
      print(summary(x))
    } else {
      u <- unique(x[!is.na(x)])
      if (length(u) < 30) {
        print(table(x, useNA = "ifany"))
      } else {
        print(utils::head(sort(table(x), decreasing = TRUE), 20))
      }
    }
    cat("\n")
  }
}

print_var_diagnostics <- function(df, vars, label) {
  if (length(vars) == 0) {
    cat("  No matches for", label, "\n")
    return(invisible(NULL))
  }

  cat("  Matches for", label, ":\n")
  for (v in vars) {
    cat("   -", v, "\n")
  }
  cat("\n")

  for (v in vars) {
    cat("  >>>", v, "\n")
    x <- df[[v]]
    x_no_na <- x[!is.na(x)]

    if (is.numeric(x)) {
      print(summary(x))
    } else {
      u <- unique(x_no_na)
      if (length(u) < 30) {
        print(table(x, useNA = "ifany"))
      } else {
        print(utils::head(sort(table(x), decreasing = TRUE), 20))
      }
    }
    cat("\n")
  }
}

print_file_inventory <- function(country_name, country_dir) {
  files <- sort(list.files(country_dir, full.names = TRUE))
  cat("===", toupper(country_name), "PHIA DIAGNOSTICS ===\n\n")
  cat("Files in", file.path("data", "phia", country_name), "/:\n")
  cat(paste(" ", basename(files), collapse = "\n"), "\n\n")
  invisible(files)
}

print_phia_file_summary <- function(file, label) {
  df <- read_csv(file, show_col_types = FALSE, progress = FALSE)
  cat("---", toupper(label), "FILE ---\n")
  cat("File:", basename(file), "\n")
  cat("Dimensions:", nrow(df), "x", ncol(df), "\n")
  cat("All column names:\n")
  for (i in seq_along(names(df))) {
    cat(sprintf("  [%03d] %s\n", i, names(df)[i]))
  }
  cat("\nFirst 3 rows:\n")
  print(utils::head(df, 3))
  cat("\n")
  invisible(df)
}

print_centroid_summary <- function(file) {
  if (!is.na(file) && file.exists(file)) {
    cent <- read_csv(file, show_col_types = FALSE, progress = FALSE)
    cat("--- CENTROIDS FILE ---\n")
    cat("File:", basename(file), "\n")
    cat("Dimensions:", nrow(cent), "x", ncol(cent), "\n")
    cat("Columns:\n")
    for (i in seq_along(names(cent))) {
      cat(sprintf("  [%03d] %s\n", i, names(cent)[i]))
    }
    cat("\n")
  } else {
    cat("--- CENTROIDS FILE ---\n")
    cat("Not found.\n\n")
  }
}

print_variable_inventory <- function(df, file_label) {
  key_patterns <- list(
    hiv_status = "hiv|status",
    awareness = "aware|knowledge|know",
    art = "arv|art|treat",
    testing = "tst|test",
    strata_psu = "strat|unit|psu|clust",
    age_sex = "age|gender|sex",
    geography = "region|province|district|zone|area"
  )

  cat("=== VARIABLE SEARCH:", toupper(file_label), "FILE ===\n\n")
  for (nm in names(key_patterns)) {
    vars <- grep(key_patterns[[nm]], names(df), ignore.case = TRUE, value = TRUE)
    print_var_diagnostics(df, vars, nm)
  }
  summarise_weight_vars(df, "weights")
}

print_survey_validation <- function(df, country_label, published_prev = NULL, published_awareness = NULL, published_art_aware = NULL) {
  des <- svydesign(
    ids = ~varunit,
    strata = ~varstrat,
    weights = ~btwt0,
    data = df,
    nest = TRUE
  )
  prev <- svymean(~hivstatus, des, na.rm = TRUE)
  des_pos <- subset(des, hivstatus == 1)
  des_pos <- update(des_pos, aware_num = ifelse(is.na(aware_bin), NA_real_, ifelse(aware_bin, 1, 0)))
  aware_est <- svymean(~aware_num, des_pos, na.rm = TRUE)

  cat("\n", country_label, "survey-weighted validation:\n", sep = "")
  cat("  Weighted HIV prevalence:", sprintf("%.1f%%", 100 * coef(prev)[1]), "\n")
  cat("  Weighted awareness:", sprintf("%.1f%%", 100 * coef(aware_est)[1]), "\n")
  if (!is.null(published_prev) || !is.null(published_awareness) || !is.null(published_art_aware)) {
    cat("\n", country_label, "published values for comparison:\n", sep = "")
    if (!is.null(published_prev)) {
      cat("  Prevalence:", published_prev, "\n")
    }
    if (!is.null(published_awareness)) {
      cat("  Awareness:", published_awareness, "\n")
    }
    if (!is.null(published_art_aware)) {
      cat("  ART among aware:", published_art_aware, "\n")
    }
  }
}

standardise_lesotho_phia <- function() {
  les_dir <- file.path(phia_dir, "lesotho")
  files <- print_file_inventory("lesotho", les_dir)

  bio_file <- files[str_detect(basename(files), regex("adultbio.*\\.csv$", ignore_case = TRUE))][1]
  ind_file <- files[str_detect(basename(files), regex("adultind.*\\.csv$", ignore_case = TRUE))][1]
  cent_file <- files[str_detect(basename(files), regex("centroids.*\\.csv$", ignore_case = TRUE))][1]

  if (is.na(bio_file) || !file.exists(bio_file)) stop("Could not find Lesotho adultbio CSV.")
  if (is.na(ind_file) || !file.exists(ind_file)) stop("Could not find Lesotho adultind CSV.")

  bio <- print_phia_file_summary(bio_file, "biomarker")
  ind <- print_phia_file_summary(ind_file, "interview")
  print_centroid_summary(cent_file)
  print_variable_inventory(bio, "biomarker")
  print_variable_inventory(ind, "interview")

  cat("=== LESOTHO MERGE + CODING DIAGNOSTICS ===\n\n")
  bio_base <- read.csv(bio_file, stringsAsFactors = FALSE)
  ind_base <- read.csv(ind_file, stringsAsFactors = FALSE)
  cent_file_csv <- file.path(les_dir, "lephia2020centroids.csv")
  hh_file <- file.path(les_dir, "lephia2020hh.csv")
  roster_file <- file.path(les_dir, "lephia2020roster.csv")

  phia <- merge(bio_base, ind_base, by = "personid", suffixes = c(".bio", ".ind"))

  cat("Merged rows:", nrow(phia), "\n")
  cat("Bio rows:", nrow(bio_base), "  Ind rows:", nrow(ind_base), "\n")

  dup_core <- names(phia)[grepl("\\.(bio|ind)$", names(phia))]
  cat("\nDuplicate-suffixed columns after merge:\n")
  if (length(dup_core) == 0) {
    cat("  None\n")
  } else {
    cat(paste(" ", dup_core, collapse = "\n"), "\n")
  }

  cat("\n--- hivstatusfinal ---\n")
  if ("hivstatusfinal.bio" %in% names(phia)) {
    print(table(phia$`hivstatusfinal.bio`, useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- aware ---\n")
  if ("aware" %in% names(phia)) {
    print(table(phia$aware, useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- gender ---\n")
  gender_col <- grep("^gender(\\.|$)", names(phia), value = TRUE)[1]
  if (!is.na(gender_col)) {
    print(table(phia[[gender_col]], useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- district (bio) ---\n")
  dist_col <- grep("^district(\\.bio|$)", names(phia), value = TRUE)[1]
  if (!is.na(dist_col)) {
    print(table(phia[[dist_col]], useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- age range ---\n")
  age_col <- grep("^age(\\.bio|$)", names(phia), value = TRUE)[1]
  if (!is.na(age_col)) {
    cat("Range:", paste(range(phia[[age_col]], na.rm = TRUE), collapse = " to "), "\n")
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- btwt0 (biomarker weight) ---\n")
  btwt_col <- grep("^btwt0$", names(phia), value = TRUE)[1]
  if (!is.na(btwt_col)) {
    cat("Non-zero:", sum(phia[[btwt_col]] > 0, na.rm = TRUE), "\n")
    cat("Zero:", sum(phia[[btwt_col]] == 0, na.rm = TRUE), "\n")
    cat("NA:", sum(is.na(phia[[btwt_col]])), "\n")
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- arvscurrent (from ind file) ---\n")
  arv_col <- grep("arvscurrent", names(phia), value = TRUE)[1]
  if (!is.null(arv_col) && !is.na(arv_col) && nzchar(arv_col)) {
    print(table(phia[[arv_col]], useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- awareselfreported ---\n")
  asr_col <- grep("awareselfreported", names(phia), value = TRUE)[1]
  if (!is.null(asr_col) && !is.na(asr_col) && nzchar(asr_col)) {
    print(table(phia[[asr_col]], useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- centroids file preview ---\n")
  if (file.exists(cent_file_csv)) {
    cent <- read.csv(cent_file_csv, stringsAsFactors = FALSE)
    print(head(cent))
    print(names(cent))
  } else {
    cat("Centroid file not found.\n")
  }

  cat("\n--- household file district counts ---\n")
  if (file.exists(hh_file)) {
    hh <- read.csv(hh_file, stringsAsFactors = FALSE)
    print(grep("district|region|province|zone|area", names(hh), ignore.case = TRUE, value = TRUE))
    if ("district" %in% names(hh)) {
      print(table(hh$district, useNA = "always"))
    }
  } else {
    cat("Household file not found.\n")
  }

  cat("\n--- roster file district counts ---\n")
  if (file.exists(roster_file)) {
    roster <- read.csv(roster_file, stringsAsFactors = FALSE)
    print(grep("district|region|province|zone|area", names(roster), ignore.case = TRUE, value = TRUE))
    if ("district" %in% names(roster)) {
      print(table(roster$district, useNA = "always"))
    }
  } else {
    cat("Roster file not found.\n")
  }

  district_map <- c(
    "1" = "Berea",
    "2" = "Butha Buthe",
    "3" = "Leribe",
    "4" = "Mafeteng",
    "5" = "Maseru",
    "6" = "Mohale's Hoek",
    "7" = "Mokhotlong",
    "8" = "Qacha's Nek",
    "9" = "Quthing",
    "10" = "Thaba Tseka"
  )

  # hivstatusfinal: 1 = HIV positive, 2 = HIV negative; aware: 1 = aware, 2 = unaware
  phia_les <- phia %>%
    transmute(
      country = "Lesotho",
      personid = personid,
      province = recode(as.character(`district.bio`), !!!district_map, .default = NA_character_),
      district = province,
      age = `age.bio`,
      sex = case_when(
        `gender.bio` == 1 ~ "Male",
        `gender.bio` == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      age_group = case_when(
        age < 25 ~ "15-24",
        age < 50 ~ "25-49",
        age >= 50 ~ "50+",
        TRUE ~ NA_character_
      ),
      hivstatus = ifelse(`hivstatusfinal.bio` == 1, 1L,
                  ifelse(`hivstatusfinal.bio` == 2, 0L, NA_integer_)),
      aware_bin = ifelse(aware == 1, TRUE,
                  ifelse(aware == 2, FALSE, NA)),
      awareselfreported = awareselfreported,
      arvscurrent = arvscurrent,
      indweight = intwt0,
      btwt0 = btwt0,
      varstrat = `varstrat.bio`,
      varunit = `varunit.bio`
    )

  cat("\n--- Processed diagnostics ---\n")
  cat("Total rows:", nrow(phia_les), "\n")
  cat("HIV+ count:", sum(phia_les$hivstatus == 1, na.rm = TRUE), "\n")
  hiv_pos <- phia_les %>% filter(hivstatus == 1)
  cat("Among HIV+:\n")
  cat("  aware = TRUE:", sum(hiv_pos$aware_bin == TRUE, na.rm = TRUE), "\n")
  cat("  aware = FALSE:", sum(hiv_pos$aware_bin == FALSE, na.rm = TRUE), "\n")
  cat("  aware = NA:", sum(is.na(hiv_pos$aware_bin)), "\n")

  prev_dat <- phia_les %>% filter(!is.na(hivstatus), !is.na(btwt0), btwt0 > 0)
  weighted_prev <- with(prev_dat, sum(btwt0 * hivstatus, na.rm = TRUE) / sum(btwt0, na.rm = TRUE))
  cat(sprintf("Weighted HIV prevalence (btwt0): %.2f%%\n", 100 * weighted_prev))

  by_province <- phia_les %>%
    filter(!is.na(province)) %>%
    group_by(province) %>%
    summarise(
      n = n(),
      hiv_positive = sum(hivstatus == 1, na.rm = TRUE),
      weighted_prev = 100 * sum(btwt0 * hivstatus, na.rm = TRUE) / sum(btwt0[!is.na(hivstatus)], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
  cat("\nBy province:\n")
  print(as.data.frame(by_province))

  rhis_les_file <- file.path(project_dir, "data", "rhis", "standardised", "rhis_lesotho.rds")
  if (file.exists(rhis_les_file)) {
    rhis_les <- readRDS(rhis_les_file)
    phia_names <- sort(unique(phia_les$province))
    rhis_names <- sort(unique(rhis_les$province))
    cat("\nPHIA province names:\n")
    print(phia_names)
    cat("\nRHIS province names:\n")
    print(rhis_names)
    cat("\nNames in PHIA not in RHIS:\n")
    print(setdiff(phia_names, rhis_names))
    cat("\nNames in RHIS not in PHIA:\n")
    print(setdiff(rhis_names, phia_names))
  }

  out_dir <- file.path(project_dir, "data", "phia", "standardised")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(phia_les, file.path(out_dir, "phia_lesotho.rds"))

  cat("========================================\n")
  cat("LESOTHO PHIA: Review codings above.\n")
  cat("Saved phia_lesotho.rds. Stop before Zimbabwe.\n")
  cat("========================================\n")

  invisible(list(bio = bio, ind = ind, merged = phia, phia_les = phia_les))
}

standardise_zimbabwe_phia <- function() {
  zim_dir <- file.path(phia_dir, "zimbabwe")
  files <- print_file_inventory("zimbabwe", zim_dir)

  bio_file <- files[str_detect(basename(files), regex("adultbio.*\\.csv$", ignore_case = TRUE))][1]
  ind_file <- files[str_detect(basename(files), regex("adultind.*\\.csv$", ignore_case = TRUE))][1]
  cent_file <- files[str_detect(basename(files), regex("centroids.*\\.csv$", ignore_case = TRUE))][1]

  if (is.na(bio_file) || !file.exists(bio_file)) stop("Could not find Zimbabwe adultbio CSV.")
  if (is.na(ind_file) || !file.exists(ind_file)) stop("Could not find Zimbabwe adultind CSV.")

  bio <- print_phia_file_summary(bio_file, "biomarker")
  ind <- print_phia_file_summary(ind_file, "interview")
  print_centroid_summary(cent_file)
  print_variable_inventory(bio, "biomarker")
  print_variable_inventory(ind, "interview")

  cat("=== ZIMBABWE MERGE + CODING DIAGNOSTICS ===\n\n")
  bio_base <- read.csv(bio_file, stringsAsFactors = FALSE)
  ind_base <- read.csv(ind_file, stringsAsFactors = FALSE)
  phia <- merge(bio_base, ind_base, by = "personid", suffixes = c(".bio", ".ind"))

  cat("Merged rows:", nrow(phia), "\n")
  cat("Bio rows:", nrow(bio_base), "  Ind rows:", nrow(ind_base), "\n")

  dup_core <- names(phia)[grepl("\\.(bio|ind)$", names(phia))]
  cat("\nDuplicate-suffixed columns after merge:\n")
  if (length(dup_core) == 0) {
    cat("  None\n")
  } else {
    cat(paste(" ", dup_core, collapse = "\n"), "\n")
  }

  cat("\n--- hivstatusfinal ---\n")
  hiv_col <- grep("^hivstatusfinal(\\.bio|$)", names(phia), value = TRUE)[1]
  if (!is.na(hiv_col)) {
    print(table(phia[[hiv_col]], useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- aware ---\n")
  aware_col <- grep("^aware(\\.bio|$)|^aware$", names(phia), value = TRUE)[1]
  if (!is.na(aware_col)) {
    print(table(phia[[aware_col]], useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- gender ---\n")
  gender_col <- grep("^gender(\\.(bio|ind)|$)|^sex(\\.(bio|ind)|$)", names(phia), value = TRUE)[1]
  if (!is.na(gender_col)) {
    print(table(phia[[gender_col]], useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  geo_cols <- grep("province|district|region|zone|area", names(phia), ignore.case = TRUE, value = TRUE)
  cat("\n--- geography variables found ---\n")
  print(geo_cols)
  if (length(geo_cols) > 0) {
    for (gc in geo_cols) {
      cat("\n", gc, "\n", sep = "")
      x <- phia[[gc]]
      if (length(unique(x[!is.na(x)])) < 30) {
        print(table(x, useNA = "always"))
      } else {
        print(utils::head(sort(table(x), decreasing = TRUE), 20))
      }
    }
  }

  cat("\n--- age ---\n")
  age_col <- grep("^age(\\.(bio|ind)|$)", names(phia), value = TRUE)[1]
  if (!is.na(age_col)) {
    cat("Range:", paste(range(phia[[age_col]], na.rm = TRUE), collapse = " to "), "\n")
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- btwt0 ---\n")
  btwt_col <- grep("^btwt0$", names(phia), value = TRUE)[1]
  if (!is.na(btwt_col)) {
    cat("Non-zero:", sum(phia[[btwt_col]] > 0, na.rm = TRUE), "\n")
    cat("Zero:", sum(phia[[btwt_col]] == 0, na.rm = TRUE), "\n")
    cat("NA:", sum(is.na(phia[[btwt_col]])), "\n")
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- arvscurrent ---\n")
  arv_col <- grep("arvscurrent", names(phia), value = TRUE)[1]
  if (!is.na(arv_col)) {
    print(table(phia[[arv_col]], useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  cat("\n--- awareselfreported ---\n")
  asr_col <- grep("awareselfreported", names(phia), value = TRUE)[1]
  if (!is.na(asr_col)) {
    print(table(phia[[asr_col]], useNA = "always"))
  } else {
    cat("Variable not found.\n")
  }

  province_map <- c(
    "0" = "Bulawayo",
    "1" = "Manicaland",
    "2" = "Mashonaland Central",
    "3" = "Mashonaland East",
    "4" = "Mashonaland West",
    "5" = "Matabeleland North",
    "6" = "Matabeleland South",
    "7" = "Midlands",
    "8" = "Masvingo",
    "9" = "Harare"
  )

  phia_zim <- phia %>%
    transmute(
      country = "Zimbabwe",
      personid = personid,
      province = recode(as.character(`province.bio`), !!!province_map, .default = NA_character_),
      district = NA_character_,
      age = `age.bio`,
      sex = case_when(
        `gender.bio` == 1 ~ "Male",
        `gender.bio` == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      age_group = case_when(
        age < 25 ~ "15-24",
        age < 50 ~ "25-49",
        age >= 50 ~ "50+",
        TRUE ~ NA_character_
      ),
      hivstatus = ifelse(`hivstatusfinal.bio` == 1, 1L,
                  ifelse(`hivstatusfinal.bio` == 2, 0L, NA_integer_)),
      aware_bin = ifelse(aware == 1, TRUE,
                  ifelse(aware == 2, FALSE, NA)),
      awareselfreported = awareselfreported,
      arvscurrent = arvscurrent,
      indweight = intwt0,
      btwt0 = btwt0,
      varstrat = `varstrat.bio`,
      varunit = `varunit.bio`
    )

  cat("\n--- Processed diagnostics ---\n")
  cat("Total rows:", nrow(phia_zim), "\n")
  cat("HIV+ count:", sum(phia_zim$hivstatus == 1, na.rm = TRUE), "\n")
  hiv_pos <- phia_zim %>% filter(hivstatus == 1)
  cat("Among HIV+:\n")
  cat("  aware = TRUE:", sum(hiv_pos$aware_bin == TRUE, na.rm = TRUE), "\n")
  cat("  aware = FALSE:", sum(hiv_pos$aware_bin == FALSE, na.rm = TRUE), "\n")
  cat("  aware = NA:", sum(is.na(hiv_pos$aware_bin)), "\n")
  cat(sprintf("Awareness rate among HIV+: %.2f%%\n",
              100 * mean(hiv_pos$aware_bin, na.rm = TRUE)))

  prev_dat <- phia_zim %>% filter(!is.na(hivstatus), !is.na(btwt0), btwt0 > 0)
  weighted_prev <- with(prev_dat, sum(btwt0 * hivstatus, na.rm = TRUE) / sum(btwt0, na.rm = TRUE))
  cat(sprintf("Weighted HIV prevalence (btwt0): %.2f%%\n", 100 * weighted_prev))

  by_province <- phia_zim %>%
    filter(!is.na(province)) %>%
    group_by(province) %>%
    summarise(
      n = n(),
      hiv_positive = sum(hivstatus == 1, na.rm = TRUE),
      weighted_prev = 100 * sum(btwt0 * hivstatus, na.rm = TRUE) / sum(btwt0[!is.na(hivstatus)], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
  cat("\nBy province:\n")
  print(as.data.frame(by_province))

  cat("\nZimbabwe province mapping used:\n")
  print(as.data.frame(province_map))

  print_survey_validation(
    phia_zim,
    country_label = "ZIMPHIA 2020",
    published_prev = "12.9%",
    published_awareness = "86.8%",
    published_art_aware = "97.4%"
  )

  les_std_file <- file.path(project_dir, "data", "phia", "standardised", "phia_lesotho.rds")
  if (file.exists(les_std_file)) {
    phia_les <- readRDS(les_std_file)
    print_survey_validation(
      phia_les,
      country_label = "LEPHIA 2020",
      published_prev = "~21-23%",
      published_awareness = "~90%"
    )
  }

  rhis_zim_file <- file.path(project_dir, "data", "rhis", "standardised", "rhis_zimbabwe.rds")
  if (file.exists(rhis_zim_file)) {
    rhis_zim <- readRDS(rhis_zim_file)
    phia_names <- sort(unique(trimws(phia_zim$province)))
    rhis_names <- sort(unique(trimws(rhis_zim$province)))
    cat("\nPHIA province names:\n")
    print(phia_names)
    cat("\nRHIS province names (trimmed):\n")
    print(rhis_names)
    cat("\nNames in PHIA not in RHIS:\n")
    print(setdiff(phia_names, rhis_names))
    cat("\nNames in RHIS not in PHIA:\n")
    print(setdiff(rhis_names, phia_names))
  }

  out_dir <- file.path(project_dir, "data", "phia", "standardised")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(phia_zim, file.path(out_dir, "phia_zimbabwe.rds"))

  cat("\n========================================\n")
  cat("ZIMBABWE PHIA: Processed and saved phia_zimbabwe.rds.\n")
  cat("Stop before Malawi.\n")
  cat("========================================\n")

  invisible(list(bio = bio, ind = ind, merged = phia, phia_zim = phia_zim))
}

standardise_zimbabwe_phia()

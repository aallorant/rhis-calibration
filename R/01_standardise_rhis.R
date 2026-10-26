# ============================================================
# 01_standardise_rhis.R
# Standardise the routine HIV-testing registers (RHIS) of Zimbabwe, Lesotho
# and Malawi to one district x period x age-group x sex schema (n_tests,
# n_positive); Mozambique is taken from the estimation results of 03.
# Reads:  data/rhis/{zimbabwe,lesotho,malawi}/*.xls(x) (not in the release)
#         and output/mozambique/ (district estimation results from 03)
# Writes: data/rhis/standardised/rhis_{country}.rds and rhis_all_countries.rds
# Run:    Rscript R/01_standardise_rhis.R  (from the repository root)
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
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

  marker_rel <- file.path("data", "rhis", "zimbabwe", "HTS data Jan-Dec 2020 Final.xlsx")

  for (cand in candidates) {
    cand_norm <- tryCatch(normalizePath(cand, mustWork = FALSE), error = function(e) cand)
    if (file.exists(file.path(cand_norm, marker_rel))) {
      return(cand_norm)
    }
  }

  for (cand in candidates) {
    cand_norm <- tryCatch(normalizePath(cand, mustWork = FALSE), error = function(e) cand)
    if (dir.exists(file.path(cand_norm, "data", "rhis"))) {
      return(cand_norm)
    }
  }

  stop(
    "Could not locate the project root containing data/rhis. ",
    "Tried: ", paste(candidates, collapse = " | ")
  )
}

project_dir <- find_project_dir()

rhis_dir <- file.path(project_dir, "data", "rhis")
out_dir <- file.path(rhis_dir, "standardised")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

parse_zim_col <- function(cn) {
  cn_clean <- str_squish(cn)
  sex <- case_when(
    str_detect(cn_clean, regex("Female", ignore_case = TRUE)) ~ "Female",
    str_detect(cn_clean, regex("Male", ignore_case = TRUE)) ~ "Male",
    TRUE ~ NA_character_
  )

  nums <- str_extract_all(cn_clean, "\\d+")[[1]]
  nums <- suppressWarnings(as.integer(nums))

  age_lower <- case_when(
    str_detect(cn_clean, regex("months|<=", ignore_case = TRUE)) ~ 0L,
    length(nums) > 0 ~ nums[1],
    TRUE ~ NA_integer_
  )

  tibble(
    col_name = cn_clean,
    sex = sex,
    age_lower = age_lower,
    age_group = case_when(
      is.na(age_lower) ~ NA_character_,
      age_lower < 15 ~ "0-14",
      age_lower >= 15 & age_lower < 25 ~ "15-24",
      age_lower >= 25 & age_lower < 50 ~ "25-49",
      age_lower >= 50 ~ "50+"
    )
  )
}

standardise_zimbabwe <- function() {
  zim_files <- c(
    file.path(rhis_dir, "zimbabwe", "HTS data Jan-Dec 2020 Final.xlsx"),
    file.path(rhis_dir, "zimbabwe", "Jan - December HTS 2021 - Final.xlsx"),
    file.path(rhis_dir, "zimbabwe", "Copy of HTS Data Jan-De 2022 Final.xlsx")
  )

  cat("=== ZIMBABWE ===\n\n")
  cat("All workbook sheet names:\n")
  for (f in zim_files) {
    cat(" ", basename(f), ":\n", sep = "")
    cat("   ", paste(excel_sheets(f), collapse = ", "), "\n", sep = "")
  }
  cat("\n")

  test_df <- read_excel(zim_files[1], sheet = "Received")
  names(test_df) <- str_squish(names(test_df))

  cat("\nAll column names (Received sheet):\n")
  for (i in seq_along(names(test_df))) {
    cat(sprintf("  [%02d] '%s'\n", i, names(test_df)[i]))
  }

  id_cols <- c("Province", "District", "Facility", "Period")
  agg_patterns <- c(
    "^\\d+-\\d+\\s*Yrs$",
    "^\\d+\\+\\s*Yrs$",
    "^Male$",
    "^Female$",
    "^Total$"
  )

  all_cols <- names(test_df)
  is_agg <- vapply(all_cols, function(x) {
    any(vapply(agg_patterns, function(p) grepl(p, str_trim(x)), logical(1)))
  }, logical(1))
  agg_cols <- all_cols[is_agg]

  cat("\nAggregate columns to exclude:", paste(agg_cols, collapse = ", "), "\n")

  age_sex_cols <- setdiff(all_cols, c(id_cols, agg_cols))
  cat("Individual age-sex columns:", length(age_sex_cols), "\n")

  col_map <- bind_rows(lapply(age_sex_cols, parse_zim_col))

  cat("\nColumn map:\n")
  print(col_map %>% select(col_name, sex, age_lower, age_group))

  unparsed <- col_map %>% filter(is.na(sex) | is.na(age_group))
  if (nrow(unparsed) > 0) {
    cat("\n*** WARNING: unparsed columns ***\n")
    print(unparsed)
  }

  adult_cols <- col_map %>%
    filter(age_group %in% c("15-24", "25-49", "50+"), !is.na(sex))

  cat("\nAdult age-sex columns:", nrow(adult_cols), "\n")

  pivot_sheet <- function(filepath, sheet_name, value_name) {
    df <- read_excel(filepath, sheet = sheet_name)
    names(df) <- str_squish(names(df))

    missing_cols <- setdiff(adult_cols$col_name, names(df))
    if (length(missing_cols) > 0) {
      cat(
        "  WARNING: missing columns in", basename(filepath),
        "sheet", sheet_name, ":", paste(missing_cols, collapse = ", "), "\n"
      )
    }

    use_cols <- intersect(adult_cols$col_name, names(df))

    df %>%
      select(all_of(c(id_cols, use_cols))) %>%
      mutate(across(all_of(use_cols), ~ replace_na(suppressWarnings(as.numeric(.)), 0))) %>%
      pivot_longer(
        cols = all_of(use_cols),
        names_to = "col_name",
        values_to = value_name
      ) %>%
      left_join(adult_cols %>% select(col_name, sex, age_group), by = "col_name") %>%
      group_by(Province, District, Period, sex, age_group) %>%
      summarise(!!value_name := sum(.data[[value_name]], na.rm = TRUE), .groups = "drop")
  }

  all_years <- vector("list", length(zim_files))

  for (i in seq_along(zim_files)) {
    f <- zim_files[i]
    sheets <- excel_sheets(f)
    cat("\nProcessing:", basename(f), "\n")
    cat("  Sheets:", paste(sheets, collapse = ", "), "\n")

    test_sheet <- sheets[str_detect(sheets, regex("^Received$", ignore_case = TRUE))]
    pos_sheet <- sheets[str_detect(sheets, regex("^PosTot$", ignore_case = TRUE))]
    if (length(pos_sheet) == 0) {
      pos_sheet <- sheets[str_detect(sheets, regex("^Positives$", ignore_case = TRUE))]
    }
    if (length(pos_sheet) == 0) {
      pos_sheet <- sheets[str_detect(sheets, regex("^Pos", ignore_case = TRUE))]
    }

    if (length(test_sheet) == 0) test_sheet <- sheets[1]
    if (length(pos_sheet) == 0) pos_sheet <- sheets[2]

    test_sheet <- test_sheet[1]
    pos_sheet <- pos_sheet[1]

    cat("  Using sheet", test_sheet, "for tests and", pos_sheet, "for positives\n")

    tests <- pivot_sheet(f, test_sheet, "n_tests")
    pos <- pivot_sheet(f, pos_sheet, "n_positive")

    joined <- tests %>%
      select(Province, District, Period, sex, age_group, n_tests) %>%
      left_join(
        pos %>% select(Province, District, Period, sex, age_group, n_positive),
        by = c("Province", "District", "Period", "sex", "age_group")
      ) %>%
      mutate(n_positive = replace_na(n_positive, 0))

    total_tests <- sum(joined$n_tests, na.rm = TRUE)
    total_pos <- sum(joined$n_positive, na.rm = TRUE)
    rate <- if (total_tests > 0) 100 * total_pos / total_tests else NA_real_

    cat(
      "  Rows:", nrow(joined),
      " Total tests:", format(total_tests, big.mark = ","),
      " Total pos:", format(total_pos, big.mark = ","),
      " Rate:", sprintf("%.1f%%", rate), "\n"
    )

    all_years[[i]] <- joined
  }

  zim_raw <- bind_rows(all_years) %>%
    mutate(
      Period = str_squish(Period),
      period_date = suppressWarnings(parse_date_time(Period, orders = c("B Y", "b Y", "my", "m/Y", "Ym"))),
      year = year(period_date),
      quarter = quarter(period_date),
      period = paste(year, paste0("Q", quarter)),
      country = "Zimbabwe",
      province_raw = str_squish(Province),
      province = case_when(
        province_raw == "Bulawayo Central Hospitals" ~ "Bulawayo",
        province_raw == "Harare Central Hospitals" ~ "Harare",
        province_raw == "Chitungwiza City" ~ "Harare",
        TRUE ~ province_raw
      ),
      district = str_squish(District),
      testing_point = NA_character_
    ) %>%
    filter(!is.na(year), !is.na(quarter)) %>%
    select(country, province_raw, province, district, year, quarter, period, age_group, sex, testing_point, n_tests, n_positive)

  expected_provinces <- c(
    "Bulawayo", "Harare", "Manicaland", "Mashonaland Central",
    "Mashonaland East", "Mashonaland West", "Masvingo",
    "Matabeleland North", "Matabeleland South", "Midlands"
  )

  extra_provinces <- zim_raw %>%
    count(province_raw, name = "rows") %>%
    filter(!province_raw %in% expected_provinces)

  cat("\n--- Province values before harmonisation ---\n")
  print(as.data.frame(zim_raw %>% count(province_raw, name = "rows") %>% arrange(desc(rows))))

  if (nrow(extra_provinces) > 0) {
    cat("\n--- Extra province values to merge ---\n")
    print(as.data.frame(extra_provinces))
  }

  zim <- zim_raw %>%
    group_by(country, province, district, year, quarter, period, age_group, sex, testing_point) %>%
    summarise(
      n_tests = sum(n_tests, na.rm = TRUE),
      n_positive = sum(n_positive, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_tests > 0)

  over_pos <- zim %>%
    filter(n_positive > n_tests) %>%
    arrange(desc(n_positive - n_tests), province, district, period, age_group, sex)

  cat("\n--- Rows where n_positive > n_tests ---\n")
  if (nrow(over_pos) == 0) {
    cat("None\n")
  } else {
    print(as.data.frame(over_pos %>%
      select(province, district, period, age_group, sex, n_tests, n_positive)))
  }

  zim$n_positive <- pmin(zim$n_positive, zim$n_tests)

  cat("\n========================================\n")
  cat("ZIMBABWE DIAGNOSTICS\n")
  cat("========================================\n\n")

  cat("Provinces (", n_distinct(zim$province), "):\n", sep = "")
  cat(paste(" ", sort(unique(zim$province)), collapse = "\n"), "\n\n")

  cat("Districts:", n_distinct(zim$district), "\n\n")
  cat("Years:", paste(sort(unique(zim$year)), collapse = ", "), "\n\n")
  cat("Age groups:", paste(sort(unique(zim$age_group)), collapse = ", "), "\n")
  cat("Sexes:", paste(sort(unique(zim$sex)), collapse = ", "), "\n\n")

  cat("--- Annual totals ---\n")
  annual <- zim %>%
    group_by(year) %>%
    summarise(
      n_tests = sum(n_tests),
      n_positive = sum(n_positive),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    )
  print(as.data.frame(annual))

  cat("\n--- By province (all years pooled) ---\n")
  by_prov <- zim %>%
    group_by(province) %>%
    summarise(
      districts = n_distinct(district),
      n_tests = sum(n_tests),
      n_positive = sum(n_positive),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    ) %>%
    arrange(desc(n_tests))
  print(as.data.frame(by_prov))

  cat("\n--- By age group x sex (all years pooled) ---\n")
  by_demo <- zim %>%
    group_by(age_group, sex) %>%
    summarise(
      n_tests = sum(n_tests),
      n_positive = sum(n_positive),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    )
  print(as.data.frame(by_demo))

  cat("\n--- Province x age group positivity (2020 only) ---\n")
  prov_age_2020 <- zim %>%
    filter(year == 2020) %>%
    group_by(province, age_group) %>%
    summarise(
      positivity = sprintf("%.1f%%", 100 * sum(n_positive) / sum(n_tests)),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = age_group, values_from = positivity)
  print(as.data.frame(prov_age_2020))

  cat("\n--- Any NAs? ---\n")
  cat("NA province:", sum(is.na(zim$province)), "\n")
  cat("NA district:", sum(is.na(zim$district)), "\n")
  cat("NA age_group:", sum(is.na(zim$age_group)), "\n")
  cat("NA sex:", sum(is.na(zim$sex)), "\n")
  cat("Zero-test rows:", sum(zim$n_tests == 0), "\n")
  cat("Negative test rows:", sum(zim$n_tests < 0), "\n")
  cat("n_positive > n_tests:", sum(zim$n_positive > zim$n_tests), "\n")

  cat("\n========================================\n")
  cat("STOP HERE. Review the diagnostics above.\n")
  cat("Confirm before proceeding to Lesotho.\n")
  cat("========================================\n")

  zim
}

standardise_lesotho <- function() {
  les_file <- file.path(rhis_dir, "lesotho", "HTS_TST, POS & HTS_SELF-v2.xls")

  cat("\n\n=== LESOTHO ===\n\n")

  sheets <- excel_sheets(les_file)
  cat("Sheets:", paste(sheets, collapse = ", "), "\n\n")

  for (i in seq_along(sheets)) {
    cat("--- Sheet", i, ":", sheets[i], "---\n")
    df <- read_excel(les_file, sheet = i, n_max = 5)
    cat("  Columns:", paste(names(df), collapse = ", "), "\n")
    cat("  Rows (total):", nrow(read_excel(les_file, sheet = i)), "\n\n")
  }

  les_raw <- read_excel(les_file, sheet = 1)

  cat("Columns:", paste(names(les_raw), collapse = ", "), "\n\n")
  cat("Districts:", paste(sort(unique(les_raw$District)), collapse = ", "), "\n\n")
  cat("Sex values:", paste(unique(les_raw$Sex), collapse = ", "), "\n\n")
  cat("Age values:\n")
  cat(paste(" ", sort(unique(les_raw$Age)), collapse = "\n"), "\n\n")
  cat("Period IDs:", paste(sort(unique(les_raw$periodid)), collapse = ", "), "\n\n")

  les <- les_raw %>%
    mutate(
      age_lower = case_when(
        grepl("^<1", Age) ~ 0L,
        grepl("^1-4", Age) ~ 1L,
        TRUE ~ suppressWarnings(as.integer(str_extract(Age, "^\\d+")))
      ),
      age_group = case_when(
        is.na(age_lower) ~ NA_character_,
        age_lower < 15 ~ "0-14",
        age_lower >= 15 & age_lower < 25 ~ "15-24",
        age_lower >= 25 & age_lower < 50 ~ "25-49",
        age_lower >= 50 ~ "50+"
      ),
      sex = case_when(
        grepl("Female", Sex, ignore.case = TRUE) ~ "Female",
        grepl("Male", Sex, ignore.case = TRUE) ~ "Male",
        TRUE ~ NA_character_
      ),
      n_positive = suppressWarnings(as.numeric(`HIV Positive (Inclusive)`)),
      n_negative = suppressWarnings(as.numeric(`HIV Negative (Inclusive)`)),
      n_positive = replace_na(n_positive, 0),
      n_negative = replace_na(n_negative, 0),
      n_tests = n_positive + n_negative,
      year = suppressWarnings(as.integer(str_extract(periodid, "^\\d{4}"))),
      quarter = suppressWarnings(as.integer(str_extract(periodid, "\\d$"))),
      period = paste(year, paste0("Q", quarter)),
      country = "Lesotho",
      province = str_squish(District),
      district = str_squish(District),
      testing_point = NA_character_
    ) %>%
    filter(age_group %in% c("15-24", "25-49", "50+"), !is.na(sex)) %>%
    group_by(country, province, district, year, quarter, period, age_group, sex, testing_point) %>%
    summarise(
      n_tests = sum(n_tests, na.rm = TRUE),
      n_positive = sum(n_positive, na.rm = TRUE),
      .groups = "drop"
    )

  les$n_positive <- pmin(les$n_positive, les$n_tests)

  cat("\n========================================\n")
  cat("LESOTHO DIAGNOSTICS\n")
  cat("========================================\n\n")

  cat("Districts/provinces (", n_distinct(les$province), "):\n", sep = "")
  cat(paste(" ", sort(unique(les$province)), collapse = "\n"), "\n\n")

  cat("Years:", paste(sort(unique(les$year)), collapse = ", "), "\n\n")
  cat("Age groups:", paste(sort(unique(les$age_group)), collapse = ", "), "\n")
  cat("Sexes:", paste(sort(unique(les$sex)), collapse = ", "), "\n\n")

  cat("--- Annual totals ---\n")
  annual <- les %>%
    group_by(year) %>%
    summarise(
      n_tests = sum(n_tests),
      n_positive = sum(n_positive),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    )
  print(as.data.frame(annual))

  cat("\n--- By district (all years pooled) ---\n")
  by_dist <- les %>%
    group_by(province) %>%
    summarise(
      n_tests = sum(n_tests),
      n_positive = sum(n_positive),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    ) %>%
    arrange(desc(n_tests))
  print(as.data.frame(by_dist))

  cat("\n--- By age group x sex (all years pooled) ---\n")
  by_demo <- les %>%
    group_by(age_group, sex) %>%
    summarise(
      n_tests = sum(n_tests),
      n_positive = sum(n_positive),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    )
  print(as.data.frame(by_demo))

  cat("\n--- District x age group positivity (2020 only) ---\n")
  dist_age_2020 <- les %>%
    filter(year == 2020) %>%
    group_by(province, age_group) %>%
    summarise(
      positivity = sprintf("%.1f%%", 100 * sum(n_positive) / sum(n_tests)),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = age_group, values_from = positivity)
  print(as.data.frame(dist_age_2020))

  cat("\n--- Any NAs? ---\n")
  cat("NA province:", sum(is.na(les$province)), "\n")
  cat("NA age_group:", sum(is.na(les$age_group)), "\n")
  cat("NA sex:", sum(is.na(les$sex)), "\n")
  cat("Zero-test rows:", sum(les$n_tests == 0), "\n")
  cat("n_positive > n_tests:", sum(les$n_positive > les$n_tests), "\n")

  cat("\n========================================\n")
  cat("STOP HERE. Review the diagnostics above.\n")
  cat("Confirm before proceeding to Malawi.\n")
  cat("========================================\n")

  les
}

standardise_malawi <- function() {
  mwi_file <- file.path(rhis_dir, "malawi", "Malawi HTS HST ART program data 2019-2022.xlsx")

  cat("\n\n=== MALAWI ===\n\n")

  sheets <- excel_sheets(mwi_file)
  cat("Sheets:", paste(sheets, collapse = ", "), "\n\n")

  for (i in seq_along(sheets)) {
    cat("--- Sheet", i, ":", sheets[i], "---\n")
    df <- read_excel(mwi_file, sheet = i, n_max = 3)
    cat("  Columns:", paste(names(df), collapse = ", "), "\n\n")
  }

  mwi_raw <- read_excel(mwi_file, sheet = "HTS")

  cat("All columns:\n")
  for (i in seq_along(names(mwi_raw))) {
    cat(sprintf("  [%02d] '%s'\n", i, names(mwi_raw)[i]))
  }

  cat("\nRegions:", paste(sort(unique(mwi_raw$Region)), collapse = ", "), "\n")
  cat("Zones:", paste(sort(unique(mwi_raw$Zone)), collapse = ", "), "\n")
  cat("Districts (n):", n_distinct(mwi_raw$district), "\n")
  cat("Periods:", paste(sort(unique(mwi_raw$Period)), collapse = ", "), "\n")
  cat("Sub groups:", paste(sort(unique(mwi_raw$`Sub groups`)), collapse = ", "), "\n")

  mwi_raw <- mwi_raw %>%
    mutate(
      testing_point = case_when(
        grepl("health facility|static", `Sub groups`, ignore.case = TRUE) ~ "facility",
        grepl("community", `Sub groups`, ignore.case = TRUE) ~ "community",
        grepl("stand-alone|HTC sites", `Sub groups`, ignore.case = TRUE) ~ "htc_standalone",
        grepl("self-test|self test", `Sub groups`, ignore.case = TRUE) ~ "hivst_return",
        TRUE ~ "other"
      )
    )

  cat("\nTesting point mapping:\n")
  print(as.data.frame(table(Original = mwi_raw$`Sub groups`, Mapped = mwi_raw$testing_point)))

  mwi <- mwi_raw %>%
    transmute(
      Region,
      Zone,
      district = str_squish(district),
      Period = str_squish(Period),
      testing_point,
      males_tested = replace_na(suppressWarnings(as.numeric(`Males tested`)), 0),
      females_nonpreg = replace_na(suppressWarnings(as.numeric(`Females non-pregnant`)), 0),
      females_preg = replace_na(suppressWarnings(as.numeric(`Females pregnant`)), 0),
      pos_male = replace_na(suppressWarnings(as.numeric(`New positive male`)), 0),
      pos_female = replace_na(suppressWarnings(as.numeric(`New positive female`)), 0)
    ) %>%
    mutate(
      year = suppressWarnings(as.integer(str_extract(Period, "\\d{4}"))),
      quarter = suppressWarnings(as.integer(str_extract(Period, "(?<=Q)\\d"))),
      period = paste(year, paste0("Q", quarter)),
      country = "Malawi",
      province = str_squish(Zone)
    ) %>%
    select(country, province, district, year, quarter, period, testing_point, males_tested, females_nonpreg, females_preg, pos_male, pos_female) %>%
    pivot_longer(
      cols = c(males_tested, females_nonpreg, females_preg, pos_male, pos_female),
      names_to = "measure",
      values_to = "value"
    ) %>%
    mutate(
      sex = case_when(
        measure %in% c("males_tested", "pos_male") ~ "Male",
        TRUE ~ "Female"
      ),
      metric = case_when(
        measure %in% c("males_tested", "females_nonpreg", "females_preg") ~ "n_tests",
        TRUE ~ "n_positive"
      )
    ) %>%
    group_by(country, province, district, year, quarter, period, sex, testing_point, metric) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = metric, values_from = value, values_fill = 0) %>%
    mutate(age_group = "all_ages") %>%
    filter(n_tests > 0) %>%
    select(country, province, district, year, quarter, period, age_group, sex, testing_point, n_tests, n_positive)

  mwi$n_positive <- pmin(mwi$n_positive, mwi$n_tests)

  cat("\n========================================\n")
  cat("MALAWI DIAGNOSTICS\n")
  cat("========================================\n\n")

  cat("Zones/provinces (", n_distinct(mwi$province), "):\n", sep = "")
  cat(paste(" ", sort(unique(mwi$province)), collapse = "\n"), "\n\n")

  cat("Districts:", n_distinct(mwi$district), "\n\n")
  cat("Districts per zone:\n")
  print(as.data.frame(mwi %>% distinct(province, district) %>% count(province, name = "n_districts")))

  cat("\nYears:", paste(sort(unique(mwi$year)), collapse = ", "), "\n\n")
  cat("Testing points:", paste(sort(unique(mwi$testing_point)), collapse = ", "), "\n\n")

  cat("--- Annual totals ---\n")
  annual <- mwi %>%
    group_by(year) %>%
    summarise(
      n_tests = sum(n_tests),
      n_positive = sum(n_positive),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    )
  print(as.data.frame(annual))

  cat("\n--- By zone (all years pooled) ---\n")
  by_zone <- mwi %>%
    group_by(province) %>%
    summarise(
      districts = n_distinct(district),
      n_tests = sum(n_tests),
      n_positive = sum(n_positive),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    ) %>%
    arrange(desc(n_tests))
  print(as.data.frame(by_zone))

  cat("\n--- By sex (all years pooled) ---\n")
  by_sex <- mwi %>%
    group_by(sex) %>%
    summarise(
      n_tests = sum(n_tests),
      n_positive = sum(n_positive),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    )
  print(as.data.frame(by_sex))

  cat("\n--- By testing point (2020 only) ---\n")
  by_tp <- mwi %>%
    filter(year == 2020) %>%
    group_by(testing_point) %>%
    summarise(
      n_tests = sum(n_tests),
      n_positive = sum(n_positive),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    ) %>%
    arrange(desc(n_tests))
  print(as.data.frame(by_tp))

  cat("\n--- Any NAs? ---\n")
  cat("NA province (zone):", sum(is.na(mwi$province)), "\n")
  cat("NA district:", sum(is.na(mwi$district)), "\n")
  cat("NA sex:", sum(is.na(mwi$sex)), "\n")
  cat("Zero-test rows:", sum(mwi$n_tests == 0), "\n")
  cat("n_positive > n_tests:", sum(mwi$n_positive > mwi$n_tests), "\n")

  cat("\n========================================\n")
  cat("STOP HERE. Review the diagnostics above.\n")
  cat("Confirm before proceeding to Mozambique.\n")
  cat("========================================\n")

  mwi
}

standardise_mozambique <- function() {
  est_dat <- readRDS(file.path(project_dir, "output", "mozambique", "estimation_results.rds"))

  cat("\n\n=== MOZAMBIQUE ===\n\n")

  moz <- est_dat %>%
    filter(has_rhis == TRUE) %>%
    transmute(
      country = "Mozambique",
      province = str_squish(province),
      district = str_squish(area_name),
      year = 2021L,
      quarter = NA_integer_,
      period = "2021",
      age_group = "all_ages",
      sex = "both",
      n_tests = as.numeric(n_d),
      n_positive = as.numeric(ny_d),
      testing_point = NA_character_
    )

  cat("Provinces (", n_distinct(moz$province), "):\n", sep = "")
  cat(paste(" ", sort(unique(moz$province)), collapse = "\n"), "\n\n")

  cat("Districts:", n_distinct(moz$district), "\n\n")

  cat("--- Annual totals ---\n")
  annual <- moz %>%
    group_by(year) %>%
    summarise(
      n_tests = sum(n_tests, na.rm = TRUE),
      n_positive = sum(n_positive, na.rm = TRUE),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    )
  print(as.data.frame(annual))

  cat("\n--- By province ---\n")
  by_prov <- moz %>%
    group_by(province) %>%
    summarise(
      districts = n_distinct(district),
      n_tests = sum(n_tests, na.rm = TRUE),
      n_positive = sum(n_positive, na.rm = TRUE),
      positivity = sprintf("%.2f%%", 100 * n_positive / n_tests),
      .groups = "drop"
    ) %>%
    arrange(desc(n_tests))
  print(as.data.frame(by_prov))

  cat("\n--- Any NAs? ---\n")
  cat("NA province:", sum(is.na(moz$province)), "\n")
  cat("NA district:", sum(is.na(moz$district)), "\n")
  cat("NA year:", sum(is.na(moz$year)), "\n")
  cat("NA age_group:", sum(is.na(moz$age_group)), "\n")
  cat("NA sex:", sum(is.na(moz$sex)), "\n")
  cat("Zero-test rows:", sum(moz$n_tests == 0, na.rm = TRUE), "\n")
  cat("n_positive > n_tests:", sum(moz$n_positive > moz$n_tests, na.rm = TRUE), "\n")

  moz
}

zim <- standardise_zimbabwe()
saveRDS(zim, file.path(out_dir, "rhis_zimbabwe.rds"))

les <- standardise_lesotho()
saveRDS(les, file.path(out_dir, "rhis_lesotho.rds"))

mwi <- standardise_malawi()
saveRDS(mwi, file.path(out_dir, "rhis_malawi.rds"))

moz <- standardise_mozambique()
saveRDS(moz, file.path(out_dir, "rhis_mozambique.rds"))

all_rhis <- bind_rows(
  readRDS(file.path(out_dir, "rhis_zimbabwe.rds")),
  readRDS(file.path(out_dir, "rhis_lesotho.rds")),
  readRDS(file.path(out_dir, "rhis_malawi.rds")),
  moz
)

saveRDS(all_rhis, file.path(out_dir, "rhis_all_countries.rds"))

cat("\n\n=== ALL COUNTRIES SUMMARY ===\n\n")
summary_tbl <- all_rhis %>%
  group_by(country) %>%
  summarise(
    provinces = n_distinct(province),
    districts = n_distinct(district),
    years = paste(sort(unique(year)), collapse = "-"),
    total_tests = format(sum(n_tests, na.rm = TRUE), big.mark = ","),
    total_positive = format(sum(n_positive, na.rm = TRUE), big.mark = ","),
    positivity = sprintf("%.2f%%", 100 * sum(n_positive, na.rm = TRUE) / sum(n_tests, na.rm = TRUE)),
    .groups = "drop"
  )
print(as.data.frame(summary_tbl))

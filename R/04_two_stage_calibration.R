#!/usr/bin/env Rscript
# 04_two_stage_calibration.R: two-stage RHIS-first calibration, one pipeline for
# every country. Stage 1: province calibration (TMB intercept-only national alpha
# via tmb/full_model.cpp; province lambda_p from the raw register-survey gaps).
# Stage 2: district prediction with the province lambda_p and propagated uncertainty.
# Reads:  data/rhis/standardised/rhis_<country>.rds, data/phia/standardised/phia_<country>.rds
# Writes: output/<country>/{province_calibration,district_predictions,tmb_province_fit}.rds,
#         output/two_stage_results.rds, output/two_stage_all_districts.rds
# Run:    Rscript R/04_two_stage_calibration.R   (from the repository root)

suppressPackageStartupMessages({
  library(TMB)
  library(dplyr)
  library(tidyr)
  library(survey)
  library(stringr)
})

options(survey.lonely.psu = "adjust")

# ── Project root ─────────────────────────────────────────────
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(path) > 0) return(dirname(normalizePath(path)))
  normalizePath(getwd())
}

find_project_dir <- function() {
  candidates <- unique(c(
    get_script_dir(), dirname(get_script_dir()),
    getwd(), file.path(getwd(), "rhis-calibration")
  ))
  for (d in candidates) {
    d <- tryCatch(normalizePath(d, mustWork = FALSE), error = function(e) d)
    if (dir.exists(file.path(d, "data", "rhis", "standardised"))) return(d)
  }
  stop("Cannot locate project root.")
}

proj <- find_project_dir()

# ── Helpers ──────────────────────────────────────────────────
logit    <- function(p) log(p / (1 - p))
se_logit <- function(p, se_p) se_p / (p * (1 - p))

safe_logodds <- function(pos, neg) {
  log((pos + 0.5) / (neg + 0.5))
}

safe_var_logodds <- function(pos, neg) {
  1 / (pos + 0.5) + 1 / (neg + 0.5)
}

# ── RHIS year closest to each country's PHIA ─────────────────
rhis_year <- c(
  Lesotho    = 2020L,
  Malawi     = 2020L,
  Mozambique = 2021L,
  Zimbabwe   = 2020L
)


# ---- PART 1: Province-level data assembly ----

# ── 1a. RHIS: province AND district level summaries ──────────

rhis_summaries <- function(country, year) {
  path <- file.path(proj, "data", "rhis", "standardised",
                    paste0("rhis_", tolower(country), ".rds"))
  if (!file.exists(path)) stop("Missing RHIS: ", path)
  rhis <- readRDS(path)

  rhis <- rhis %>% filter(year == !!year)
  if (nrow(rhis) == 0) stop("No RHIS data for ", country, " year ", year)

  # Province-level
  prov <- rhis %>%
    group_by(province) %>%
    summarise(
      n_tests          = sum(n_tests,    na.rm = TRUE),
      n_positive       = sum(n_positive, na.rm = TRUE),
      n_districts_rhis = n_distinct(district),
      .groups = "drop"
    ) %>%
    mutate(
      n_negative = n_tests - n_positive,
      z_p = safe_logodds(n_positive, n_negative),
      v_p = safe_var_logodds(n_positive, n_negative)
    )

  # District-level
  dist <- rhis %>%
    group_by(province, district) %>%
    summarise(
      n_tests    = sum(n_tests,    na.rm = TRUE),
      n_positive = sum(n_positive, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      n_negative = n_tests - n_positive,
      z_d = safe_logodds(n_positive, n_negative),
      v_d = safe_var_logodds(n_positive, n_negative)
    )

  list(province = prov, district = dist)
}


# ── 1b. PHIA: province-level survey estimates ────────────────
# btwt0: biomarker weight (aware_bin is defined for biomarker respondents only)

phia_province_survey <- function(country) {
  path <- file.path(proj, "data", "phia", "standardised",
                    paste0("phia_", tolower(country), ".rds"))
  if (!file.exists(path)) stop("Missing PHIA: ", path)
  phia <- readRDS(path) %>% filter(age >= 15, age <= 64)

  phia <- phia %>%
    mutate(
      hiv_pos     = as.integer(hivstatus == 1L),
      hiv_unaware = as.integer(hivstatus == 1L & aware_bin == FALSE)
    )

  des <- svydesign(
    id      = ~varunit,
    strata  = ~varstrat,
    weights = ~btwt0,
    data    = phia,
    nest    = TRUE
  )

  # P(HIV+) by province
  theta_raw <- svyby(~hiv_pos, ~province, des, svymean, na.rm = TRUE)
  theta_tbl <- theta_raw %>%
    transmute(province, theta_N = hiv_pos, se_theta = se)

  # P(aware | HIV+) by province
  des_pos <- subset(des, hivstatus == 1L)
  kappa_raw <- svyby(~aware_bin, ~province, des_pos, svymean, na.rm = TRUE)
  kappa_tbl <- kappa_raw %>%
    transmute(province, kappa = aware_binTRUE, se_kappa = se.aware_binTRUE)

  # P(HIV+ & unaware) by province: the estimand
  zhat_raw <- svyby(~hiv_unaware, ~province, des, svymean, na.rm = TRUE)
  zhat_tbl <- zhat_raw %>%
    transmute(province, Z_hat_prob = hiv_unaware, se_Z_hat = se)

  # Combine and compute logit(Z_hat) + delta-method variance
  theta_tbl %>%
    left_join(kappa_tbl, by = "province") %>%
    left_join(zhat_tbl,  by = "province") %>%
    mutate(
      Z_hat_p    = logit(Z_hat_prob),
      se_logit_Z = se_logit(Z_hat_prob, se_Z_hat),
      psi_p      = se_logit_Z^2
    )
}


# ── 1c. Assemble one country ─────────────────────────────────

assemble_country <- function(country) {
  yr <- rhis_year[country]
  cat(sprintf("\n== %s (RHIS year: %d) ==\n", country, yr))

  # RHIS
  rhis <- rhis_summaries(country, yr)

  # PHIA
  phia <- phia_province_survey(country)

  # Harmonise province names for join
  rhis$province$province <- str_squish(rhis$province$province)
  phia$province <- str_squish(phia$province)

  # Join province-level
  prov_cal <- rhis$province %>%
    inner_join(phia, by = "province") %>%
    mutate(
      lambda_p     = z_p - Z_hat_p,         # raw gap = log(lambda_p)
      var_lambda_p = v_p + psi_p,            # uncertainty in lambda
      se_lambda_p  = sqrt(var_lambda_p),
      OR_p         = exp(lambda_p),
      country      = country
    )

  # Diagnostics
  rhis_provs <- sort(unique(rhis$province$province))
  phia_provs <- sort(unique(phia$province))
  only_rhis  <- setdiff(rhis_provs, phia_provs)
  only_phia  <- setdiff(phia_provs, rhis_provs)

  cat("Provinces: RHIS=", length(rhis_provs),
      " PHIA=", length(phia_provs),
      " Joined=", nrow(prov_cal), "\n")
  if (length(only_rhis) > 0)
    cat("  WARNING in RHIS not PHIA:", paste(only_rhis, collapse = ", "), "\n")
  if (length(only_phia) > 0)
    cat("  WARNING in PHIA not RHIS:", paste(only_phia, collapse = ", "), "\n")

  cat("\nProvince-level calibration:\n")
  prov_display <- prov_cal %>%
    transmute(
      province,
      n_tests,
      positivity = round(n_positive / n_tests, 4),
      theta_N    = round(theta_N, 4),
      kappa      = round(kappa, 3),
      Z_hat_prob = round(Z_hat_prob, 4),
      Z_hat_p    = round(Z_hat_p, 3),
      psi_p      = round(psi_p, 4),
      lambda_p   = round(lambda_p, 3),
      OR_p       = round(OR_p, 2)
    )
  print(as.data.frame(prov_display), row.names = FALSE)

  # District-level RHIS with province lookup
  dist_rhis <- rhis$district %>%
    mutate(country = country) %>%
    left_join(
      prov_cal %>% select(province, lambda_p, var_lambda_p, OR_p,
                           Z_hat_p, psi_p, theta_N, kappa),
      by = "province"
    )

  n_orphan <- sum(is.na(dist_rhis$lambda_p))
  if (n_orphan > 0)
    cat("\n  WARNING:", n_orphan, "districts with no province calibration\n")

  list(
    province_cal  = prov_cal,
    district_rhis = dist_rhis,
    country       = country,
    rhis_year     = yr
  )
}


# ---- PART 2: TMB province-level fit (Stage 1a) ----

compile_tmb <- function() {
  tmb_path <- file.path(proj, "tmb", "full_model.cpp")
  if (!file.exists(tmb_path)) stop("Cannot find tmb/full_model.cpp")
  if (!is.loaded("full_model")) {
    compile(tmb_path)
    dyn.load(dynlib(sub("\\.cpp$", "", tmb_path)))
  }
}

fit_tmb_province <- function(prov_cal) {
  # Intercept-only TMB: same specification for all countries
  n_p <- nrow(prov_cal)

  data_tmb <- list(
    z        = prov_cal$z_p,
    v        = prov_cal$v_p,
    Z_hat    = prov_cal$Z_hat_p,
    psi      = prov_cal$psi_p,
    d_rhis   = as.integer(0:(n_p - 1)),
    d_survey = as.integer(0:(n_p - 1)),
    X_lambda = matrix(0, nrow = n_p, ncol = 0)
  )

  params_tmb <- list(
    alpha       = 0.5,
    beta_lambda = numeric(0),
    Z_d         = rep(-3, n_p),
    log_sigma_Z = log(1.5)
  )

  obj <- MakeADFun(
    data       = data_tmb,
    parameters = params_tmb,
    random     = "Z_d",
    DLL        = "full_model",
    silent     = TRUE
  )

  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 2000))

  rep <- sdreport(obj)

  list(obj = obj, opt = opt, rep = rep)
}


# ---- PART 3: District-level prediction (Stage 2) ----

predict_districts <- function(country_data) {
  dist <- country_data$district_rhis %>%
    filter(!is.na(lambda_p)) %>%
    mutate(
      # Calibrated undiagnosed prevalence
      calibrated_logit = z_d - lambda_p,
      calibrated_p     = plogis(calibrated_logit),

      # Uncertainty: district RHIS noise + province calibration noise
      var_calibrated = v_d + var_lambda_p,
      se_calibrated  = sqrt(var_calibrated),

      # Approximate CI on probability scale
      calibrated_lower = plogis(calibrated_logit - 1.96 * se_calibrated),
      calibrated_upper = plogis(calibrated_logit + 1.96 * se_calibrated),

      # CI width as audit priority metric
      ci_width = calibrated_upper - calibrated_lower,

      # Variance decomposition: RHIS volume vs calibration
      pct_var_rhis  = 100 * v_d / var_calibrated,
      pct_var_calib = 100 * var_lambda_p / var_calibrated
    )

  # Flag: admin2 = admin1 (Lesotho)
  n_admin2 <- n_distinct(dist$district)
  n_admin1 <- n_distinct(dist$province)
  if (n_admin2 == n_admin1)
    cat("  NOTE: admin2 = admin1 for", country_data$country,
        "- district predictions do not add sub-province resolution.\n")

  dist
}


# ---- PART 4: Main execution ----

countries <- c("Lesotho", "Zimbabwe", "Malawi", "Mozambique")

compile_tmb()

results <- list()

for (cty in countries) {
  cat("\n", strrep("=", 60), "\n")
  cat(" ", cty, "\n")
  cat(strrep("=", 60), "\n")

  # Assemble province + district data
  cty_data <- tryCatch(
    assemble_country(cty),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
  )
  if (is.null(cty_data)) next

  # ── Stage 1a: TMB fit ──
  cat("\n--- TMB province-level fit ---\n")
  tmb_res <- tryCatch(
    fit_tmb_province(cty_data$province_cal),
    error = function(e) { cat("  TMB ERROR:", e$message, "\n"); NULL }
  )

  alpha_tmb <- NA_real_; se_tmb <- NA_real_; sigma_Z <- NA_real_
  if (!is.null(tmb_res)) {
    tmb_fixed <- summary(tmb_res$rep, "fixed")
    alpha_tmb <- tmb_fixed["alpha", "Estimate"]
    se_tmb    <- tmb_fixed["alpha", "Std. Error"]
    sigma_Z   <- exp(tmb_fixed["log_sigma_Z", "Estimate"])

    cat(sprintf("  alpha_TMB = %.4f (SE = %.4f), OR = %.3f\n",
                alpha_tmb, se_tmb, exp(alpha_tmb)))
    cat(sprintf("  sigma_Z = %.2f\n", sigma_Z))
    cat("  Converged:", tmb_res$opt$convergence == 0, "\n")
  }

  # ── Stage 1b: Province-specific lambda (raw gaps) ──
  prov_cal <- cty_data$province_cal
  wt_alpha <- weighted.mean(prov_cal$lambda_p, 1 / prov_cal$var_lambda_p)
  cat(sprintf("\n--- Province-specific lambda ---\n"))
  cat(sprintf("  IV-weighted mean lambda: %.3f (OR = %.2f)\n",
              wt_alpha, exp(wt_alpha)))
  cat(sprintf("  Province OR range: %.2f - %.2f\n",
              min(prov_cal$OR_p), max(prov_cal$OR_p)))

  # ── Stage 2: District predictions ──
  cat("\n--- District-level predictions ---\n")
  dist_pred <- predict_districts(cty_data)
  cat(sprintf("  Districts predicted: %d\n", nrow(dist_pred)))

  if (nrow(dist_pred) > 0) {
    cat(sprintf("  Calibrated p: median %.1f%%, range %.1f%%-%.1f%%\n",
                100 * median(dist_pred$calibrated_p, na.rm = TRUE),
                100 * min(dist_pred$calibrated_p, na.rm = TRUE),
                100 * max(dist_pred$calibrated_p, na.rm = TRUE)))
    cat(sprintf("  Median CI width: %.1f pp\n",
                100 * median(dist_pred$ci_width, na.rm = TRUE)))
    cat(sprintf("  Variance from RHIS: median %.0f%%, from calibration: %.0f%%\n",
                median(dist_pred$pct_var_rhis, na.rm = TRUE),
                median(dist_pred$pct_var_calib, na.rm = TRUE)))

    # Top 10 audit priorities (widest CI)
    cat("\n--- Top 10 audit-priority districts ---\n")
    top10 <- dist_pred %>%
      arrange(desc(ci_width)) %>%
      slice_head(n = 10) %>%
      transmute(
        Province   = province,
        District   = district,
        Tests      = format(n_tests, big.mark = ","),
        `RHIS pos` = sprintf("%.1f%%", 100 * n_positive / n_tests),
        `Calib p`  = sprintf("%.1f%%", 100 * calibrated_p),
        `95% CI`   = sprintf("[%.1f, %.1f]%%",
                              100 * calibrated_lower,
                              100 * calibrated_upper),
        `CI width`  = sprintf("%.1f pp", 100 * ci_width),
        `%var RHIS` = sprintf("%.0f%%", pct_var_rhis)
      )
    print(as.data.frame(top10), row.names = FALSE)
  }

  # ── Save ──
  out_cty <- file.path(proj, "output", tolower(cty))
  dir.create(out_cty, showWarnings = FALSE, recursive = TRUE)

  saveRDS(prov_cal,  file.path(out_cty, "province_calibration.rds"))
  saveRDS(dist_pred, file.path(out_cty, "district_predictions.rds"))
  if (!is.null(tmb_res))
    saveRDS(tmb_res, file.path(out_cty, "tmb_province_fit.rds"))

  results[[cty]] <- list(
    country        = cty,
    n_provinces    = nrow(prov_cal),
    n_districts    = nrow(dist_pred),
    alpha_shortcut = wt_alpha,
    or_shortcut    = exp(wt_alpha),
    alpha_tmb      = alpha_tmb,
    se_tmb         = se_tmb,
    or_tmb         = exp(alpha_tmb),
    sigma_Z        = sigma_Z,
    or_range_prov  = sprintf("%.2f-%.2f", min(prov_cal$OR_p), max(prov_cal$OR_p)),
    province_cal   = prov_cal,
    district_pred  = dist_pred
  )
}


# ---- PART 5: Cross-country summary tables ----

cat("\n\n")
cat(strrep("=", 75), "\n")
cat("CROSS-COUNTRY SUMMARY\n")
cat(strrep("=", 75), "\n\n")

# ── Table 1: National alpha ──
cat("Table 1: National-level selection bias\n\n")
national_tbl <- bind_rows(lapply(results, function(r) {
  data.frame(
    Country   = r$country,
    Provinces = r$n_provinces,
    Districts = r$n_districts,
    alpha_IV  = round(r$alpha_shortcut, 3),
    OR_IV     = round(r$or_shortcut, 2),
    alpha_TMB = round(r$alpha_tmb, 3),
    OR_TMB    = round(r$or_tmb, 2),
    sigma_Z   = round(r$sigma_Z, 2),
    OR_range  = r$or_range_prov,
    stringsAsFactors = FALSE
  )
}))
print(national_tbl, row.names = FALSE)

# ── Table 2: All provinces ──
cat("\nTable 2: Province-level selection bias (all countries)\n\n")
all_provinces <- bind_rows(lapply(results, function(r) {
  r$province_cal %>%
    transmute(
      Country  = r$country,
      Province = province,
      Tests    = format(n_tests, big.mark = ","),
      `RHIS z` = round(z_p, 3),
      `PHIA Z` = round(Z_hat_p, 3),
      psi      = round(psi_p, 4),
      lambda   = round(lambda_p, 3),
      OR       = round(OR_p, 2),
      SE       = round(se_lambda_p, 3)
    )
}))
print(as.data.frame(all_provinces), row.names = FALSE)

# ── Table 3: District summary ──
cat("\nTable 3: District-level prediction summary\n\n")
dist_summary <- bind_rows(lapply(results, function(r) {
  d <- r$district_pred
  if (nrow(d) == 0) return(NULL)
  data.frame(
    Country       = r$country,
    n_districts   = nrow(d),
    median_p      = sprintf("%.1f%%", 100 * median(d$calibrated_p, na.rm = TRUE)),
    range_p       = sprintf("%.1f-%.1f%%",
                             100 * min(d$calibrated_p, na.rm = TRUE),
                             100 * max(d$calibrated_p, na.rm = TRUE)),
    median_CI     = sprintf("%.1f pp", 100 * median(d$ci_width, na.rm = TRUE)),
    pct_var_RHIS  = sprintf("%.0f%%", median(d$pct_var_rhis, na.rm = TRUE)),
    pct_var_calib = sprintf("%.0f%%", median(d$pct_var_calib, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}))
print(dist_summary, row.names = FALSE)

# ── Save ──
saveRDS(results, file.path(proj, "output", "two_stage_results.rds"))

all_districts <- bind_rows(lapply(results, function(r) {
  r$district_pred %>% mutate(country = r$country)
}))
saveRDS(all_districts, file.path(proj, "output", "two_stage_all_districts.rds"))

cat("\nAll outputs saved.\n")
cat("=== Done ===\n")

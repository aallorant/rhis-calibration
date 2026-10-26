## RUN_canonical.R: canonical fit pipeline. For each country separately, at two
## resolutions (district-marginal and district x age-sex), fits the 2x2 grid
## estimator {TMB, EBLUP} x covariates {none, district-level covariates}.
## Reads:  output/<country>/fits/ (district predictions, conditional survey anchors,
##         covariate matrices) and the TMB templates in tmb/.
## Writes: output/canonical/estimates/, output/canonical/figures/, output/canonical/MANIFEST.csv.
## Run:    Rscript R/RUN_canonical.R [marginal|agesex|figures|all]  (from the repository root)
##
## Estimand: undiagnosed (conditional) HIV positivity on the logit scale.
##   register: z = logit(theta) + lambda + e,  e ~ N(0, v)    (known v)
##   survey:   Zhat = logit(theta) + r,         r ~ N(0, psi)  (known psi)
##   lambda: log relative test propensity, the selection offset
## TMB templates: tmb/province_anchored.cpp (baseline, lambda_p per province or per
##   province x stratum), tmb/mu_free.cpp (marginal +cov, lambda_d = [1, x_d] . beta),
##   tmb/full_model_covariates.cpp (age-sex +cov, district eta + x_d' beta_g).
## EBLUP counterparts: measurement-error Fay-Herriot fits of the same models, by ML.

suppressPackageStartupMessages({library(TMB); library(dplyr); library(tidyr)})
source("R/utils.R")   # harmonise_register_units()
options(stringsAsFactors = FALSE)

PROJ    <- normalizePath(".")
OUT     <- file.path(PROJ, "output", "canonical")
EST_DIR <- file.path(OUT, "estimates")
FIG_DIR <- file.path(OUT, "figures")
for (d in c(OUT, EST_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

STAGE <- {a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "all"}

COUNTRIES <- c("Mozambique", "Zimbabwe", "Malawi", "Lesotho")
short <- function(c) tolower(c)

## per-country covariate spec (self-test + PEPFAR tier)
COV_SPEC <- list(
  Mozambique = c("pepfar_priority_scaleup", "pepfar_priority_notset", "pepfar_selftest_log_raw"),
  Zimbabwe   = c("pepfar_priority_scaleup", "pepfar_priority_notset", "selftest_log_norm", "pepfar_covered"),
  Malawi     = c("pepfar_priority_scaleup", "pepfar_priority_notset", "selftest_log_norm"),
  Lesotho    = c("pepfar_priority_scaleup", "pepfar_priority_notset", "selftest_binary")
)
## age-sex stratum layout per country: base_vars = stratum of the baseline fits
## (age x sex where the RHIS carries it; Malawi sex-only); var/map = collapsed
## stratum index for the covariate model (tmb/full_model_covariates.cpp).
STRATUM_CFG <- list(
  Mozambique = list(base_vars = c("age_group_harmonised","sex"),
                    var = "age_group_harmonised", map = c("15-24"=0L,"25-49"=1L,"50+"=2L)),
  Zimbabwe   = list(base_vars = c("age_group_harmonised","sex"),
                    var = "age_group_harmonised", map = c("15-24"=0L,"25-49"=1L,"50+"=2L)),
  Malawi     = list(base_vars = "sex",
                    var = "sex",                  map = c("Female"=0L,"Male"=1L)),
  Lesotho    = list(base_vars = c("age_group_harmonised","sex"),
                    var = "age_group_harmonised", map = c("15-24"=0L,"25-49"=1L,"50+"=2L))
)

## ---- compile / load TMB templates ----------------------------------------
load_dll <- function(cpp) {
  so <- sub("\\.cpp$", ".so", cpp)
  if (!file.exists(so) || file.info(so)$mtime < file.info(cpp)$mtime) TMB::compile(cpp)
  dyn.load(TMB::dynlib(sub("\\.cpp$", "", cpp)))
}
load_dll(file.path(PROJ, "tmb", "province_anchored.cpp"))
load_dll(file.path(PROJ, "tmb", "mu_free.cpp"))
load_dll(file.path(PROJ, "tmb", "full_model_covariates.cpp"))

## inverse-variance weighted mean (national grand mean of a per-unit offset)
ivw <- function(est, se) {w <- 1/se^2; c(m = sum(est*w)/sum(w), se = sqrt(1/sum(w)))}

## ML profile of sigma_u^2 in the differenced Fay-Herriot measurement-error model
##   y_i = mu_{group(i)} + u_i + e_i - r_{group(i)},  u_i ~ N(0, s2), e_i ~ N(0, v_i), r_g ~ N(0, psi_g),
##   block covariance V_g = s2 I + diag(v_i) + psi_g J, handled by Sherman-Morrison.
fh_profile <- function(y, v, psi, grp) {
  grp <- factor(grp); gl <- levels(grp)
  pg  <- tapply(psi, grp, function(z) z[1])[gl]
  nll <- function(ls2) {
    s2 <- exp(ls2); A <- s2 + v
    S  <- tapply(1/A, grp, sum)[gl]; T <- tapply(y/A, grp, sum)[gl]
    m  <- setNames(as.numeric(T/S), gl); r <- y - m[as.character(grp)]
    0.5 * ( sum(log(2*pi*A)) + sum(log(1 + pg*S)) + sum(r^2/A) )
  }
  o  <- optimize(nll, c(-12, 4)); s2 <- exp(o$minimum); A <- s2 + v
  S  <- tapply(1/A, grp, sum)[gl]; T <- tapply(y/A, grp, sum)[gl]
  list(sigma_u2 = s2, mean = T/S, se_mean = sqrt(1/S + pg), nll = o$objective)
}

## ---- STAGE 1: DISTRICT-MARGINAL fits (no age-sex) ----

fit_marginal <- function(cty) {
  dp <- harmonise_register_units(as.data.frame(readRDS(file.path("output", short(cty), "fits", "district_predictions.rds"))), short(cty)) |>
    filter(is.finite(z_d), is.finite(v_d), v_d > 0)
  pc <- as.data.frame(readRDS(file.path("output", short(cty), "fits", "province_calibration_conditional.rds")))
  provs <- sort(unique(dp$province)); pc <- pc |> filter(province %in% provs) |> arrange(match(province, provs))
  dp$pi <- match(dp$province, provs) - 1L
  npr <- length(provs)

  ## --- TMB baseline: province-anchored ---
  dat <- list(z = dp$z_d, v = dp$v_d, prov_d = as.integer(dp$pi), Zhat = pc$Z_hat_p, psi = pc$psi_p)
  par <- list(lambda = rep(0.5, npr), eta_p = pc$Z_hat_p, u = rep(0, nrow(dp)), log_sigma_u = log(0.3))
  obj <- MakeADFun(dat, par, random = "u", DLL = "province_anchored", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 3000, eval.max = 5000))
  rep <- sdreport(obj); ss <- summary(rep, "report"); lam <- ss[rownames(ss) == "lambda", , drop = FALSE]
  tmb_lam <- data.frame(province = provs, log_lambda = lam[, "Estimate"], se = lam[, "Std. Error"])
  tn <- ivw(tmb_lam$log_lambda, tmb_lam$se)
  sigma_u_tmb <- obj$report()$sigma_u

  ## --- EBLUP baseline: differenced FH-ME ---
  y  <- dp$z_d - pc$Z_hat_p[dp$pi + 1]
  fh <- fh_profile(y, dp$v_d, pc$psi_p[dp$pi + 1], factor(dp$province, levels = provs))
  eb_lam <- data.frame(province = provs, log_lambda = as.numeric(fh$mean), se = as.numeric(fh$se_mean))
  en <- ivw(eb_lam$log_lambda, eb_lam$se)

  ## district predictors: (i) reporting equation theta_hat = expit(z - lambda_hat_p),
  ##   Var(logit) = v_d + se(lambda_hat_p)^2; (ii) shrinkage EBLUP predictor
  se_lam_d <- eb_lam$se[dp$pi + 1]
  eta_rep  <- dp$z_d - eb_lam$log_lambda[dp$pi + 1]
  var_rep  <- dp$v_d + se_lam_d^2
  share_calib <- se_lam_d^2 / var_rep
  ci_lo <- plogis(eta_rep - 1.96*sqrt(var_rep)); ci_hi <- plogis(eta_rep + 1.96*sqrt(var_rep))
  su2 <- fh$sigma_u2
  ## shrinkage factor gamma = sigma_u^2 / (sigma_u^2 + v_d), no psi_p; diagnostic column
  gamma  <- su2 / (su2 + dp$v_d)
  eta_eb  <- gamma * (dp$z_d - eb_lam$log_lambda[dp$pi + 1]) + (1 - gamma) * pc$Z_hat_p[dp$pi + 1]

  ## province-aggregate two-stage lambda_p from the province calibration file
  or_prov_twostage <- if (all(c("lambda_p","se_lambda_p") %in% names(pc))) {
    tp <- ivw(pc$lambda_p, pc$se_lambda_p); exp(tp["m"])
  } else NA_real_

  ## --- +covariates: mu_free (TMB) and FH regression (EBLUP) ---
  cov_out <- tryCatch({
    cov <- readRDS(file.path("output", short(cty), "fits",
                             if (cty == "Mozambique") "moz_covariate_matrix.rds" else "country_covariate_matrix.rds"))
    ctab <- cov$table
    key_d <- if ("canonical_district" %in% names(ctab)) toupper(ctab$canonical_district) else toupper(ctab$district)
    ctab$mk <- key_d
    dp$mk <- toupper(dp$district)
    specs <- intersect(COV_SPEC[[cty]], names(ctab))
    Xc <- ctab[match(dp$mk, ctab$mk), specs, drop = FALSE]
    if (length(specs) == 0 || any(is.na(Xc[[1]])) || nrow(na.omit(Xc)) != nrow(dp))
      stop("covariates unavailable at this resolution")

    X <- as.matrix(cbind(1, Xc)); storage.mode(X) <- "double"
    w <- ave(rep(1, nrow(dp)), dp$province, FUN = function(z) z/length(z))
    datc <- list(z = dp$z_d, v = dp$v_d, prov_d = as.integer(dp$pi), Xlam = X,
                 Zhat = pc$Z_hat_p, psi = pc$psi_p, w = w, agg = 0L)
    parc <- list(eta_p = pc$Z_hat_p, beta = rep(0, ncol(X)), log_sigma_u = log(0.3), u = rep(0, nrow(dp)))
    objc <- MakeADFun(datc, parc, random = "u", DLL = "mu_free", silent = TRUE)
    optc <- nlminb(objc$par, objc$fn, objc$gr, control = list(iter.max = 3000, eval.max = 5000))
    repc <- sdreport(objc); su2c_tmb <- objc$report()$sigma_u^2
    bt <- summary(repc, "report"); bt <- bt[rownames(bt) == "beta", , drop = FALSE]

    ## EBLUP covariate: FH-ME regression y = province-FE + X.beta + u, GLS with ML sigma_u.
    ## QR drops any aliased province dummy (e.g. singleton provinces) so small-N countries fit.
    Dp <- model.matrix(~ 0 + factor(dp$province, levels = provs))
    Xr0 <- cbind(Dp, as.matrix(Xc)); storage.mode(Xr0) <- "double"
    keep <- {qr0 <- qr(Xr0); qr0$pivot[seq_len(qr0$rank)]}
    Xr <- Xr0[, sort(keep), drop = FALSE]
    cov_names <- colnames(Xr0)[sort(keep)]
    fh_reg <- function(ls2) {
      s2 <- exp(ls2); V <- ve + s2; Wv <- 1/V
      b <- solve(t(Xr*Wv) %*% Xr, t(Xr*Wv) %*% y); r <- y - Xr %*% b
      0.5 * sum(log(2*pi*V) + r^2/V)
    }
    oc <- optimize(fh_reg, c(-12, 4)); s2c <- exp(oc$minimum); V <- ve + s2c; Wv <- 1/V
    XtWX <- t(Xr*Wv) %*% Xr; b_eb <- as.numeric(solve(XtWX, t(Xr*Wv) %*% y))
    se_eb <- sqrt(diag(solve(XtWX)))
    names(b_eb) <- names(se_eb) <- cov_names
    sel <- cov_names %in% specs
    list(spec = specs, tmb_beta = bt, tmb_sigma_u2 = su2c_tmb,
         eblup_beta = b_eb[sel], eblup_se = se_eb[sel], eblup_sigma_u2 = s2c)
  }, error = function(e) {cat("   [", cty, "cov skipped:", conditionMessage(e), "]\n"); NULL})

  list(country = cty, n_units = nrow(dp), n_prov = npr, provinces = provs,
       tmb = list(lambda = tmb_lam, OR = exp(tn["m"]), OR_lo = exp(tn["m"]-1.96*tn["se"]),
                  OR_hi = exp(tn["m"]+1.96*tn["se"]), sigma_u2 = sigma_u_tmb^2,
                  converged = opt$convergence == 0, pdHess = rep$pdHess),
       eblup = list(lambda = eb_lam, OR = exp(en["m"]), OR_lo = exp(en["m"]-1.96*en["se"]),
                    OR_hi = exp(en["m"]+1.96*en["se"]), sigma_u2 = fh$sigma_u2),
       OR_prov_twostage = or_prov_twostage, cov = cov_out,
       district = data.frame(country = cty, province = dp$province, district = dp$district,
                             theta_hat = plogis(eta_rep), ci_lo = ci_lo, ci_hi = ci_hi,
                             v_d = dp$v_d, se_lambda = se_lam_d, share_calib = share_calib,
                             theta_hat_shrunk = plogis(eta_eb), z_d = dp$z_d))
}

run_marginal <- function() {
  cat("\n================  STAGE: DISTRICT-MARGINAL  ================\n")
  res <- lapply(COUNTRIES, function(c) tryCatch(fit_marginal(c),
                error = function(e) {cat(c, "ERR:", conditionMessage(e), "\n"); NULL}))
  names(res) <- COUNTRIES
  saveRDS(res, file.path(EST_DIR, "collapsed_admin2_fits.rds"))

  tab <- bind_rows(lapply(res, function(r) if (is.null(r)) NULL else data.frame(
    country = r$country, n_units = r$n_units, n_prov = r$n_prov,
    OR_tmb = round(r$tmb$OR, 3), OR_tmb_lo = round(r$tmb$OR_lo, 3), OR_tmb_hi = round(r$tmb$OR_hi, 3),
    OR_eblup_diff = round(r$eblup$OR, 3),
    OR_eblup_lo = round(r$eblup$OR_lo, 3), OR_eblup_hi = round(r$eblup$OR_hi, 3),
    OR_prov_twostage_legacy = round(r$OR_prov_twostage, 3),
    sigma_u2_tmb = round(r$tmb$sigma_u2, 4), sigma_u2_eblup = round(r$eblup$sigma_u2, 4),
    sigma_u2_cov_tmb = if (is.null(r$cov)) NA else round(r$cov$tmb_sigma_u2, 4),
    sigma_u2_cov_eblup = if (is.null(r$cov)) NA else round(r$cov$eblup_sigma_u2, 4))))
  write.csv(tab, file.path(EST_DIR, "collapsed_admin2_summary.csv"), row.names = FALSE)
  cat("\n--- collapsed_admin2_summary.csv ---\n"); print(tab)

  ## per-province offset, both estimators
  pl <- bind_rows(lapply(res, function(r) {
    if (is.null(r)) return(NULL)
    bind_rows(
      r$tmb$lambda   |> transmute(country = r$country, province, estimator = "TMB",
                                  log_lambda, se, OR = exp(log_lambda),
                                  lo = exp(log_lambda - 1.96*se), hi = exp(log_lambda + 1.96*se)),
      r$eblup$lambda |> transmute(country = r$country, province, estimator = "EBLUP",
                                  log_lambda, se, OR = exp(log_lambda),
                                  lo = exp(log_lambda - 1.96*se), hi = exp(log_lambda + 1.96*se)))
  }))
  pl <- pl |> mutate(across(c(log_lambda, se, OR, lo, hi), ~round(.x, 3)))
  write.csv(pl, file.path(EST_DIR, "collapsed_admin2_province_lambda.csv"), row.names = FALSE)
  cat("\n--- province_lambda.csv (", nrow(pl), "rows ) ---\n")

  ## district predictions (reporting equation) + per-country summary table
  dd <- bind_rows(lapply(res, function(r) if (is.null(r)) NULL else r$district))
  dd <- dd |> mutate(across(c(theta_hat, ci_lo, ci_hi, se_lambda, share_calib, theta_hat_shrunk), ~round(.x, 5)))
  write.csv(dd, file.path(EST_DIR, "district_predictions.csv"), row.names = FALSE)
  dt <- dd |> group_by(country) |> summarise(
    n_areas = n(),
    median_theta_pct = round(100*median(theta_hat), 1),
    min_theta_pct = round(100*min(theta_hat), 1),
    max_theta_pct = round(100*max(theta_hat), 1),
    median_ci_width_pp = round(100*median(ci_hi - ci_lo), 1),
    median_share_calib_pct = round(100*median(share_calib)),
    .groups = "drop") |>
    arrange(match(country, COUNTRIES))
  ## printed only; the district table file is written by R/regen_table_districts.R
  cat("\n--- district_table.csv ---\n"); print(as.data.frame(dt))

  ## covariate beta table (EBLUP FH-regression betas that survived QR; TMB betas alongside)
  cb <- bind_rows(lapply(res, function(r) {
    if (is.null(r) || is.null(r$cov)) return(NULL)
    nm <- names(r$cov$eblup_beta)
    tb <- r$cov$tmb_beta                         # mu_free ADREPORT beta = [intercept, x1, x2, ...]
    data.frame(country = r$country, covariate = nm,
               beta_eblup = round(as.numeric(r$cov$eblup_beta), 4),
               se_eblup   = round(as.numeric(r$cov$eblup_se), 4),
               beta_tmb   = if (!is.null(tb) && nrow(tb) >= length(nm) + 1)
                              round(tb[1 + seq_along(nm), "Estimate"], 4) else NA)
  }))
  if (nrow(cb)) {write.csv(cb, file.path(EST_DIR, "marginal_covariate_beta.csv"), row.names = FALSE)
    cat("\n--- marginal_covariate_beta.csv ---\n"); print(cb)}
  invisible(res)
}

if (STAGE %in% c("marginal", "all")) run_marginal()

## ---- STAGE 2: AGE-SEX fits (district x stratum) ----
## Baseline: province-anchored, unit = province x stratum, aggregated to national-by-stratum
##   (TMB and EBLUP). +covariates: stratum-aware district model (tmb/full_model_covariates.cpp)
##   with stratum-specific effects beta_{k,g} (TMB only).

fit_agesex <- function(cty) {
  cfg <- STRATUM_CFG[[cty]]; svar <- cfg$var
  rh <- harmonise_register_units(as.data.frame(readRDS(file.path("output", short(cty), "fits", "district_predictions_demo.rds"))), short(cty)) |>
    filter(is.finite(z), is.finite(v), v > 0)
  sv <- as.data.frame(readRDS(file.path("output", short(cty), "fits", "province_calibration_demo_conditional.rds"))) |>
    filter(is.finite(Z_hat), is.finite(psi), psi > 0)
  ## BASELINE stratum = full interaction of base_vars (age x sex where carried)
  rh$stratum <- do.call(paste, c(rh[cfg$base_vars], sep = " | "))
  sv$stratum <- do.call(paste, c(sv[cfg$base_vars], sep = " | "))
  rh$combo <- paste(rh$province, rh$stratum); sv$combo <- paste(sv$province, sv$stratum)
  combos <- intersect(unique(rh$combo), unique(sv$combo))
  rh <- rh |> filter(combo %in% combos)
  sv <- sv |> filter(combo %in% combos) |> distinct(combo, .keep_all = TRUE) |> arrange(match(combo, combos))
  rh$idx <- match(rh$combo, combos) - 1L

  ## --- TMB baseline (province x stratum unit) ---
  dat <- list(z = rh$z, v = rh$v, prov_d = as.integer(rh$idx), Zhat = sv$Z_hat, psi = sv$psi)
  par <- list(lambda = rep(0.5, length(combos)), eta_p = sv$Z_hat, u = rep(0, nrow(rh)), log_sigma_u = log(0.3))
  obj <- MakeADFun(dat, par, random = "u", DLL = "province_anchored", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 3000, eval.max = 5000))
  rep <- sdreport(obj); ss <- summary(rep, "report"); lam <- ss[rownames(ss) == "lambda", , drop = FALSE]
  cb <- data.frame(combo = combos, log_lambda = lam[, "Estimate"], se = lam[, "Std. Error"],
                   stratum = sv$stratum[match(combos, sv$combo)])
  tmb_ns <- cb |> group_by(stratum) |>
    summarise(w = sum(1/se^2), log_lambda = sum(log_lambda/se^2)/w, se = sqrt(1/w), .groups = "drop") |>
    transmute(stratum, estimator = "TMB", log_lambda, se, OR = exp(log_lambda))

  ## --- EBLUP baseline (differenced FH per province x stratum) ---
  y  <- rh$z - sv$Z_hat[rh$idx + 1]
  fh <- fh_profile(y, rh$v, sv$psi[rh$idx + 1], factor(rh$combo, levels = combos))
  eb <- data.frame(combo = combos, log_lambda = as.numeric(fh$mean), se = as.numeric(fh$se_mean),
                   stratum = sv$stratum[match(combos, sv$combo)])
  eb_ns <- eb |> group_by(stratum) |>
    summarise(w = sum(1/se^2), log_lambda = sum(log_lambda/se^2)/w, se = sqrt(1/w), .groups = "drop") |>
    transmute(stratum, estimator = "EBLUP", log_lambda, se, OR = exp(log_lambda))

  ## --- TMB +covariates (stratum-aware district model) ---
  cov_tmb <- tryCatch({
    cov <- readRDS(file.path("output", short(cty), "fits",
                             if (cty == "Mozambique") "moz_covariate_matrix.rds" else "country_covariate_matrix.rds"))
    ctab <- cov$table |> arrange(if ("canonical_province" %in% names(cov$table)) canonical_province else province,
                                 if ("canonical_district" %in% names(cov$table)) canonical_district else district)
    provcol <- if ("canonical_province" %in% names(ctab)) "canonical_province" else "province"
    distcol <- if ("canonical_district" %in% names(ctab)) "canonical_district" else "district"
    ctab$district_idx0 <- seq_len(nrow(ctab)) - 1L
    provs <- sort(unique(ctab[[provcol]])); ctab$province_idx0 <- match(ctab[[provcol]], provs) - 1L
    specs <- intersect(COV_SPEC[[cty]], names(ctab)); if (!length(specs)) stop("no covariates")
    rh2 <- rh; rh2$stratum_idx0 <- cfg$map[as.character(rh2[[cfg$var]])]
    ck <- paste(toupper(ctab[[provcol]]), toupper(ctab[[distcol]])); dk <- paste(toupper(rh2$province), toupper(rh2$district))
    rh2$district_idx0 <- match(dk, ck) - 1L; rh2 <- rh2 |> filter(!is.na(district_idx0), !is.na(stratum_idx0))
    sv2 <- sv; sv2$stratum_idx0 <- cfg$map[as.character(sv2[[cfg$var]])]
    sv2$province_idx0 <- match(sv2$province, provs) - 1L
    ph <- ctab |> select(province = all_of(provcol), district_idx0) |>
      inner_join(sv2 |> transmute(province, stratum_idx0, Z_hat, psi), by = "province", relationship = "many-to-many")
    X <- as.matrix(ctab[, specs, drop = FALSE]); storage.mode(X) <- "double"
    ns <- length(cfg$map)
    datc <- list(z = rh2$z, v = rh2$v, d_rhis = as.integer(rh2$district_idx0), g_rhis = as.integer(rh2$stratum_idx0),
                 Z_hat = ph$Z_hat, psi = ph$psi, d_phia = as.integer(ph$district_idx0), g_phia = as.integer(ph$stratum_idx0),
                 prov_of_district = as.integer(ctab$province_idx0), X = X)
    parc <- list(alpha = matrix(0, length(provs), ns), delta = rep(0, ns - 1), beta = matrix(0, ncol(X), ns),
                 eta = rep(-3, nrow(ctab)), u = matrix(0, nrow(ctab), ns), log_sigma_eta = log(1), log_sigma_u = log(0.5))
    objc <- MakeADFun(datc, parc, random = c("eta", "u"), DLL = "full_model_covariates", silent = TRUE)
    optc <- nlminb(objc$par, objc$fn, objc$gr, control = list(iter.max = 1500, eval.max = 3000, rel.tol = 1e-8))
    repc <- sdreport(objc); fx <- summary(repc, "fixed")
    brow <- grepl("^beta", rownames(fx)); bt <- fx[brow, , drop = FALSE]
    list(spec = specs, sigma_u2 = exp(optc$par["log_sigma_u"])^2, n_cov = ncol(X), n_strata = ns,
         beta = if (nrow(bt)) matrix(bt[, "Estimate"], ncol(X), ns) else NULL,
         beta_se = if (nrow(bt)) matrix(bt[, "Std. Error"], ncol(X), ns) else NULL,
         pdHess = repc$pdHess, converged = optc$convergence == 0, strata = names(cfg$map))
  }, error = function(e) {cat("   [", cty, "agesex cov skipped:", conditionMessage(e), "]\n"); NULL})

  list(country = cty, strata = names(cfg$map), national_by_stratum = bind_rows(tmb_ns, eb_ns),
       sigma_u2_tmb = obj$report()$sigma_u^2, sigma_u2_eblup = fh$sigma_u2,
       converged = opt$convergence == 0, cov = cov_tmb)
}

run_agesex <- function() {
  cat("\n================  STAGE: AGE-SEX  ================\n")
  res <- lapply(COUNTRIES, function(c) tryCatch(fit_agesex(c),
                error = function(e) {cat(c, "ERR:", conditionMessage(e), "\n"); NULL}))
  names(res) <- COUNTRIES
  saveRDS(res, file.path(EST_DIR, "mainfit_fits.rds"))

  ns <- bind_rows(lapply(res, function(r) if (is.null(r)) NULL else
    cbind(country = r$country, r$national_by_stratum)))
  parts <- strsplit(ns$stratum, " | ", fixed = TRUE)
  ns$age_group <- vapply(parts, function(p) if (length(p) == 2) p[1] else NA_character_, "")
  ns$sex       <- vapply(parts, function(p) if (length(p) == 2) p[2] else p[1], "")
  ns$OR <- round(ns$OR, 3); ns$lo <- round(exp(ns$log_lambda - 1.96*ns$se), 3)
  ns$hi <- round(exp(ns$log_lambda + 1.96*ns$se), 3)
  write.csv(ns, file.path(EST_DIR, "mainfit_national_by_stratum.csv"), row.names = FALSE)
  cat("\n--- mainfit_national_by_stratum.csv (OR by age x sex cell, both estimators) ---\n"); print(ns)

  ## national aggregate over strata (IVW of cell log-lambda), both estimators
  nat <- ns |> group_by(country, estimator) |>
    summarise(w = sum(1/se^2), log_lambda = sum(log_lambda/se^2)/w, se = sqrt(1/w), .groups = "drop") |>
    mutate(OR = round(exp(log_lambda), 3),
           lo = round(exp(log_lambda - 1.96*se), 3), hi = round(exp(log_lambda + 1.96*se), 3)) |>
    select(-w) |> arrange(match(country, COUNTRIES), estimator)
  write.csv(nat, file.path(EST_DIR, "mainfit_national_aggregate.csv"), row.names = FALSE)
  cat("\n--- mainfit_national_aggregate.csv (national OR, IVW over cells) ---\n"); print(as.data.frame(nat))

  covb <- bind_rows(lapply(res, function(r) {
    if (is.null(r) || is.null(r$cov) || is.null(r$cov$beta)) return(NULL)
    cv <- r$cov; out <- expand.grid(k = seq_len(cv$n_cov), g = seq_len(cv$n_strata))
    data.frame(country = r$country, covariate = cv$spec[out$k], stratum = cv$strata[out$g],
               beta = round(cv$beta[cbind(out$k, out$g)], 4), se = round(cv$beta_se[cbind(out$k, out$g)], 4))
  }))
  if (!is.null(covb) && nrow(covb)) {
    write.csv(covb, file.path(EST_DIR, "UNRELIABLE_agesex_covariate_beta.csv"), row.names = FALSE)
    cat("\n--- UNRELIABLE_agesex_covariate_beta.csv (", nrow(covb), "rows ) ---\n")
  }
  su <- bind_rows(lapply(res, function(r) if (is.null(r)) NULL else data.frame(
    country = r$country, sigma_u2_baseline = round(r$sigma_u2_tmb, 4),
    sigma_u2_cov = if (is.null(r$cov)) NA else round(r$cov$sigma_u2, 4),
    cov_pdHess = if (is.null(r$cov)) NA else r$cov$pdHess)))
  write.csv(su, file.path(EST_DIR, "mainfit_sigma_u.csv"), row.names = FALSE)
  cat("\n--- mainfit_sigma_u.csv ---\n"); print(su)
  invisible(res)
}

if (STAGE %in% c("agesex", "all")) run_agesex()

## ---- FIGURES + MANIFEST ----
make_figures <- function() {
  cat("\n================  FIGURES  ================\n")
  ## Fig 1: national selection OR by country, three estimators
  m <- tryCatch(read.csv(file.path(EST_DIR, "collapsed_admin2_summary.csv")), error = function(e) NULL)
  if (!is.null(m)) {
    m <- m[match(COUNTRIES, m$country), ]
    pdf(file.path(FIG_DIR, "fig_national_or.pdf"), width = 8, height = 5)
    cols <- c(TMB = "#1b7837", `EBLUP (consistent)` = "#2166ac", `two-stage (legacy)` = "#b2182b")
    x <- barplot(rbind(m$OR_tmb, m$OR_eblup_diff, m$OR_prov_twostage_legacy),
                 beside = TRUE, names.arg = m$country, col = cols, border = NA,
                 ylab = "national selection odds ratio", ylim = c(0, max(m$OR_prov_twostage_legacy, na.rm=TRUE)*1.15),
                 main = "National selection OR by estimator (district-marginal)")
    abline(h = 1, lty = 3, col = "grey40")
    legend("topright", legend = names(cols), fill = cols, border = NA, bty = "n", cex = 0.9)
    dev.off(); cat("wrote fig_national_or.pdf\n")
  }
  ## Fig 2: age-sex OR by country x stratum, TMB vs EBLUP
  a <- tryCatch(read.csv(file.path(EST_DIR, "mainfit_national_by_stratum.csv")), error = function(e) NULL)
  if (!is.null(a)) {
    pdf(file.path(FIG_DIR, "fig_agesex_or.pdf"), width = 9, height = 5.5)
    op <- par(mfrow = c(2, 2), mar = c(3, 3.5, 2.5, 0.6), mgp = c(2.1, 0.7, 0))
    for (cty in COUNTRIES) {
      ac <- a[a$country == cty, ]; if (!nrow(ac)) {plot.new(); title(cty); next}
      st <- unique(ac$stratum); xt <- seq_along(st)
      tm <- ac[ac$estimator == "TMB", ]; eb <- ac[ac$estimator == "EBLUP", ]
      tm <- tm[match(st, tm$stratum), ]; eb <- eb[match(st, eb$stratum), ]
      plot(NA, xlim = c(0.5, length(st)+0.5), ylim = c(0.7, max(ac$hi, na.rm=TRUE)*1.05),
           xaxt = "n", xlab = "", ylab = "selection OR", main = cty, log = "y")
      axis(1, xt, st, cex.axis = 0.85); abline(h = 1, lty = 3, col = "grey60")
      arrows(xt-0.08, tm$lo, xt-0.08, tm$hi, angle = 90, code = 3, length = 0.03, col = "#1b7837")
      points(xt-0.08, tm$OR, pch = 19, col = "#1b7837")
      arrows(xt+0.08, eb$lo, xt+0.08, eb$hi, angle = 90, code = 3, length = 0.03, col = "#2166ac")
      points(xt+0.08, eb$OR, pch = 17, col = "#2166ac")
    }
    par(op); mtext("Age-sex selection OR: TMB (green) vs EBLUP (blue)", outer = TRUE, line = -1.2, cex = 0.95)
    dev.off(); cat("wrote fig_agesex_or.pdf\n")
  }
  ## Fig 3: within-area residual sigma_u^2 with vs without covariates (age-sex model)
  s <- tryCatch(read.csv(file.path(EST_DIR, "mainfit_sigma_u.csv")), error = function(e) NULL)
  if (!is.null(s)) {
    s <- s[match(COUNTRIES, s$country), ]
    pdf(file.path(FIG_DIR, "fig_cov_sigmau.pdf"), width = 7.5, height = 4.6)
    x <- barplot(rbind(s$sigma_u2_baseline, s$sigma_u2_cov), beside = TRUE, names.arg = s$country,
                 col = c("#7fbf7b", "#af8dc3"), border = NA, ylab = expression(sigma[u]^2~"(within-area residual)"),
                 main = "Covariates absorb within-area selection heterogeneity (age-sex model)")
    legend("topright", c("baseline", "+ district covariates"), fill = c("#7fbf7b", "#af8dc3"), border = NA, bty = "n")
    dev.off(); cat("wrote fig_cov_sigmau.pdf\n")
  }
}

if (STAGE %in% c("figures", "all")) make_figures()

## Manifest of everything written
write_manifest <- function() {
  files <- list.files(OUT, recursive = TRUE, full.names = FALSE)
  info <- file.info(file.path(OUT, files))
  man <- data.frame(file = files, bytes = info$size, modified = as.character(info$mtime))
  write.csv(man, file.path(OUT, "MANIFEST.csv"), row.names = FALSE)
  cat("\nwrote MANIFEST.csv (", nrow(man), "files )\n")
}
if (STAGE %in% c("all", "figures")) write_manifest()

cat("\n[RUN_canonical] stage '", STAGE, "' done. Outputs in output/canonical/\n", sep = "")

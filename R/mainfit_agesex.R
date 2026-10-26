## =====================================================================
## mainfit_agesex.R -- the main fit: block-covariance EBLUP at area x age x sex, one fit
## per country (admin-2 for Mozambique, Zimbabwe and Malawi; admin-1 for Lesotho).
## Groups are province x stratum with V = sigma_u^2 I + diag(v) + psi J, so lambda_{p,g} = T/S
## and Var(lambda_{p,g}) = 1/S + psi; sigma_u^2 is profiled by ML per country.
## Reads:  output/<country>/fits/{district_predictions_demo,province_calibration_demo_conditional}.rds,
##         output/canonical/estimates/mainfit_fits.rds and mainfit_national_by_stratum.csv (self-checks)
## Writes: output/canonical/estimates/mainfit_{province_stratum_lambda,province_offset,sigma_u_profile}.csv
## Usage:  Rscript R/mainfit_agesex.R   (AGESEX_OUT=<dir> redirects the outputs)
## =====================================================================

suppressMessages({library(dplyr)})
source("R/utils.R")   # harmonise_register_units()

CTS <- c("Mozambique", "Zimbabwe", "Malawi", "Lesotho")
OUT <- Sys.getenv("AGESEX_OUT", unset = "output/canonical/estimates")

ivw <- function(est, se) {w <- 1/se^2; c(m = sum(est*w)/sum(w), se = sqrt(1/sum(w)))}

## ---- block-covariance profile, group = province x stratum ----------------
fh_block <- function(y, v, psi, grp) {
  gs <- unique(grp)
  nll <- function(ls2) {
    s2 <- exp(ls2); tot <- 0
    for (g in gs) {
      ix <- which(grp == g); A <- s2 + v[ix]; S <- sum(1/A); p <- psi[[g]]
      m  <- sum(y[ix]/A) / S
      tot <- tot + 0.5 * (sum(log(A)) + log1p(p * S) + sum((y[ix] - m)^2 / A))
    }
    tot
  }
  ls2 <- optimize(nll, c(-12, 3))$minimum; s2 <- exp(ls2)

  ## profile-likelihood 95% interval for sigma_u^2 (boundary-robust)
  tgt <- nll(ls2) + qchisq(0.95, 1)/2
  fr  <- function(l) nll(l) - tgt
  lo  <- tryCatch(exp(uniroot(fr, c(-20, ls2))$root), error = function(e) NA_real_)
  hi  <- tryCatch(exp(uniroot(fr, c(ls2,  8))$root), error = function(e) NA_real_)

  est <- lapply(gs, function(g) {
    ix <- which(grp == g); A <- s2 + v[ix]; S <- sum(1/A); p <- psi[[g]]
    data.frame(group = g, n = length(ix),
               log_lambda = sum(y[ix]/A)/S, se = sqrt(1/S + p))
  })
  list(sigma_u2 = s2, sigma_u2_lo = lo, sigma_u2_hi = hi, lambda = bind_rows(est))
}

rows_pg <- list(); rows_p <- list(); rows_s <- list(); chk <- list()
ref_fits <- readRDS("output/canonical/estimates/mainfit_fits.rds")

for (cty in CTS) {
  s  <- tolower(cty)
  rh <- harmonise_register_units(as.data.frame(readRDS(file.path("output", s, "fits", "district_predictions_demo.rds"))), s) |>
    filter(is.finite(z), is.finite(v), v > 0) |>
    mutate(stratum = trimws(paste(age_group_harmonised, sex)))
  sv <- as.data.frame(readRDS(file.path("output", s, "fits", "province_calibration_demo_conditional.rds"))) |>
    filter(is.finite(Z_hat), is.finite(psi), psi > 0) |>
    mutate(stratum = trimws(paste(age_group_harmonised, sex)))

  ## Malawi's register is sex-only: collapse the survey over age within sex (inverse-variance)
  if (all(grepl("^all ", unique(rh$stratum)))) {
    sv <- sv |> mutate(sex2 = ifelse(grepl("Female", stratum), "Female", "Male")) |>
      group_by(province, sex2) |>
      summarise(Z_hat = sum(Z_hat/psi)/sum(1/psi), psi = 1/sum(1/psi), .groups = "drop") |>
      mutate(stratum = paste("all", sex2)) |> select(province, stratum, Z_hat, psi)
  }

  d <- rh |> inner_join(select(sv, province, stratum, Z_hat, psi),
                        by = c("province", "stratum")) |>
    mutate(grp = paste(province, stratum, sep = " | "), y = z - Z_hat)

  psi_by <- d |> distinct(grp, .keep_all = TRUE) |> select(grp, psi)
  fit <- fh_block(d$y, d$v, setNames(as.list(psi_by$psi), psi_by$grp), d$grp)

  ## Lesotho: 1 area per province x stratum, so sigma_u^2 is NOT ESTIMABLE.
  estimable <- max(fit$lambda$n) > 1
  rows_s[[cty]] <- data.frame(
    country = cty, estimable = estimable, n_areas = length(unique(d$district)),
    sigma_u2    = if (estimable) fit$sigma_u2    else NA_real_,
    sigma_u2_lo = if (estimable) fit$sigma_u2_lo else NA_real_,
    sigma_u2_hi = if (estimable) fit$sigma_u2_hi else NA_real_)

  lam <- fit$lambda |>
    tidyr::separate(group, c("province", "stratum"), sep = " \\| ", remove = FALSE) |>
    mutate(country = cty, OR = exp(log_lambda),
           lo = exp(log_lambda - 1.96*se), hi = exp(log_lambda + 1.96*se))
  rows_pg[[cty]] <- lam

  ## ---- province offset: direct standardisation over strata (equal weights) ---
  ## strata are independent blocks, so Var(sum_g w_g lambda_g) = sum_g w_g^2 se_g^2;
  ## the IVW columns come first because summarise() overwrites log_lambda/se in sequence
  rows_p[[cty]] <- lam |> group_by(country, province) |>
    summarise(n_strata   = n(),
              log_lambda_ivw = ivw(log_lambda, se)["m"],           # kept for comparison
              se_ivw         = ivw(log_lambda, se)["se"],
              log_lambda = mean(log_lambda),                       # w_g = 1/n
              se         = sqrt(sum((se / n())^2)),
              .groups = "drop") |>
    mutate(OR = exp(log_lambda), lo = exp(log_lambda - 1.96*se), hi = exp(log_lambda + 1.96*se),
           OR_ivw = exp(log_lambda_ivw))

  ## ---- self-checks --------------------------------------------------------
  st <- read.csv("output/canonical/estimates/mainfit_national_by_stratum.csv",
                 stringsAsFactors = FALSE) |>
    filter(country == cty, estimator == "EBLUP") |>
    mutate(stratum = sub("^all ", "", trimws(gsub(" \\| ", " ", stratum))))
  mine <- lam |> mutate(stratum = sub("^all ", "", stratum)) |> group_by(stratum) |>
    summarise(mine = ivw(log_lambda, se)["m"], .groups = "drop")
  cmp <- inner_join(mine, select(st, stratum, canon = log_lambda), by = "stratum")
  chk[[cty]] <- data.frame(country = cty,
    sigma_u2 = round(fit$sigma_u2, 4),
    ref_sigma_u2 = round(ref_fits[[cty]]$sigma_u2_eblup, 4),
    strata_matched = nrow(cmp),
    max_abs_diff = if (nrow(cmp)) max(abs(cmp$mine - cmp$canon)) else NA_real_)
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(bind_rows(rows_pg) |> select(country, province, stratum, n, log_lambda, se, OR, lo, hi),
          file.path(OUT, "mainfit_province_stratum_lambda.csv"), row.names = FALSE)
write.csv(bind_rows(rows_p),  file.path(OUT, "mainfit_province_offset.csv"),  row.names = FALSE)
write.csv(bind_rows(rows_s),  file.path(OUT, "mainfit_sigma_u_profile.csv"),  row.names = FALSE)

cat("\n=== SELF-CHECK 1+2 (sigma_u^2 vs mainfit_fits.rds; lambda aggregation vs canonical) ===\n")
print(bind_rows(chk), row.names = FALSE)
cat("\n=== sigma_u^2 with 95% profile interval (MAIN FIT) ===\n")
print(as.data.frame(bind_rows(rows_s) |> mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
cat("\n=== province offsets, DIRECT-STANDARDISED over strata (equal weights) ===\n")
print(as.data.frame(bind_rows(rows_p) |> group_by(country) |>
  summarise(n_prov = n(), OR_min = round(min(OR),2), OR_max = round(max(OR),2),
            n_excl_1 = sum(lo > 1 | hi < 1), .groups = "drop")), row.names = FALSE)

cat("\n=== how far the standardisation moves each province, against the old IVW ===\n")
print(as.data.frame(bind_rows(rows_p) |>
  mutate(pct = 100*(OR/OR_ivw - 1)) |> group_by(country) |>
  summarise(median_pct = round(median(pct),1), max_abs_pct = round(max(abs(pct)),1),
            .groups = "drop")), row.names = FALSE)
cat("\nmoved most:\n")
print(as.data.frame(bind_rows(rows_p) |> mutate(pct = 100*(OR/OR_ivw - 1)) |>
  arrange(desc(abs(pct))) |> head(5) |>
  transmute(country, province, n_strata, OR_std = round(OR,2),
            OR_ivw = round(OR_ivw,2), pct = round(pct,1))), row.names = FALSE)
cat("\nwrote mainfit_province_stratum_lambda.csv, mainfit_province_offset.csv, mainfit_sigma_u_profile.csv to", OUT, "\n")

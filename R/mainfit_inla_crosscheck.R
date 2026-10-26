## =====================================================================
## mainfit_inla_crosscheck.R -- INLA cross-check of the main fit (block-covariance EBLUP
## at area x age x sex) as a fully Bayesian model with the same block structure:
##   y_{d,g} = lambda_{p,g} + b_{p,g} + u_{d,g} + e_{d,g},  lambda_{p,g} fixed (vague prior),
##   b_{p,g} ~ N(0, psi_{p,g}) known, u_{d,g} ~ N(0, sigma_u^2), e_{d,g} ~ N(0, v_{d,g}) known,
## so Var(lambda_{p,g}) = psi_{p,g} + 1/S_{p,g}. Lesotho (one area per cell) does not identify sigma_u^2.
## Reads:  output/<country>/fits/*.rds, output/canonical/estimates/mainfit_{sigma_u_profile,national_aggregate}.csv
## Writes: output/canonical/estimates/mainfit_inla_crosscheck.csv, mainfit_inla_lambda.csv
## Usage:  Rscript R/mainfit_inla_crosscheck.R
## =====================================================================

suppressMessages({library(INLA); library(dplyr); library(Matrix)})
source("R/utils.R")   # harmonise_register_units()

CTS <- c("Mozambique", "Zimbabwe", "Malawi", "Lesotho")
EST <- "output/canonical/estimates"

ivw <- function(est, se) {w <- 1/se^2; c(m = sum(est*w)/sum(w), se = sqrt(1/sum(w)))}

## ---- data: same construction as R/mainfit_agesex.R, incl. the Malawi sex-only collapse ----
load_mainfit <- function(cty) {
  s  <- tolower(cty)
  rh <- harmonise_register_units(as.data.frame(readRDS(file.path("output", s, "fits", "district_predictions_demo.rds"))), s) |>
    filter(is.finite(z), is.finite(v), v > 0) |>
    mutate(stratum = trimws(paste(age_group_harmonised, sex)))
  sv <- as.data.frame(readRDS(file.path("output", s, "fits", "province_calibration_demo_conditional.rds"))) |>
    filter(is.finite(Z_hat), is.finite(psi), psi > 0) |>
    mutate(stratum = trimws(paste(age_group_harmonised, sex)))

  if (all(grepl("^all ", unique(rh$stratum)))) {
    sv <- sv |> mutate(sex2 = ifelse(grepl("Female", stratum), "Female", "Male")) |>
      group_by(province, sex2) |>
      summarise(Z_hat = sum(Z_hat/psi)/sum(1/psi), psi = 1/sum(1/psi), .groups = "drop") |>
      mutate(stratum = paste("all", sex2)) |> select(province, stratum, Z_hat, psi)
  }

  d <- rh |> inner_join(select(sv, province, stratum, Z_hat, psi),
                        by = c("province", "stratum")) |>
    mutate(grp = paste(province, stratum, sep = " | "), y = z - Z_hat)
  d$grp_f   <- factor(d$grp)
  d$grp_idx <- as.integer(d$grp_f)
  d$unit_idx <- seq_len(nrow(d))
  d
}

## --------------------------------------------------------------- model ----
fit_inla_mainfit <- function(d) {
  psi_g <- d |> group_by(grp_idx) |> summarise(psi = first(psi), .groups = "drop") |>
    arrange(grp_idx) |> pull(psi)
  Cpsi <- Diagonal(x = 1/psi_g)

  ff <- y ~ 0 + grp_f +
    f(grp_idx, model = "generic0", Cmatrix = Cpsi,
      hyper = list(prec = list(initial = 0, fixed = TRUE)), constr = FALSE) +
    f(unit_idx, model = "iid",
      hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
  environment(ff) <- environment()

  inla(ff, family = "gaussian", data = d, scale = 1/d$v,
       control.family    = list(hyper = list(prec = list(initial = 0, fixed = TRUE))),
       control.compute   = list(dic = TRUE, waic = TRUE),
       control.predictor = list(compute = FALSE),
       control.inla      = list(strategy = "simplified.laplace", int.strategy = "ccd"))
}

sigma_u2_post <- function(res) {
  nm <- grep("^Precision for unit_idx", rownames(res$summary.hyperpar), value = TRUE)[1]
  if (is.na(nm)) return(c(mean = NA, q025 = NA, q975 = NA))
  m <- INLA::inla.tmarginal(function(x) 1/x, res$marginals.hyperpar[[nm]])
  z <- INLA::inla.zmarginal(m, silent = TRUE)
  c(mean = z$mean, q025 = z$quant0.025, q975 = z$quant0.975)
}

## Sherman-Morrison target from the MAIN FIT's sigma_u^2:
##   Var(lambda_{p,g}) = psi_{p,g} + 1/S_{p,g}
validate_block <- function(d, lam, su2) {
  tgt <- d |> group_by(grp) |>
    summarise(block = sqrt(first(psi) + 1/sum(1/(su2 + v))),
              indep = sqrt(1/sum(1/(su2 + v + first(psi)))), .groups = "drop")
  m <- inner_join(lam, tgt, by = "grp")
  c(block = max(abs(m$sd/m$block - 1)), indep = max(abs(m$sd/m$indep - 1)))
}

## ---------------------------------------------------------------- main ----
su_main  <- read.csv(file.path(EST, "mainfit_sigma_u_profile.csv"), stringsAsFactors = FALSE)
nat_main <- read.csv(file.path(EST, "mainfit_national_aggregate.csv"), stringsAsFactors = FALSE)

rows <- list(); lamrows <- list()

for (cty in CTS) {
  d   <- load_mainfit(cty)
  su2 <- su_main$sigma_u2[su_main$country == cty]
  ident <- su_main$estimable[su_main$country == cty]

  res <- fit_inla_mainfit(d)

  fx  <- res$summary.fixed
  lam <- data.frame(grp = sub("^grp_f", "", rownames(fx)),
                    mean = fx$mean, sd = fx$sd, row.names = NULL) |>
    mutate(country = cty)
  agg <- ivw(lam$mean, lam$sd)
  s2  <- sigma_u2_post(res)
  vb  <- if (ident) validate_block(d, lam, su2) else c(block = NA, indep = NA)

  eb <- nat_main |> filter(country == cty, estimator == "EBLUP")
  tm <- nat_main |> filter(country == cty, estimator == "TMB")

  rows[[cty]] <- data.frame(
    country = cty, n_units = nrow(d), n_groups = nlevels(d$grp_f),
    OR_inla  = exp(agg["m"]),
    OR_lo    = exp(agg["m"] - 1.96*agg["se"]),
    OR_hi    = exp(agg["m"] + 1.96*agg["se"]),
    OR_eblup = eb$OR, OR_tmb = tm$OR,
    se_inla  = agg["se"], se_eblup = eb$se, se_tmb = tm$se,
    sigma_u2_identified = ident,
    sigma_u2_inla = if (ident) s2["mean"] else NA_real_,
    sigma_u2_lo   = if (ident) s2["q025"] else NA_real_,
    sigma_u2_hi   = if (ident) s2["q975"] else NA_real_,
    sigma_u2_main = if (ident) su2 else NA_real_,
    max_rel_err_block = vb["block"], max_rel_err_indep = vb["indep"],
    dic = res$dic$dic, row.names = NULL)
  lamrows[[cty]] <- lam
  message(sprintf("  %-11s done (%d units, %d province x stratum groups)",
                  cty, nrow(d), nlevels(d$grp_f)))
}

out <- bind_rows(rows)
write.csv(out, file.path(EST, "mainfit_inla_crosscheck.csv"), row.names = FALSE)
write.csv(bind_rows(lamrows), file.path(EST, "mainfit_inla_lambda.csv"), row.names = FALSE)

cat("\n=== national selection OR: three routes, ONE model (area x age x sex) ===\n")
print(as.data.frame(out |> transmute(country, n_units, n_groups,
  EBLUP = round(OR_eblup,3), TMB = round(OR_tmb,3), INLA = round(OR_inla,3),
  `INLA 95%` = sprintf("(%.2f, %.2f)", OR_lo, OR_hi),
  se_ratio_inla_eblup = round(se_inla/se_eblup, 3))), row.names = FALSE)

cat("\n=== sigma_u^2 ===\n")
print(as.data.frame(out |> transmute(country, identified = sigma_u2_identified,
  main_fit_ML = round(sigma_u2_main,4), INLA_post_mean = round(sigma_u2_inla,4),
  `INLA 95%` = ifelse(is.na(sigma_u2_lo), "-", sprintf("(%.3f, %.3f)", sigma_u2_lo, sigma_u2_hi)))),
  row.names = FALSE)

cat("\n=== does the INLA posterior sd match the block or the independent form? ===\n")
print(as.data.frame(out |> filter(sigma_u2_identified) |> transmute(country,
  vs_block = sprintf("%.2f%%", 100*max_rel_err_block),
  vs_indep = sprintf("%.2f%%", 100*max_rel_err_indep))), row.names = FALSE)

cat("\nwrote", file.path(EST, "mainfit_inla_crosscheck.csv"), "\n")

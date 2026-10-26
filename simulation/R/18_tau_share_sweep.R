## ============================================================================
## 18_tau_share_sweep.R -- the resolution ladder against the transport-failure
##   share of sigma_u^2, with sigma_u^2 held at TOT = 0.20: the offset's district
##   sd is set to give the target tau variance and the truth's district sds are
##   rescaled by a common factor so the total stays at TOT.
## Reads : simulation/R/00_config.R to 03_designs.R (sourced)
## Writes: tau_share_sweep.csv under simulation/output/
## Usage : Rscript simulation/R/18_tau_share_sweep.R [REPS]   (default 150)
## ============================================================================
suppressWarnings(suppressMessages({ library(data.table) }))
here <- "simulation"
for (f in c("00_config","01_dgp","02_estimators","03_designs")) source(file.path(here,"R",paste0(f,".R")))
REPS <- { a <- commandArgs(TRUE); if (length(a)) as.integer(a[1]) else 150 }
K_AUDIT <- 8; SHARES <- c(0.20, 0.34, 0.44, 0.64, 0.80); TOT <- 0.20
comp <- function(kw, lds, reps = 60) rowMeans(vapply(seq_len(reps), function(r) {
  c2 <- cfg; c2$dist_sd <- cfg$dist_sd * kw; c2$aware_dist_sd <- cfg$aware_dist_sd * kw
  p <- simulate_population(c2, "constant", seed = 9000 + r, lambda_dist_sd = lds)
  key <- paste(p$idx$prov, p$idx$s)
  w <- p$theta[, 1] - ave(p$theta[, 1], key); tau <- p$lambda[, 1] - ave(p$lambda[, 1], key)
  c(var(w + tau), var(tau))
}, numeric(2)))
rmse_at <- function(th, p, g) { agg <- function(pr) { den <- rowsum(p$pop, g)
  sapply(seq_len(ncol(pr)), function(t) rowsum(p$pop * pr[, t], g)[, 1] / den[, 1]) }
  sqrt(mean((agg(plogis(th)) - agg(plogis(p$theta)))^2)) * 100 }
null_surv <- function(p) list(Z = matrix(NA, p$U, p$T), psi = matrix(NA, p$U, p$T))
rows <- list()
for (sh in SHARES) {
  ## tau variance scales with lds^2: one calibration point, then solve kw for the total
  base <- comp(1, cfg$lambda_dist_sd)[2]
  lds <- cfg$lambda_dist_sd * sqrt(sh * TOT / base)
  kw <- uniroot(function(k) comp(k, lds)[1] - TOT, c(0.05, 4), tol = 1e-3)$root
  ach <- comp(kw, lds, reps = 200)
  c2 <- cfg; c2$dist_sd <- cfg$dist_sd * kw; c2$aware_dist_sd <- cfg$aware_dist_sd * kw
  cat(sprintf("share %.2f: lambda sd %.3f, truth factor %.3f -> sigma_u^2 %.3f, share %.2f\n",
              sh, lds, kw, ach[1], ach[2] / ach[1]))
  for (scn in c("constant", "drift")) for (rep in seq_len(REPS)) {
    p <- simulate_population(c2, scn, seed = 5000 + rep, lambda_dist_sd = lds)
    dr <- run_design(p, c2, "random", k = K_AUDIT, seed = rep)
    es <- list(register = p$z, survey = est_survey_only(p, dr$surv, c2)$theta,
               calibrated = est_static_calib(p, dr$surv, c2)$theta)
    lv <- list(national = rep(1, p$U), `admin-1` = p$idx$prov, district = p$idx$d, `district x age-sex` = seq_len(p$U))
    for (m in names(es)) for (L in names(lv))
      rows[[length(rows) + 1]] <- data.table(share_target = sh, share = ach[2] / ach[1], sigma_u2 = ach[1],
                                             scenario = scn, method = m, level = L, rmse = rmse_at(es[[m]], p, lv[[L]]))
  }
}
res <- rbindlist(rows)[, .(rmse = mean(rmse)), by = .(share_target, share, sigma_u2, scenario, method, level)]
w <- dcast(res, share_target + share + sigma_u2 + scenario + level ~ method, value.var = "rmse")
w[, margin_over_survey := survey - calibrated]
fwrite(w, file.path(OUT, "tau_share_sweep.csv"))
print(w[, .(share_target, scenario, level, register = round(register, 2), survey = round(survey, 2),
            calibrated = round(calibrated, 2), margin = round(margin_over_survey, 2))], nrows = 200)

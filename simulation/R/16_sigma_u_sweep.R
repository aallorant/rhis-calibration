## ============================================================================
## 16_sigma_u_sweep.R -- the resolution ladder across a grid of sigma_u^2, the
##   within-province variance of the district residual (truth + offset). For
##   each target the three district sds (dist_sd, aware_dist_sd, lambda_dist_sd)
##   are scaled by one common factor, so the offset share of sigma_u^2 is held
##   fixed while its level moves; the achieved sigma_u^2 is reported.
## Reads : simulation/R/00_config.R to 03_designs.R (sourced)
## Writes: sigma_u_sweep.csv, sigma_u_sweep_wide.csv, sigma_u_sweep_by_scenario.csv
##         and sigma_u_sweep_by_scenario_wide.csv under simulation/output/
## Usage : Rscript simulation/R/16_sigma_u_sweep.R [REPS]   (default 150)
## ============================================================================
suppressWarnings(suppressMessages({ library(data.table) }))
here <- "simulation"
for (f in c("00_config","01_dgp","02_estimators","03_designs")) source(file.path(here,"R",paste0(f,".R")))
REPS    <- { a <- commandArgs(TRUE); if (length(a)) as.integer(a[1]) else 150 }
K_AUDIT <- 8
TARGETS <- c(0.05, 0.10, 0.20, 0.30, 0.50)
COUNTRY <- c(Lesotho = 0.108, Malawi = 0.112, Zimbabwe = 0.126, Mozambique = 0.202)   # age-collapsed; Zimbabwe pools Harare and Bulawayo

scaled_cfg <- function(k) { c2 <- cfg; c2$dist_sd <- cfg$dist_sd*k; c2$aware_dist_sd <- cfg$aware_dist_sd*k; c2 }

## realised sigma_u^2 = var(w + tau) within province x stratum, at wave 0
realised <- function(k, reps = 60) {
  mean(vapply(seq_len(reps), function(r) {
    p <- simulate_population(scaled_cfg(k), "constant", seed = 9000 + r,
                             lambda_dist_sd = cfg$lambda_dist_sd * k)
    key <- paste(p$idx$prov, p$idx$s)
    var((p$theta[,1] - ave(p$theta[,1], key)) + (p$lambda[,1] - ave(p$lambda[,1], key)))
  }, numeric(1)))
}

cat("solving for the scaling factor at each target sigma_u^2\n")
kk <- vapply(TARGETS, function(tg)
  uniroot(function(k) realised(k) - tg, c(0.4, 6), tol = 1e-3)$root, numeric(1))
ach <- vapply(kk, realised, numeric(1), reps = 200)
for (i in seq_along(TARGETS))
  cat(sprintf("  target %.2f  ->  k = %.3f  ->  achieved %.3f\n", TARGETS[i], kk[i], ach[i]))

agg <- function(prob, group, pop) {
  den <- rowsum(pop, group)
  sapply(seq_len(ncol(prob)), function(t) rowsum(pop*prob[,t], group)[,1]/den[,1])
}
rmse_at <- function(th, p, group)
  sqrt(mean((agg(plogis(th), group, p$pop) - agg(plogis(p$theta), group, p$pop))^2))*100
null_surv  <- function(p) list(Z = matrix(NA,p$U,p$T), psi = matrix(NA,p$U,p$T))
levels_of  <- function(p) list(national = rep(1,p$U), `admin-1 (province)` = p$idx$prov,
                               district = p$idx$d, `district x age-sex` = seq_len(p$U))

rows <- list(); ri <- 0
for (i in seq_along(TARGETS)) {
  c2 <- scaled_cfg(kk[i]); lds <- cfg$lambda_dist_sd * kk[i]
  for (scn in names(cfg$scenarios)) for (rep in seq_len(REPS)) {
    p  <- simulate_population(c2, scn, seed = 5000 + rep, lambda_dist_sd = lds)
    dr <- run_design(p, c2, "random", k = K_AUDIT, seed = rep)
    es <- list(`routine only`      = est_routine_only(p, null_surv(p), c2),
               `survey audit only` = est_survey_only(p, dr$surv, c2),
               `calibrated`        = est_static_calib(p, dr$surv, c2))
    lv <- levels_of(p)
    for (m in names(es)) for (L in names(lv)) {
      ri <- ri + 1
      rows[[ri]] <- data.table(sigma_u2_target = TARGETS[i], sigma_u2 = ach[i],
                               scenario = scn, method = m, level = L,
                               rmse = rmse_at(es[[m]]$theta, p, lv[[L]]))
    }
  }
  cat(sprintf("  done sigma_u^2 = %.2f\n", TARGETS[i]))
}
all <- rbindlist(rows)
## per-scenario results
bs <- all[, .(rmse = mean(rmse)), by = .(sigma_u2_target, sigma_u2, scenario, method, level)]
fwrite(bs, file.path(OUT, "sigma_u_sweep_by_scenario.csv"))
wb <- dcast(bs, sigma_u2_target + sigma_u2 + scenario + level ~ method, value.var = "rmse")
wb[, margin_over_survey := `survey audit only` - calibrated]
wb[, margin_over_register := `routine only` - calibrated]
fwrite(wb, file.path(OUT, "sigma_u_sweep_by_scenario_wide.csv"))
res <- all[, .(rmse = mean(rmse)), by = .(sigma_u2_target, sigma_u2, method, level)]
fwrite(res, file.path(OUT, "sigma_u_sweep.csv"))

w <- dcast(res, sigma_u2_target + sigma_u2 + level ~ method, value.var = "rmse")
w[, margin_over_survey := `survey audit only` - calibrated]
w[, margin_over_register := `routine only` - calibrated]
fwrite(w, file.path(OUT, "sigma_u_sweep_wide.csv"))

cat("\n=== resolution ladder across sigma_u^2 (RMSE, percentage points) ===\n\n")
for (L in c("national","admin-1 (province)","district","district x age-sex")) {
  cat(" ", L, "\n")
  x <- w[level == L]
  print(as.data.frame(x[, .(`sigma_u^2` = round(sigma_u2,3),
                            register = round(`routine only`,2),
                            survey = round(`survey audit only`,2),
                            calibrated = round(calibrated,2),
                            `vs survey` = round(margin_over_survey,2))]), row.names = FALSE)
  cat("\n")
}
cat("countries sit at: ")
cat(paste(sprintf("%s %.3f", names(COUNTRY), COUNTRY), collapse = ", "), "\n")
cat("written: ", file.path(OUT,"sigma_u_sweep.csv"), ", _wide.csv and _by_scenario(_wide).csv\n", sep = "")

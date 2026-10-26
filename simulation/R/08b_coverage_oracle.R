## ============================================================================
## 08b_coverage_oracle.R -- coverage of the static calibrated estimator when
##   sigma_tau^2 is fixed at its generating value (lambda_dist_sd^2) instead of
##   the REML estimate of 02_estimators.R. Same populations and audits as
##   08_coverage.R (population seeds 6000 + r, design seed r).
## Reads : simulation/R/00_config.R to 03_designs.R (sourced)
## Writes: coverage_oracle.csv under simulation/output/
## Usage : Rscript simulation/R/08b_coverage_oracle.R
## ============================================================================
suppressWarnings(suppressMessages(library(data.table)))
here <- "simulation"
for (f in c("00_config","01_dgp","02_estimators","03_designs")) source(file.path(here,"R",paste0(f,".R")))
REPS <- 150; K <- 8; s2 <- cfg$lambda_dist_sd^2
rows <- list()
for (scn in names(cfg$scenarios)) for (r in seq_len(REPS)) {
  p <- simulate_population(cfg, scn, seed = 6000 + r); d <- run_design(p, cfg, "random", k = K, seed = r)
  aud <- !is.na(d$surv$Z); th <- v <- matrix(NA, p$U, p$T)
  for (t in 1:p$T) { lo <- p$z[, t] - d$surv$Z[, t]; rr <- d$surv$psi[, t] + cfg$v_rhis
    for (g in unique(p$idx$s)) { rows_g <- which(p$idx$s == g); sm <- rows_g[aud[rows_g, t]]
      w <- 1 / (s2 + rr[sm]); th[rows_g, t] <- p$z[rows_g, t] - sum(w * lo[sm]) / sum(w)
      v[rows_g, t] <- cfg$v_rhis + 1 / sum(w) + s2 } }
  rows[[length(rows) + 1]] <- data.table(scenario = scn, coverage = mean(abs(th - p$theta) <= 1.96 * sqrt(v)))
}
res <- rbindlist(rows)[, .(coverage_oracle = mean(coverage)), by = scenario]
fwrite(res, file.path(OUT, "coverage_oracle.csv")); print(res)

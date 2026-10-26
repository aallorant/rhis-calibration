## ============================================================================
## 17_sensitivity.R -- three sensitivity sweeps on simulate_population():
##   A  offset drift         how fast may lambda move before carrying the offset
##                           forward stops beating freezing the last estimate
##   B  pooling age and sex  what aggregating over strata costs, against the
##                           between-stratum spread of the offset
##   C  minimal covariates   whether the POOLED calibration still improves on
##                           both its sources when the register has (d, p) but no g
## Reads : simulation/R/00_config.R to 03_designs.R (sourced)
## Writes: sens_drift.csv, sens_stratum.csv under simulation/output/
## Usage : Rscript simulation/R/17_sensitivity.R [REPS]   (default 150)
## ============================================================================
suppressWarnings(suppressMessages({ library(data.table) }))
here <- "simulation"
for (f in c("00_config","01_dgp","02_estimators","03_designs")) source(file.path(here,"R",paste0(f,".R")))
REPS <- { a <- commandArgs(TRUE); if (length(a)) as.integer(a[1]) else 150 }
K_AUDIT <- 8

SD_CROSS   <- cfg$lambda_dist_sd            # offset spread across areas within a year
STR_SD0    <- sd(cfg$lambda0)               # the DGP's own between-stratum spread
rmse_pp <- function(th, p) sqrt(mean((plogis(th) - plogis(p$theta))^2))*100
null_surv <- function(p) list(Z = matrix(NA,p$U,p$T), psi = matrix(NA,p$U,p$T))
cat(sprintf("anchors: cross-sectional offset spread %.3f, between-stratum spread %.3f\n",
            SD_CROSS, STR_SD0))

## ---- A. offset drift -------------------------------------------------------
## axis: drift sd as a fraction of SD_CROSS, the offset's cross-sectional spread.
## One-step-ahead: calibrate at wave t, judge the carried estimate at wave t+1.
## freeze  carries the calibrated ESTIMATE forward; update carries the OFFSET
## forward and applies it to the live register.
s_grid <- seq(0, 1, by = 0.1)
A <- rbindlist(lapply(s_grid, function(s) {
  ef <- eu <- numeric(REPS)
  for (r in seq_len(REPS)) {
    scn <- list(type = "drift", slope_sd = s * SD_CROSS)
    p <- simulate_population(cfg, "constant", seed = 5000 + r)
    p$lambda <- p$lambda + lambda_dynamics(scn, p$U, p$T)
    p$z <- p$theta + p$lambda + matrix(rnorm(p$U*p$T, 0, sqrt(cfg$v_rhis)), p$U)
    d <- run_design(p, cfg, "random", k = K_AUDIT, seed = r)
    e <- est_static_calib(p, d$surv, cfg)
    tt <- 1:(p$T - 1)
    ef[r] <- sqrt(mean((plogis(e$theta[, tt]) - plogis(p$theta[, tt + 1]))^2))*100
    lam_hat <- p$z[, tt] - e$theta[, tt]                      # offset implied at t
    eu[r] <- sqrt(mean((plogis(p$z[, tt + 1] - lam_hat) - plogis(p$theta[, tt + 1]))^2))*100
  }
  data.table(s = s, rmse_freeze = mean(ef), rmse_update = mean(eu))
}))
fwrite(A, file.path(OUT, "sens_drift.csv"))

## ---- B and C. between-stratum spread ---------------------------------------
## axis: between-stratum sd of the offset, set by scaling cfg$lambda0 about its mean
sp_grid <- seq(0, 0.6, by = 0.06)
BC <- rbindlist(lapply(sp_grid, function(sp) {
  k <- sp / STR_SD0
  er <- ep <- eg <- es <- numeric(REPS)
  for (r in seq_len(REPS)) {
    p <- simulate_population(cfg, "constant", seed = 5000 + r, lambda_stratum_scale = k)
    d <- run_design(p, cfg, "random", k = K_AUDIT, seed = r)
    er[r] <- rmse_pp(est_static_calib(p, d$surv, cfg)$theta, p)          # stratum-resolved
    ep[r] <- rmse_pp(est_static_calib_pooled(p, d$surv, cfg)$theta, p)   # pooled over strata
    eg[r] <- rmse_pp(est_routine_only(p, null_surv(p), cfg)$theta, p)    # register alone
    es[r] <- rmse_pp(est_survey_only(p, d$surv, cfg)$theta, p)           # survey alone
  }
  data.table(str_spread = sp, rmse_resolved = mean(er), rmse_pooled = mean(ep),
             rmse_register = mean(eg), rmse_survey = mean(es))
}))
fwrite(BC, file.path(OUT, "sens_stratum.csv"))

## ---- report ----------------------------------------------------------------
cat("\n=== A. offset drift (RMSE pp, one step ahead) ===\n")
print(as.data.frame(A[, .(`drift / spread` = s, freeze = round(rmse_freeze,2),
                          update = round(rmse_update,2),
                          `update gain` = round(rmse_freeze - rmse_update,2))]), row.names = FALSE)
cr <- A[rmse_update > rmse_freeze]
cat(if (nrow(cr)) sprintf("  updating stops winning beyond drift = %.1f of the spread\n", min(cr$s))
    else "  updating wins across the whole range\n")

cat("\n=== B/C. between-stratum spread ===\n")
print(as.data.frame(BC[, .(spread = str_spread, resolved = round(rmse_resolved,2),
                           pooled = round(rmse_pooled,2), register = round(rmse_register,2),
                           survey = round(rmse_survey,2),
                           `pool gain` = round(rmse_resolved - rmse_pooled,2),
                           `pooled vs survey` = round(rmse_survey - rmse_pooled,2))]), row.names = FALSE)
ipl <- function(x, y) approx(BC$str_spread, y, x)$y
for (k in list(c("Mozambique",0.16), c("Zimbabwe",0.28), c("Malawi",0.30), c("Lesotho",0.43))) {
  s <- as.numeric(k[2])
  cat(sprintf("  %-11s spread %.2f | resolved %.2f pooled %.2f | pool gain %+.2f | pooled vs register %+.2f vs survey %+.2f\n",
      k[1], s, ipl(s, BC$rmse_resolved), ipl(s, BC$rmse_pooled),
      ipl(s, BC$rmse_resolved) - ipl(s, BC$rmse_pooled),
      ipl(s, BC$rmse_register) - ipl(s, BC$rmse_pooled),
      ipl(s, BC$rmse_survey)   - ipl(s, BC$rmse_pooled)))
}
cat("\nwritten: sens_drift.csv, sens_stratum.csv\n")

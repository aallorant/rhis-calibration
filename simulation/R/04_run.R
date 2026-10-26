## ============================================================================
## 04_run.R -- run the Monte-Carlo experiment and save results + example data.
## Reads : simulation/R/00_config.R to 03_designs.R (sourced)
## Writes (under simulation/output/):
##   * frontier.csv      : cost vs accuracy for every method/design/budget/rep
##   * example_pops.rds  : one fixed population per lambda scenario (for figures)
##   * example_fit.rds   : truth vs estimators at a mid budget (calibration demo)
##   * epi_summary_stratum.csv, epi_summary_province.csv : age-sex & province epidemiology
## Usage:  Rscript simulation/R/04_run.R [REPS]      (REPS default 200)
## ============================================================================

suppressWarnings(suppressMessages({ library(data.table) }))
here <- "simulation"
source(file.path(here, "R/00_config.R"))
source(file.path(here, "R/01_dgp.R"))
source(file.path(here, "R/02_estimators.R"))
source(file.path(here, "R/03_designs.R"))

args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1) as.integer(args[1]) else 200
SCEN <- names(cfg$scenarios)            # constant, drift, shock, reshuffle
KGRID <- c(2, 4, 6, 8, 12, 20, 40)
S <- cfg$S
cat(sprintf("Running %d reps x %d scenarios x %d budgets ...\n", REPS, length(SCEN), length(KGRID)))

null_surv <- function(pop) list(Z = matrix(NA, pop$U, pop$T),
                                psi = matrix(NA, pop$U, pop$T))

rows <- list(); ri <- 0
add <- function(scn, method, design, k, cost, sc) {
  ri <<- ri + 1
  rows[[ri]] <<- data.table(scenario = scn, method = method, design = design,
                            k = k, cost = cost,
                            rmse_logit = sc$rmse_logit, rmse_pp = sc$rmse_pp,
                            bias_logit = sc$bias_logit)
}

t0 <- Sys.time()
for (scn in SCEN) {
  for (rep in seq_len(REPS)) {
    pop <- simulate_population(cfg, scenario = scn, seed = 1000 * match(scn, SCEN) + rep)

    ## routine-only (free, biased)
    add(scn, "routine_only", "none", 0, 0, score(pop, est_routine_only(pop, null_surv(pop), cfg)))

    ## survey-only at full budget (expensive ceiling)
    df <- run_design(pop, cfg, "full", k = pop$D, seed = rep)
    add(scn, "survey_only", "full", pop$D, cost_of(df$nmat, S),
        score(pop, est_survey_only(pop, df$surv, cfg)))

    for (k in KGRID) {
      ## random sparse design
      dr <- run_design(pop, cfg, "random", k = k, seed = 7 * rep + k)
      cst <- cost_of(dr$nmat, S)
      add(scn, "survey_only",   "random", k, cst, score(pop, est_survey_only(pop, dr$surv, cfg)))
      add(scn, "static_calib",  "random", k, cst, score(pop, est_static_calib(pop, dr$surv, cfg)))
      add(scn, "dynamic_calib", "random", k, cst, score(pop, est_dynamic_calib(pop, dr$surv, cfg)))

      ## adaptive design (dynamic calibration)
      da <- run_design(pop, cfg, "adaptive", k = k, seed = 13 * rep + k)
      add(scn, "dynamic_calib", "adaptive", k, cost_of(da$nmat, S),
          score(pop, est_dynamic_calib(pop, da$surv, cfg)))
    }
  }
  cat(sprintf("  scenario %-9s done (%.1fs)\n", scn, as.numeric(Sys.time() - t0, units = "secs")))
}

res <- rbindlist(rows)
fwrite(res, file.path(OUT, "frontier.csv"))
cat("wrote frontier.csv:", nrow(res), "rows\n")

## ---- example populations + one fitted example (for figures) ---------------
example_pops <- lapply(SCEN, function(s) simulate_population(cfg, s, seed = 101))
names(example_pops) <- SCEN
saveRDS(example_pops, file.path(OUT, "example_pops.rds"))

popD <- example_pops[["drift"]]
## Monte-Carlo distribution of the dynamic estimator over survey draws
## (fixed population, random design k = 8)
MCEX <- 120; k_ex <- 8
acc_dyn <- acc_dyn2 <- array(0, dim = c(popD$U, popD$T))
qs_dyn  <- vector("list", MCEX); one_survey <- NULL
for (r in seq_len(MCEX)) {
  dr <- run_design(popD, cfg, "random", k = k_ex, seed = 1000 + r)
  ed <- est_dynamic_calib(popD, dr$surv, cfg)
  acc_dyn  <- acc_dyn  + plogis(ed$theta)
  acc_dyn2 <- acc_dyn2 + plogis(ed$theta)^2
  qs_dyn[[r]] <- plogis(ed$theta)
  if (r == 1) one_survey <- est_survey_only(popD, dr$surv, cfg)  # one realisation, for points
}
arr <- simplify2array(qs_dyn)                                   # U x T x MCEX
## summarised by the median (transformation-invariant) and the 10-90% band
dyn_med <- apply(arr, c(1, 2), median)
dyn_lo  <- apply(arr, c(1, 2), quantile, 0.10)
dyn_hi  <- apply(arr, c(1, 2), quantile, 0.90)

## also keep one adaptive vs random allocation for the design-map figure
set.seed(42)
da6 <- run_design(popD, cfg, "adaptive", k = 6, seed = 1)
dr6 <- run_design(popD, cfg, "random",   k = 6, seed = 1)
example_fit <- list(
  pop = popD, k_ex = k_ex, MCEX = MCEX,
  routine = plogis(popD$z),                       # biased, no survey dependence
  truth   = plogis(popD$theta),
  dyn_med = dyn_med, dyn_lo = dyn_lo, dyn_hi = dyn_hi,
  survey_one = one_survey,
  alloc_adaptive = da6$nmat, alloc_random = dr6$nmat
)
saveRDS(example_fit, file.path(OUT, "example_fit.rds"))

## ---- epidemiology realism summary (from example drift population) ---------
pop <- popD; idx <- pop$idx
epi_strat <- data.table(
  stratum = cfg$strata$stratum,
  prevalence_pct  = sapply(cfg$strata$stratum, function(s)
                      mean(pop$prev[idx$stratum == s, 1]) * 100),
  undiagnosed_pct = sapply(cfg$strata$stratum, function(s)
                      mean(pop$undiag[idx$stratum == s, 1]) * 100),
  awareness_pct   = sapply(cfg$strata$stratum, function(s)
                      mean(pop$aware[idx$stratum == s, 1]) * 100))
epi_strat[, undiagnosed_share_pct := undiagnosed_pct / prevalence_pct * 100]
fwrite(epi_strat, file.path(OUT, "epi_summary_stratum.csv"))

epi_prov <- data.table(
  province = seq_len(pop$G),
  prevalence_pct = sapply(seq_len(pop$G), function(g)
    mean(pop$prev[idx$prov == g, 1]) * 100))[order(prevalence_pct)]
fwrite(epi_prov, file.path(OUT, "epi_summary_province.csv"))

cat("wrote example_pops.rds, example_fit.rds, epi summaries\n")
print(epi_strat)

## ============================================================================
## 07_resolution.R -- the resolution ladder: RMSE of each method as the target
##   gets finer (national -> admin-1 -> district -> district x age-sex), with a
##   fixed periodic audit of K_AUDIT districts per wave, averaged over the four
##   offset scenarios.
## Reads : simulation/R/00_config.R to 03_designs.R (sourced)
## Writes: resolution_ladder.csv, resolution_ladder_by_scenario.csv and
##         figures/fig7_resolution_ladder.{pdf,png} under simulation/output/
## Usage : Rscript simulation/R/07_resolution.R [REPS]   (default 150)
## ============================================================================
suppressWarnings(suppressMessages({ library(data.table); library(ggplot2) }))
here <- "simulation"
for (f in c("00_config","01_dgp","02_estimators","03_designs")) source(file.path(here, "R", paste0(f, ".R")))
REPS <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 150 }
K_AUDIT <- 8                                   # fixed periodic audit: 8 districts/wave

## population-weighted prevalence aggregated to a grouping, per wave
agg_level <- function(prob, group, pop) {
  den <- rowsum(pop, group)
  sapply(seq_len(ncol(prob)), function(t) rowsum(pop * prob[, t], group)[, 1] / den[, 1])
}
rmse_at <- function(est_theta, p, group)
  sqrt(mean((agg_level(plogis(est_theta), group, p$pop) -
             agg_level(plogis(p$theta), group, p$pop))^2)) * 100

null_surv <- function(p) list(Z = matrix(NA, p$U, p$T), psi = matrix(NA, p$U, p$T))
levels_of <- function(p) list(national = rep(1, p$U), `admin-1 (province)` = p$idx$prov,
                              district = p$idx$d, `district x age-sex` = seq_len(p$U))

rows <- list(); ri <- 0
## REPS replications per offset scenario; the ladder averages over all four
for (scn in names(cfg$scenarios)) for (rep in seq_len(REPS)) {
  p <- simulate_population(cfg, scn, seed = 5000 + rep)
  dr <- run_design(p, cfg, "random", k = K_AUDIT, seed = rep)
  ests <- list(`routine only`         = est_routine_only(p, null_surv(p), cfg),
               `survey audit only`    = est_survey_only(p, dr$surv, cfg),
               `calibrated (static)`  = est_static_calib(p, dr$surv, cfg),
               `calibrated (dynamic)` = est_dynamic_calib(p, dr$surv, cfg))
  lv <- levels_of(p)
  for (m in names(ests)) for (L in names(lv)) {
    ri <- ri + 1
    rows[[ri]] <- data.table(scenario = scn, method = m, level = L,
                             rmse = rmse_at(ests[[m]]$theta, p, lv[[L]]))
  }
}
all_rows <- rbindlist(rows)
## per-scenario breakdown
by_scn <- all_rows[, .(rmse = mean(rmse)), by = .(scenario, method, level)]
fwrite(by_scn, file.path(OUT, "resolution_ladder_by_scenario.csv"))
res <- all_rows[, .(rmse = mean(rmse)), by = .(method, level)]
lev <- c("national","admin-1 (province)","district","district x age-sex")
mth <- c("routine only","survey audit only","calibrated (static)","calibrated (dynamic)")
res[, level := factor(level, levels = lev)]; res[, method := factor(method, levels = mth)]
fwrite(res, file.path(OUT, "resolution_ladder.csv"))

th <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "bottom")
## figure: the three core methods (calibrated = static)
res3 <- res[method %in% c("routine only","survey audit only","calibrated (static)")]
res3[, method := factor(as.character(method),
     levels = c("routine only","survey audit only","calibrated (static)"),
     labels = c("routine only (biased)","survey audit only","register + audit (calibrated)"))]
p7 <- ggplot(res3, aes(level, rmse, colour = method, group = method)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  scale_colour_manual(values = c("routine only (biased)" = "#d1495b",
                                 "survey audit only" = "#8d99ae",
                                 "register + audit (calibrated)" = "#00798c"), name = NULL) +
  labs(title = "The resolution ladder (average over the four offset scenarios, audit of 8 districts per wave)",
       subtitle = "Surveys suffice at national level; the register's value grows as the target gets finer. Honest ceiling: the rarest cells.",
       x = "Target resolution (coarse -> fine)", y = "RMSE of undiagnosed positivity (pts)") + th
ggsave(file.path(OUT, "figures/fig7_resolution_ladder.pdf"), p7, width = 9, height = 5.5)
ggsave(file.path(OUT, "figures/fig7_resolution_ladder.png"), p7, width = 9, height = 5.5, dpi = 150)
cat("resolution ladder (RMSE pts):\n"); print(dcast(res, method ~ level, value.var = "rmse"))

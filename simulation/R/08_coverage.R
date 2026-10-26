## ============================================================================
## 08_coverage.R -- empirical coverage and median width of the 95% intervals of
##   the calibrated estimators (static with and without the sigma_tau^2 term,
##   dynamic), fixed periodic audit of K_AUDIT districts per wave; reported for
##   all cells and split by audited / unaudited cells.
## Reads : simulation/R/00_config.R to 03_designs.R (sourced)
## Writes: coverage.csv, coverage_by_cells.csv and figures/fig8_coverage.{pdf,png}
##         under simulation/output/
## Usage : Rscript simulation/R/08_coverage.R [REPS]   (default 150)
## ============================================================================
suppressWarnings(suppressMessages({ library(data.table); library(ggplot2) }))
here <- "simulation"
for (f in c("00_config","01_dgp","02_estimators","03_designs")) source(file.path(here, "R", paste0(f, ".R")))
REPS <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 150 }
K_AUDIT <- 8

## coverage and median interval width (prevalence points) over the cells in `mask`
cw1 <- function(est, p, mask) {
  sd <- sqrt(pmax(est$var, 0)); lo <- est$theta - 1.96 * sd; hi <- est$theta + 1.96 * sd
  c(coverage = mean((p$theta >= lo & p$theta <= hi)[mask]),
    width    = median(((plogis(hi) - plogis(lo)) * 100)[mask]))   # prevalence pts
}
rows <- list(); ri <- 0
for (scn in names(cfg$scenarios)) for (rep in seq_len(REPS)) {
  p  <- simulate_population(cfg, scn, seed = 6000 + rep)
  dr <- run_design(p, cfg, "random", k = K_AUDIT, seed = rep)
  ests <- list(static        = est_static_calib(p, dr$surv, cfg, tau = TRUE),
               static_no_tau = est_static_calib(p, dr$surv, cfg, tau = FALSE),
               dynamic       = est_dynamic_calib(p, dr$surv, cfg))
  aud <- !is.na(dr$surv$Z)
  for (m in names(ests)) for (cell in c("all", "audited", "unaudited")) {
    mask <- switch(cell, all = matrix(TRUE, p$U, p$T), audited = aud, unaudited = !aud)
    cw <- cw1(ests[[m]], p, mask); ri <- ri + 1
    rows[[ri]] <- data.table(scenario = scn, method = m, cells = cell,
                             coverage = cw[["coverage"]], width = cw[["width"]])
  }
}
res_all <- rbindlist(rows)[, .(coverage = mean(coverage), width = mean(width)),
                           by = .(scenario, method, cells)]
fwrite(res_all, file.path(OUT, "coverage_by_cells.csv"))
res <- res_all[cells == "all", .(scenario, method, coverage, width)]
fwrite(res, file.path(OUT, "coverage.csv"))

th <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "bottom")
res[, scenario := factor(scenario, levels = names(cfg$scenarios))]
p8 <- ggplot(res, aes(scenario, coverage, fill = method)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "#d1495b") +
  annotate("text", x = 0.7, y = 0.965, label = "nominal 95%", size = 3, colour = "#d1495b", hjust = 0) +
  scale_fill_manual(values = c(static = "#edae49", static_no_tau = "#d1495b", dynamic = "#00798c"), name = NULL) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(title = "Interval coverage of the calibrated estimator (95% target)",
       subtitle = "Fixed audit, 8 districts/wave. Closer to the dashed line is better-calibrated uncertainty.",
       x = NULL, y = "empirical coverage of 95% interval") + th
ggsave(file.path(OUT, "figures/fig8_coverage.pdf"), p8, width = 9, height = 5)
ggsave(file.path(OUT, "figures/fig8_coverage.png"), p8, width = 9, height = 5, dpi = 150)
cat("coverage & median interval width (pts):\n"); print(dcast(res_all, scenario + method ~ cells, value.var = c("coverage","width")))

## ============================================================================
## 13_scenario_paths.R -- trajectories of the offset lambda and the estimand eta
##   under the four offset scenarios, from one seed so that eta is identical
##   across scenarios and only lambda differs (checked with stopifnot).
## Reads : simulation/R/00_config.R and 01_dgp.R (sourced)
## Writes: figures/fig_scenario_paths.{pdf,png}         (both rows)
##         figures/fig_scenario_paths_lambda.{pdf,png}  (offset row only)
##         under simulation/output/
## Usage : Rscript simulation/R/13_scenario_paths.R
## ============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2)
}))
here <- "simulation"
for (f in c("00_config", "01_dgp")) source(file.path(here, "R", paste0(f, ".R")))

SEED    <- 4242L      # one seed for all four, so eta is common by construction
N_SHOW  <- 40L        # units drawn as background lines
LAB <- c(constant = "Constant", drift = "Drift",
         shock = "Shock", reshuffle = "Reshuffle")

set.seed(1)
rows <- list()
for (scn in names(cfg$scenarios)) {
  p <- simulate_population(cfg, scn, seed = SEED)
  show <- sort(sample(nrow(p$lambda), N_SHOW))
  for (q in c("lambda", "eta")) {
    M <- if (q == "lambda") p$lambda else p$theta   # p$theta is already logit(theta)
    rows[[length(rows) + 1L]] <- data.table(
      scenario = LAB[[scn]], quantity = q,
      unit = rep(show, times = ncol(M)),
      wave = rep(seq_len(ncol(M)) - 1L, each = length(show)),
      value = as.vector(M[show, ]))
    rows[[length(rows) + 1L]] <- data.table(
      scenario = LAB[[scn]], quantity = q, unit = -1L,
      wave = seq_len(ncol(M)) - 1L, value = colMeans(M))
  }
}
d <- rbindlist(rows)
d[, scenario := factor(scenario, levels = unname(LAB))]
d[, quantity := factor(quantity, levels = c("lambda", "eta"),
     labels = c("Offset~lambda", "Estimand~eta==logit(theta)"))]

## check that the mean eta path is identical across the four scenarios
e <- dcast(d[grepl("Estimand", quantity) & unit == -1L],
           wave ~ scenario, value.var = "value")
stopifnot(max(abs(e[[2]] - e[[3]])) < 1e-12,
          max(abs(e[[2]] - e[[4]])) < 1e-12,
          max(abs(e[[2]] - e[[5]])) < 1e-12)
message("verified: mean eta path identical across all four scenarios")

g <- ggplot(d[unit != -1L], aes(wave, value, group = unit)) +
  geom_line(colour = "grey72", linewidth = 0.25) +
  geom_line(data = d[unit == -1L], aes(group = 1), colour = "#b2182b", linewidth = 1) +
  facet_grid(quantity ~ scenario, scales = "free_y",
             labeller = labeller(quantity = label_parsed)) +
  scale_x_continuous(breaks = 0:(cfg$T - 1)) +
  labs(x = "Survey wave", y = "Log-odds scale") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.9, "lines"))

FIG <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIG, "fig_scenario_paths.pdf"), g, width = 10, height = 5)
ggsave(file.path(FIG, "fig_scenario_paths.png"), g, width = 10, height = 5, dpi = 150)
cat("wrote", file.path(FIG, "fig_scenario_paths.pdf"), "\n")

## the offset row alone
gl <- ggplot(d[unit != -1L & grepl("Offset", quantity)], aes(wave, value, group = unit)) +
  geom_line(colour = "grey72", linewidth = 0.25) +
  geom_line(data = d[unit == -1L & grepl("Offset", quantity)], aes(group = 1),
            colour = "#b2182b", linewidth = 1) +
  facet_grid(. ~ scenario) +
  scale_x_continuous(breaks = 0:(cfg$T - 1)) +
  labs(x = "Survey wave", y = "Offset (log-odds scale)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.9, "lines"))
ggsave(file.path(FIG, "fig_scenario_paths_lambda.pdf"), gl, width = 10, height = 3)
ggsave(file.path(FIG, "fig_scenario_paths_lambda.png"), gl, width = 10, height = 3, dpi = 150)
cat("wrote", file.path(FIG, "fig_scenario_paths_lambda.pdf"), "\n")

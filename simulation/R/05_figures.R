## ============================================================================
## 05_figures.R -- pedagogical and results figures.
## Reads : example_pops.rds, example_fit.rds, epi_summary_stratum.csv and
##         frontier.csv from simulation/output/ (written by 04_run.R)
## Writes: PDF + PNG of each figure to simulation/output/figures/
##   fig1_epidemiology           : admin-1 gradient + age-sex prevalence vs undiagnosed
##   fig2_lambda_scenarios       : the four selection-offset scenarios over time
##   fig3_calibration_in_action  : truth vs routine/calibrated for example units
##   fig4_cost_frontier_drift, fig4b_cost_frontier_all : cost vs accuracy by method
##   fig5_design_allocation      : where the survey budget lands (random vs adaptive)
## Usage : Rscript simulation/R/05_figures.R
## ============================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
}))
here <- "simulation"
source(file.path(here, "R/00_config.R"))
FIG <- file.path(OUT, "figures")
sv <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
}
th <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")
pal <- c(routine_only = "#d1495b", survey_only = "#8d99ae",
         static_calib = "#edae49", dynamic_calib = "#00798c",
         "dynamic_calib (adaptive)" = "#003d5b")

pops <- readRDS(file.path(OUT, "example_pops.rds"))
fit  <- readRDS(file.path(OUT, "example_fit.rds"))
pop  <- fit$pop; idx <- pop$idx

## ---- FIG 1: epidemiology realism ------------------------------------------
epi_s <- fread(file.path(OUT, "epi_summary_stratum.csv"))
mlt <- melt(epi_s[, .(stratum, Prevalence = prevalence_pct,
                      `Undiagnosed positivity` = undiagnosed_pct)],
            id.vars = "stratum")
p1a <- ggplot(mlt, aes(stratum, value, fill = variable)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Prevalence" = "#8d99ae",
                               "Undiagnosed positivity" = "#d1495b"), name = NULL) +
  labs(title = "A. Age-sex structure", subtitle =
       "Undiagnosed burden shifts toward men",
       x = NULL, y = "Percent") + th
prov <- data.table(province = idx$prov, prev = pop$prev[, 1] * 100)[
  , .(prev = mean(prev)), by = province][order(prev)]
prov[, rank := .I]
p1b <- ggplot(prov, aes(rank, prev)) +
  geom_col(fill = "#00798c") +
  ## integer breaks: the x axis is a rank
  scale_x_continuous(breaks = seq_len(nrow(prov))) +
  labs(title = "B. Admin-1 gradient", subtitle =
       "Realistic province range",
       x = "Province (ranked)", y = "HIV prevalence (%)") + th
sv(p1a + p1b, "fig1_epidemiology", 11, 4.5)

## ---- FIG 2: lambda scenarios ----------------------------------------------
sc_df <- rbindlist(lapply(names(pops), function(s) {
  P <- pops[[s]]; lam <- P$lambda
  data.table(scenario = s, wave = rep(1:P$T, each = P$U),
             unit = rep(1:P$U, P$T), lambda = as.vector(lam))
}))
sc_mean <- sc_df[, .(lambda = mean(lambda)), by = .(scenario, wave)]
set.seed(1); samp_u <- sample(pop$U, 12)
p2 <- ggplot(sc_df[unit %in% samp_u], aes(wave, lambda, group = unit)) +
  geom_line(alpha = 0.25, colour = "#00798c") +
  geom_line(data = sc_mean, aes(wave, lambda, group = 1),
            colour = "#d1495b", linewidth = 1.1) +
  facet_wrap(~ scenario, nrow = 1) +
  labs(title = "Selection-offset scenarios (lambda over time)",
       subtitle = "Grey: example units. Red: mean. The offset a single survey must learn and transport.",
       x = "Survey wave", y = expression(lambda~"= logit(routine) - logit(undiagnosed)")) + th
sv(p2, "fig2_lambda_scenarios", 12, 3.8)

## ---- FIG 3: calibration in action (Monte-Carlo behaviour) ------------------
## three illustrative units on the undiagnosed-positivity (%) scale; calibrated
## line = median over survey draws, band = 10-90%
pick <- c(which(idx$stratum == "M.25-49")[1],
          which(idx$stratum == "F.25-49")[1],
          which(idx$stratum == "M.15-24")[1])
mk <- function(u) data.table(
  unit = u, wave = 1:pop$T,
  truth   = fit$truth[u, ] * 100,
  routine = fit$routine[u, ] * 100,
  cal     = fit$dyn_med[u, ] * 100,
  lo      = fit$dyn_lo[u, ] * 100,
  hi      = fit$dyn_hi[u, ] * 100)
d3 <- rbindlist(lapply(pick, mk))
d3[, label := factor(unit, levels = pick,
       labels = paste0(idx$stratum[pick], " (district ", idx$d[pick], ")"))]
p3 <- ggplot(d3, aes(wave)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#00798c", alpha = 0.18) +
  geom_line(aes(y = routine, colour = "Routine only (biased)"), linewidth = 1) +
  geom_line(aes(y = cal, colour = "Calibrated (dynamic), median"), linewidth = 1.1) +
  geom_line(aes(y = truth, colour = "Truth"), linewidth = 1, linetype = "22") +
  facet_wrap(~ label, scales = "free_y") +
  scale_colour_manual(values = c("Truth" = "black",
      "Routine only (biased)" = "#d1495b",
      "Calibrated (dynamic), median" = "#003d5b"), name = NULL) +
  labs(title = "Calibration in action (drift scenario)",
       subtitle = paste0("Median over ", fit$MCEX, " survey realisations (random design, k = ",
                         fit$k_ex, "); shaded band = 10-90%. Routine alone is biased by lambda; calibration tracks truth."),
       x = "Survey wave", y = "Undiagnosed prevalence (%)") + th
sv(p3, "fig3_calibration_in_action", 12, 4.5)

## ---- FIG 4: cost-accuracy frontier ----------------------------------------
d <- fread(file.path(OUT, "frontier.csv"))
agg <- d[, .(rmse = mean(rmse_pp), se = sd(rmse_pp) / sqrt(.N),
             bias = mean(bias_logit), cost = mean(cost)),
         by = .(scenario, method, design, k)]
agg[, lab := method]
## same-budget comparison: survey audit alone vs register + audit (calibrated);
## routine-only is the biased free floor
fpal <- c(`survey audit only` = "#8d99ae", `register + audit (calibrated)` = "#00798c")
rout <- agg[method == "routine_only"]

mkfront <- function(scn) {
  dd <- agg[scenario == scn & design == "random" & method %in% c("survey_only", "static_calib")]
  dd[, lab := ifelse(method == "survey_only", "survey audit only", "register + audit (calibrated)")]
  ggplot(dd, aes(cost / 1000, rmse, colour = lab)) +
    geom_hline(data = rout[scenario == scn], aes(yintercept = rmse),
               linetype = "dotted", colour = "#d1495b") +
    annotate("text", x = max(dd$cost)/1000*0.5, y = rout[scenario == scn]$rmse,
             label = "routine register only (biased, but free)", vjust = -0.6, size = 3, colour = "#d1495b") +
    geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
    scale_colour_manual(values = fpal, name = NULL) +
    labs(title = paste0("Same budget, with vs without the register -- ", scn, " scenario"),
         subtitle = "At a fixed audit budget, calibrating the register beats spending the same survey effort alone",
         x = "Survey audit cost (thousands of interviews)",
         y = "RMSE of undiagnosed prevalence (pts)") + th
}
sv(mkfront("drift"), "fig4_cost_frontier_drift", 8.5, 5.5)
sv((mkfront("constant") | mkfront("shock")) /
   (mkfront("drift") | mkfront("reshuffle")), "fig4b_cost_frontier_all", 13, 10)

## ---- FIG 5: where the budget lands (random vs adaptive) -------------------
alloc <- function(mat, name) data.table(
  district = rep(1:nrow(mat), ncol(mat)), wave = rep(1:ncol(mat), each = nrow(mat)),
  sampled = as.vector(mat) > 0, design = name)
ad <- rbind(alloc(fit$alloc_random, "random"), alloc(fit$alloc_adaptive, "adaptive"))
p5 <- ggplot(ad, aes(wave, district, fill = sampled)) +
  geom_tile(colour = "white") + facet_wrap(~ design) +
  scale_fill_manual(values = c(`TRUE` = "#00798c", `FALSE` = "grey90"), name = "surveyed") +
  labs(title = "Where the survey budget lands (k = 6 districts / wave)",
       subtitle = "Simple adaptive rule maximises coverage but rarely revisits to track drift",
       x = "Survey wave", y = "District") + th
sv(p5, "fig5_design_allocation", 10, 6)

cat("figures written to", FIG, "\n")

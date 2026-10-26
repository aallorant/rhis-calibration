## ============================================================================
## 14_epi_realism.R -- epidemiology-realism summary of the simulated populations
##   (national undiagnosed positivity, prevalence and awareness by wave, selection
##   odds ratio exp(lambda), cell-level spread), on the same populations as
##   07_resolution.R (REPS replications per scenario, seeds 5000 + rep). National
##   figures are population-weighted on the probability scale with cfg$pop_share.
## Reads : simulation/R/00_config.R and 01_dgp.R (sourced)
## Writes: epi_realism.csv under simulation/output/
## Usage : Rscript simulation/R/14_epi_realism.R [REPS]   (default 150)
## ============================================================================
suppressWarnings(suppressMessages(library(data.table)))
here <- "simulation"
for (f in c("00_config", "01_dgp")) source(file.path(here, "R", paste0(f, ".R")))
REPS <- { a <- commandArgs(TRUE); if (length(a)) as.integer(a[1]) else 150 }

one <- function(p) {
  w <- p$pop / sum(p$pop)
  or0 <- exp(p$lambda[, 1])
  data.table(
    undiag_w0 = sum(p$undiag[, 1] * w),  undiag_w5 = sum(p$undiag[, p$T] * w),
    undiag_mean = mean(colSums(p$undiag * w)),
    prev_w0 = sum(p$prev[, 1] * w),      prev_w5 = sum(p$prev[, p$T] * w),
    aware_w0 = sum(p$aware[, 1] * w),    aware_w5 = sum(p$aware[, p$T] * w),
    or_mean = mean(or0), or_min = min(or0), or_max = max(or0),
    cell_min = min(p$undiag), cell_max = max(p$undiag))
}

rows <- list(); i <- 0
for (scn in names(cfg$scenarios)) for (rep in seq_len(REPS)) {
  i <- i + 1; rows[[i]] <- cbind(scenario = scn, one(simulate_population(cfg, scn, seed = 5000 + rep)))
}
all <- rbindlist(rows)

## range rows: mean over replications of each population's own minimum and maximum
out <- data.table(
  quantity = c("undiagnosed positivity, wave 0", "undiagnosed positivity, wave 5",
               "undiagnosed positivity, mean over waves",
               "adult prevalence, wave 0", "adult prevalence, wave 5",
               "awareness of status, wave 0", "awareness of status, wave 5",
               "selection odds ratio at wave 0, mean",
               "selection odds ratio at wave 0, mean minimum cell",
               "selection odds ratio at wave 0, mean maximum cell",
               "cell undiagnosed positivity, mean minimum over waves",
               "cell undiagnosed positivity, mean maximum over waves"),
  value = c(mean(all$undiag_w0) * 100, mean(all$undiag_w5) * 100,
            mean(all$undiag_mean) * 100,
            mean(all$prev_w0) * 100, mean(all$prev_w5) * 100,
            mean(all$aware_w0) * 100, mean(all$aware_w5) * 100,
            mean(all$or_mean), mean(all$or_min), mean(all$or_max),
            mean(all$cell_min) * 100, mean(all$cell_max) * 100),
  unit = c(rep("per cent", 7), rep("odds ratio", 3), rep("per cent", 2)))

fwrite(out, file.path(OUT, "epi_realism.csv"))
cat("\n=== epidemiology realism, ", REPS, " reps x 4 scenarios, seeds 5000+rep ===\n\n", sep = "")
print(out[, .(quantity, value = round(value, 2), unit)], row.names = FALSE)
cat("\nwritten: ", file.path(OUT, "epi_realism.csv"), "\n", sep = "")

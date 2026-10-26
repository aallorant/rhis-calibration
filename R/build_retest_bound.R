## =====================================================================
## build_retest_bound.R -- the re-testing bound on the national selection odds ratio. A share
## f of routine tests on already-diagnosed adults that is common to a country cannot exceed the
## positivity of its least positive area, rho = min_d y_d/n_d; removing it from the national
## positivity theta_bar gives the log-odds shift Delta = logit(theta_bar) - logit((theta_bar - rho)/(1 - rho))
## and the adjusted odds ratio exp(log OR_pop - Delta). Units: admin-2 register units of the
## age-collapsed fit after harmonise_register_units() (Lesotho: its 72 community councils).
## Reads:  output/<country>/fits/district_predictions.rds, output/canonical/estimates/mainfit_national_popweighted.csv
## Writes: output/canonical/estimates/table_retest_bound.csv.  Usage: Rscript R/build_retest_bound.R
## =====================================================================
suppressMessages(library(dplyr))
source("R/utils.R")
lg <- function(p) log(p / (1 - p))
orp <- read.csv("output/canonical/estimates/mainfit_national_popweighted.csv")
out <- bind_rows(lapply(c("Mozambique", "Zimbabwe", "Malawi", "Lesotho"), function(c) {
  d <- harmonise_register_units(as.data.frame(readRDS(file.path("output", tolower(c), "fits", "district_predictions.rds"))), c) |>
         filter(n_tests > 0)
  pos <- d$n_positive / d$n_tests; rho <- min(pos); tb <- sum(d$n_positive) / sum(d$n_tests)
  Delta <- lg(tb) - lg((tb - rho) / (1 - rho))
  OR <- orp$OR_pop[orp$country == c]
  data.frame(country = c, units = nrow(d), least_positive_unit = d$district[which.min(pos)],
             rho = rho, theta_bar = tb, Delta = Delta, OR_pop = OR, OR_adj = exp(log(OR) - Delta))
}))
write.csv(out, "output/canonical/estimates/table_retest_bound.csv", row.names = FALSE)
print(out |> mutate(across(where(is.numeric), ~ signif(.x, 4))), row.names = FALSE)

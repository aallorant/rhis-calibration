## regen_table_districts.R -- per-country summary of the calibrated district estimates of the
## age-collapsed admin-2 fit (conditional province anchors, EBLUP): median and range of
## positivity, median 95% interval width, the calibration share of the prediction variance,
## and the width with the transport variance sigma_tau^2 added.
## Reads:  output/<country>/fits/{district_predictions,province_calibration_conditional}.rds,
##         output/canonical/estimates/transport_diagnostic.csv
## Writes: output/canonical/estimates/table_districts.csv
## Usage:  Rscript R/regen_table_districts.R
suppressMessages(library(dplyr))
source("R/utils.R")   # harmonise_register_units()
cts <- c("Mozambique","Zimbabwe","Malawi","Lesotho")

## sigma_tau^2 from the transport diagnostic, added to the prediction variance (negative estimates
## set to zero; Lesotho has no diagnostic, one area per province, so its column is NA)
td <- read.csv("output/canonical/estimates/transport_diagnostic.csv")
td <- td[!(td$country == "Mozambique" & !td$boundary_dropped), ]
s2tau <- setNames(pmax(0, td$sigma_tau2), td$country)

out <- bind_rows(lapply(cts, function(c){
  dp <- harmonise_register_units(as.data.frame(readRDS(file.path("output", tolower(c), "fits", "district_predictions.rds"))), tolower(c)) |>
          select(province, district, z_d, v_d)
  pc <- as.data.frame(readRDS(file.path("output", tolower(c), "fits", "province_calibration_conditional.rds"))) |>
          select(province, lambda_p, var_lambda_p)
  d <- dp |> inner_join(pc, by="province") |>
    filter(is.finite(z_d), is.finite(v_d), v_d > 0, is.finite(var_lambda_p), var_lambda_p > 0) |>
    mutate(lgt   = z_d - lambda_p,
           vtot  = v_d + var_lambda_p,
           se    = sqrt(vtot),
           theta = plogis(lgt),
           lo    = plogis(lgt - 1.96*se),
           hi    = plogis(lgt + 1.96*se),
           ciw   = (hi - lo) * 100,
           pcal  = 100 * var_lambda_p / vtot,
           st2   = if (c %in% names(s2tau)) s2tau[[c]] else NA_real_,
           ciw_tau = (plogis(lgt + 1.96*sqrt(vtot + st2)) - plogis(lgt - 1.96*sqrt(vtot + st2))) * 100)
  data.frame(country = c,
             areas   = nrow(d),
             med     = 100*median(d$theta),
             lo      = 100*min(d$theta),
             hi      = 100*max(d$theta),
             ciw     = median(d$ciw),
             pcal    = median(d$pcal),
             sigma_tau2 = d$st2[1],
             ciw_tau = median(d$ciw_tau))
}))

write.csv(out, "output/canonical/estimates/table_districts.csv", row.names = FALSE)

cat("\n=== table_districts regenerated, conditional anchors (EBLUP) ===\n")
for(i in seq_len(nrow(out))) with(out[i,],
  cat(sprintf("%-11s & %3d & %.1f\\%% & %.1f--%.1f\\%% & %.1f pp & %.0f\\%% \\\\\n",
              country, areas, med, lo, hi, ciw, pcal)))
cat(sprintf("\ncalibration-share range across countries: %.0f--%.0f%%\n",
            min(out$pcal), max(out$pcal)))

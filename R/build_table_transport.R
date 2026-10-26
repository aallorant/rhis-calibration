## =====================================================================
## build_table_transport.R -- per-country table of the transport diagnostic: the test
## statistics and the district accounting (districts surveyed, and those dropped for no
## register record, no undiagnosed case in the sample, a single sampled cluster, or being
## the only usable district of their province). Mozambique is the boundary-dropped run on
## the design-based inputs of R/build_moz_district_survey_design.R. Run after R/transport_diagnostic.R.
## Reads:  output/canonical/estimates/transport_diagnostic{,_by_province}.csv, output/<country>/fits/*.rds
## Writes: output/canonical/estimates/table_transport.csv.  Usage: Rscript R/build_table_transport.R
## =====================================================================
suppressMessages(library(dplyr))
EST <- "output/canonical/estimates"
nat <- read.csv(file.path(EST, "transport_diagnostic.csv"))
nat <- nat[!(nat$country == "Mozambique" & !nat$boundary_dropped), ]
per <- read.csv(file.path(EST, "transport_diagnostic_by_province.csv"))
per <- per[!(per$country == "Mozambique" & !per$boundary_dropped), ]

acct <- function(cty) {
  if (cty == "Mozambique") {
    ph <- readRDS("output/mozambique/fits/phia_district_estimates_design.rds")
    pr <- table(ph$province[ph$informative])
    return(data.frame(provinces_total = length(unique(ph$province)), districts_surveyed = nrow(ph),
      drop_no_register = 0L, drop_zero_count = sum(ph$p_hat_raw_is_boundary),
      drop_single_cluster_only = sum(ph$single_psu & !ph$p_hat_raw_is_boundary),
      drop_singleton_province = sum(pr == 1)))
  }
  f  <- tolower(cty); pc <- if (cty == "Malawi") "zone" else "province"
  e  <- as.data.frame(readRDS(sprintf("output/%s/fits/phia_district_estimates.rds", f)))
  d  <- as.data.frame(readRDS(sprintf("output/%s/fits/district_estimation_dataset.rds", f)))
  u  <- d[d$has_both, ]; pr <- table(u[[pc]])
  data.frame(provinces_total = length(unique(d[[pc]])), districts_surveyed = nrow(e),
    drop_no_register = sum(d$has_phia & !d$has_rhis), drop_zero_count = sum(e$is_boundary),
    drop_single_cluster_only = 0L, drop_singleton_province = sum(pr == 1),
    singleton_provinces = paste(names(pr)[pr == 1], collapse = "; "))
}
tab <- bind_rows(lapply(c("Mozambique", "Zimbabwe", "Malawi"), function(cty) {
  n <- nat[nat$country == cty, ]; p <- per[per$country == cty, ]
  cbind(country = cty, acct(cty),
        provinces_tested = n$n_prov, districts_tested = n$n_dist,
        t = n$t, z = n$zstat, p_one_sided = n$p_one_sided,
        sigma_tau2 = n$sigma_tau2, sigma_tau2_lo = n$sigma_tau2_lo, sigma_tau2_hi = n$sigma_tau2_hi,
        sigma_u2_collapsed = n$sigma_u2,
        province_tests = nrow(p), province_rejections_05 = sum(p$p_one_sided < 0.05),
        provinces_rejecting = paste(p$prov[p$p_one_sided < 0.05], collapse = "; "))
}))
## the accounting must close
stopifnot(with(tab, districts_surveyed - drop_no_register - drop_zero_count -
                    drop_single_cluster_only - drop_singleton_province == districts_tested))
print(tab, row.names = FALSE)
write.csv(tab, file.path(EST, "table_transport.csv"), row.names = FALSE)
cat("written:", file.path(EST, "table_transport.csv"), "\n")

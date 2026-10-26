## =====================================================================
## build_table_weights.R -- shrinkage weights of the two fits, per country:
##   gamma = sigma_u^2 / (sigma_u^2 + v), median over units;
##   separation = median absolute gap, on the logit scale, between the shrinkage predictor and
##   the calibrated register, i.e. |(1 - gamma) (survey estimate - calibrated register)|.
## Both resolutions are rebuilt from the fit inputs (age-sex as in R/mainfit_inla_crosscheck.R,
## incl. the Malawi sex-only collapse); Lesotho's age-sex fit is saturated, so its weights are NA.
## Reads:  output/<country>/fits/*.rds, output/canonical/estimates/{collapsed_admin2_summary,mainfit_province_stratum_lambda,mainfit_sigma_u}.csv
## Writes: output/canonical/estimates/table_weights.csv.  Usage: Rscript R/build_table_weights.R
## =====================================================================
suppressMessages({library(dplyr)})
source("R/utils.R")   # harmonise_register_units()
CTS <- c("Mozambique", "Zimbabwe", "Malawi", "Lesotho")
EST <- "output/canonical/estimates"

weights_of <- function(g, Z_hat, cal) {
  s <- abs((1 - g) * (Z_hat - cal))
  c(gamma = median(g), sep = median(s), gamma_min = min(g),
    gamma_ratio = max(g) / min(g), sep_p90 = unname(quantile(s, 0.9)))
}

## ---- admin-2, the age-collapsed fit --------------------------------
a2 <- read.csv(file.path(EST, "collapsed_admin2_summary.csv"))
adm2 <- bind_rows(lapply(CTS, function(c) {
  s  <- tolower(c)
  d  <- harmonise_register_units(as.data.frame(readRDS(file.path("output", s, "fits", "district_predictions.rds"))), s) |>
    select(province, district, z_d, v_d)
  pc <- as.data.frame(readRDS(file.path("output", s, "fits", "province_calibration_conditional.rds"))) |>
    select(province, lambda_p, Z_hat_p)
  d  <- d |> inner_join(pc, by = "province") |>
    filter(is.finite(z_d), is.finite(v_d), v_d > 0, is.finite(Z_hat_p))
  s2 <- a2$sigma_u2_eblup[a2$country == c]
  w  <- weights_of(s2 / (s2 + d$v_d), d$Z_hat_p, d$z_d - d$lambda_p)
  data.frame(country = c, gamma_admin2 = w[["gamma"]], sep_admin2 = w[["sep"]])
}))

## ---- age-sex, rebuilding the main fit's cells ------------------------
lam <- read.csv(file.path(EST, "mainfit_province_stratum_lambda.csv"))
load_mainfit <- function(cty) {
  s  <- tolower(cty)
  rh <- harmonise_register_units(as.data.frame(readRDS(file.path("output", s, "fits", "district_predictions_demo.rds"))), s) |>
    filter(is.finite(z), is.finite(v), v > 0) |>
    mutate(stratum = trimws(paste(age_group_harmonised, sex)))
  sv <- as.data.frame(readRDS(file.path("output", s, "fits", "province_calibration_demo_conditional.rds"))) |>
    filter(is.finite(Z_hat), is.finite(psi), psi > 0) |>
    mutate(stratum = trimws(paste(age_group_harmonised, sex)))
  if (all(grepl("^all ", unique(rh$stratum)))) {          # Malawi: sex only
    sv <- sv |> mutate(sex2 = ifelse(grepl("Female", stratum), "Female", "Male")) |>
      group_by(province, sex2) |>
      summarise(Z_hat = sum(Z_hat / psi) / sum(1 / psi), psi = 1 / sum(1 / psi), .groups = "drop") |>
      mutate(stratum = paste("all", sex2)) |> select(province, stratum, Z_hat, psi)
  }
  rh |> inner_join(select(sv, province, stratum, Z_hat, psi), by = c("province", "stratum"))
}
su <- read.csv(file.path(EST, "mainfit_sigma_u.csv"))
ages <- bind_rows(lapply(CTS, function(c) {
  d  <- load_mainfit(c)
  s2 <- su$sigma_u2_baseline[su$country == c]
  ## saturated where every province x stratum cell holds one area
  n_per_cell <- d |> count(province, stratum) |> pull(n)
  if (s2 <= 0 || all(n_per_cell <= 1))
    return(data.frame(country = c, gamma_agesex = NA_real_, sep_agesex = NA_real_,
                      gamma_min_agesex = NA_real_, gamma_ratio_agesex = NA_real_,
                      sep_p90_agesex = NA_real_))
  lam_c <- lam |> filter(country == c) |> select(province, stratum, log_lambda)
  d <- d |> left_join(lam_c, by = c("province", "stratum")) |>
    filter(is.finite(log_lambda))
  w <- weights_of(s2 / (s2 + d$v), d$Z_hat, d$z - d$log_lambda)
  data.frame(country = c, gamma_agesex = w[["gamma"]], sep_agesex = w[["sep"]],
             gamma_min_agesex = w[["gamma_min"]], gamma_ratio_agesex = w[["gamma_ratio"]],
             sep_p90_agesex = w[["sep_p90"]])
}))

out <- adm2 |> left_join(ages, by = "country")
write.csv(out, file.path(EST, "table_weights.csv"), row.names = FALSE)

cat("\n=== table_weights regenerated ===\n\n")
print(as.data.frame(out |> transmute(country,
  `admin-2 gamma` = round(gamma_admin2, 3), `admin-2 sep` = round(sep_admin2, 3),
  `age-sex gamma` = round(gamma_agesex, 2), `age-sex sep` = round(sep_agesex, 3),
  `min gamma` = round(gamma_min_agesex, 3), `gamma ratio` = round(gamma_ratio_agesex, 1),
  `sep p90` = round(sep_p90_agesex, 3))), row.names = FALSE)
cat("\nwritten: ", file.path(EST, "table_weights.csv"), "\n", sep = "")

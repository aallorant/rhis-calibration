## =====================================================================
## build_table_predvar_h.R -- 95% interval width of a calibrated register estimate h years
## after the last survey, at each country's median district (percentage points, interval
## built on the logit scale):  Var_h = v_d + Var(lambda_hat_p) + sigma_tau^2 + h sigma_Delta^2 + h^2 mu^2.
## sigma_tau^2: transport diagnostic (negative set to zero; Lesotho has none and is omitted);
## sigma_Delta^2 = sigma2_lambda and mu = sqrt(drift2_lambda) per year, pooled 3-country row of volatility_ratio.csv.
## Reads:  output/<country>/fits/*.rds, output/canonical/estimates/transport_diagnostic.csv, output/volatility_ratio.csv
## Writes: output/canonical/estimates/table_predvar_h.csv.  Usage: Rscript R/build_table_predvar_h.R
## =====================================================================
suppressMessages(library(dplyr))
source("R/utils.R")
cts <- c("Mozambique", "Zimbabwe", "Malawi")
H   <- c(0, 1, 3, 5)
vr  <- read.csv("output/volatility_ratio.csv"); pr <- vr[vr$group == "Pooled (3 countries)", ]
s2D <- pr$sigma2_lambda; mu <- sqrt(pr$drift2_lambda)
td  <- read.csv("output/canonical/estimates/transport_diagnostic.csv")
td  <- td[!(td$country == "Mozambique" & !td$boundary_dropped), ]
wid <- function(l, v) 100 * (plogis(l + 1.96 * sqrt(v)) - plogis(l - 1.96 * sqrt(v)))
out <- bind_rows(lapply(cts, function(c) {
  dp <- harmonise_register_units(as.data.frame(readRDS(file.path("output", tolower(c), "fits", "district_predictions.rds"))), c) |>
          select(province, district, z_d, v_d)
  pc <- as.data.frame(readRDS(file.path("output", tolower(c), "fits", "province_calibration_conditional.rds"))) |>
          select(province, lambda_p, var_lambda_p)
  d <- dp |> inner_join(pc, by = "province") |>
    filter(is.finite(z_d), is.finite(v_d), v_d > 0, var_lambda_p > 0) |>
    mutate(lgt = z_d - lambda_p, th = plogis(lgt))
  x <- d[which.min(abs(d$th - median(d$th))), ]
  st2 <- max(0, td$sigma_tau2[td$country == c])
  bind_rows(lapply(H, function(h) data.frame(
    country = c, h = h, theta = 100 * x$th, v_d = x$v_d, var_lambda = x$var_lambda_p,
    sigma_tau2 = st2, sigma_Delta2 = s2D, mu = mu,
    w_current = wid(x$lgt, x$v_d + x$var_lambda_p),
    w_tau     = wid(x$lgt, x$v_d + x$var_lambda_p + st2),
    w_diff    = wid(x$lgt, x$v_d + x$var_lambda_p + st2 + h * s2D),
    w_drift   = wid(x$lgt, x$v_d + x$var_lambda_p + st2 + h * s2D + (h * mu)^2))))
}))
write.csv(out, "output/canonical/estimates/table_predvar_h.csv", row.names = FALSE)
cat(sprintf("sigma_Delta^2 = %.4f / yr, |mu| = %.3f / yr\n", s2D, mu))
print(out |> mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)

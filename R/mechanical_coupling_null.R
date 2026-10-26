## =====================================================================
## mechanical_coupling_null.R -- the mechanical value of the correlation between the change in
## the offset lambda = z - Zhat and the change in prevalence: the part that shared survey noise
## alone produces. With u = P(HIV+ and unaware) and p = P(HIV+) from one sample of effective size n,
## Cov(logit u_hat, p_hat) = (1 - p) / (n (1 - u)), n = 1 / (se_lambda^2 u (1 - u)); rounds are
## independent, so the noise covariance of (d_lambda, d_prev) is -(c_round0 + c_round1), and the
## mechanical correlation is its mean over areas divided by sd(d_lambda) sd(d_prev).
## Reads:  output/temporal_lambda_paired.csv (R/04d_temporal_lambda.R)
## Writes: output/mechanical_coupling_null.csv.  Usage: Rscript R/mechanical_coupling_null.R
## =====================================================================
p <- read.csv("output/temporal_lambda_paired.csv")
one <- function(d) {
  u0 <- plogis(d$Z_hat_round0); u1 <- plogis(d$Z_hat_round1)
  p0 <- d$prev_round0;          p1 <- d$prev_round1
  n0 <- 1 / (d$se_lambda_round0^2 * u0 * (1 - u0))
  n1 <- 1 / (d$se_lambda_round1^2 * u1 * (1 - u1))
  cov_noise <- -((1 - p0) / (n0 * (1 - u0)) + (1 - p1) / (n1 * (1 - u1)))
  data.frame(country = d$country[1], n_areas = nrow(d),
             r_observed = cor(d$d_lambda, d$d_prev),
             r_mechanical = mean(cov_noise) / (sd(d$d_lambda) * sd(d$d_prev)))
}
res <- do.call(rbind, lapply(split(p, p$country), one))
print(res, row.names = FALSE, digits = 3)
write.csv(res, "output/mechanical_coupling_null.csv", row.names = FALSE)
cat("written: output/mechanical_coupling_null.csv\n")

## refit_block_eblup_canonical.R -- refit the EBLUP rows of the canonical estimates under
## the block covariance V_g = sigma_u^2 I + diag(v_i) + psi_g J. With A_g = diag(sigma_u^2 + v_i),
## S_g = 1'A^-1 1 and T_g = 1'A^-1 y: mean T_g/S_g, Var = 1/S_g + psi_g, logdet V_g = logdet A_g + log(1 + psi_g S_g).
## Reads:  output/<country>/fits/*.rds, output/canonical/estimates/collapsed_admin2_summary.csv
## Rewrites the estimator == "EBLUP" rows of output/canonical/estimates/collapsed_admin2_province_lambda.csv
##   and mainfit_national_by_stratum.csv (TMB rows untouched) and prints old/new/TMB for validation.
## Usage:  Rscript R/refit_block_eblup_canonical.R
suppressMessages({library(dplyr)})
source("R/utils.R")   # harmonise_register_units()
suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
cts <- c("Mozambique","Zimbabwe","Malawi","Lesotho")

## profile ML for sigma_u^2 under the block covariance; returns per-group mean and SE
fh_block <- function(y, v, psi_g, grp) {
  grp <- factor(grp); gl <- levels(grp)
  psi <- psi_g[gl]                                   # one psi per group
  nll <- function(ls2) {
    s2 <- exp(ls2); A <- s2 + v
    S  <- tapply(1/A, grp, sum); T <- tapply(y/A, grp, sum)
    m  <- T/S; r <- y - m[grp]
    0.5*( sum(log(2*pi*A)) + sum(log(1 + psi[gl]*S[gl])) + sum(r^2/A) )
  }
  o  <- optimize(nll, c(-12, 4)); s2 <- exp(o$minimum); A <- s2 + v
  S  <- tapply(1/A, grp, sum); T <- tapply(y/A, grp, sum)
  m  <- T/S
  se <- sqrt(1/S[gl] + psi[gl])                      # block: psi enters at full strength
  list(sigma_u2 = s2, mean = m[gl], se = as.numeric(se), groups = gl)
}
ivw <- function(est, se){ w <- 1/se^2; c(m = sum(est*w)/sum(w), se = sqrt(1/sum(w))) }

## ---------------------------------------------------------------- marginal (admin-2)
prov_new <- bind_rows(lapply(cts, function(c){
  dp <- harmonise_register_units(as.data.frame(readRDS(file.path("output", tolower(c), "fits", "district_predictions.rds"))), tolower(c))
  pc <- as.data.frame(readRDS(file.path("output",tolower(c),"fits","province_calibration_conditional.rds")))
  d  <- dp |> select(province, z_d, v_d) |>
        inner_join(pc |> select(province, Z_hat_p, psi_p), by="province") |>
        filter(is.finite(z_d), is.finite(v_d), v_d>0, is.finite(psi_p), psi_p>0)
  psi <- tapply(d$psi_p, d$province, function(x) x[1])
  f   <- fh_block(d$z_d - d$Z_hat_p, d$v_d, psi, d$province)
  n   <- ivw(f$mean, f$se)
  data.frame(country=c, province=f$groups, log_lambda=f$mean, se=f$se,
             nat_m=n["m"], nat_se=n["se"], sigma_u2=f$sigma_u2, row.names=NULL)
}))

## ---------------------------------------------------------------- age x sex x province
strat_new <- bind_rows(lapply(cts, function(c){
  rh <- harmonise_register_units(as.data.frame(readRDS(file.path("output", tolower(c), "fits", "district_predictions_demo.rds"))), tolower(c)) |>
        filter(is.finite(z), is.finite(v), v>0) |>
        mutate(stratum=trimws(paste(age_group_harmonised, sex)), unit=paste(province, stratum))
  sv <- as.data.frame(readRDS(file.path("output",tolower(c),"fits","province_calibration_demo_conditional.rds"))) |>
        mutate(stratum=trimws(paste(age_group_harmonised, sex)), unit=paste(province, stratum)) |>
        filter(is.finite(Z_hat), psi>0) |> distinct(unit, .keep_all=TRUE)
  if (all(grepl("^all ", unique(rh$stratum)))) {          # Malawi: RHIS sex-only
    sv <- sv |> mutate(sex2=ifelse(grepl("Female", stratum), "Female","Male")) |>
          group_by(province, sex2) |>
          summarise(Z_hat=sum(Z_hat/psi)/sum(1/psi), psi=1/sum(1/psi), .groups="drop") |>
          mutate(stratum=paste("all", sex2), unit=paste(province, stratum)) |> select(unit, Z_hat, psi)
  }
  d <- rh |> inner_join(sv |> select(unit, Z_hat, psi), by="unit")
  psi <- tapply(d$psi, d$unit, function(x) x[1])
  f <- fh_block(d$z - d$Z_hat, d$v, psi, d$unit)
  st <- sub("^[^ ]+ ", "", f$groups); st[!grepl(" ", f$groups)] <- f$groups[!grepl(" ", f$groups)]
  # stratum label = everything after the province name; rebuild from the unit strings
  lut <- d |> distinct(unit, stratum); st <- lut$stratum[match(f$groups, lut$unit)]
  data.frame(country=c, stratum=st, log_lambda=f$mean, se=f$se, row.names=NULL) |>
    group_by(country, stratum) |>
    summarise(w=sum(1/se^2), log_lambda=sum(log_lambda/se^2)/w, se=sqrt(1/w), .groups="drop") |>
    select(-w)
}))

## ---------------------------------------------------------------- validate + write
pl <- read.csv("output/canonical/estimates/collapsed_admin2_province_lambda.csv", stringsAsFactors=FALSE)
sn <- read.csv("output/canonical/estimates/mainfit_national_by_stratum.csv", stringsAsFactors=FALSE)

cat("\n=== national selection OR, admin-2 marginal fit: EBLUP old (indep) vs new (block) vs TMB ===\n")
ms <- read.csv("output/canonical/estimates/collapsed_admin2_summary.csv", stringsAsFactors=FALSE)
se_ci <- function(lo,hi)(log(hi)-log(lo))/(2*1.96)
for (c in cts) {
  nn <- prov_new |> filter(country==c) |> slice(1)
  cat(sprintf("  %-11s OR %.3f (%.3f, %.3f) | se old %.4f  new %.4f  TMB %.4f  new/old %.2f  new/TMB %.2f\n",
      c, exp(nn$nat_m), exp(nn$nat_m-1.96*nn$nat_se), exp(nn$nat_m+1.96*nn$nat_se),
      se_ci(ms$OR_eblup_lo[ms$country==c], ms$OR_eblup_hi[ms$country==c]), nn$nat_se,
      se_ci(ms$OR_tmb_lo[ms$country==c], ms$OR_tmb_hi[ms$country==c]),
      nn$nat_se/se_ci(ms$OR_eblup_lo[ms$country==c], ms$OR_eblup_hi[ms$country==c]),
      nn$nat_se/se_ci(ms$OR_tmb_lo[ms$country==c], ms$OR_tmb_hi[ms$country==c])))
  cat(sprintf("              sigma_u^2 new %.4f (canonical EBLUP %.4f)\n", nn$sigma_u2, ms$sigma_u2_eblup[ms$country==c]))
}

## province_lambda.csv: replace EBLUP rows
pl_new <- pl |> filter(estimator!="EBLUP") |>
  bind_rows(prov_new |> transmute(country, province, estimator="EBLUP",
              log_lambda=round(log_lambda,3), se=round(se,3), OR=round(exp(log_lambda),3),
              lo=round(exp(log_lambda-1.96*se),3), hi=round(exp(log_lambda+1.96*se),3))) |>
  arrange(country, estimator, province)
write.csv(pl_new, "output/canonical/estimates/collapsed_admin2_province_lambda.csv", row.names=FALSE)

## mainfit_national_by_stratum.csv: replace EBLUP rows, keep the age_group/sex columns
key <- sn |> filter(estimator=="TMB") |> select(country, stratum, age_group, sex)
sn_new <- sn |> filter(estimator!="EBLUP") |>
  bind_rows(strat_new |> left_join(key, by=c("country","stratum")) |>
              transmute(country, stratum, estimator="EBLUP",
                        log_lambda, se, OR=round(exp(log_lambda),3), age_group, sex,
                        lo=round(exp(log_lambda-1.96*se),3), hi=round(exp(log_lambda+1.96*se),3))) |>
  arrange(country, estimator, stratum)
write.csv(sn_new, "output/canonical/estimates/mainfit_national_by_stratum.csv", row.names=FALSE)

cat("\n=== age-sex cells: 95% interval excludes 1 ===\n")
for (e in c("EBLUP","TMB")) {
  s <- sn_new |> filter(estimator==e)
  cat(sprintf("  %-6s %d / %d cells\n", e, sum(s$lo>1 | s$hi<1), nrow(s)))
}
cat("\nwrote province_lambda.csv and mainfit_national_by_stratum.csv (EBLUP rows = block covariance)\n")

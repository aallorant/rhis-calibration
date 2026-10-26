## build_conditional_anchors.R: conditional province survey anchors,
## Z_hat_p = P(HIV+ | not previously diagnosed), estimated on the survey subpopulation
## not known-positive (svydesign on the full data, subset, svyby HIV+ by province).
## Reads:  data/processed/phia/phia_<country>.rds, output/<country>/fits/province_calibration.rds
## Writes: output/<country>/fits/province_calibration_conditional.rds (with conditional lambda_p);
##         prints the joint vs conditional national selection OR.
## Run:    Rscript R/build_conditional_anchors.R   (from the repository root)
suppressMessages({library(survey); library(dplyr)})
options(survey.lonely.psu="adjust")
logit <- function(p) log(p/(1-p)); se_logit <- function(p,se) se/(p*(1-p))
cts <- c("Mozambique","Zimbabwe","Malawi","Lesotho")

report <- list()
for (cty in cts) {
  ph <- as.data.frame(readRDS(file.path("data/processed/phia", paste0("phia_", tolower(cty), ".rds")))) %>%
    filter(age>=15, age<=64)
  ph$province <- trimws(ph$province)
  des <- svydesign(id=~varunit, strata=~varstrat, weights=~btwt0, data=ph, nest=TRUE)
  # conditional: subpopulation = NOT known-positive
  sub <- subset(des, !(hivstatus==1 & aware_bin==TRUE))
  zc <- svyby(~hivstatus, ~province, sub, svymean, na.rm=TRUE) %>%
    transmute(province, Zc_prob=hivstatus, se_Zc=se)
  # joint anchor + RHIS province values from province_calibration.rds
  pc <- as.data.frame(readRDS(file.path("output", tolower(cty), "fits", "province_calibration.rds")))
  m <- pc %>% select(province, z_p, v_p, theta_N, kappa, Z_hat_prob, Z_hat_p_joint=Z_hat_p, psi_joint=psi_p) %>%
    inner_join(zc, by="province") %>%
    mutate(Z_hat_p   = logit(Zc_prob),                       # conditional anchor (logit)
           psi_p     = se_logit(Zc_prob, se_Zc)^2,
           lambda_p  = z_p - Z_hat_p,                         # conditional selection log-OR
           var_lambda_p = v_p + psi_p,
           se_lambda_p  = sqrt(var_lambda_p),
           OR_p      = exp(lambda_p))
  saveRDS(m, file.path("output", tolower(cty), "fits", "province_calibration_conditional.rds"))
  # national grand means, joint vs conditional, IV-weighted
  oldll <- pc$lambda_p; oldse <- pc$se_lambda_p
  ivw <- function(ll,se){w<-1/se^2; exp(sum(ll*w)/sum(w))}
  report[[cty]] <- data.frame(country=cty, n_prov=nrow(m),
    OR_joint=round(ivw(oldll,oldse),2), OR_cond=round(ivw(m$lambda_p,m$se_lambda_p),2),
    mean_Zhat_joint=round(mean(plogis(pc$Z_hat_p)),4), mean_Zhat_cond=round(mean(m$Zc_prob),4),
    min_OR_cond=round(min(m$OR_p),2))
}
res <- bind_rows(report)
cat("=== survey anchor: joint vs conditional, national selection OR ===\n"); print(as.data.frame(res), row.names=FALSE)
cat("\nselection positive everywhere (min province OR_cond > 1)?", all(res$min_OR_cond>1), "\n")

## build_conditional_stratum_anchors.R: conditional province x (age x sex) survey anchors,
## the subpopulation svymean of HIV+ over respondents not known-positive.
## Reads:  data/processed/phia/phia_<country>.rds
## Writes: output/<country>/fits/province_calibration_demo_conditional.rds
##         (the province_calibration_demo columns with conditional Z_hat and psi)
## Run:    Rscript R/build_conditional_stratum_anchors.R   (from the repository root)
suppressMessages({library(survey); library(dplyr)})
options(survey.lonely.psu="adjust")
logit <- function(p) log(p/(1-p)); se_logit <- function(p,se) se/(p*(1-p))
agg <- "age_group_harmonised"
## Malawi's routine data carry sex but not age, so its anchor is estimated at province x sex.
SEX_ONLY <- c("Malawi")
for (cty in c("Mozambique","Zimbabwe","Malawi","Lesotho")) {
  ph <- as.data.frame(readRDS(file.path("data/processed/phia", paste0("phia_", tolower(cty), ".rds")))) %>%
    filter(age>=15, age<=64)
  ph$province <- trimws(ph$province)
  # harmonised age group + sex; match the demo file's strata
  ph$age_group_harmonised <- ph$age_group
  des <- svydesign(id=~varunit, strata=~varstrat, weights=~btwt0, data=ph, nest=TRUE)
  sub <- subset(des, !(hivstatus==1 & aware_bin==TRUE))
  sex_only <- cty %in% SEX_ONLY
  forml <- if (sex_only) as.formula("~province + sex") else as.formula(paste("~province +", agg, "+ sex"))
  zc <- tryCatch(svyby(~hivstatus, forml, sub, svymean, na.rm=TRUE), error=function(e) NULL)
  if (is.null(zc)) { cat(cty, "stratum svyby failed; skipping\n"); next }
  zc <- as.data.frame(zc)
  if (sex_only) zc[[agg]] <- "all"          # register carries "all" as the age label
  zc <- zc %>%
    transmute(province, age_group_harmonised = .data[[agg]], sex,
              Zc_prob = hivstatus, se_Zc = se) %>%
    filter(is.finite(Zc_prob), Zc_prob>0, Zc_prob<1, is.finite(se_Zc), se_Zc>0) %>%
    mutate(Z_hat = logit(Zc_prob), psi = se_logit(Zc_prob, se_Zc)^2)
  saveRDS(zc, file.path("output", tolower(cty), "fits", "province_calibration_demo_conditional.rds"))
  cat(sprintf("%-11s conditional stratum anchors: %d cells\n", cty, nrow(zc)))
}

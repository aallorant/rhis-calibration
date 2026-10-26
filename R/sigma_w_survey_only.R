## =====================================================================
## sigma_w_survey_only.R -- survey-only estimate of sigma_w^2, the within-province variance of
## the true undiagnosed positivity, against the deployment threshold sigma_u^2/2.
## Binomial mixed model on the not-yet-diagnosed subset:  y_d ~ Bin(n_d^eff, theta_d),
## logit theta_d = mu_{p(d)} + w_d, w_d ~ N(0, sigma_w^2), province fixed effects, Kish effective
## sizes n_eff = (sum w)^2 / sum w^2; EA nested in district for Zimbabwe and Malawi (tagged microdata),
## two-level from district effective counts for Mozambique; Lesotho is skipped (no within-province data).
## Reads:  output/<country>/fits/{phia_tagged_admin2,phia_district_estimates}.rds
## Writes: output/canonical/estimates/mainfit_sigma_w_survey.csv.  Usage: Rscript R/sigma_w_survey_only.R
## =====================================================================

suppressMessages({library(lme4); library(dplyr)})

THRESH <- c(Mozambique = 0.2019, Zimbabwe = 0.2351, Malawi = 0.1116) / 2

kish <- function(w) sum(w)^2 / sum(w^2)

## ---- Zimbabwe and Malawi: three-level from tagged microdata ---------------
fit_micro <- function(cty) {
  d <- as.data.frame(readRDS(file.path("output", tolower(cty), "fits", "phia_tagged_admin2.rds")))
  d$w <- if (any(is.finite(d$btwt0))) d$btwt0 else d$indweight
  dist_col <- if ("shp_district" %in% names(d)) "shp_district" else "district"
  d$dist <- d[[dist_col]]
  d$prov <- if ("shp_province" %in% names(d)) d$shp_province else d$province
  d <- d |> filter(is.finite(w), w > 0, !is.na(hivstatus), !is.na(dist), !is.na(prov)) |>
    filter(hivstatus == 0 | (hivstatus == 1 & (is.na(aware_bin) | aware_bin == 0))) |>
    mutate(ea = paste(varstrat, varunit, sep = "_"))

  ## aggregate to EA, with Kish effective sizes
  a <- d |> group_by(prov, dist, ea) |>
    summarise(n_eff = kish(w), p = sum(w * (hivstatus == 1))/sum(w), .groups = "drop") |>
    filter(is.finite(n_eff), n_eff >= 1) |>
    mutate(n = round(n_eff), y = round(n_eff * p)) |>
    filter(n > 0, y >= 0, y <= n)
  list(dat = a, three = TRUE)
}

## ---- Mozambique: two-level from district-level effective counts -----------
fit_agg <- function(cty) {
  ph <- as.data.frame(readRDS(file.path("output", tolower(cty), "fits", "phia_district_estimates.rds")))
  a <- ph |> filter(is.finite(n_eff), n_eff >= 1, is.finite(p_hat_raw)) |>
    transmute(prov = province, dist = area_name, ea = area_name,
              n = round(n_eff), y = round(n_eff * p_hat_raw)) |>
    filter(n > 0, y >= 0, y <= n)
  list(dat = a, three = FALSE)
}

rows <- list()
for (cty in c("Mozambique", "Zimbabwe", "Malawi")) {
  o <- if (cty == "Mozambique") fit_agg(cty) else fit_micro(cty)
  a <- o$dat
  form <- if (o$three) cbind(y, n - y) ~ prov + (1|dist) + (1|ea)
          else         cbind(y, n - y) ~ prov + (1|dist)
  m <- tryCatch(glmer(form, data = a, family = binomial,
                      control = glmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 3e5))),
                error = function(e) {message(cty, ": ", conditionMessage(e)); NULL})
  if (is.null(m)) next

  vc <- as.data.frame(VarCorr(m))
  s2 <- vc$vcov[vc$grp == "dist"]
  ## Wald interval on the log-variance scale via a profile of the SD if it converges
  ci <- tryCatch({
    p <- confint(m, parm = "theta_", method = "profile", oldNames = FALSE, quiet = TRUE)
    r <- p[grep("dist", rownames(p)), , drop = FALSE]
    c(r[1,1]^2, r[1,2]^2)
  }, error = function(e) c(NA, NA))

  rows[[cty]] <- data.frame(country = cty, units = nrow(a),
    districts = length(unique(a$dist)), three_level = o$three,
    sigma_w2 = s2, lo = ci[1], hi = ci[2],
    threshold = THRESH[[cty]],
    verdict = ifelse(is.na(ci[2]), "no interval",
              ifelse(ci[2] < THRESH[[cty]], "RULED OUT",
              ifelse(ci[1] > THRESH[[cty]], "criterion met", "inconclusive"))))
}

out <- bind_rows(rows)
cat("\n=== survey-only sigma_w^2 against the deployment threshold sigma_u^2/2 ===\n\n")
print(as.data.frame(out |> transmute(country, units, districts,
  levels = ifelse(three_level, "3", "2"),
  `sigma_w^2` = round(sigma_w2, 4),
  `95% CI` = ifelse(is.na(lo), "-", sprintf("(%.3f, %.3f)", lo, hi)),
  threshold = round(threshold, 4), verdict)), row.names = FALSE)
cat("\n  RULED OUT needs the UPPER limit below the threshold.\n")
write.csv(out, "output/canonical/estimates/mainfit_sigma_w_survey.csv", row.names = FALSE)

## =====================================================================
## Design-based district survey estimates for Mozambique (INSIDA 2021) on the
## method used for Zimbabwe and Malawi: Taylor-linearised svydesign on the
## survey strata, PSUs and biomarker weights, domain mean of HIV+ among adults
## not known positive, logit by the delta method; district flags defined below.
## Reads:  output/mozambique/fits/phia_insida2021_with_awareness.rds and
##         output/mozambique/fits/phia_district_estimates.rds (comparison values)
## Writes: output/mozambique/fits/phia_district_estimates_design.{rds,csv}
## Run:    Rscript R/build_moz_district_survey_design.R  (from the repository root)
## =====================================================================
suppressMessages({library(survey); library(dplyr)})
options(survey.lonely.psu = "adjust")

ph <- as.data.frame(readRDS("output/mozambique/fits/phia_insida2021_with_awareness.rds"))
ph <- ph[!is.na(ph$hivstatus) & !is.na(ph$area_id) & !is.na(ph$btwt0) & ph$btwt0 > 0, ]
ph$hivpos <- as.numeric(ph$hivstatus == 1)
ph$known  <- ph$hivpos == 1 & ph$aware_bin %in% TRUE
ph$psu    <- paste(ph$varstrat, ph$varunit)
des <- svydesign(id = ~varunit, strata = ~varstrat, weights = ~btwt0, data = ph, nest = TRUE)

safe_qlogis <- function(p, n) { eps <- pmin(0.5 / pmax(n, 1), 1e-3); qlogis(pmin(pmax(p, eps), 1 - eps)) }

rows <- lapply(sort(unique(ph$area_id)), function(a) {
  k   <- ph[ph$area_id == a & !ph$known, ]
  sub <- des[des$variables$area_id == a & !des$variables$known, ]
  m   <- svymean(~hivpos, sub, na.rm = TRUE)
  p <- as.numeric(coef(m)); se <- as.numeric(SE(m))
  data.frame(area_id = a, area_name = k$area_name[1], province = k$province_name[1],
             n = nrow(k), n_psu = length(unique(k$psu)), n_pos = sum(k$hivpos),
             w_total = sum(k$btwt0), w_pos = sum(k$btwt0 * k$hivpos),
             p_hat_raw = p, se_p = se)
})
## single_psu: all sampled respondents in one cluster (linearised domain variance
## exactly zero); boundary: no undiagnosed positive; informative: neither
d <- bind_rows(rows) |>
  mutate(p_hat_raw_is_boundary = n_pos == 0 | n_pos == n,
         single_psu = n_psu < 2,
         Z_hat = safe_qlogis(p_hat_raw, n),
         psi   = ifelse(p_hat_raw_is_boundary, NA_real_, (se_p / (p_hat_raw * (1 - p_hat_raw)))^2),
         informative = !p_hat_raw_is_boundary & !single_psu & is.finite(psi) & psi > 0)

## psi_deff2, p_deff2: the binomial deff = 2 values from phia_district_estimates.rds, for comparison
old <- as.data.frame(readRDS("output/mozambique/fits/phia_district_estimates.rds"))
d$psi_deff2 <- old$psi[match(d$area_id, old$area_id)]
d$p_deff2   <- old$p_hat_raw[match(d$area_id, old$area_id)]

cat(sprintf("districts %d | boundary (no undiagnosed positive) %d | single cluster %d | both %d | informative %d\n",
            nrow(d), sum(d$p_hat_raw_is_boundary), sum(d$single_psu),
            sum(d$p_hat_raw_is_boundary & d$single_psu), sum(d$informative)))
inf <- d[d$informative, ]
cat(sprintf("informative districts: median n %d, median clusters %d; psi design median %.3f (IQR %.3f-%.3f); psi deff2 median %.3f; ratio deff2/design median %.2f\n",
            median(inf$n), median(inf$n_psu), median(inf$psi), quantile(inf$psi, .25), quantile(inf$psi, .75),
            median(inf$psi_deff2), median(inf$psi_deff2 / inf$psi)))
cat(sprintf("point estimates: max |p design - p deff2| = %.4f\n", max(abs(d$p_hat_raw - d$p_deff2), na.rm = TRUE)))
dir.create("output/mozambique/fits", showWarnings = FALSE, recursive = TRUE)
saveRDS(d, "output/mozambique/fits/phia_district_estimates_design.rds")
write.csv(d, "output/mozambique/fits/phia_district_estimates_design.csv", row.names = FALSE)
cat("written: output/mozambique/fits/phia_district_estimates_design.{rds,csv}\n")

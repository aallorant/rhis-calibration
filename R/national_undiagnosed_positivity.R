## national_undiagnosed_positivity.R -- national undiagnosed positivity
## theta = P(HIV+ | not previously diagnosed), design-weighted from the PHIA microdata
## (ages 15-64, no province split), with prevalence q, awareness a among people with HIV and
## undiagnosed prevalence U = q(1-a), each with a Wald 95% CI; the direct theta is
## cross-checked against the algebraic form q(1-a)/(1-aq).
## Reads:  data/processed/phia/phia_<country>.rds
## Writes: output/national_undiagnosed_positivity.csv.  Usage: Rscript R/national_undiagnosed_positivity.R
suppressMessages({library(survey); library(dplyr)})
options(survey.lonely.psu = "adjust")

cts <- c("lesotho", "malawi", "mozambique", "zimbabwe")
yrs <- c(lesotho = "2020", malawi = "2020-21", mozambique = "2021", zimbabwe = "2020")
pct <- function(m) { e <- as.numeric(coef(m)); s <- as.numeric(SE(m))
  round(100 * c(e, e - 1.96 * s, e + 1.96 * s), 1) }
rows <- lapply(cts, function(c) {
  ph <- as.data.frame(readRDS(file.path("data/processed/phia", paste0("phia_", c, ".rds"))))
  ph <- subset(ph, age >= 15 & age <= 64 & !is.na(hivstatus))
  ph$hivpos   <- as.numeric(ph$hivstatus == 1)
  ph$awarepos <- as.numeric(ph$aware_bin == TRUE)          # NA for HIV-negatives
  ph$undx     <- as.numeric(ph$hivpos == 1 & ph$awarepos == 0) # HIV+ and not yet diagnosed
  des <- svydesign(id = ~varunit, strata = ~varstrat, weights = ~btwt0, data = ph, nest = TRUE)

  qm  <- svymean(~hivpos, des, na.rm = TRUE)                                          # prevalence q
  am  <- svymean(~awarepos, subset(des, hivpos == 1), na.rm = TRUE)                   # first 95, a
  um  <- svymean(~undx, des, na.rm = TRUE)                                            # undiagnosed prevalence U
  thm <- svymean(~hivpos, subset(des, !(hivpos == 1 & awarepos == 1)), na.rm = TRUE)  # theta
  q <- pct(qm); a <- pct(am); u <- pct(um); th <- pct(thm)
  p <- as.numeric(coef(qm)); aa <- as.numeric(coef(am))

  data.frame(country = tools::toTitleCase(c), year = yrs[[c]], n_1564 = nrow(ph),
             prev_pct = q[1], prev_lo = q[2], prev_hi = q[3],
             aware_pct = a[1], aware_lo = a[2], aware_hi = a[3],
             undx_prev_pct = u[1], undx_lo = u[2], undx_hi = u[3],
             theta_pct = th[1], theta_lo = th[2], theta_hi = th[3],
             theta_formula_pct = round(100 * (p * (1 - aa) / (1 - aa * p)), 1),
             undx_formula_pct  = round(100 * p * (1 - aa), 1))
})
res <- do.call(rbind, rows)
print(res, row.names = FALSE)
dir.create("output", showWarnings = FALSE)
write.csv(res, "output/national_undiagnosed_positivity.csv", row.names = FALSE)

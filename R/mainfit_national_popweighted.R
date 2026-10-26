## =====================================================================
## mainfit_national_popweighted.R -- the national selection odds ratio, standardised to
## the not-yet-diagnosed population by province x stratum: weights are the sums of PHIA
## biomarker weights over adults who are HIV-negative, or HIV-positive and unaware.
## The inverse-variance and equal-weight aggregates of the same cells are reported beside it.
## Reads:  output/canonical/estimates/mainfit_province_stratum_lambda.csv, data/processed/phia/phia_<country>.rds
## Writes: output/canonical/estimates/mainfit_national_popweighted.csv
## Usage:  Rscript R/mainfit_national_popweighted.R
## =====================================================================

suppressMessages(library(dplyr))
CTS <- c("Mozambique", "Zimbabwe", "Malawi", "Lesotho")
EST <- "output/canonical/estimates"

## decode the literal "<U+00E9>" escapes of the canonical CSVs so province names join the PHIA microdata
fix_uesc <- function(x) {
  m <- gregexpr("<U\\+([0-9A-Fa-f]{4})>", x)
  regmatches(x, m) <- lapply(regmatches(x, m), function(s)
    if (length(s)) intToUtf8(strtoi(sub("<U\\+([0-9A-Fa-f]{4})>", "\\1", s), 16L), multiple = TRUE) else s)
  x
}

## ---- not-yet-diagnosed population by province x stratum -------------------
pop_weights <- function(cty, sex_only) {
  d <- as.data.frame(readRDS(file.path("data", "processed", "phia",
                                       paste0("phia_", tolower(cty), ".rds"))))
  w <- if ("btwt0" %in% names(d) && any(is.finite(d$btwt0))) d$btwt0 else d$indweight
  d$w <- w
  ## the estimand's denominator: HIV-negative, or positive and not yet aware
  d <- d |> filter(is.finite(w), w > 0, !is.na(hivstatus)) |>
    filter(hivstatus == 0 | (hivstatus == 1 & (is.na(aware_bin) | aware_bin == 0)))
  d$stratum <- if (sex_only) paste("all", d$sex) else trimws(paste(d$age_group, d$sex))
  d |> group_by(province, stratum) |> summarise(W = sum(w), .groups = "drop")
}

lam <- read.csv(file.path(EST, "mainfit_province_stratum_lambda.csv"), stringsAsFactors = FALSE) |>
  mutate(province = fix_uesc(province))
ivw <- function(e, s) {w <- 1/s^2; sum(e*w)/sum(w)}

rows <- list()
for (cty in CTS) {
  L <- lam |> filter(country == cty)
  sex_only <- all(grepl("^all ", L$stratum))
  W <- pop_weights(cty, sex_only)

  j <- L |> inner_join(W, by = c("province", "stratum"))
  if (nrow(j) < nrow(L))
    message(sprintf("  %s: %d of %d cells matched a population weight", cty, nrow(j), nrow(L)))

  pw   <- sum(j$W * j$log_lambda) / sum(j$W)
  pwse <- sqrt(sum((j$W/sum(j$W))^2 * j$se^2))

  rows[[cty]] <- data.frame(
    country = cty, cells = nrow(j),
    OR_ivw  = exp(ivw(L$log_lambda, L$se)),
    OR_equal = exp(mean(L$log_lambda)),
    OR_pop  = exp(pw),
    lo = exp(pw - 1.96*pwse), hi = exp(pw + 1.96*pwse),
    pct_vs_ivw = 100 * (exp(pw)/exp(ivw(L$log_lambda, L$se)) - 1))
}

out <- bind_rows(rows)
write.csv(out, file.path(EST, "mainfit_national_popweighted.csv"), row.names = FALSE)

cat("\n=== national selection odds ratio, three weightings of the same cells ===\n")
cat("    IVW   = minimum-variance estimator of a COMMON offset (superseded)\n")
cat("    equal = arbitrary\n")
cat("    pop   = standardised to the not-yet-diagnosed population (reported)\n\n")
print(as.data.frame(out |> transmute(country, cells,
  IVW = round(OR_ivw, 2), equal = round(OR_equal, 2),
  population = round(OR_pop, 2),
  `95% CI` = sprintf("(%.2f, %.2f)", lo, hi),
  `pop vs IVW` = sprintf("%+.1f%%", pct_vs_ivw))), row.names = FALSE)

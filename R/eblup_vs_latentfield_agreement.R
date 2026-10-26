## ---------------------------------------------------------------------------
## eblup_vs_latentfield_agreement.R -- agreement between the block-covariance EBLUP and the
## latent-field (TMB) district predictor at the age-collapsed admin-2 resolution: maximum
## absolute difference on the logit scale, correlation, and in percentage points. Units are
## those the fit uses (Zimbabwe's Harare and Bulawayo pooled to one area each by R/utils.R).
## Reuses fit_country() from R/fig_eblup_vs_tmb_scatter.R, sourced up to where the figure is built.
## Reads:  R/fig_eblup_vs_tmb_scatter.R, output/<country>/fits/*.rds
## Writes: output/canonical/estimates/eblup_vs_latentfield_agreement.csv
## Usage:  Rscript R/eblup_vs_latentfield_agreement.R
## ---------------------------------------------------------------------------

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

SRC <- "R/fig_eblup_vs_tmb_scatter.R"
OUT <- "output/canonical/estimates/eblup_vs_latentfield_agreement.csv"

## reuse the fit, stopping before the figure is built
src <- readLines(SRC)
stop_at <- grep("^df <- bind_rows", src)[1]
stopifnot(!is.na(stop_at))
eval(parse(text = paste(src[seq_len(stop_at - 1L)], collapse = "\n")),
     envir = globalenv())

rows <- lapply(c("Mozambique", "Zimbabwe", "Malawi", "Lesotho"), function(k) {
  r <- tryCatch(fit_country(k), error = function(e) {
    message(k, " failed: ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  le <- qlogis(r$eblup); lt <- qlogis(r$tmb)
  ok <- is.finite(le) & is.finite(lt)
  data.frame(country      = k,
             n_units      = sum(ok),
             max_abs_logit = max(abs(le[ok] - lt[ok])),
             correlation   = cor(le[ok], lt[ok]),
             max_abs_pp    = 100 * max(abs(r$eblup[ok] - r$tmb[ok])),
             stringsAsFactors = FALSE)
})
res <- do.call(rbind, rows)

cat("\nBlock-covariance EBLUP against the latent-field predictor,\n")
cat("age-collapsed resolution (the one whose predictions the paper prints).\n\n")
print(transform(res,
                max_abs_logit = round(max_abs_logit, 4),
                correlation   = round(correlation, 6),
                max_abs_pp    = round(max_abs_pp, 3)), row.names = FALSE)

cat("\nSummary:\n")
cat(sprintf("  max absolute difference on the logit scale = %.3f (largest, %s)\n",
            max(res$max_abs_logit), res$country[which.max(res$max_abs_logit)]))
cat(sprintf("  correlations above %.4f\n", floor(min(res$correlation) * 1e4) / 1e4))

write.csv(res, OUT, row.names = FALSE)
cat("\nwrote ", OUT, "\n", sep = "")

# utils.R: shared helpers, loaded with source("R/utils.R") from the repository root.
# Reads:  output/qa_zero_positive_district_units.csv (register QA flags). Writes: nothing.

# harmonise_register_units(df, country, path): every unit of a province that contains a
# unit flagged zero_positive_artifact in the QA file is pooled into one area per province
# (and per stratum where df carries age_group_harmonised and sex); tests and positives
# are summed and the logit z and its variance v = 1/y + 1/(n - y) recomputed.
#   df:      data frame with province, district, n_tests, n_positive and (z_d, v_d) or (z, v)
#   country: country name in any case
#   path:    QA flag file, relative to the repository root
# Returns df with the units of flagged provinces pooled; other rows unchanged.
harmonise_register_units <- function(df, country,
                                     path = "output/qa_zero_positive_district_units.csv") {
  qa <- utils::read.csv(path, stringsAsFactors = FALSE)
  qa <- qa[qa$flag == "zero_positive_artifact" & tolower(qa$country) == tolower(country), ]
  provs <- unique(qa$province)
  if (!length(provs)) return(df)
  zc <- if ("z_d" %in% names(df)) c("z_d", "v_d") else c("z", "v")
  keys <- intersect(c("province", "age_group_harmonised", "sex"), names(df))
  inp <- df$province %in% provs
  pooled <- lapply(split(df[inp, , drop = FALSE], df[inp, keys, drop = FALSE], drop = TRUE), function(g) {
    r <- g[1, , drop = FALSE]
    r$district <- r$province
    r$n_tests <- sum(g$n_tests); r$n_positive <- sum(g$n_positive)
    if ("n_negative" %in% names(r)) r$n_negative <- sum(g$n_negative)
    y <- r$n_positive; n <- r$n_tests
    r[[zc[1]]] <- if (y > 0 && y < n) qlogis(y / n) else NA_real_
    r[[zc[2]]] <- if (y > 0 && y < n) 1 / y + 1 / (n - y) else NA_real_
    r
  })
  out <- rbind(df[!inp, , drop = FALSE], do.call(rbind, pooled))
  message(sprintf("[register QA] %s: pooled %d rows of %s into %d",
                  country, sum(inp), paste(provs, collapse = ", "), length(pooled)))
  out
}

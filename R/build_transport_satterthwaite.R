## build_transport_satterthwaite.R -- Satterthwaite sensitivity check on the standard-normal
## reference of the per-province transport test: S^2 ~ c chi^2_k with moments matched to
## (tau, Var(S^2)), k = 2 tau^2 / Var, c = Var / (2 tau); the chi-square p-value is reported
## beside the normal one, which remains the reported test.
## Reads:  output/canonical/estimates/transport_diagnostic_by_province.csv (Mozambique: boundary-dropped run)
## Writes: output/canonical/estimates/transport_satterthwaite.csv.  Usage: Rscript R/build_transport_satterthwaite.R
tp <- read.csv("output/canonical/estimates/transport_diagnostic_by_province.csv")
tp <- tp[!(tp$country == "Mozambique" & !tp$boundary_dropped), ]
tp$k <- 2 * tp$tau^2 / tp$sd_S2^2; tp$c <- tp$sd_S2^2 / (2 * tp$tau)
tp$p_satt <- pchisq(tp$S2 / tp$c, tp$k, lower.tail = FALSE)
out <- tp[, c("country", "prov", "n_dist", "zstat", "p_one_sided", "k", "c", "p_satt")]
write.csv(out, "output/canonical/estimates/transport_satterthwaite.csv", row.names = FALSE)
cat(sprintf("provinces %d; k %.2f to %.2f; rejections normal %d, Satterthwaite %d; min p normal %.3f (%s), min p satt %.3f\n",
    nrow(out), min(out$k), max(out$k), sum(out$p_one_sided < 0.05), sum(out$p_satt < 0.05),
    min(out$p_one_sided), out$prov[which.min(out$p_one_sided)], min(out$p_satt)))

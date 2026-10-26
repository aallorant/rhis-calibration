## =====================================================================
## transport_diagnostic.R -- double-difference diagnostic of the transport assumption
## lambda_{d,g} = lambda_{p(d),g}, at district-marginal (age-collapsed) resolution, for
## Mozambique, Zimbabwe and Malawi (Lesotho's survey resolves only its own province level).
## s_d = (z_d - zbar_p) - (y_d - ybar_p) removes the true positivity; under transport
## E[S_p^2] = tau_p, with S_p^2 = sum_d s_d^2 and tau_p = (n_p-1)/n_p sum_d (v_d + psi_d), and
## sigma_tau^2 = sum_p (S_p^2 - tau_p) / sum_p (n_p - 1) estimates the within-province offset variance.
## Reads:  output/<country>/fits/*.rds, output/canonical/estimates/{collapsed_admin2_summary,mainfit_sigma_w_survey}.csv
## Writes: output/canonical/estimates/transport_diagnostic{,_by_province,_power}.csv.  Usage: Rscript R/transport_diagnostic.R [--power]
## =====================================================================

suppressMessages({library(dplyr)})

OUT <- "output/canonical/estimates"

## ---- name handling for the Mozambique join ---------------------------------
fix_uesc <- function(x) {
  m <- gregexpr("<U\\+([0-9A-Fa-f]{4})>", x)
  regmatches(x, m) <- lapply(regmatches(x, m), function(s)
    if (length(s)) intToUtf8(strtoi(sub("<U\\+([0-9A-Fa-f]{4})>", "\\1", s), 16L),
                             multiple = TRUE) else s)
  x
}
norm_name <- function(x) {
  x <- fix_uesc(x); x <- iconv(x, "UTF-8", "ASCII//TRANSLIT")
  gsub("[^a-z]", "", tolower(x))
}
alias_name <- function(x) {
  k <- norm_name(x)
  k <- sub("^cidaded[ae]", "", k)
  k <- sub("^cidade", "", k)
  k <- sub("city$", "urban", k)
  k
}
## PHIA district spellings not carried to the register's by the rules above, applied to PHIA keys
## only (the two "Distrito De" entries are the city districts; PHIA's "Nacala" is Nacala-Porto)
phia_to_register <- c(vilankulo = "vilankulos", sussendenga = "sussundenga",
  zumbu = "zumbo", namaroi = "namarroi", nampula = "distritodenampula",
  morrupula = "murrupula", nacala = "nacalaporto", lichinga = "distritodelichinga",
  chimbunila = "chimbonila", kampfumu = "kampfumo", nhlamankulu = "nlhamankulu",
  kamaxakeni = "kamaxakene")
phia_key <- function(x) { k <- alias_name(x); ifelse(k %in% names(phia_to_register), phia_to_register[k], k) }

## ---- assemble (province, district, z_d, v_d, y_d, psi_d) per country -------
get_zwe <- function() {
  d <- as.data.frame(readRDS("output/zimbabwe/fits/district_estimation_dataset.rds"))
  d |> filter(has_both) |>
    transmute(prov = province, dist = district,
              z = z_d, v = v_d, y = Z_hat_d, psi = psi_d, boundary = FALSE)
}
get_mwi <- function() {
  d <- as.data.frame(readRDS("output/malawi/fits/district_estimation_dataset.rds"))
  d |> filter(has_both) |>
    transmute(prov = zone, dist = district,
              z = z_d, v = v_d, y = Z_hat_d, psi = psi_d, boundary = FALSE)
}
get_moz <- function() {
  rh <- as.data.frame(readRDS("output/mozambique/fits/rhis_district_2021.rds"))
  ## design-based district survey estimates from R/build_moz_district_survey_design.R;
  ## boundary = unusable: no undiagnosed positive in the sample, or a single sampled cluster
  ph <- as.data.frame(readRDS("output/mozambique/fits/phia_district_estimates_design.rds"))
  rh$key <- alias_name(rh$district)
  ph$key <- phia_key(ph$area_name)
  j <- inner_join(
    rh |> transmute(key, dist = fix_uesc(district), z = z_d, v = v_d),
    ph |> transmute(key, prov = fix_uesc(province), y = Z_hat,
                    psi = ifelse(is.finite(psi), psi, 0),
                    boundary = !informative,
                    zero_count = p_hat_raw_is_boundary, single_psu = single_psu),
    by = "key")
  cat(sprintf("  Mozambique join: %d of %d PHIA districts matched to the register (%d unusable: %d with no undiagnosed positive, %d in a single cluster)\n",
              nrow(j), nrow(ph), sum(j$boundary), sum(j$zero_count), sum(j$single_psu)))
  j |> select(prov, dist, z, v, y, psi, boundary)
}

## ---- the diagnostic, province by province ---------------------------------
province_stats <- function(g) {
  n <- nrow(g)
  s <- (g$z - mean(g$z)) - (g$y - mean(g$y))
  S2 <- sum(s^2)
  sig2 <- g$v + g$psi
  tot <- sum(sig2)
  tau <- (n - 1) / n * tot
  V_dd <- 2 * ((n - 2) / n * sig2 + tot / n^2)^2
  ## Var(S_p^2) under normal errors: diagonal terms V_dd plus ordered pairs d != l
  A <- outer(sig2, sig2, function(a, b) 2 * (-a / n - b / n + tot / n^2)^2)
  diag(A) <- 0
  varS2 <- sum(V_dd) + sum(A)
  data.frame(n_dist = n, S2 = S2, tau = tau, t = S2 / tau - 1,
             sd_S2 = sqrt(varS2),
             zstat = (S2 - tau) / sqrt(varS2))
}

run_country <- function(cty, dat, drop_boundary) {
  d <- if (drop_boundary) dat |> filter(!boundary) else dat
  d <- d |> filter(is.finite(z), is.finite(v), is.finite(y), is.finite(psi))
  kept <- d |> group_by(prov) |> filter(n() >= 2) |> ungroup()
  dropped <- nrow(d) - nrow(kept)
  per <- kept |> group_by(prov) |> group_modify(~ province_stats(.x)) |> ungroup() |>
    mutate(country = cty, boundary_dropped = drop_boundary,
           p_one_sided = pnorm(zstat, lower.tail = FALSE))

  S2t <- sum(per$S2); taut <- sum(per$tau); sdt <- sqrt(sum(per$sd_S2^2))
  df <- sum(per$n_dist - 1)
  data.frame(country = cty,
             boundary_dropped = drop_boundary,
             n_prov = nrow(per), n_dist = sum(per$n_dist),
             singleton_dropped = dropped,
             S2 = S2t, tau = taut,
             t = S2t / taut - 1,
             zstat = (S2t - taut) / sdt,
             p_one_sided = pnorm((S2t - taut) / sdt, lower.tail = FALSE),
             sigma_tau2 = (S2t - taut) / df,
             sigma_tau2_se = sdt / df,
             sigma_tau2_lo = (S2t - taut - 1.96 * sdt) / df,
             sigma_tau2_hi = (S2t - taut + 1.96 * sdt) / df) |>
    list(national = _, per_province = per)
}

## ---- power of the pooled test at a given country's configuration ----------
power_curve <- function(dat, grid = c(0, 0.02, 0.05, 0.10, 0.20), B = 2000) {
  d <- dat |> filter(!boundary, is.finite(z), is.finite(v), is.finite(y), is.finite(psi)) |>
    group_by(prov) |> filter(n() >= 2) |> ungroup()
  parts <- split(d$v + d$psi, d$prov)
  mom <- lapply(parts, function(s2) {
    n <- length(s2); tot <- sum(s2)
    V_dd <- 2 * ((n - 2) / n * s2 + tot / n^2)^2
    A <- outer(s2, s2, function(a, b) 2 * (-a / n - b / n + tot / n^2)^2); diag(A) <- 0
    c(tau = (n - 1) / n * tot, var = sum(V_dd) + sum(A))
  })
  taut <- sum(sapply(mom, `[`, "tau")); sdt <- sqrt(sum(sapply(mom, `[`, "var")))
  do.call(rbind, lapply(grid, function(st2) {
    zs <- replicate(B, {
      S2 <- sum(vapply(parts, function(s2) {
        x <- rnorm(length(s2), 0, sqrt(st2)) + rnorm(length(s2), 0, sqrt(s2))
        sum((x - mean(x))^2)
      }, numeric(1)))
      (S2 - taut) / sdt
    })
    data.frame(sigma_tau2 = st2, mean_z = mean(zs),
               power_05 = mean(zs > qnorm(0.95)))
  }))
}

## ---- run -------------------------------------------------------------------
cat("\n=== transport diagnostic, district-marginal resolution ===\n\n")
dats <- list(Mozambique = get_moz(), Zimbabwe = get_zwe(), Malawi = get_mwi())

res <- list()
res[["Zimbabwe"]]         <- run_country("Zimbabwe", dats$Zimbabwe, FALSE)
res[["Malawi"]]           <- run_country("Malawi",   dats$Malawi,   FALSE)
res[["Mozambique_all"]]   <- run_country("Mozambique", dats$Mozambique, FALSE)
res[["Mozambique_clean"]] <- run_country("Mozambique", dats$Mozambique, TRUE)

nat <- bind_rows(lapply(res, `[[`, "national"))
per <- bind_rows(lapply(res, `[[`, "per_province"))

## sigma_u^2 from the collapsed admin-2 fit, the matching resolution
su <- read.csv(file.path(OUT, "collapsed_admin2_summary.csv"))
nat$sigma_u2 <- su$sigma_u2_eblup[match(nat$country, su$country)]
nat$tau_share <- nat$sigma_tau2 / nat$sigma_u2

## survey-only sigma_w^2, the same split attacked from the other side
sw <- tryCatch(read.csv(file.path(OUT, "mainfit_sigma_w_survey.csv")),
               error = function(e) NULL)
if (!is.null(sw)) {
  nat$sigma_w2_hi <- sw$hi[match(nat$country, sw$country)]
  nat$unexplained <- nat$sigma_u2 - nat$sigma_tau2 - sw$sigma_w2[match(nat$country, sw$country)]
}

cat("\n--- national (provinces pooled) ---\n\n")
print(as.data.frame(nat |> transmute(
  country, boundary_dropped, n_prov, n_dist,
  S2 = round(S2, 2), tau = round(tau, 2), t = round(t, 3),
  z = round(zstat, 2), p = signif(p_one_sided, 3),
  `sigma_tau^2` = round(sigma_tau2, 4),
  `95% CI` = sprintf("(%.3f, %.3f)", sigma_tau2_lo, sigma_tau2_hi),
  `sigma_u^2` = sigma_u2,
  `tau share` = round(tau_share, 3))), row.names = FALSE)

cat("\n--- per province ---\n\n")
print(as.data.frame(per |> transmute(country, prov, n_dist,
  S2 = round(S2, 2), tau = round(tau, 2), t = round(t, 2),
  z = round(zstat, 2), p = signif(p_one_sided, 3))), row.names = FALSE)

write.csv(nat, file.path(OUT, "transport_diagnostic.csv"), row.names = FALSE)
write.csv(per, file.path(OUT, "transport_diagnostic_by_province.csv"), row.names = FALSE)
cat("\nwritten: ", file.path(OUT, "transport_diagnostic.csv"), " and _by_province.csv\n", sep = "")

if ("--power" %in% commandArgs(TRUE)) {
  set.seed(20260812)
  cat("\n--- power of the pooled one-sided 0.05 test ---\n\n")
  pw <- bind_rows(lapply(c("Zimbabwe", "Malawi"), function(k)
    cbind(country = k, power_curve(dats[[k]]))))
  print(as.data.frame(pw |> transmute(country, `sigma_tau^2` = sigma_tau2,
    `mean z` = round(mean_z, 2), power = round(power_05, 3))), row.names = FALSE)
  write.csv(pw, file.path(OUT, "transport_diagnostic_power.csv"), row.names = FALSE)
  cat("\nSize at sigma_tau^2 = 0 exceeds 0.05 because S^2 is right-skewed.\n")
}
cat("\nt <= 0 or small supports transport; large positive does not.\n")
cat("sigma_tau^2 is the moment estimate of the within-province offset variance,\n")
cat("to be read against sigma_w^2 in mainfit_sigma_w_survey.csv (same split, other side).\n")

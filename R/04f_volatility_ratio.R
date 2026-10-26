## ---------------------------------------------------------------------------
## 04f_volatility_ratio.R -- the volatility ratio r = M2_lambda / M2_eta per year from the paired
## two-round data. A K-year change decomposes as K mu + random walk, so drift mu = mean(delta)/K and
## diffusion sigma^2 = [var(delta) - mean(se^2)]/K are annualised separately (sampling error subtracted
## by method of moments) and recombined as M2 = mu^2 + sigma^2. Zimbabwe, Malawi and Lesotho (one
## awareness definition at both rounds) form the pooled reading; Mozambique is reported separately with
## the IMASIDA proxy displacement removed from round 0. Also runs the retrospective freeze-or-update test.
## Reads:  output/temporal_lambda_{paired,by_area_wave,summary}.csv, output/imasida_proxy_displacement.csv
## Writes: output/volatility_ratio.csv, output/freeze_vs_update.csv.  Usage: Rscript R/04f_volatility_ratio.R
## ---------------------------------------------------------------------------

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

PAIR <- "output/temporal_lambda_paired.csv"
WAVE <- "output/temporal_lambda_by_area_wave.csv"
SUMM <- "output/temporal_lambda_summary.csv"
DISP <- "output/imasida_proxy_displacement.csv"
OUT  <- "output/volatility_ratio.csv"

BOOT <- 4000L
SEED <- 20260804L

## countries measured on one awareness definition at both rounds
SAME_DEF <- c("Zimbabwe", "Malawi", "Lesotho")

pair <- read.csv(PAIR, stringsAsFactors = FALSE)
wave <- read.csv(WAVE, stringsAsFactors = FALSE)
summ <- read.csv(SUMM, stringsAsFactors = FALSE)

## ---- 1. attach se_Z at each round, and the inter-round gap in years ----

w0 <- wave[wave$wave == "round0", c("country", "area", "se_Z", "year")]
w1 <- wave[wave$wave == "round1", c("country", "area", "se_Z", "year")]
names(w0) <- c("country", "area", "se_Z_round0", "year0")
names(w1) <- c("country", "area", "se_Z_round1", "year1")

d <- merge(pair, w0, by = c("country", "area"), all.x = TRUE)
d <- merge(d,    w1, by = c("country", "area"), all.x = TRUE)

stopifnot(!any(is.na(d$se_Z_round0)), !any(is.na(d$se_Z_round1)))
stopifnot(nrow(d) == nrow(pair))

d$gap <- d$year1 - d$year0
stopifnot(all(d$gap > 0))

## ---- 2. consistency check against temporal_lambda_summary.csv ----

chk <- do.call(rbind, lapply(split(d, d$country), function(g) {
  s <- summ[summ$country == g$country[1], ]
  data.frame(country = g$country[1],
             sd_dl_here  = sd(g$d_lambda),
             sd_dl_summ  = s$sd_dlambda,
             mean_dl_here = mean(g$d_lambda),
             mean_dl_summ = s$mean_dlambda,
             stringsAsFactors = FALSE)
}))
if (max(abs(chk$sd_dl_here - chk$sd_dl_summ),
        abs(chk$mean_dl_here - chk$mean_dl_summ)) > 1e-8) {
  print(chk)
  stop("paired file no longer reproduces temporal_lambda_summary.csv")
}
cat("regression check against temporal_lambda_summary.csv: PASS\n\n")

## ---- 3. Mozambique: remove the IMASIDA proxy displacement from round 0 ----

disp <- read.csv(DISP, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
moz  <- d[d$country == "Mozambique", ]
moz  <- merge(moz, disp[, c("province", "D")],
              by.x = "area", by.y = "province", all.x = TRUE)
if (any(is.na(moz$D)))
  stop("unmatched Mozambique province in displacement file: ",
       paste(moz$area[is.na(moz$D)], collapse = ", "))

## lambda_0 - D  ==  Z_hat_0 + D, so the correction moves both sides
moz$d_lambda <- moz$d_lambda + moz$D
moz$d_Zhat   <- moz$d_Zhat   - moz$D
moz$country  <- "Mozambique (corrected)"

## ---- 4. per-area annualised innovations and their sampling variances ----

prep <- function(g) {
  data.frame(country = g$country,
             area    = g$area,
             K       = g$gap,
             ## observed K-year changes
             dl = g$d_lambda,
             de = g$d_Zhat,
             ## sampling variances: lambda carries both z and Z error, eta only Z
             v_dl = g$se_lambda_round0^2 + g$se_lambda_round1^2,
             v_de = g$se_Z_round0^2      + g$se_Z_round1^2,
             d_prev = g$d_prev,
             ## kept for the within-year cross-area spread of lambda
             lambda_round1    = g$lambda_round1,
             se_lambda_round1 = g$se_lambda_round1,
             stringsAsFactors = FALSE)
}

A         <- prep(d[d$country %in% SAME_DEF, ])
A_moz     <- prep(moz)
A_moz_raw <- prep(d[d$country == "Mozambique", ])

## ---- 5. the estimator, per year, pooled over countries c:
##   sigma^2 = [sum_i (delta_i - mean_c(delta))^2 / K_c] / (N - C) - mean_i(v_i / K_c)
##   mu_c = mean_c(delta) / K_c, with mu_c^2 debiased by Var(mu_c) = (K_c sigma^2 + mean_c(v)) / (n_c K_c^2)
##   M2 = weighted mean_c(mu_c^2 debiased) + sigma^2;  r = M2_lambda / M2_eta ----

## sigma^2 and the debiased drift term for one quantity over a set of areas
components <- function(delta, v, K, ctry) {
  cs <- split(seq_along(delta), ctry)
  C  <- length(cs); N <- length(delta)
  if (N <= C) return(c(sigma2 = NA_real_, drift2 = NA_real_, M2 = NA_real_))

  ## diffusion: centre within country, annualise, correct for sampling error
  ss <- sum(unlist(lapply(cs, function(i)
    sum((delta[i] - mean(delta[i]))^2) / K[i][1])))
  sigma2 <- ss / (N - C) - mean(v / K)

  ## drift: per-country mean slope, squared, debiased
  dr <- vapply(cs, function(i) {
    n <- length(i); k <- K[i][1]
    mu <- mean(delta[i]) / k
    mu^2 - (k * max(sigma2, 0) + mean(v[i])) / (n * k^2)
  }, numeric(1))
  nc <- vapply(cs, length, integer(1))
  drift2 <- sum(dr * nc) / N

  c(sigma2 = sigma2, drift2 = drift2, M2 = drift2 + sigma2)
}

both <- function(a) {
  L <- components(a$dl, a$v_dl, a$K, a$country)
  E <- components(a$de, a$v_de, a$K, a$country)
  c(sigma2_lambda = L[["sigma2"]], drift2_lambda = L[["drift2"]],
    M2_lambda     = L[["M2"]],
    sigma2_eta    = E[["sigma2"]], drift2_eta    = E[["drift2"]],
    M2_eta        = E[["M2"]],
    r_M2          = L[["M2"]] / E[["M2"]],
    r_sigma       = L[["sigma2"]] / E[["sigma2"]])
}

## cross-area spread of lambda within one year (round 1), centred within country, corrected for sampling error
spread_lambda <- function(g) {
  cs <- split(seq_len(nrow(g)), g$country)
  N <- nrow(g); C <- length(cs)
  if (N <= C) return(NA_real_)
  ss <- sum(unlist(lapply(cs, function(i)
    sum((g$lambda_round1[i] - mean(g$lambda_round1[i]))^2))))
  sqrt(max(ss / (N - C) - mean(g$se_lambda_round1^2), 0))
}

## bootstrap over areas within country: 2.5% and 97.5% quantiles of r and the share of resamples with r < 1
boot_r <- function(a, B = BOOT) {
  if (nrow(a) < 4) return(c(NA_real_, NA_real_, NA_real_))
  cs <- split(seq_len(nrow(a)), a$country)
  set.seed(SEED)
  r <- replicate(B, {
    i <- unlist(lapply(cs, function(j) sample(j, replace = TRUE)))
    both(a[i, ])[["r_M2"]]
  })
  r <- r[is.finite(r)]
  if (!length(r)) return(c(NA_real_, NA_real_, NA_real_))
  c(unname(quantile(r, c(0.025, 0.975))), mean(r < 1))
}

## r after residualising dLambda on the change in prevalence (illustrative only)
r_resid <- function(a) {
  if (nrow(a) < 5 || all(is.na(a$d_prev))) return(NA_real_)
  f <- if (length(unique(a$country)) > 1) dl ~ d_prev + country else dl ~ d_prev
  b <- a
  b$dl <- residuals(lm(f, data = a))
  both(b)[["r_M2"]]
}

rows <- list()
add <- function(label, a, n_note = "") {
  q  <- boot_r(a)
  v  <- both(a)
  sl <- sqrt(max(v[["M2_lambda"]], 0))
  sp <- spread_lambda(a)
  rows[[length(rows) + 1L]] <<- data.frame(
    group   = label,
    n_areas = nrow(a),
    sd_lambda = sl,
    sd_eta    = sqrt(max(v[["M2_eta"]], 0)),
    t(v),
    r_M2_lo = q[1], r_M2_hi = q[2], p_r_lt_1 = q[3],
    ## annual movement of lambda as a fraction of its cross-area spread within one year
    spread_lambda    = sp,
    drift_over_spread = sl / sp,
    r_resid_prev = r_resid(a),
    note = n_note,
    stringsAsFactors = FALSE)
}

for (cn in SAME_DEF) add(cn, A[A$country == cn, ], "one awareness definition")
add("Pooled (3 countries)", A, "headline reading")
add("Pooled, Malawi dropped", A[A$country != "Malawi", ],
    "sensitivity: Malawi has 3 areas and no identified diffusion")
add("Mozambique (corrected)", A_moz, "round-0 proxy displacement removed")
add("Mozambique (uncorrected)", A_moz_raw, "rounds differ in definition; not evidence")

res <- do.call(rbind, rows)
rownames(res) <- NULL

## ---- 6. report ----

cat("Annualised movement (log-odds per year), noise-corrected.\n")
cat("sd_lambda = selection, sd_eta = the estimand. r = (sd_lambda/sd_eta)^2.\n\n")
pr <- res[, c("group", "n_areas", "sd_lambda", "sd_eta",
              "r_M2", "r_M2_lo", "r_M2_hi", "p_r_lt_1")]
pr[, 3:8] <- round(pr[, 3:8], 3)
print(pr, row.names = FALSE)

cat("\nr < 1 means updating the offset beats freezing the survey estimate.\n")
cat("p_r_lt_1 is the share of bootstrap resamples below that crossover, which\n")
cat("is the one-sided statement at issue; a negative lower quantile is a\n")
cat("small-sample artefact of the noise correction and carries no meaning.\n\n")

cat("Annual movement as a fraction of the cross-area spread (reversal at about two thirds):\n")
sp <- res[, c("group", "sd_lambda", "spread_lambda", "drift_over_spread")]
sp[, 2:4] <- round(sp[, 2:4], 3)
print(sp, row.names = FALSE)
cat("\n")
cat("Where each M2 comes from (drift^2 + diffusion, per year):\n")
dc <- res[, c("group", "drift2_lambda", "sigma2_lambda",
              "drift2_eta", "sigma2_eta", "r_sigma")]
dc[, 2:6] <- round(dc[, 2:6], 4)
print(dc, row.names = FALSE)
cat("\nr_sigma is the diffusion-only ratio, which drops each country's common\n")
cat("level shift. It answers a different question from the one at issue, and\n")
cat("its denominator is near zero because the epidemic moved coherently across\n")
cat("areas within a country, so it is unstable and is not the headline.\n")

cat("\nCovariate-residualised r (dLambda on dPrev), illustrative only:\n")
print(data.frame(group = res$group, r_resid_prev = round(res$r_resid_prev, 3)),
      row.names = FALSE)

neg <- which(res$M2_lambda < 0 | res$M2_eta < 0 | res$sigma2_eta < 0)
if (length(neg))
  cat("\nNOTE: a noise-corrected component fell below zero for: ",
      paste(res$group[neg], collapse = ", "),
      "\n  (sampling error exceeds observed spread; too few areas to separate)\n",
      sep = "")

write.csv(res, OUT, row.names = FALSE)
cat("\nwrote ", OUT, "\n", sep = "")

## ---- 7. retrospective test over the real gap: carry round 0 forward to round 1 by freezing
##   (expit(Z_hat_0)) or updating (expit(z_1 - lambda_0)), judged against truth = expit(Z_hat_1),
##   with the uncalibrated register expit(z_1) alongside; RMSE and bias in percentage points ----

OUT2 <- "output/freeze_vs_update.csv"

retro <- function(g, label) {
  truth <- plogis(g$Z_hat_round1)
  fr    <- plogis(g$Z_hat_round0)
  up    <- plogis(g$z_round1 - g$lambda_round0)
  rg    <- plogis(g$z_round1)
  rmse  <- function(p) sqrt(mean((100 * (p - truth))^2))
  bias  <- function(p) mean(100 * (p - truth))
  data.frame(group = label, n_areas = nrow(g),
             gap_years = if (length(unique(g$gap)) == 1L) g$gap[1] else NA_integer_,
             rmse_freeze = rmse(fr), rmse_update = rmse(up),
             rmse_register = rmse(rg),
             bias_freeze = bias(fr), bias_update = bias(up),
             gain_from_updating = rmse(fr) - rmse(up),
             stringsAsFactors = FALSE)
}

## Mozambique corrected: lambda_0 - D is the same as Z_hat_0 + D
mozc <- d[d$country == "Mozambique", ]
mozc <- merge(mozc, disp[, c("province", "D", "D_sens")],
              by.x = "area", by.y = "province", all.x = TRUE)
mozs <- mozc
mozc$lambda_round0 <- mozc$lambda_round0 - mozc$D
mozc$Z_hat_round0  <- mozc$Z_hat_round0  + mozc$D

## sensitivity: the narrower reading of the round-0 proxy (D_sens), and no correction at all
mozs$lambda_round0 <- mozs$lambda_round0 - mozs$D_sens
mozs$Z_hat_round0  <- mozs$Z_hat_round0  + mozs$D_sens

r2 <- rbind(
  do.call(rbind, lapply(SAME_DEF, function(cn) retro(d[d$country == cn, ], cn))),
  retro(mozc, "Mozambique (corrected)"),
  retro(mozs, "Mozambique (corrected, narrow proxy)"),
  retro(d[d$country == "Mozambique", ], "Mozambique (uncorrected)"),
  retro(d[d$country %in% SAME_DEF, ], "Pooled (3 countries)"))

cat("\n\nRetrospective test: round 0 carried to round 1, judged against round 1.\n")
cat("RMSE of undiagnosed positivity in percentage points.\n\n")
p2 <- r2[, c("group", "n_areas", "gap_years", "rmse_register",
             "rmse_freeze", "rmse_update", "gain_from_updating")]
p2[, 4:7] <- round(p2[, 4:7], 2)
print(p2, row.names = FALSE)
cat("\nMean signed error (percentage points), freeze then update:\n")
print(data.frame(group = r2$group,
                 freeze = round(r2$bias_freeze, 2),
                 update = round(r2$bias_update, 2)), row.names = FALSE)

write.csv(r2, OUT2, row.names = FALSE)
cat("\nwrote ", OUT2, "\n", sep = "")

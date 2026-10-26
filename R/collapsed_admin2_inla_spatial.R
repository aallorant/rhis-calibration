## =====================================================================
## collapsed_admin2_inla_spatial.R -- fully Bayesian (INLA) cross-check of the age-collapsed
## admin-2 block-covariance fit, with an iid, besag or bym2 within-province field and the
## DIC comparison between them:
##   y_d = z_d - Zhat_{p(d)} = lambda_p + b_p + u_d + e_d,  b_p ~ N(0, psi_p) known,
##   u_d ~ field with sigma_u^2 estimated, e_d ~ N(0, v_d) known;  Var(lambda_p) = psi_p + 1/S_p.
## Reads:  output/<country>/fits/*.rds, output/canonical/estimates/collapsed_admin2_summary.csv
## Writes: output/canonical/estimates/collapsed_admin2_inla_spatial.csv, collapsed_admin2_inla_lambda.csv
## Usage:  Rscript R/collapsed_admin2_inla_spatial.R [field ...]   field in {iid, besag, bym2}, default iid
## =====================================================================

suppressMessages({library(INLA); library(dplyr); library(Matrix)})

CTS   <- c("Mozambique", "Zimbabwe", "Malawi", "Lesotho")
OUTD  <- "output/canonical/estimates"
GRAPH <- function(cn) file.path("output", tolower(cn), "fits", "district_adj_canonical.graph")

## countries whose fitted units are areal, so that an ICAR/BYM2 prior is meaningful; Zimbabwe
## is excluded because 19 of its 80 units are hospitals or city health districts with no footprint
SPATIAL_OK <- c("Mozambique", "Malawi", "Lesotho")

## the one aggregation rule (identical to ivw() in RUN_canonical.R)
ivw <- function(est, se) {w <- 1/se^2; c(m = sum(est * w) / sum(w), se = sqrt(1 / sum(w)))}

## ---------------------------------------------------------------- data ----
load_country <- function(cn) {
  s  <- tolower(cn)
  dp <- as.data.frame(readRDS(file.path("output", s, "fits", "district_predictions.rds")))
  pc <- as.data.frame(readRDS(file.path("output", s, "fits", "province_calibration_conditional.rds")))
  d  <- dp |>
    select(province, district, z_d, v_d) |>
    inner_join(select(pc, province, Z_hat_p, psi_p), by = "province") |>
    filter(is.finite(z_d), is.finite(v_d), v_d > 0,
           is.finite(Z_hat_p), is.finite(psi_p), psi_p > 0) |>
    mutate(y = z_d - Z_hat_p)
  d$prov_f   <- factor(d$province)
  d$prov_idx <- as.integer(d$prov_f)
  d$dist_idx <- seq_len(nrow(d))
  d
}

## --------------------------------------------------------------- model ----
fit_inla <- function(d, field = "iid", graph = NULL) {
  psi_p <- d |> group_by(prov_idx) |> summarise(psi = first(psi_p), .groups = "drop") |>
    arrange(prov_idx) |> pull(psi)
  Cpsi <- Diagonal(x = 1 / psi_p)

  pc_prec <- "hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)))"
  fld <- switch(field,
    iid   = sprintf("f(dist_idx, model='iid', %s)", pc_prec),
    besag = sprintf("f(dist_idx, model='besag', graph=graph, scale.model=TRUE, %s)", pc_prec),
    bym2  = sprintf(paste0("f(dist_idx, model='bym2', graph=graph, scale.model=TRUE, ",
                           "hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)), ",
                           "phi=list(prior='pc', param=c(0.5,0.5))))")),
    stop("unknown field: ", field))

  ff <- as.formula(paste0(
    "y ~ 0 + prov_f",
    " + f(prov_idx, model='generic0', Cmatrix=Cpsi,",
    "     hyper=list(prec=list(initial=0, fixed=TRUE)), constr=FALSE)",
    " + ", fld))
  environment(ff) <- environment()          # so Cpsi / graph resolve

  inla(ff, family = "gaussian", data = d, scale = 1 / d$v_d,
       control.family    = list(hyper = list(prec = list(initial = 0, fixed = TRUE))),
       control.compute   = list(dic = TRUE, waic = TRUE),
       control.predictor = list(compute = FALSE),
       control.inla      = list(strategy = "simplified.laplace", int.strategy = "ccd"))
}

## ----------------------------------------------------- summarise a fit ----
## sigma_u^2 posterior: transform the precision marginal.  For bym2 the
## relevant hyperparameter is the total precision of the field.
sigma_u2_post <- function(res) {
  nm <- grep("^Precision for dist_idx", rownames(res$summary.hyperpar), value = TRUE)[1]
  if (is.na(nm)) return(c(mean = NA, sd = NA, q025 = NA, q975 = NA))
  m <- INLA::inla.tmarginal(function(x) 1 / x, res$marginals.hyperpar[[nm]])
  z <- INLA::inla.zmarginal(m, silent = TRUE)
  c(mean = z$mean, sd = z$sd, q025 = z$quant0.025, q975 = z$quant0.975)
}

lambda_post <- function(res) {
  fx <- res$summary.fixed
  data.frame(province = sub("^prov_f", "", rownames(fx)),
             mean = fx$mean, sd = fx$sd, row.names = NULL)
}

## ---- block-structure check: INLA sd of lambda_p against the block target
## psi_p + 1/S_p, S_p = sum_d 1/(sigma_u^2 + v_d), and the independent 1/sum_d 1/(sigma_u^2+v_d+psi_p) ----
validate_block <- function(d, lam, su2) {
  tgt <- d |> group_by(prov_idx, province) |>
    summarise(psi   = first(psi_p),
              block = sqrt(first(psi_p) + 1 / sum(1 / (su2 + v_d))),
              indep = sqrt(1 / sum(1 / (su2 + v_d + first(psi_p)))),
              .groups = "drop")
  m <- inner_join(lam, tgt, by = "province")
  c(max_rel_err_block = max(abs(m$sd / m$block - 1)),
    max_rel_err_indep = max(abs(m$sd / m$indep - 1)))
}

## ---------------------------------------------------------------- main ----
fields <- commandArgs(trailingOnly = TRUE)
if (!length(fields)) fields <- "iid"

ms  <- read.csv(file.path(OUTD, "collapsed_admin2_summary.csv"), stringsAsFactors = FALSE)
rows <- list(); lamrows <- list()

for (cn in CTS) {
  d   <- load_country(cn)
  su2 <- ms$sigma_u2_eblup[ms$country == cn]

  for (fd in fields) {
    gfile <- GRAPH(cn)
    if (fd != "iid" && !(cn %in% SPATIAL_OK)) {
      message(sprintf("  %-11s %-5s SKIPPED (units are not areal; see SPATIAL_OK)", cn, fd)); next
    }
    if (fd != "iid" && !file.exists(gfile)) {
      message(sprintf("  %-11s %-5s SKIPPED (no adjacency graph at %s)", cn, fd, gfile)); next
    }
    graph <- if (fd == "iid") NULL else gfile

    res <- try(fit_inla(d, fd, graph), silent = TRUE)
    if (inherits(res, "try-error")) {
      message(sprintf("  %-11s %-5s FAILED: %s", cn, fd, conditionMessage(attr(res, "condition")))); next
    }

    lam <- lambda_post(res)
    agg <- ivw(lam$mean, lam$sd)
    s2  <- sigma_u2_post(res)
    vb  <- if (fd == "iid") validate_block(d, lam, su2) else c(NA, NA)

    rows[[length(rows) + 1]] <- data.frame(
      country = cn, field = fd, n_units = nrow(d), n_prov = nlevels(d$prov_f),
      OR_inla = exp(agg["m"]), OR_lo = exp(agg["m"] - 1.96 * agg["se"]),
      OR_hi = exp(agg["m"] + 1.96 * agg["se"]), se_log_inla = agg["se"],
      sigma_u2_mean = s2["mean"], sigma_u2_q025 = s2["q025"], sigma_u2_q975 = s2["q975"],
      dic = res$dic$dic, waic = res$waic$waic,
      max_rel_err_block = vb[1], max_rel_err_indep = vb[2],
      row.names = NULL)

    lam$country <- cn; lam$field <- fd
    lamrows[[length(lamrows) + 1]] <- lam
    message(sprintf("  %-11s %-5s done (DIC %.1f)", cn, fd, res$dic$dic))
  }
}

out <- bind_rows(rows)
write.csv(out, file.path(OUTD, "collapsed_admin2_inla_spatial.csv"), row.names = FALSE)
write.csv(bind_rows(lamrows), file.path(OUTD, "collapsed_admin2_inla_lambda.csv"), row.names = FALSE)

cat("\n=== INLA block-covariance cross-check ===\n")
print(as.data.frame(out |> transmute(country, field, n_units,
  OR_inla = round(OR_inla, 3), OR_lo = round(OR_lo, 3), OR_hi = round(OR_hi, 3),
  sigma_u2 = round(sigma_u2_mean, 4), dic = round(dic, 1))), row.names = FALSE)

if ("iid" %in% out$field) {
  cat("\n=== does the INLA posterior sd match the block or the independent form? ===\n")
  cat("    (max relative error of INLA sd against each analytic target)\n")
  print(as.data.frame(out |> filter(field == "iid") |> transmute(country,
    vs_block = sprintf("%.2f%%", 100 * max_rel_err_block),
    vs_indep = sprintf("%.2f%%", 100 * max_rel_err_indep))), row.names = FALSE)
}
cat("\nwrote", file.path(OUTD, "collapsed_admin2_inla_spatial.csv"), "\n")

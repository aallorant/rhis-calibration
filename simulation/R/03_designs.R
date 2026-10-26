## ============================================================================
## 03_designs.R -- sampling designs (which DISTRICTS to survey each wave).
##                 Sourced by the run scripts; not run directly.
##
##   full     : every district every wave              (gold-standard ceiling)
##   random   : k random districts each wave            (naive sparse design)
##   adaptive : k districts each wave chosen ONLINE by largest predicted lambda
##              uncertainty, aggregated over a district's strata. Wave 1 random;
##              thereafter an online random-walk filter tracks each unit's lambda
##              and we revisit the districts where calibration is least certain.
##
## A surveyed district observes all S strata. Cost = total interviews conducted.
## Returns list(nmat = D x T district counts, surv = unit-level observations).
## ============================================================================

run_design <- function(pop, cfg, design = c("full", "random", "adaptive"),
                       k, q = 0.05, seed = NULL) {
  design <- match.arg(design)
  if (!is.null(seed)) set.seed(seed)
  D <- pop$D; T <- pop$T; U <- pop$U
  nmat <- matrix(0, D, T)

  if (design == "full") {
    nmat[, ] <- cfg$n_per
    return(list(nmat = nmat, surv = observe_survey(pop, nmat, cfg)))
  }
  if (design == "random") {
    for (t in 1:T) nmat[sample(D, min(k, D)), t] <- cfg$n_per
    return(list(nmat = nmat, surv = observe_survey(pop, nmat, cfg)))
  }

  ## ---- adaptive: interleave sampling with an online per-unit lambda filter --
  Z <- psi <- matrix(NA, U, T)
  a <- rep(0, U); P <- rep(10, U)               # filter state per unit
  du <- pop$idx$d                                # district of each unit
  for (t in 1:T) {
    if (t == 1) {
      seld <- sample(D, min(k, D))
    } else {
      Ppred_u <- P + q
      d_unc <- tapply(Ppred_u, du, mean)          # district uncertainty = mean over strata
      seld  <- as.integer(names(sort(d_unc, decreasing = TRUE)))[1:min(k, D)]
    }
    nmat[seld, t] <- cfg$n_per
    sel_u <- which(du %in% seld)

    pic  <- pmin(pmax(pop$undiag[, t], 1e-4), 1 - 1e-4)
    psis <- cfg$deff / (cfg$n_per * pic * (1 - pic))
    Z[sel_u, t]   <- pop$theta[sel_u, t] + rnorm(length(sel_u), 0, sqrt(psis[sel_u]))
    psi[sel_u, t] <- psis[sel_u]

    for (u in 1:U) {                              # predict, update if sampled
      Ppred <- P[u] + q; apred <- a[u]
      if (!is.na(Z[u, t])) {
        robs <- psi[u, t] + cfg$v_rhis; yo <- pop$z[u, t] - Z[u, t]
        K <- Ppred / (Ppred + robs)
        a[u] <- apred + K * (yo - apred); P[u] <- (1 - K) * Ppred
      } else { a[u] <- apred; P[u] <- Ppred }
    }
  }
  list(nmat = nmat, surv = list(Z = Z, psi = psi, n_d = nmat,
                                n_u = nmat[du, , drop = FALSE]))
}

## ---- scoring ---------------------------------------------------------------
score <- function(pop, est) {
  err <- est$theta - pop$theta
  ## also report on the natural prevalence scale (percentage points)
  pp  <- (plogis(est$theta) - plogis(pop$theta)) * 100
  data.frame(rmse_logit = sqrt(mean(err^2)),
             bias_logit = mean(err),
             rmse_pp    = sqrt(mean(pp^2)),
             mae_pp     = mean(abs(pp)))
}
cost_of <- function(nmat, S) sum(nmat) * S        # total interviews (all strata)

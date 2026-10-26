## ============================================================================
## 02_estimators.R -- estimators of theta_{u,t} = logit(undiagnosed positivity).
##                    Sourced by the run scripts; not run directly.
##   routine_only  : theta_hat = z (ignores lambda; biased, free)
##   survey_only   : theta_hat = Zhat where surveyed, same-stratum wave mean else
##   static_calib  : per-wave calibration theta_hat = z - lambda_hat, with
##                   lambda_hat transported from the audited cells of the wave
##   dynamic_calib : random-walk smoothing of each unit's lambda over waves
## All return list(theta = U x T, var = U x T); `var` feeds the adaptive design.
## ============================================================================

## ---- per-series random-walk Kalman filter + RTS smoother (NA-aware) --------
kalman_rw <- function(y, r, q, P0 = 10) {
  n <- length(y); af <- Pf <- ppred <- numeric(n)
  a_prev <- mean(y, na.rm = TRUE); if (!is.finite(a_prev)) a_prev <- 0
  P_prev <- P0
  for (t in 1:n) {
    apred <- a_prev; Ppred <- P_prev + q; ppred[t] <- Ppred
    if (!is.na(y[t])) {
      K <- Ppred / (Ppred + r[t]); af[t] <- apred + K * (y[t] - apred)
      Pf[t] <- (1 - K) * Ppred
    } else { af[t] <- apred; Pf[t] <- Ppred }
    a_prev <- af[t]; P_prev <- Pf[t]
  }
  as <- af; Ps <- Pf
  if (n >= 2) for (t in (n - 1):1) {
    Pp <- Pf[t] + q; C <- Pf[t] / Pp
    as[t] <- af[t] + C * (as[t + 1] - af[t]); Ps[t] <- Pf[t] + C^2 * (Ps[t + 1] - Pp)
  }
  list(m = as, v = Ps, fpred = ppred)
}

est_routine_only <- function(pop, surv, cfg)
  list(theta = pop$z, var = matrix(Inf, pop$U, pop$T))

est_survey_only <- function(pop, surv, cfg) {
  U <- pop$U; T <- pop$T; s <- pop$idx$s
  th <- surv$Z; v <- surv$psi
  for (t in 1:T) for (ss in unique(s)) {              # fill gaps by stratum-wave mean
    rows <- which(s == ss); obs <- rows[!is.na(surv$Z[rows, t])]
    m <- if (length(obs)) mean(surv$Z[obs, t]) else 0
    miss <- rows[is.na(th[rows, t])]; th[miss, t] <- m; v[miss, t] <- 5
  }
  list(theta = th, var = v)
}

## ---- calibrated register: theta_hat = z - lambda_hat -----------------------
## lambda_hat: GLS mean of lambda-tilde = z - Zhat over the audited cells of the
##   broad-area cell (stratum-wave; one cell per wave if pooled), with weights
##   1/(sigma_tau^2 + psi + v)
## sigma_tau^2: transport-failure variance, REML over all audited cells of the
##   wave with one fixed mean per broad-area cell (reml_s2)
## var: v + Var(lambda_hat) [+ sigma_tau^2 if tau = TRUE]
## Returns theta, var, aud (logical U x T, audited cells) and s2tau (per wave).
reml_s2 <- function(y, r, g, upper = 3) {
  nll <- function(s2) { V <- s2 + r; out <- 0
    for (k in unique(g)) { i <- g == k; w <- 1 / V[i]; m <- sum(w * y[i]) / sum(w)
      out <- out + 0.5 * (sum(log(V[i])) + log(sum(w)) + sum(w * (y[i] - m)^2)) }
    out }
  optimize(nll, c(0, upper))$minimum
}
est_static_calib <- function(pop, surv, cfg, tau = TRUE, pooled = FALSE) {
  U <- pop$U; T <- pop$T
  grp <- if (pooled) rep(1L, U) else pop$idx$s
  th <- v <- matrix(NA, U, T); s2tau <- numeric(T)
  aud <- !is.na(surv$Z)
  for (t in 1:T) {
    lo <- pop$z[, t] - surv$Z[, t]; r <- surv$psi[, t] + cfg$v_rhis
    ## REML estimate of sigma_tau^2 over the audited cells of the wave
    sm_all <- which(aud[, t])
    s2 <- if (length(sm_all) > length(unique(grp[sm_all])))
            reml_s2(lo[sm_all], r[sm_all], grp[sm_all]) else 0
    s2tau[t] <- s2
    for (g in unique(grp)) {
      rows <- which(grp == g); sm <- rows[aud[rows, t]]
      if (length(sm)) {
        w <- 1 / (s2 + r[sm]); lam <- sum(w * lo[sm]) / sum(w); vlam <- 1 / sum(w)
      } else { lam <- 0; vlam <- 1 }
      th[rows, t] <- pop$z[rows, t] - lam
      v[rows, t]  <- cfg$v_rhis + vlam + if (tau) s2 else 0
    }
  }
  list(theta = th, var = v, aud = aud, s2tau = s2tau)
}

## Pooled variant: ONE offset per wave for every stratum, which is what a register
## carrying (d, p) but no g can support. Same estimator, one broad-area cell.
est_static_calib_pooled <- function(pop, surv, cfg, tau = TRUE)
  est_static_calib(pop, surv, cfg, tau = tau, pooled = TRUE)

est_dynamic_calib <- function(pop, surv, cfg, q = 0.05) {
  U <- pop$U; T <- pop$T; s <- pop$idx$s
  lam_obs <- pop$z - surv$Z; r_obs <- surv$psi + cfg$v_rhis
  lam_sm <- lam_v <- matrix(NA, U, T)
  for (u in 1:U) {
    kf <- kalman_rw(lam_obs[u, ], r_obs[u, ], q = q)
    lam_sm[u, ] <- kf$m; lam_v[u, ] <- kf$v
  }
  for (t in 1:T) for (ss in unique(s)) {              # never-sampled units -> stratum mean
    rows <- which(s == ss); miss <- rows[!is.finite(lam_sm[rows, t])]
    if (length(miss)) {
      gm <- mean(lam_sm[rows[is.finite(lam_sm[rows, t])], t]); if (!is.finite(gm)) gm <- 0
      lam_sm[miss, t] <- gm; lam_v[miss, t] <- 1
    }
  }
  cal <- pop$z - lam_sm; v <- lam_v + cfg$v_rhis; th <- cal
  samp <- !is.na(surv$Z)
  wv <- 1 / surv$psi[samp]; wc <- 1 / v[samp]
  th[samp] <- (wv * surv$Z[samp] + wc * cal[samp]) / (wv + wc)
  v[samp]  <- 1 / (wv + wc)
  list(theta = th, var = v)
}

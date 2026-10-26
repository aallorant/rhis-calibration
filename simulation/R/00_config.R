## ============================================================================
## 00_config.R -- parameters for the temporal calibration / adaptive-design
##                simulation. Sourced by every other simulation script; sets
##                the parameter list `cfg` and the output directory OUT.
##                Not run directly.
## ============================================================================

## Estimand: theta_{u,t} = logit( undiagnosed positivity ) for unit u
##   (unit = district x age-sex stratum) at wave t, where
##   undiagnosed positivity  P(HIV+ | not yet diagnosed) = p(1-a) / (1 - p a)
##                           for prevalence p and awareness of status a
##   routine positivity      logit(r) = theta + lambda + noise   (register/RHIS)
##   survey estimate         Zhat = theta + sampling noise        (PHIA/DHS)

cfg <- list(

  ## ---- panel structure ------------------------------------------------------
  G = 10,          # provinces (admin-1)
  d_per_prov = 4,  # districts per province  -> D = 40 districts
  T = 6,           # survey waves / years

  ## age-sex strata (2 sex x 2 age = 4); order fixed and reused everywhere
  strata = data.frame(
    stratum = c("F.15-24", "F.25-49", "M.15-24", "M.25-49"),
    sex     = c("F", "F", "M", "M"),
    age     = c("Y", "O", "Y", "O"),
    stringsAsFactors = FALSE
  ),

  ## population share by stratum (for aggregating to coarser resolutions);
  ## district populations assumed equal
  pop_share = c("F.15-24" = 0.27, "F.25-49" = 0.26,
                "M.15-24" = 0.26, "M.25-49" = 0.21),

  ## ---- epidemiology: stratum baseline prevalence & awareness ----------------
  ## high-burden setting: adult prevalence ~12-13%, female > male, rising with age
  prev0 = c("F.15-24" = 0.060, "F.25-49" = 0.240,   # baseline prevalence by stratum
            "M.15-24" = 0.025, "M.25-49" = 0.150),
  aware0 = c("F.15-24" = 0.62, "F.25-49" = 0.82,    # baseline awareness by stratum
             "M.15-24" = 0.55, "M.25-49" = 0.72),   #   women/older -> higher
  prev_decline = 0.035,   # annual prevalence decline (exp decay), incidence falling
  aware_gain   = 0.12,    # annual rise in awareness (logit scale), treat-scale-up
  prov_spread  = 0.85,    # sd of province effect on logit-prevalence (admin-1 gradient)
  dist_sd      = 0.18,    # district-within-province sd on logit-prevalence
  aware_dist_sd= 0.15,    # district sd on logit-awareness

  ## ---- selection offset lambda: stratum baselines ---------------------------
  ## lambda > 0: routine testing selects higher-yield contacts; men (tested more
  ## selectively) carry a larger offset than women (reached via ANC), youth slightly higher
  lambda0 = c("F.15-24" = 0.55, "F.25-49" = 0.35,
              "M.15-24" = 1.00, "M.25-49" = 0.80),
  lambda_dist_sd = 0.15,  # district variation in the offset
  v_rhis  = 0.010,        # routine positivity logit-noise variance (near-census)

  ## ---- survey design --------------------------------------------------------
  deff   = 2.0,    # survey design effect
  n_per  = 120,    # individuals sampled PER STRATUM in a surveyed district-wave

  ## ---- lambda scenarios -----------------------------------------------------
  ## Each defines how the offset evolves over T waves (added to the baseline):
  ##   constant : no change
  ##   drift    : smooth per-unit linear trend (slope sd = par)
  ##   shock    : structural break at wave t0 of size `par`
  ##   reshuffle: AR(1) area drift, ranks reshuffle (phi, sd in par list)
  scenarios = list(
    constant  = list(type = "constant"),
    drift     = list(type = "drift",     slope_sd = 0.10),
    shock     = list(type = "shock",     t0 = 4, size = 0.50),
    reshuffle = list(type = "reshuffle", phi = 0.30, sd_omega = 0.25)
  )
)

## ---- sigma_u^2 anchoring ---------------------------------------------------
## sigma_u^2: within-province variance of the district residual (truth + offset)
## sigma_u2_k: common scaling of the three district sds that puts the simulated
##   population at sigma_u2_target; obtained with 16_sigma_u_sweep.R and to be
##   re-solved after any change to prev0, aware0, dist_sd, aware_dist_sd or
##   lambda_dist_sd
cfg$sigma_u2_target <- 0.20
cfg$sigma_u2_k      <- 2.050
cfg$dist_sd        <- cfg$dist_sd        * cfg$sigma_u2_k
cfg$aware_dist_sd  <- cfg$aware_dist_sd  * cfg$sigma_u2_k
cfg$lambda_dist_sd <- cfg$lambda_dist_sd * cfg$sigma_u2_k

cfg$D <- cfg$G * cfg$d_per_prov
cfg$S <- nrow(cfg$strata)
OUT <- "simulation/output"

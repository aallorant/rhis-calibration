// full_model_covariates.cpp
//
// Stratum-aware RHIS-PHIA joint model with district-level covariates.
// One latent log-odds prevalence per district (eta_d, iid normal).
// Province x stratum selection bias alpha[p, g].
// District x stratum residual within-province heterogeneity u[d, g] ~ N(0, sigma_u).
// Stratum-specific covariate effects beta[k, g] (k covariates, g strata).
// Stratum prevalence shift delta[g] with delta[0] = 0 by convention.
//
// Model:
//   z_{d,g}   ~ N(eta_d + delta_g + alpha_{p(d),g} + x_d' beta_g + u_{d,g}, v_{d,g})  RHIS
//   Z_hat_{d,g} ~ N(eta_d + delta_g, psi_{d,g})                                       PHIA
//   eta_d   ~ N(0, sigma_eta^2)
//   u_{d,g} ~ N(0, sigma_u^2)
//
// For the baseline (no covariates), pass an X with zero columns -> the beta
// contribution drops to zero.
//
// Inputs (0-indexed integer vectors):
//   d_rhis[i], g_rhis[i]    district and stratum of RHIS obs i
//   d_phia[i], g_phia[i]    district and stratum of survey obs i
//   prov_of_district[d]     province index per district

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {

  // Data
  DATA_VECTOR(z);                  // RHIS log-odds (n_rhis)
  DATA_VECTOR(v);                  // RHIS variance (n_rhis)
  DATA_IVECTOR(d_rhis);            // RHIS district idx
  DATA_IVECTOR(g_rhis);            // RHIS stratum idx

  DATA_VECTOR(Z_hat);              // PHIA log-odds (n_phia)
  DATA_VECTOR(psi);                // PHIA variance (n_phia)
  DATA_IVECTOR(d_phia);            // PHIA district idx
  DATA_IVECTOR(g_phia);            // PHIA stratum idx

  DATA_IVECTOR(prov_of_district);  // province per district (n_districts)
  DATA_MATRIX(X);                  // covariates (n_districts x n_cov; n_cov may be 0)

  // Parameters
  PARAMETER_MATRIX(alpha);         // alpha[p, g]: n_provinces x n_strata
  PARAMETER_VECTOR(delta);         // delta[g - 1] for g >= 1; length n_strata - 1
  PARAMETER_MATRIX(beta);          // beta[k, g]: n_cov x n_strata (n_cov may be 0)
  PARAMETER_VECTOR(eta);           // eta[d]: n_districts (random)
  PARAMETER_MATRIX(u);             // u[d, g]: n_districts x n_strata (random)
  PARAMETER(log_sigma_eta);
  PARAMETER(log_sigma_u);

  Type sigma_eta = exp(log_sigma_eta);
  Type sigma_u   = exp(log_sigma_u);
  Type nll = 0;

  int n_districts = eta.size();
  int n_strata    = alpha.cols();
  int n_cov       = beta.rows();

  // Prior on eta (iid normal)
  for (int d = 0; d < n_districts; d++) {
    nll -= dnorm(eta(d), Type(0), sigma_eta, true);
  }

  // Prior on u (iid normal, district x stratum)
  for (int d = 0; d < n_districts; d++) {
    for (int g = 0; g < n_strata; g++) {
      nll -= dnorm(u(d, g), Type(0), sigma_u, true);
    }
  }

  // RHIS likelihood
  for (int i = 0; i < z.size(); i++) {
    int d = d_rhis(i);
    int g = g_rhis(i);
    int p = prov_of_district(d);
    Type delta_g = (g == 0) ? Type(0) : delta(g - 1);
    Type xbeta = Type(0);
    for (int k = 0; k < n_cov; k++) {
      xbeta += X(d, k) * beta(k, g);
    }
    Type mu = eta(d) + delta_g + alpha(p, g) + xbeta + u(d, g);
    nll -= dnorm(z(i), mu, sqrt(v(i)), true);
  }

  // PHIA likelihood
  for (int i = 0; i < Z_hat.size(); i++) {
    int d = d_phia(i);
    int g = g_phia(i);
    Type delta_g = (g == 0) ? Type(0) : delta(g - 1);
    Type mu = eta(d) + delta_g;
    nll -= dnorm(Z_hat(i), mu, sqrt(psi(i)), true);
  }

  ADREPORT(alpha);
  ADREPORT(beta);
  ADREPORT(delta);
  ADREPORT(sigma_eta);
  ADREPORT(sigma_u);
  REPORT(eta);
  REPORT(u);

  return nll;
}

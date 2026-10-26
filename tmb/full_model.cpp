// full_model.cpp: joint RHIS-survey model with a latent district log-odds Z_d, a global
// selection intercept alpha and optional province effects beta_lambda on lambda through
// X_lambda (zero columns for the intercept-only fit). Used by R/04_two_stage_calibration.R.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {

  // Data
  DATA_VECTOR(z);          // RHIS log-odds (n_rhis)
  DATA_VECTOR(v);          // RHIS variance (n_rhis)
  DATA_VECTOR(Z_hat);      // Survey log-odds (n_survey)
  DATA_VECTOR(psi);        // Survey variance (n_survey)
  DATA_IVECTOR(d_rhis);    // District index for RHIS obs (0-based)
  DATA_IVECTOR(d_survey);  // District index for survey obs (0-based)
  DATA_MATRIX(X_lambda);  // Province design matrix for lambda (n_rhis x p)

  // Parameters
  PARAMETER(alpha);                  // Global lambda intercept
  PARAMETER_VECTOR(beta_lambda);     // Province effects on lambda
  PARAMETER_VECTOR(Z_d);             // Latent Z_d (n_districts)
  PARAMETER(log_sigma_Z);            // SD of Z_d prior

  Type sigma_Z = exp(log_sigma_Z);
  Type nll = 0;

  // Prior on Z_d (iid normal)
  for (int d = 0; d < Z_d.size(); d++) {
    nll -= dnorm(Z_d(d), Type(0), sigma_Z, true);
  }

  // RHIS likelihood: z_d ~ N(Z_d + alpha + X'beta, v_d)
  for (int i = 0; i < z.size(); i++) {
    Type mu_rhis = Z_d(d_rhis(i)) + alpha;
    for (int j = 0; j < beta_lambda.size(); j++) {
      mu_rhis += X_lambda(i, j) * beta_lambda(j);
    }
    nll -= dnorm(z(i), mu_rhis, sqrt(v(i)), true);
  }

  // Survey likelihood: Z_hat_d ~ N(Z_d, psi_d)
  for (int i = 0; i < Z_hat.size(); i++) {
    nll -= dnorm(Z_hat(i), Z_d(d_survey(i)), sqrt(psi(i)), true);
  }

  ADREPORT(Z_d);

  return nll;
}

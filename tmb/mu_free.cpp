// mu_free.cpp: mu_p-free latent-field calibration.
//   z_d   ~ N(eta_p[d] + u_d + Xlam_d . beta, v_d)     RHIS, district
//   Zhat_p~ N(eta_p, psi_p)                            survey, province (noisy observation of mu_p)
//   u_d   ~ N(0, sigma_u^2)                            within-province prevalence deviation
// eta_p (= mu_p) is a FREE province level informed by BOTH arms; bias lambda_d = Xlam_d . beta.
// Xlam = column of ones -> global bias (no covariates); [1, x_d] -> covariate bias.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(z);
  DATA_VECTOR(v);
  DATA_IVECTOR(prov_d);
  DATA_MATRIX(Xlam);
  DATA_VECTOR(Zhat);
  DATA_VECTOR(psi);
  DATA_VECTOR(w);        // district population share within its province
  DATA_INTEGER(agg);     // 1 = include sum_{d in p} w_d u_d in survey arm, 0 = omit
  PARAMETER_VECTOR(eta_p);
  PARAMETER_VECTOR(beta);
  PARAMETER(log_sigma_u);
  PARAMETER_VECTOR(u);
  Type sigma_u = exp(log_sigma_u);
  Type nll = 0.0;
  vector<Type> lambda = Xlam * beta;
  for (int d = 0; d < z.size(); d++) {
    Type eta_d = eta_p(prov_d(d)) + u(d);
    nll -= dnorm(z(d), eta_d + lambda(d), sqrt(v(d)), true);
    nll -= dnorm(u(d), Type(0), sigma_u, true);
  }
  vector<Type> agg_p(Zhat.size()); agg_p.setZero();
  if (agg == 1) for (int d = 0; d < z.size(); d++) agg_p(prov_d(d)) += w(d) * u(d);
  for (int p = 0; p < Zhat.size(); p++)
    nll -= dnorm(Zhat(p), eta_p(p) + agg_p(p), sqrt(psi(p)), true);
  vector<Type> eta_d_all(z.size());
  for (int d = 0; d < z.size(); d++) eta_d_all(d) = eta_p(prov_d(d)) + u(d);
  ADREPORT(beta);
  REPORT(eta_p); REPORT(beta); REPORT(sigma_u); REPORT(eta_d_all); REPORT(lambda);
  return nll;
}

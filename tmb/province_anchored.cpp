// province_anchored.cpp: province-anchored joint calibration model, fitted at admin-2
// (district) with the admin-1 (province) survey estimate as the calibration plug-in,
// lambda homogeneous within province.
//
//   district d in province p:
//     eta_d   = eta_p + u_d,     u_d ~ N(0, sigma_u^2)        [latent prevalence]
//     z_d     ~ N(eta_d + lambda_p, v_d)                       [RHIS, logit]
//   province p:
//     Zhat_p  ~ N(eta_p, psi_p)                                [survey plug-in, admin-1]
//
// lambda_p (province selection log-RTP) and eta_p (province latent prevalence) are
// fixed effects; u_d are random effects integrated out by Laplace. ADREPORT(lambda)
// gives the per-province selection log-RTP with standard errors.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(z);          // district RHIS logit positivity
  DATA_VECTOR(v);          // district RHIS delta-method variance (known)
  DATA_IVECTOR(prov_d);    // province index per district (0-based)
  DATA_VECTOR(Zhat);       // province survey logit undiagnosed prevalence
  DATA_VECTOR(psi);        // province survey design variance (known)

  PARAMETER_VECTOR(lambda);   // province selection log-RTP (length nprov)
  PARAMETER_VECTOR(eta_p);    // province latent prevalence (length nprov)
  PARAMETER_VECTOR(u);        // district random effect (length ndist)
  PARAMETER(log_sigma_u);

  Type sigma_u = exp(log_sigma_u);
  Type nll = 0.0;

  // survey arm: province plug-in anchors eta_p
  for (int p = 0; p < Zhat.size(); p++)
    nll -= dnorm(Zhat(p), eta_p(p), sqrt(psi(p)), true);

  // RHIS arm + latent random effect
  for (int d = 0; d < z.size(); d++) {
    nll -= dnorm(u(d), Type(0.0), sigma_u, true);
    Type eta_d = eta_p(prov_d(d)) + u(d);
    nll -= dnorm(z(d), eta_d + lambda(prov_d(d)), sqrt(v(d)), true);
  }

  ADREPORT(lambda);
  REPORT(sigma_u);
  return nll;
}

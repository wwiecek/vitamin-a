data {
  int<lower=1> K;
  vector[K] y;      // study estimates
  vector<lower=0>[K] se;    // their s.e.
  real<lower=0> c;          // z cutoff (e.g., 1.96)
}
parameters {
  real tau;
  real<lower=0> sigma_tau;
  vector[K] eta;            // non-centered
  // real logomega;   // rel. pub prob for |z|<c
  real<lower=0, upper=1> omega;   // rel. pub prob for |z|<c
}
transformed parameters {
  vector[K] theta_k = tau + sigma_tau * eta;
  // real omega = exp(logomega);
}
model {
  // vague (as in your sketch; tighten if needed)
  tau ~ normal(0, 10);
  sigma_tau ~ normal(0, 10);
  eta ~ normal(0, 1);
  // logomega ~ normal(0, 2);
  omega ~ beta(1, 1);
  
  // likelihood + selection
  for (i in 1:K) {
    // baseline
    target += normal_lpdf(y[i] | theta_k[i], se[i]);
    
    // region weight (observed z)
    if (fabs(y[i] / se[i]) < c)
      target += log(omega);
    
    // selection normaliser (μ on z-scale)
    {
      real sigmatilde = sqrt(square(se[i]) + square(sigma_tau));
      // P(|Z_i|<c | tau, sigma_tau) where Z_i = Y_i / se[i]
      real upp = ( c * se[i] - tau) / sigmatilde;
      real low = (-c * se[i] - tau) / sigmatilde;
      real p_non_marg = normal_cdf(upp, 0, 1) - normal_cdf(low, 0, 1);
      // S_marg = omega * p_non_marg + 1 * (1 - p_non_marg)
      target += -log1m( (1 - omega) * p_non_marg );
    }
  }
}

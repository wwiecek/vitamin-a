// publication_bias_k1_collapsed.stan  — stable version of your second model

data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  real<lower=0> c;     // e.g., 1.96
}

transformed data {
  // data-only indicator for whether each study fell below the cutoff
  array[N] int below;   // 1 if |y[i]/se[i]| < c
  for (i in 1:N) below[i] = fabs(y[i] / se[i]) < c;
  real se_bar = mean(se);
}

parameters {
  real mu;
  real<lower=0> tau_sd;
  real alpha;                 // logit-omega
}

transformed parameters {
  real<lower=0, upper=1> omega = inv_logit(alpha);  // selection weight for |Z|<c
}

model {
  // Priors: tighten scales to kill the funnel; center mu weakly near data scale
  mu     ~ normal(0, 5 * se_bar);
  tau_sd ~ student_t(3, 0, 2 * se_bar);   // half-t(3, 2*se_bar)
  alpha  ~ normal(-1.58, 2);                // implies omega usually in (0.18, 0.82)

  // Selection-adjusted likelihood (collapsed over theta_i)
  for (i in 1:N) {
    real s   = sqrt(se[i]^2 + tau_sd^2);              // marginal sd
    // P(|X_i| < c * se[i]) under X_i ~ N(mu, s^2)
    real upp = ( c * se[i] - mu) / s;
    real low = (-c * se[i] - mu) / s;
    real p   = normal_cdf(upp, 0, 1) - normal_cdf(low, 0, 1);   // in (0,1)

    // log normalizer: log( omega * p + (1 - p) )
    // numerically stable when p ~ 1 and omega ~ 1
    real logZ = log1p( (omega - 1) * p );

    // baseline density
    target += normal_lpdf(y[i] | mu, s);

    // selection weight for this observed bin
    if (below[i] == 1) target += log(omega);  // else +0

    // normalization
    target += -logZ;
  }
}

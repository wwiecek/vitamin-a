  data {
    int<lower=0> K; // number of sites 
    real tau_hat_k[K]; // estimated treatment effects
    real<lower=0> se_k[K]; // s.e. of effect estimates 
  }
  parameters {
    real tau; 
    real<lower=0> sigma_tau;
    real tau_k[K];
  }
  transformed parameters {
  
  }
  model {
    sigma_tau ~ uniform(0,1000);//
    tau ~ normal(0,1000);//
    tau_k ~ normal(tau, sigma_tau); // second level normal
    tau_hat_k ~ normal(tau_k, se_k); //third level normal
  }

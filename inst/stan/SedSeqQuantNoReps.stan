data {
  // Number of RNAs
  int<lower=1> NRNA;     
  
  // Note: These are all integers
  // columns t, s, p
  int<lower=0> tot_obs_counts[NRNA];
  int<lower=0> sup_obs_counts[NRNA];
  int<lower=0> pel_obs_counts[NRNA];
}
parameters {
  // Unnormalized mixing proportions
  // real<lower=0> mixing_t;
  real<lower=0> mixing_factor_sup;
  real<lower=0> mixing_factor_pellet;
  
  // dispersion parameter for counts
  real<lower=0> phi;
}
model{
  // mixing ratios
  mixing_factor_sup ~ gamma(1,1);
  mixing_factor_pellet ~ gamma(1,1);
  // Cauchy prior for negbin dispersion parameter
  phi ~ cauchy(0,3);
  
  for(idx in 1:NRNA){ 
    // count distn negative binomial with specified means
    // Total
    tot_obs_counts[idx] ~ neg_binomial_2(mixing_factor_sup * sup_obs_counts[idx] +
                                          mixing_factor_pellet * pel_obs_counts[idx], phi);
  }
}

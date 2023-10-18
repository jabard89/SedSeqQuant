data {
  // Number of RNAs
  int<lower=1> NRNA;

  // Note: These are all integers
  // columns t, s, p
  int<lower=0> tot_obs_counts[NRNA];
  int<lower=0> sup_obs_counts[NRNA];
  int<lower=0> pel_obs_counts[NRNA];

  // Mixing factor priors
  real mixing_factor_sup_guess_mean;
  real mixing_factor_sup_guess_sd;
  real mixing_factor_pellet_guess_mean;
  real mixing_factor_pellet_guess_sd;
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
  mixing_factor_sup ~ normal(mixing_factor_sup_guess_mean, mixing_factor_sup_guess_sd);
  mixing_factor_pellet ~ normal(mixing_factor_pellet_guess_mean, mixing_factor_pellet_guess_sd);

  // Cauchy prior for negbin dispersion parameter
  phi ~ cauchy(0,3);

  for(idx in 1:NRNA){
    // count distn negative binomial with specified means
    tot_obs_counts[idx] ~ neg_binomial_2(sup_obs_counts[idx] / mixing_factor_sup +
                                          pel_obs_counts[idx] / mixing_factor_pellet, phi);
  }
}

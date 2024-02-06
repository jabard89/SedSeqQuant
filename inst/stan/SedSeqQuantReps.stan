data {
  // Number of RNAs
  int<lower=1> NRNA;

  // Number of Reps
  int<lower=1> NREP;

  // Note: These are all integers
  // columns t, s, p
  int<lower=0> tot_obs_counts[NRNA,NREP];
  int<lower=0> sup_obs_counts[NRNA,NREP];
  int<lower=0> pel_obs_counts[NRNA,NREP];

  real mean_log_tot_obs[NRNA]; // mean of log-transformed tot_obs_counts across all reps
  real sd_log_tot_obs[NRNA];   // standard deviation of log-transformed tot_obs_counts across all reps


  // Mixing factor priors
  real mixing_factor_total_guess_mean;
  real mixing_factor_total_guess_sd;
  real mixing_factor_sup_guess_mean;
  real mixing_factor_sup_guess_sd;
  real mixing_factor_pellet_guess_mean;
  real mixing_factor_pellet_guess_sd;
}

parameters {
  // Unnormalized mixing proportions
  real<lower=0> mixing_factor_total[NREP];
  real<lower=0> mixing_factor_sup[NREP];
  real<lower=0> mixing_factor_pellet[NREP];

  // dispersion parameter for counts
  // replicate and fraction specific but shared across genes
  real<lower=0> phi[3];

  // adjusted total counts
  real<lower=0> tot_latent_counts[NRNA];

  // log-odds of pSup
  real lopSup[NRNA];
}

transformed parameters {
  real pSup[NRNA];

  // pSup
  for (i in 1:NRNA) {
    pSup[i] = inv_logit( lopSup[i] );
  }
}

model{
  // Cauchy prior for negbin dispersion parameter ()
  phi ~ cauchy(0,3);

  // Gaussian prior for log-odds pSup
  lopSup ~ normal(0,3);

  // Priors for tot_latent_counts and mixing_factor_total
  for (i in 1:NRNA) {
    log(tot_latent_counts[i]) ~ normal(mean_log_tot_obs[i], sd_log_tot_obs[i]);
  }
  mixing_factor_total ~ normal(mixing_factor_total_guess_mean,mixing_factor_total_guess_sd);
  mixing_factor_sup ~ normal(mixing_factor_sup_guess_mean,mixing_factor_sup_guess_sd);
  mixing_factor_pellet ~ normal(mixing_factor_pellet_guess_mean,mixing_factor_pellet_guess_sd);

  for (rep in 1:NREP) {
    for (i in 1:NRNA) {
      sup_obs_counts[i, rep] ~ neg_binomial_2(tot_latent_counts[i] * pSup[i] * mixing_factor_sup[rep], phi[1]);
      pel_obs_counts[i, rep] ~ neg_binomial_2(tot_latent_counts[i] * (1 - pSup[i]) * mixing_factor_pellet[rep], phi[2]);
      tot_obs_counts[i, rep] ~ neg_binomial_2(mixing_factor_total[rep] * tot_latent_counts[i], phi[3]);
    }
  }
}

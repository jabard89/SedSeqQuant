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
}

parameters {
  // Unnormalized mixing proportions
  // only estimate scaling parameter for everything but the first rep
  real<lower=0> scaling_factor_total_raw[NREP > 1 ? NREP - 1 : 0];
  real<lower=0> mixing_factor_sup[NREP];
  real<lower=0> mixing_factor_pellet[NREP];
  
  // dispersion parameter for counts
  // replicate and fraction specific but shared across genes
  real<lower=0> phi[3];
  
  // adjusted fraction counts
  real<lower=0> sup_latent_counts[NRNA];
  real<lower=0> pel_latent_counts[NRNA];
}

transformed parameters {
  real scaling_factor_total[NREP];
  scaling_factor_total[1] = 1; // Fix the first rep's scaling factor to 1
  if (NREP > 1) {
    for (rep in 2:NREP) {
      scaling_factor_total[rep] = scaling_factor_total_raw[rep - 1];
    }
  }
}

model{
  // Cauchy prior for negbin dispersion parameter ()
  phi ~ cauchy(0,3);

  // mixing ratio priors
  if (NREP > 1) {
    for (rep in 2:NREP) {
      scaling_factor_total[rep] ~ gamma(1,1);
    }
  }
  mixing_factor_sup ~ gamma(1,1);
  mixing_factor_pellet ~ gamma(1,1);
    for (rep in 1:NREP) {
      for (i in 1:NRNA) {
        sup_obs_counts[i, rep] ~ neg_binomial_2(sup_latent_counts[i] * mixing_factor_sup[rep] * scaling_factor_total[rep], phi[1]);
        pel_obs_counts[i, rep] ~ neg_binomial_2(pel_latent_counts[i] * mixing_factor_pellet[rep] * scaling_factor_total[rep], phi[2]);
        tot_obs_counts[i, rep] ~ neg_binomial_2(scaling_factor_total[rep] * (sup_latent_counts[i] + pel_latent_counts[i]), phi[3]);
      }
  }
}

generated quantities {
  real pSup[NRNA]; // The derived pSup for each RNA

  for (i in 1:NRNA) {
    pSup[i] = sup_latent_counts[i] / (sup_latent_counts[i] + pel_latent_counts[i]);
  }
}

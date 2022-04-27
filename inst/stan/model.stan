// Piecewise linear modelling of Cycle Thresholds
// 
// Flexible piecewise modelling of Cycle Thresholds with adjustment for 
// varying intercepts and linear scales.
//
// @licence MIT
// @author Sam Abbott
// @author Tim Russell
// @author Joel Hellewell
functions{ 
#include functions/construct_individual_variation.stan
#include functions/piecewise_ct.stan
#include functions/combine_effects.stan
#include functions/onsets_lmpf.stan
}

data {
  // Observations and indexing parameters
  int P; // number of patients
  int N; // number of tests
  array[N] int id; // id of person
  array[P] int tests_per_id; // Tests per ID
  array[P] int cum_tests_per_id; // Cumulative tests per id
  int ncensored; // Number of censored tests
  array[ncensored] int censored; // Which tests have been censored
  int nuncensored; // Number of uncensored tests
  array[nuncensored] int uncensored; // Which tests have not been censored
  array[N] int uncensored_by_test; // Uncensored by test
  vector[N] day_rel; // day of test (integer)
  vector<lower = 0>[N] ct_value; // Ct value of test
  int any_onsets; // Are there any symptom onsets
  int nonsets; // Number of onsets
  vector[P] onset_avail; // Onsets available per ID
  vector[P] onset_time; // Time of onset per ID
  vector[P] onset_window; // Window in which onsets could have occurred per ID
  array[nonsets] int ids_with_onsets; // IDs that have onsets
  // User specified population-level parameters
  int switch; //Should a secondary breakpoint in the CT curve be modelled
  real c_lod; // Ct value at limit of detection (censoring Ct value)
  array[2] real c_thres_p; // Prior on initial/final Ct value
  array[2] real sigma_p; // Prior on the observation error
  array[2] real c_p_p; // intercept of peak Ct (mean + sd)
  array[2] real t_p_p; // intercept of time at peak Ct (mean + sd)
  array[2] real c_s_p; // intercept of switch Ct (mean + sd)
  array[2] real t_s_p; // intercept of time at switch Ct (mean + sd)
  array[2] real t_clear_p; // intercept of time at which cleared (mean + sd)
  array[2] real inc_mean_p; // mean of incubation period (mean + sd)
  array[2] real inc_sd_p; // standard deviation of incubation period (mean + sd)
  // User specified covariate parameters - what should be included
  int adj_t_p; // Should time at peak be adjusted
  int adj_t_s; // Should time at switch be adjusted
  int adj_t_clear; // Should time at LOD be adjusted
  int adj_c_p; // Should CT at peak be adjusted
  int adj_c_s; // Should CT at switch be adjusted
  int adj_inc_mean; // Should incubation period mean be adjusted
  int adj_inc_sd; // Should incubation period standard deviation be adjusted
  // User specified covariate parameters - model design
  int preds; // Number of predictors
  matrix[P, preds + 1] design; //Design matrix
  // User specified covariate parameters - priors
  vector[preds] beta_c_p_m; 
  vector[preds] beta_c_p_sd;
  vector[preds] beta_t_p_m; 
  vector[preds] beta_t_p_sd;
  vector[preds] beta_c_s_m; 
  vector[preds] beta_c_s_sd; 
  vector[preds] beta_t_s_m; 
  vector[preds] beta_t_s_sd; 
  vector[preds] beta_t_clear_m; 
  vector[preds] beta_t_clear_sd; 
  vector[preds] beta_inc_mean_m; 
  vector[preds] beta_inc_mean_sd; 
  vector[preds] beta_inc_sd_m; 
  vector[preds] beta_inc_sd_sd; 
  // Individual-level parameters
  int latent_inf; // Should latent infection time be modelled?
  array[2] real t_inf_p; // Prior on infection time adjustment
  int K; //Number of parameters with individual level variation
  int ind_var_m; // Should inividual variation be modelled
  vector[K] ind_var_mean; // Mean of individual variation
  vector[K] ind_var_sd; // Standard deviation of inividual variation
  int ind_corr; // Should individual variation be modelled with correlation
  real lkj_p; // LKJ prior for individual level variation
  // Ct adjustment model parameters
  int adj_ct; // Should cts be adjusted
  int ct_preds; // Number of predictors for CT adjustment
  matrix[N, ct_preds + 1] ct_design; // Design matrix for CT adjustment
  vector[ct_preds] beta_ct_shift_m; // Mean of Ct shifts effects
  vector[ct_preds] beta_ct_shift_sd; // Standard deviation of Ct shifts effects
  vector[ct_preds] beta_ct_scale_m; // Mean of Ct shifts effects
  vector[ct_preds] beta_ct_scale_sd; // Standard deviation of Ct shifts effects
  // Model control parameters
  int likelihood; // Include log-likelihood
  int output_loglik; // Output the log-likelihood (infection aggregated)
  int pp; // Compute posterior predictions
}

transformed data {
  vector[P] t_inf_bound;
  vector[P] t_inf_mean;
  vector[61] sim_times;
  real c_thres_mean;
  for (i in 1:P) {
    t_inf_bound[i] = max({-onset_time[i], 0});
  }
  t_inf_mean = t_inf_bound + t_inf_p[1];
  c_thres_mean = c_lod + c_thres_p[1];
  for (i in 0:60) {
    sim_times[i + 1] = i;
  }
}

parameters {
  // Population-level parameters
  real<lower = c_lod> c_thres;   // Ct value before and after infection
  real c_p_int; // Intercept of Ct value of viral load at peak
  real t_p_int; // Intercept of time at peak
  array[switch] real c_s_int; // Intercept of Ct value at switch
  array[switch] real t_s_int; // Intercept at time of switch
  array[any_onsets] real inc_mean_int; //Incubation period mean intercept
  array[any_onsets] real<lower = 0> inc_sd_int; //Incubation period sd intercept
  real t_clear_int; // Intercept of the time virus is cleared
  real<lower = 0> sigma; // Variance parameter for oobservation model
  // Individual-level parameters
  // Inferred time of infection
  vector<lower = t_inf_bound>[latent_inf ? P : 0] t_inf; 
  // Cholesky_factored correlation matrix
  cholesky_factor_corr[ind_corr ? K : 0] L_Omega;
  vector<lower = 0>[ind_var_m ? K : 0] ind_var; // SD of individual variation
  matrix[ind_var_m ? K : 0, P] ind_eta; // Individual level variation
  // Covariate Coefficients
  vector[preds && adj_t_p ? preds : 0] beta_t_p;
  vector[preds && adj_t_s ? preds : 0] beta_t_s;
  vector[preds && adj_t_clear ? preds : 0] beta_t_clear;
  vector[preds && adj_c_s ? preds : 0] beta_c_s;
  vector[preds && adj_c_p ? preds : 0] beta_c_p;
  vector[preds && adj_inc_mean ? preds : 0] beta_inc_mean;
  vector[preds && adj_inc_sd ? preds : 0] beta_inc_sd;
  vector[ct_preds && adj_ct ? ct_preds : 0] beta_ct_shift;
  vector[ct_preds && adj_ct ? ct_preds : 0] beta_ct_scale;
}

transformed parameters {
  vector[P] t_p, t_s, t_clear, c_p, c_s, t_clear_abs;
  vector[N] inf_rel, exp_ct, adj_exp_ct;
  array[nonsets] real onsets_star;
  vector[nonsets ? P :0] onsets_log_lik;
{
  matrix[P, K] eta;
  eta = construct_individual_variation(
    ind_eta, ind_var, P, K, ind_corr, L_Omega, ind_var_m
  );
  // Combine effects for each CT parameter and transform to required scale
  t_p = exp(combine_effects(t_p_int, beta_t_p, design) + eta[, 1]);
  t_clear = exp(combine_effects(t_clear_int, beta_t_clear, design) + eta[, 2]);
  c_p = inv_logit(combine_effects(c_p_int, beta_c_p, design) + eta[, 3]);
  // Optional effects if a second breakpoint is used
  if (switch) {
    t_s = exp(combine_effects(t_s_int[1], beta_t_s, design) + eta[, 4]);
    c_s = c_thres * inv_logit(
      combine_effects(c_s_int[1], beta_c_s, design) + eta[, 5]
    );
    c_p = c_s .* c_p;
  }else{
    c_s = rep_vector(0.0, P);
    t_s = rep_vector(0.0, P);
    c_p = c_thres * c_p;
  }
}
  // Make times absolute
  t_clear_abs = t_p + t_s + t_clear;
  // Adjust observed times since first test to be time since infection
  if (latent_inf) {
    inf_rel = day_rel + t_inf[id];
  }else{
    inf_rel = day_rel;
  }
  // Expected ct value given viral load parameters
  exp_ct = piecewise_ct_by_id(
    inf_rel, c_thres, c_p, c_s, c_thres, 0, t_p, t_s, t_clear_abs, id,
    tests_per_id, cum_tests_per_id, switch
  );
  // Shift and scale ct values
  adj_exp_ct = combine_effects(0, beta_ct_shift, ct_design) +
    exp(combine_effects(0, beta_ct_scale, ct_design) + log(exp_ct));
  // Model symptom onset likelihood: see onsets_lmpf.stan
  if (any_onsets) {
    vector[P] onsets_ttar;
    onsets_ttar = onsets_lmpf(
      inc_mean_int[1], inc_sd_int[1], beta_inc_mean, beta_inc_sd, design,
      onset_avail, onset_time, onset_window, t_inf, ids_with_onsets
    );
    onsets_star[1] = sum(onsets_ttar);
    if (output_loglik) {
      onsets_log_lik = onsets_ttar;
    }
  }
}

model {
  // Prior over possible infection times relative to first
  // positive test or symtom onset.
  // Assumes that the first positive test is not a false positive.
  if (latent_inf) {
    t_inf ~ normal(t_inf_mean, t_inf_p[2]); 
  }
  // CT piecewise linear intercept parameters
  c_thres ~ normal(c_thres_mean, c_thres_p[2]) T[c_lod, ];
  c_p_int ~ normal(c_p_p[1], c_p_p[2]); // Mean at 50% of switch value
  t_p_int ~ normal(t_p_p[1], t_p_p[2]); // Mean at log(5)
  //mean at log(10) + peak + scale timing 
  t_clear_int ~ normal(t_clear_p[1], t_clear_p[2]); 
  if (switch) {
    c_s_int ~ normal(c_s_p[1], c_s_p[2]); //mean at 50% of maximum ct
    t_s_int ~ normal(t_s_p[1], t_s_p[2]); //mean at log(5) + peak timing
  }
  // Individual level variation
  if (ind_var_m) {
    to_vector(ind_eta) ~ std_normal();
    ind_var ~ normal(ind_var_mean, ind_var_sd);
  }
  // LKJ prior on correlation between individual level dynamics
  if (ind_corr) {
    L_Omega ~ lkj_corr_cholesky(lkj_p);
  }
  // Variation in observation model
  sigma ~ normal(sigma_p[1], sigma_p[2]) T[0,];
  // Coefficients priors for predictors
  if (preds) {
    if (adj_t_p) {
      beta_t_p ~ normal(beta_t_p_m, beta_t_p_sd);
    }
    if (adj_t_s) {
      beta_t_s ~ normal(beta_t_s_m, beta_t_s_sd);
    }
    if (adj_t_clear) {
      beta_t_clear ~ normal(beta_t_clear_m, beta_t_clear_sd);
    }
    if (adj_c_p) {
      beta_c_p ~ normal(beta_c_p_m, beta_c_p_sd);
    }
    if (adj_c_s) {
      beta_c_s ~ normal(beta_c_s_m, beta_c_s_sd);
    }
    if (adj_inc_mean) {
      beta_inc_mean ~ normal(beta_inc_mean_m, beta_inc_mean_sd);
    } 
    if (adj_inc_sd) {
      beta_inc_sd ~ normal(beta_inc_sd_m, beta_inc_sd_sd);
    } 
  }
  if (ct_preds && adj_ct) {
    beta_ct_shift ~ normal(beta_ct_shift_m, beta_ct_shift_sd);
    beta_ct_scale ~ normal(beta_ct_scale_m, beta_ct_scale_sd);
  }
  if (any_onsets) {
    // Priors on the incubation period
    inc_mean_int ~ normal(inc_mean_p[1], inc_mean_p[2]);
    inc_sd_int[1] ~ normal(inc_sd_p[1], inc_sd_p[2]) T[0, ];
    if (likelihood) {
      // Component of likelihood for symptom onsets see onsets_lpmf.stan
      target += onsets_star[1];
    }
  }
  if (likelihood) {
    // Component of likelihood for expected ct values
    // If non-censored: P(observed ct | expected ct)
    ct_value[uncensored] ~ normal(adj_exp_ct[uncensored], sigma);
    // If censored: P(expected ct >= censored ct)
    target += normal_lccdf(c_lod | adj_exp_ct[censored], sigma);
    // All CTs are truncated above 0. P(0 <= expected ct)
    target += -normal_lccdf(0 | adj_exp_ct, sigma);
  }
}

generated quantities {
  vector[N] sim_ct;
  matrix[ind_corr ? K : 0, ind_corr ? K : 0] correlation;
  vector[output_loglik ? P : 0] log_lik;
  if (ind_corr) {
    correlation = L_Omega * L_Omega';
  }
  if (pp) {
    // Posterior predictions
    sim_ct = to_vector(normal_rng(adj_exp_ct, sigma));
    sim_ct = fmin(sim_ct, c_lod);
  }
  // Output by infection log-likelihood
  if (output_loglik) {
    log_lik = rep_vector(0, P);
    if (nonsets) {
      log_lik = onsets_log_lik;
    }
    for (i in 1:P) {
      int t_start = cum_tests_per_id[i] - tests_per_id[i] + 1;
      int t_end = cum_tests_per_id[i];
      for (j in t_start:t_end) {
        if (uncensored_by_test[j] == 1) {
          log_lik[i] += normal_lpdf(ct_value[j] | adj_exp_ct[j], sigma);
        }else{
          log_lik[i] += normal_lccdf(c_lod | adj_exp_ct[j], sigma);
        }
      }
      log_lik[i] += -normal_lccdf(0 | adj_exp_ct[t_start:t_end], sigma);
    }
  }
}

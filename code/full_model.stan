
data {
  // constants ------------------------------------
  int<lower=1> nyear; // number of years
  int<lower=1> ncountry; // number of countries
  int<lower=2> ncat;  // number of main causes categories

  // for model for main-cause observations ------------------
  int<lower=1> nregion; // number of modelling regions
  
  int<lower=1> N;  // number of observations
  array[N] int<lower=1, upper=ncountry> c; // country labels for observations
  array[ncountry] int<lower=1> region; // region labels for the countries
  
  array[N, ncat] int Y;  // response array
  array[N] int tot; // rowwise sums for observations
  
  int<lower=1> K; // number of population level effects (covariates)
  matrix[N, K] X; // design matrix
  matrix[N, ncat] zero_mat; // matrix for zeroes/NAs

  vector<lower=0, upper=1>[N] not_type_1; // boolean for whether observation is highest category
  vector<lower=0, upper=1>[N] quality;
  
  real<lower=0> delta_scale;
  
  // for DIR subcauses ---------------------------
  int<lower=2> ncat_subDIR;
  int<lower=1> N_subDIR; // number of observations with DIR subcauses
  array[N_subDIR] int<lower=1> c_subDIR;
  
  array[N_subDIR, ncat_subDIR] int Y_subDIR;
  matrix[N_subDIR, K] X_subDIR;
  
  matrix[N_subDIR, ncat_subDIR] zero_mat_subDIR;
  
  int<lower=1, upper=ncat> index_DIR;
  
  // for HEM subcauses ----------------------------
  int<lower=2> ncat_subHEM;
  int<lower=1> N_subHEM; // number of observations with HEM subcauses
  array[N_subHEM] int<lower=1> c_subHEM;
  
  array[N_subHEM, ncat_subHEM] int Y_subHEM;
  matrix[N_subHEM, K] X_subHEM;

  matrix[N_subHEM, ncat_subHEM] zero_mat_subHEM;
  
  int<lower=1, upper=ncat> index_HEM;
  
  // for SEP subcauses ----------------------------
  int<lower=2> ncat_subSEP;
  int<lower=1> N_subSEP; // number of observations with SEP subcauses
  array[N_subSEP] int<lower=1> c_subSEP;
  
  array[N_subSEP, ncat_subSEP] int Y_subSEP;
  matrix[N_subSEP, K] X_subSEP;
    
  matrix[N_subSEP, ncat_subSEP] zero_mat_subSEP;
  
  int<lower=1, upper=ncat> index_SEP;
  
  // for generated quantities ------------------
  int<lower=N> N_complete; // for prediction (number of all possible country-years)
  matrix[N_complete, K] X_complete; // model matrix for prediction
  array[N_complete] int<lower=1> c_complete; // country labels for prediction matrix
  vector<lower=0, upper=1>[N_complete] w; // weights based on coverage
  
  vector<lower=0>[N_complete] matdeaths; // number of maternal deaths in each country-year
  vector<lower=0>[N_complete] aidsdeaths; // estimated AIDS deaths in each country-year
  
  int<lower=1> nregion_sdg;
  
  matrix[ncountry, N_complete] country_mat;
  matrix[nregion_sdg, N_complete] region_sdg_mat;

}

transformed data {
}

parameters {
  // global intercept
  row_vector[ncat-1] alpha; 
  
  // region level effects
  matrix[K, ncat-1] beta_raw; 
  vector<lower=0>[ncat-1] sd_beta;
  
  // for country-level effects
  array[ncountry, ncat-1] real gamma; 
  vector<lower=0>[ncat-1] tau;
  cholesky_factor_corr[ncat-1] L_Omega;
  
  // data quality terms
  matrix[N, ncat-1] delta_raw;
  vector<lower=0>[ncat-1] sd_delta;
  
  // extra poisson parameter
  array[N] real xi; 

  // for DIR subcauses ---------------------------
  row_vector[ncat_subDIR-1] alpha_subDIR;
  matrix[K, ncat_subDIR-1] beta_subDIR_raw;
  real<lower=0> sd_beta_subDIR;
  array[ncountry, ncat_subDIR-1] real gamma_subDIR; 
  vector<lower=0>[ncat_subDIR-1] tau_subDIR;
  cholesky_factor_corr[ncat_subDIR-1] L_Omega_subDIR;
  array[N_subDIR] real xi_subDIR;
  
  // for HEM subcauses ---------------------------
  row_vector[ncat_subHEM-1] alpha_subHEM;
  matrix[K, ncat_subHEM-1] beta_subHEM_raw;
  real<lower=0> sd_beta_subHEM;
  array[ncountry, ncat_subHEM-1] real gamma_subHEM; 
  vector<lower=0>[ncat_subHEM-1] tau_subHEM;
  cholesky_factor_corr[ncat_subHEM-1] L_Omega_subHEM;
  array[N_subHEM] real xi_subHEM;

  // for SEP subcauses ---------------------------
  row_vector[ncat_subSEP-1] alpha_subSEP;
  matrix[K, ncat_subSEP-1] beta_subSEP_raw;
  real<lower=0> sd_beta_subSEP;
  array[ncountry, ncat_subSEP-1] real gamma_subSEP;
  vector<lower=0>[ncat_subSEP-1] tau_subSEP;
  cholesky_factor_corr[ncat_subSEP-1] L_Omega_subSEP;
  array[N_subSEP] real xi_subSEP;
}

transformed parameters {
  ////////////////// Declarations ///////////////////////////
  // model for main causes ---------------------------------
  array[N, ncat] real mu;  // log ratios
  matrix[ncat-1, ncat-1] L; // variance-covariance matrix
  
  matrix[K, ncat-1] beta; 
  array[N, ncat-1] real gamma_data; // country-level effect populated as array for each observation-cause
  matrix[N, ncat-1] delta;
    
  // for DIR subcauses -------------------------------------------
  array[N_subDIR, ncat_subDIR] real mu_subDIR;  
  matrix[K, ncat_subDIR-1] beta_subDIR;
  array[N_subDIR, ncat_subDIR-1] real gamma_subDIR_data; 
  matrix[ncat_subDIR-1, ncat_subDIR-1] L_subDIR;
  
  // for HEM subcauses -------------------------------------------
  array[N_subHEM, ncat_subHEM] real mu_subHEM;  
  matrix[K, ncat_subHEM-1] beta_subHEM;
  array[N_subHEM, ncat_subHEM-1] real gamma_subHEM_data;
  matrix[ncat_subHEM-1, ncat_subHEM-1] L_subHEM;
  
  // for SEP subcauses -------------------------------------------
  array[N_subSEP, ncat_subSEP] real mu_subSEP;   
  matrix[K, ncat_subHEM-1] beta_subSEP;
  array[N_subSEP, ncat_subSEP-1] real gamma_subSEP_data; 
  matrix[ncat_subSEP-1, ncat_subSEP-1] L_subSEP;
  
  //////////////////// Calculations //////////////////////////
  
  // model for main causes -----------------------------------
  
  beta = diag_post_multiply(beta_raw, sd_beta);
  
  // populate country random effect matrix 
  for (n in 1:N) {
    gamma_data[n] = gamma[c[n]];
  }
  
  L =  diag_pre_multiply(tau, L_Omega);
  
  delta = diag_post_multiply(diag_pre_multiply(not_type_1 .* quality, delta_raw), sd_delta);

  // calculate log ratios
  mu = to_array_2d(append_col(
      rep_matrix(alpha, N) + 
      X * beta + 
      to_matrix(gamma_data) + 
      delta, 
    rep_vector(0, N)));
    
  // for DIR subcauses ----------------------------------------
  beta_subDIR = sd_beta_subDIR*beta_subDIR_raw;
  
  // populate country random effect matrix 
  for (n in 1:N_subDIR) {
    gamma_subDIR_data[n] = gamma_subDIR[c_subDIR[n]];
  }
  
  // variance-covariance
  L_subDIR =  diag_pre_multiply(tau_subDIR, L_Omega_subDIR);
  
  mu_subDIR = to_array_2d(append_col(
      rep_matrix(alpha_subDIR, N_subDIR) + 
      X_subDIR * beta_subDIR + 
      to_matrix(gamma_subDIR_data), 
    rep_vector(0, N_subDIR)));
  
  // for HEM subcauses ----------------------------------------
  beta_subHEM = sd_beta_subHEM*beta_subHEM_raw;
    
  // populate country random effect matrix 
  for (n in 1:N_subHEM) {
    gamma_subHEM_data[n] = gamma_subHEM[c_subHEM[n]];
  }
  
  // variance-covariance
  L_subHEM =  diag_pre_multiply(tau_subHEM, L_Omega_subHEM);
  
  mu_subHEM = to_array_2d(append_col(
      rep_matrix(alpha_subHEM, N_subHEM) + 
      X_subHEM * beta_subHEM + 
      to_matrix(gamma_subHEM_data), 
    rep_vector(0, N_subHEM)));
  
  // for SEP subcauses ----------------------------------------
  beta_subSEP = sd_beta_subSEP*beta_subSEP_raw;
  
  // population country random effect matrix
  for (n in 1:N_subSEP) {
    gamma_subSEP_data[n] = gamma_subSEP[c_subSEP[n]];
  }
  
  // variance-covariance
  L_subSEP =  diag_pre_multiply(tau_subSEP, L_Omega_subSEP);
  
  mu_subSEP = to_array_2d(append_col(
      rep_matrix(alpha_subSEP, N_subSEP) + 
      X_subSEP * beta_subSEP + 
      to_matrix(gamma_subSEP_data), 
    rep_vector(0, N_subSEP)));
}

model {
  matrix[N, ncat] lik_mat;
  matrix[N_subDIR, ncat_subDIR] lik_mat_subDIR;
  matrix[N_subHEM, ncat_subHEM] lik_mat_subHEM;
  matrix[N_subSEP, ncat_subSEP] lik_mat_subSEP;

  // model for main causes --------------------------------------
  // Multinomial likelihood: 
  for (n in 1:N) {
    for (j in 1:ncat) {
      // Store likelihoods for each cell
      lik_mat[n, j] = poisson_log_lpmf(Y[n,j] | mu[n, j] + xi[n]);
    }
  }

  // only non-zero values contribute to the likelihood
  target += sum(zero_mat .* lik_mat);

  // intercepts
  alpha ~ normal(0, 3); 

  // region effects
  for(k in 1:K) {
    beta_raw[k, ] ~ std_normal(); 
  }
  
  // priors for beta 
  sd_beta ~ normal(0, 1);

  // country-level random effects
  for (country in 1:ncountry) {
    to_vector(gamma[country]) ~ multi_normal_cholesky(rep_vector(0, ncat-1), L);
  }
  
  tau ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(1);
  
  // noise terms
  for (j in 1:(ncat-1)) {
    delta_raw[, j] ~ std_normal();
  }
  
  sd_delta ~ normal(0, delta_scale);
  
  // for DIR subcauses --------------------------------------
  // Multinomial likelihood: 
  for (n in 1:N_subDIR) {
    for (j in 1:ncat_subDIR) {
      // Store likelihoods for each cell
      lik_mat_subDIR[n, j] = poisson_log_lpmf(Y_subDIR[n,j] | mu_subDIR[n, j] + xi_subDIR[n]);
    }
  }

  // only non-zero values contribute to the likelihood
  target += sum(zero_mat_subDIR .* lik_mat_subDIR);

  // intercepts:
  alpha_subDIR ~ normal(0, 3); 

  // region effects:
  for(k in 1:K) {
    beta_subDIR_raw[k, ] ~ std_normal();  
  }

  sd_beta_subDIR ~ normal(0, 1);

  // country-level random effects:
  for (country in 1:ncountry) {
    to_vector(gamma_subDIR[country]) ~ multi_normal_cholesky(rep_vector(0, ncat_subDIR-1), L_Omega_subDIR);
  }
  
  tau_subDIR ~ normal(0, 1);
  L_Omega_subDIR ~ lkj_corr_cholesky(1);
  
  // for HEM subcauses --------------------------------------
  // Multinomial likelihood: 
  for (n in 1:N_subHEM) {
    for (j in 1:ncat_subHEM) {
      // Store likelihoods for each cell
      lik_mat_subHEM[n, j] = poisson_log_lpmf(Y_subHEM[n,j] | mu_subHEM[n, j] + xi_subHEM[n]);
    }
  }

  // only non-zero values contribute to the likelihood
  target += sum(zero_mat_subHEM .* lik_mat_subHEM);

  // intercepts
  alpha_subHEM ~ normal(0, 3); 

  // region effects
  for(k in 1:K) {
    beta_subHEM_raw[k, ] ~ std_normal(); 
  }
  
  // priors for beta 
  sd_beta_subHEM ~ normal(0, 1);

  // country-level random effects:
  for (country in 1:ncountry) {
    to_vector(gamma_subHEM[country]) ~ multi_normal_cholesky(rep_vector(0, ncat_subHEM-1), L_Omega_subHEM);
  }
  
  tau_subHEM ~ normal(0, 1);
  L_Omega_subHEM ~ lkj_corr_cholesky(1);

  // for SEP subcauses --------------------------------------
  // Multinomial likelihood: 
  for (n in 1:N_subSEP) {
    for (j in 1:ncat_subSEP) {
      // Store likelihoods for each cell
      lik_mat_subSEP[n, j] = poisson_log_lpmf(Y_subSEP[n,j] | mu_subSEP[n, j] + xi_subSEP[n]);
    }
  }

  // only non-zero values contribute to the likelihood
  target += sum(zero_mat_subSEP .* lik_mat_subSEP);

  // intercepts
  alpha_subSEP ~ normal(0, 3); 

  // region effects
  for(k in 1:K) {
    beta_subSEP_raw[k, ] ~ std_normal(); 
  }
  
  // priors for beta 
  sd_beta_subSEP ~ normal(0, 1);

  // country-level random effects:
  for (country in 1:ncountry) {
    to_vector(gamma_subSEP[country]) ~ multi_normal_cholesky(rep_vector(0, ncat_subSEP-1), L_Omega_subSEP);
  }
  
  tau_subSEP ~ normal(0, 1);
  L_Omega_subSEP ~ lkj_corr_cholesky(1);
}

generated quantities {
  ////////////////// Declarations ///////////////////////////
  // covariance matrix
  matrix[ncat-1, ncat-1] Sigma;
  
  // diagnostic stuff
  matrix[N, ncat] lik_mat;
  matrix[N, ncat] lik_mat_0;
  vector[N] log_lik;
  array[N, ncat] int<lower=0> Y_rep;
  
  // for main distribution proportion estimates
  matrix[N_complete, ncat] count_complete;
  matrix[ncountry, ncat] count_country;
  matrix[nregion_sdg, ncat] count_region_sdg;
  row_vector[ncat] count_global;
  matrix[N_complete, ncat] p_complete;
  matrix[ncountry, ncat] p_country;
  matrix[nregion_sdg, ncat] p_region_sdg;
  row_vector[ncat] p_global;

  // DIR subcauses
  matrix[N_complete, ncat_subDIR] p_subDIR_complete;
  matrix[nregion_sdg, ncat_subDIR] p_subDIR_region_sdg;
  matrix[nregion_sdg, ncat_subDIR] count_subDIR_region_sdg;
  row_vector[ncat_subDIR] p_subDIR_global;
  row_vector[ncat_subDIR] count_subDIR_global;
  matrix[ncat_subDIR-1, ncat_subDIR-1] Sigma_subDIR;
  
  // HEM subcauses
  matrix[N_complete, ncat_subHEM] p_subHEM_complete;
  matrix[nregion_sdg, ncat_subHEM] p_subHEM_region_sdg;
  matrix[nregion_sdg, ncat_subHEM] count_subHEM_region_sdg;
  row_vector[ncat_subHEM] p_subHEM_global;
  row_vector[ncat_subHEM] count_subHEM_global;
  matrix[ncat_subHEM-1, ncat_subHEM-1] Sigma_subHEM;
  
  // SEP subcauses
  matrix[N_complete, ncat_subSEP] p_subSEP_complete;
  matrix[nregion_sdg, ncat_subSEP] p_subSEP_region_sdg;
  matrix[nregion_sdg, ncat_subSEP] count_subSEP_region_sdg;
  row_vector[ncat_subSEP] p_subSEP_global;
  row_vector[ncat_subSEP] count_subSEP_global;
  matrix[ncat_subSEP-1, ncat_subSEP-1] Sigma_subSEP;

  ////////////////// Calculations ///////////////////////////
  // reconstruct variance-covariance matrices
  Sigma = multiply_lower_tri_self_transpose(L);
  Sigma_subDIR = multiply_lower_tri_self_transpose(L_subDIR);
  Sigma_subHEM = multiply_lower_tri_self_transpose(L_subHEM);
  Sigma_subSEP = multiply_lower_tri_self_transpose(L_subSEP);
  
  // for diagnostics ---------------------------------------------
  for (n in 1:N) {
    for (j in 1:ncat) {
    // store likelihoods (even zeroes)
      lik_mat[n, j] = poisson_log_lpmf(Y[n,j] | mu[n,j] + xi[n]);
    }
  }
  
  lik_mat_0 = zero_mat .* lik_mat;
  
  for (i in 1:N) {
    log_lik[i] = sum(lik_mat_0[i]);
  }
  
  for (n in 1:N) {
    Y_rep[n] =  multinomial_rng(softmax(to_vector(mu[n])), tot[n]);  
  }
  
  // main cause distribution -------
  {
    array[N_complete, ncat-1] real gamma_complete;
    array[ncountry, ncat-1] real gamma_unobs;
    array[N_complete, ncat-1] real gamma_unobs_complete;
    matrix[N_complete, ncat] mu_complete;

    
    // observed portion
    for (n in 1:N_complete) {
      gamma_complete[n] = gamma[c_complete[n]];
    }
    
    // unobserved portion
    for (country in 1:ncountry) {
      gamma_unobs[country] = to_array_1d(multi_normal_cholesky_rng(rep_vector(0, ncat-1), L));
    }
    
    for (n in 1:N_complete) {
      gamma_unobs_complete[n] = gamma_unobs[c_complete[n]];
    }

    mu_complete = append_col(
        rep_matrix(alpha, N_complete) +
        X_complete*beta +
        diag_pre_multiply(w, to_matrix(gamma_complete)) + diag_pre_multiply(1-w, to_matrix(gamma_unobs_complete)), 
      rep_vector(0, N_complete)
    );
      
    for (i in 1:N_complete) {
      p_complete[i] = to_row_vector(softmax(to_vector(mu_complete[i])));
    }

     // convert each observation into counts
    count_complete = diag_pre_multiply(matdeaths, p_complete);
    
    // add in HIV/AIDS deaths
    count_complete[, 6] = count_complete[, 6] + aidsdeaths;
    
    // re-calculate proportions after AIDS deaths
    for (n in 1:N_complete) {
      p_complete[n] = count_complete[n] / sum(count_complete[n]);
    }
    
    count_country = country_mat * count_complete;
    for (i in 1:ncountry) {
      p_country[i] = count_country[i] ./ sum(count_country[i]);
    }
    
    // sum counts by region
    count_region_sdg = region_sdg_mat * count_complete;
    
    // turn region counts into proportions
    for (r in 1:nregion_sdg) {
      p_region_sdg[r] = count_region_sdg[r] ./ sum(count_region_sdg[r]);
    }
  
    count_global = rep_row_vector(1, nregion_sdg) * count_region_sdg;
    p_global = count_global ./sum(count_global);
    
  }

  // DIR subcauses ------ 
  {
    array[N_complete, ncat_subDIR-1] real gamma_subDIR_complete;
    array[ncountry, ncat_subDIR-1] real gamma_subDIR_unobs;
    array[N_complete, ncat_subDIR-1] real gamma_subDIR_unobs_complete;
    matrix[N_complete, ncat_subDIR] mu_subDIR_complete;
    matrix[N_complete, ncat_subDIR] count_subDIR_complete;
  
        // observed portion
    for (n in 1:N_complete) {
      gamma_subDIR_complete[n] = gamma_subDIR[c_complete[n]];
    }
    
    // unobserved portion
    for (country in 1:ncountry) {
      gamma_subDIR_unobs[country] = to_array_1d(multi_normal_cholesky_rng(rep_vector(0, ncat_subDIR-1), L_subDIR));
    }
    
    for (n in 1:N_complete) {
      gamma_subDIR_unobs_complete[n] = gamma_subDIR_unobs[c_complete[n]];
    }

    mu_subDIR_complete = append_col(
        rep_matrix(alpha_subDIR, N_complete) +
        X_complete*beta_subDIR +
        diag_pre_multiply(w, to_matrix(gamma_subDIR_complete)) + diag_pre_multiply(1-w, to_matrix(gamma_subDIR_unobs_complete)), 
      rep_vector(0, N_complete)
    );
      
    for (i in 1:N_complete) {
      p_subDIR_complete[i] = to_row_vector(softmax(to_vector(mu_subDIR_complete[i])));
    }

    count_subDIR_complete = diag_pre_multiply(count_complete[, index_DIR], p_subDIR_complete);
    count_subDIR_region_sdg = region_sdg_mat * count_subDIR_complete;
    
    for (r in 1:nregion_sdg) {
      p_subDIR_region_sdg[r] = count_subDIR_region_sdg[r] ./ sum(count_region_sdg[r]);
    }
  
    count_subDIR_global = rep_row_vector(1, nregion_sdg) * count_subDIR_region_sdg;
    p_subDIR_global = count_subDIR_global ./sum(count_global);
    
  }

  // HEM subcauses ----------
  {
    array[N_complete, ncat_subHEM-1] real gamma_subHEM_complete;
    array[ncountry, ncat_subHEM-1] real gamma_subHEM_unobs;
    array[N_complete, ncat_subHEM-1] real gamma_subHEM_unobs_complete;
    matrix[N_complete, ncat_subHEM] mu_subHEM_complete;
    matrix[N_complete, ncat_subHEM] count_subHEM_complete;
    
        // observed portion
    for (n in 1:N_complete) {
      gamma_subHEM_complete[n] = gamma_subHEM[c_complete[n]];
    }
    
    // unobserved portion
    for (country in 1:ncountry) {
      gamma_subHEM_unobs[country] = to_array_1d(multi_normal_cholesky_rng(rep_vector(0, ncat_subHEM-1), L_subHEM));
    }
    
    for (n in 1:N_complete) {
      gamma_subHEM_unobs_complete[n] = gamma_subHEM_unobs[c_complete[n]];
    }

    mu_subHEM_complete = append_col(
        rep_matrix(alpha_subHEM, N_complete) +
        X_complete*beta_subHEM +
        diag_pre_multiply(w, to_matrix(gamma_subHEM_complete)) + diag_pre_multiply(1-w, to_matrix(gamma_subHEM_unobs_complete)), 
      rep_vector(0, N_complete)
    );
      
    for (i in 1:N_complete) {
      p_subHEM_complete[i] = to_row_vector(softmax(to_vector(mu_subHEM_complete[i])));
    }

    count_subHEM_complete = diag_pre_multiply(count_complete[, index_HEM], p_subHEM_complete);
    count_subHEM_region_sdg = region_sdg_mat * count_subHEM_complete;
    
    for (r in 1:nregion_sdg) {
      p_subHEM_region_sdg[r] = count_subHEM_region_sdg[r] ./ sum(count_region_sdg[r]);
    }
  
    count_subHEM_global = rep_row_vector(1, nregion_sdg) * count_subHEM_region_sdg;
    p_subHEM_global = count_subHEM_global ./sum(count_global);
    
  }
  
  // SEP subcauses ----------
  {
    array[N_complete, ncat_subSEP-1] real gamma_subSEP_complete;
    array[ncountry, ncat_subSEP-1] real gamma_subSEP_unobs;
    array[N_complete, ncat_subSEP-1] real gamma_subSEP_unobs_complete;
    matrix[N_complete, ncat_subSEP] mu_subSEP_complete;
    matrix[N_complete, ncat_subSEP] count_subSEP_complete;
    
        // observed portion
    for (n in 1:N_complete) {
      gamma_subSEP_complete[n] = gamma_subSEP[c_complete[n]];
    }
    
    // unobserved portion
    for (country in 1:ncountry) {
      gamma_subSEP_unobs[country] = to_array_1d(multi_normal_cholesky_rng(rep_vector(0, ncat_subSEP-1), L_subSEP));
    }
    
    for (n in 1:N_complete) {
      gamma_subSEP_unobs_complete[n] = gamma_subSEP_unobs[c_complete[n]];
    }

    mu_subSEP_complete = append_col(
        rep_matrix(alpha_subSEP, N_complete) +
        X_complete*beta_subSEP +
        diag_pre_multiply(w, to_matrix(gamma_subSEP_complete)) + diag_pre_multiply(1-w, to_matrix(gamma_subSEP_unobs_complete)), 
      rep_vector(0, N_complete)
    );
      
    for (i in 1:N_complete) {
      p_subSEP_complete[i] = to_row_vector(softmax(to_vector(mu_subSEP_complete[i])));
    }

    count_subSEP_complete = diag_pre_multiply(count_complete[, index_SEP], p_subSEP_complete);
    count_subSEP_region_sdg = region_sdg_mat * count_subSEP_complete;
    
    for (r in 1:nregion_sdg) {
      p_subSEP_region_sdg[r] = count_subSEP_region_sdg[r] ./ sum(count_region_sdg[r]);
    }
  
    count_subSEP_global = rep_row_vector(1, nregion_sdg) * count_subSEP_region_sdg;
    p_subSEP_global = count_subSEP_global ./sum(count_global);
    
  }
  
}

data {
  int k; // traits
  int m; // taxa
  cov_matrix[m] C; // phylo
  int ni[m]; // amostras
  int ni_max; // máximo amostra
  vector[k] X[m,ni_max]; // dados
  vector[k] priorX;
  cov_matrix[k] priorS;
}

transformed data {
  vector[k] zero_vector;
  real ldetC;
  cov_matrix[m] invC;
  
  for (i in 1:k)
    zero_vector[i] <- 0;

  ldetC <- log_determinant(C);
  invC <- inverse_spd(C);
 
}

parameters {
  vector[k] Xbar[m]; 
  vector[k] alpha;
  
  corr_matrix[k] Gamma; // Sigma = sGs
  real<lower=0> sigma[k];

  corr_matrix[k] Gamma_bm; 
  real<lower=0> sigma_bm[k];

}

transformed parameters {
  vector[k] eX[m,ni_max];
  matrix[m,k] eXbar;
  
  for(i in 1:m)
    {
      eXbar[i] <- to_row_vector ((Xbar[i] - alpha) ./ to_vector(sigma_bm));
      for(j in 1:ni[i])
	eX[i][j] <- (X[i][j] - Xbar[i]) ./ to_vector(sigma);
    }
}

model {

  real ldet_BM;
  real llik_BM; 
    
  /** priors **/
  
  for (i in 1:m)
    Xbar[i] ~ multi_normal(priorX, priorS);
  
  alpha ~ multi_normal(priorX, priorS);
  
  // lambda ~ uniform(0, 1);
  
  for (i in 1:k)
    {
      sigma[i] ~ chi_square(min(ni) - 1);
      sigma_bm[i] ~ chi_square(min(ni) - 1);
    }

  // Gamma ~ ldk_corr(1);
  // Gamma_bm ~ ldk_corr(1);

  /** brownian **/

  ldet_BM <- log_determinant(Gamma_bm) 
  ldet_BM <- ldet_BM + (2 * sum(sigma_bm));
  // jacobian of transform to correlation matrix

  llik_BM  <- -0.5 * (trace_gen_quad_form(inverse_spd(Gamma_bm), C, eXbar) +
		      k * ldetC + m * ldet_BM);

  increment_log_prob(llik_BM);
  
  /** pops **/

  for (i in 1:m)
    for (j in 1:(ni[i]))
      eX[i][j] ~ multi_normal(zero_vector, Gamma);

  increment_log_prob(- sum(sigma)); 
  // jacobian of transform to correlation matrix
}

generated quantities {
  cov_matrix[k] Sigma;
  cov_matrix[k] Sigma_bm;
  
  Sigma <- quad_form(Gamma, diag_matrix(to_vector(sigma)));
  Sigma_bm <- quad_form(Gamma_bm, diag_matrix(to_vector(sigma_bm)));
}
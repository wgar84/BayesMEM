data {
  int k; // traits
  int m; // taxa
  cov_matrix[2 * (m-1)] C; // phylo
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
  vector[k] Xbar[2 * (m-1)]; 
  vector[k] alpha;
  
  cholesky_factor_corr[k] Gamma; // Sigma = sGG's
  real<lower=0> sigma[k];

  cholesky_factor_corr[k] Gamma_bm; 
  real<lower=0> sigma_bm[k];

}

transformed parameters {
  vector[k] eX[m,ni_max];
  matrix[2 * (m-1),k] eXbar;
  
  for(i in 1:(2 * (m-1)))
    eXbar[i] <- to_row_vector ((Xbar[i] - alpha) ./ to_vector(sigma_bm));
      
  for(i in 1:m)
    for(j in 1:ni[i])
      eX[i,j] <- (X[i,j] - Xbar[i]) ./ to_vector(sigma);
}

model {

  real ldet_BM;
  real llik_BM; 
    
  /** priors **/
  
  for (i in 1:(2 * (m-1)))
    Xbar[i] ~ multi_normal(priorX, priorS);
  
  alpha ~ multi_normal(priorX, priorS);
  
  // lambda ~ uniform(0, 1);
  
  for (i in 1:k)
    {
      sigma[i] ~ chi_square(min(ni) - 1);
      sigma_bm[i] ~ chi_square(min(ni) - 1);
    }

  // Gamma ~ ldk_corr_chol(1);
  // Gamma_bm ~ ldk_corr_chol(1);

  /** brownian **/

  ldet_BM <- 2 * log_determinant(Gamma_bm); 
  // detS = det(GG) = det(G)det(G) = det(G)^2 (e tira log) 
  ldet_BM <- ldet_BM + (2 * sum(sigma_bm));
  // jacobian of transform to correlation matrix

  llik_BM  <- - 0.5 * (trace_gen_quad_form(inverse_spd(Gamma_bm * Gamma_bm'), C, eXbar) +
		      k * ldetC + m * ldet_BM);

  increment_log_prob(llik_BM);
  
  /** pops **/

  for (i in 1:m)
    for (j in 1:(ni[i]))
      eX[i,j] ~ multi_normal_cholesky(zero_vector, Gamma);

  increment_log_prob(- k * sum(sigma)); 
  // jacobian of transform to correlation matrix
}

generated quantities {
  cov_matrix[k] Sigma;
  cov_matrix[k] Sigma_bm;
  
  Sigma <- quad_form(Gamma * Gamma', diag_matrix(to_vector(sigma)));
  Sigma_bm <- quad_form(Gamma_bm * Gamma_bm', diag_matrix(to_vector(sigma_bm)));
}
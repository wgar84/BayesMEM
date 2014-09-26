data {
  int k; // traits
  int m; // taxa
  cov_matrix[2*m - 2] C; // phylo
  vector[k] Xbar[m]; // médias
  vector[k] priorX;
  cov_matrix[k] priorS;
}

transformed data {
  vector[k] zero_vector;
  real ldetC;
  cov_matrix[2*m - 2] invC;
  
  for (i in 1:k)
    zero_vector[i] <- 0;

  ldetC <- log_determinant(C);
  invC <- inverse_spd(C);
 
}

parameters {
  vector[k] ancestor[m - 2]; 
  vector[k] root;
  
  cholesky_factor_corr[k] Gamma_bm; 
  real<lower=0> sigma_bm[k];

}

transformed parameters {
  matrix[2 * m - 2,k] eXbar;
  
  for(i in 1:m)
    eXbar[i] <- to_row_vector ((Xbar[i] - root) ./ to_vector(sigma_bm));
    
  for(i in 1:(m-2))
    eXbar[i+m] <- to_row_vector ((ancestor[i] - root) ./ to_vector(sigma_bm));
}

model {

  real ldet_BM;
  real llik_BM; 
    
  /** priors **/
  
  for (i in 1:(m - 2))
    ancestor[i] ~ multi_normal(priorX, (10 ^ 3) * priorS);
  
  root ~ multi_normal(priorX, (10 ^ 3) * priorS);
  
  for (i in 1:k)
    sigma_bm[i] ~ chi_square(m - 1);
  
  // Gamma ~ ldk_corr_chol(1);
  // Gamma_bm ~ ldk_corr_chol(1);

  /** brownian **/

  ldet_BM <- log_determinant(Gamma_bm); 
  ldet_BM <- ldet_BM + sum(sigma_bm);
  // jacobian of transform to correlation matrix

  llik_BM  <- -0.5 * (trace_gen_quad_form(inverse_spd(Gamma_bm * Gamma_bm'), C, eXbar) +
		      k * ldetC) - m * ldet_BM;

  increment_log_prob(llik_BM);
}

generated quantities {
  cov_matrix[k] Sigma_bm;
  
  Sigma_bm <- quad_form(Gamma_bm * Gamma_bm', diag_matrix(to_vector(sigma_bm)));
}
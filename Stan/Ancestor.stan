data {
  int k; // traits
  int m; // taxa
  cov_matrix[2*m - 2] C; // phylo
  vector[k] Xbar[m]; // médias
  //vector[k] priorX;
  //cov_matrix[k] priorS;
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

model {

  real ldet_BM;
  real llik_BM; 
  matrix[2 * m - 2,k] eXbar;
  matrix[k,k] invGamma_bm;  

  for(i in 1:m)
    eXbar[i] <- to_row_vector ((Xbar[i] - root) ./ to_vector(sigma_bm));
    
  for(i in 1:(m-2))
    eXbar[i+m] <- to_row_vector ((ancestor[i] - root) ./ to_vector(sigma_bm));
    
  /** priors **/

  /**

  for (i in 1:(m - 2))
    ancestor[i] ~ multi_normal(priorX, (10 ^ 3) * priorS);
  
  root ~ multi_normal(priorX, (10 ^ 3) * priorS);
  
  for (i in 1:k)
    sigma_bm[i] ~ chi_square(m - 1);

  **/

  // Gamma ~ ldk_corr_chol(1);
  // Gamma_bm ~ ldk_corr_chol(1);

  /** brownian **/

  ldet_BM <- 2 * (log_determinant(Gamma_bm) + 
		  sum(sigma_bm));
  // jacobian of transform to correlation matrix

  invGamma_bm <- inverse_spd(multiply_lower_tri_self_transpose(Gamma_bm));

  llik_BM  <- - 0.5 * (trace_gen_quad_form(invC, invGamma_bm, eXbar') +
		       k * ldetC + (2 * (m - 1) * ldet_BM)); // puta que pariu 2

  increment_log_prob(llik_BM);

}


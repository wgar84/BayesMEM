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
  matrix[2 * (m-1),2 * (m-1)] invC;
  
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

model {

  real ldet_BM;
  real llik_BM; 
  
  matrix[k,k] invGamma_bm;
  vector[k] eX[m,ni_max];
  matrix[2 * (m-1),k] eXbar;

  /** priors **/
  
  /**

  for (i in 1:(2 * (m-1)))
    Xbar[i] ~ multi_normal(priorX, (10^3) * priorS);
  
  alpha ~ multi_normal(priorX, (10^3) * priorS);
  
  lambda ~ uniform(0, 1);
  
  for (i in 1:k)
    {
      sigma[i] ~ chi_square(min(ni) - 1);
      sigma_bm[i] ~ chi_square(min(ni) - 1);
    }
  
  Gamma ~ ldk_corr_chol(1);
  Gamma_bm ~ ldk_corr_chol(1);
  
  **/

  /** brownian **/

  for(i in 1:(2 * (m-1)))
    eXbar[i] <- to_row_vector ((Xbar[i] - alpha) ./ to_vector(sigma_bm));
  
  ldet_BM <- 2 * (sum (log (diagonal (Gamma_bm))) + 
		  sum(sigma_bm));
  // jacobian of transform to correlation matrix

  invGamma_bm <- inverse_spd(multiply_lower_tri_self_transpose(Gamma_bm));

  llik_BM  <- - 0.5 * (trace_gen_quad_form(invC, invGamma_bm, eXbar') +
		       k * ldetC + (2 * (m - 1) * ldet_BM)); // puta que pariu 2

  increment_log_prob(llik_BM);
  
  /** pops **/

  for(i in 1:m)
    for(j in 1:ni[i])
      {
	eX[i,j] <- (X[i,j] - Xbar[i]) ./ to_vector(sigma);
	eX[i,j] ~ multi_normal_cholesky(zero_vector, Gamma);
      }

  increment_log_prob(- m * sum (ni) * sum(sigma)); 
  // jacobian of transform to correlation matrix
}

generated quantities {
  cov_matrix[k] Sigma;
  cov_matrix[k] Sigma_bm;
  
  Sigma <- quad_form_diag(multiply_lower_tri_self_transpose(Gamma), 
			  to_vector(sigma));
  Sigma_bm <- quad_form_diag(multiply_lower_tri_self_transpose(Gamma_bm), 
			     to_vector(sigma_bm));
}
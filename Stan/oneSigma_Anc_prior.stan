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
  cov_matrix[2 * (m-1)] invC;
  
  for (i in 1:k)
    zero_vector[i] <- 0;

  ldetC <- log_determinant(C);
  invC <- inverse_spd(C);
  
  
}

parameters {
  vector[k] Xbar[2 * (m-1)]; 
  vector[k] alpha;
  
  cholesky_factor_corr[k] Gamma; // Sigma = sGG's
  vector<lower=0>[k] sigma;

  cholesky_factor_corr[k] Gamma_bm; 
  vector<lower=0>[k] sigma_bm;

}

model {

  real ldet_BM;
  real llik_BM; 
  
  matrix[k,k] invGamma_bm;
  vector[k] eX[m,ni_max];
  matrix[2 * (m-1),k] eXbar;

  /** priors **/
  
  for (i in 1:(2 * (m-1)))
    Xbar[i] ~ multi_normal(priorX, priorS);
  
  alpha ~ multi_normal(priorX, priorS);
  
  //lambda ~ uniform(0, 1);
  
  for (i in 1:k)
    {
      sigma[i] ~ chi_square(max(ni));
      sigma_bm[i] ~ chi_square(max(ni));
    }
  
  Gamma ~ ldk_corr_cholesky(2);
  Gamma_bm ~ ldk_corr_cholesky(2);
  
  /** brownian **/

  for(i in 1:(2 * (m-1)))
    eXbar[i] <- to_row_vector ((Xbar[i] - alpha) ./ sigma_bm);
  
  ldet_BM <- 2 * (sum (log (diagonal (Gamma_bm))) + 
		  sum(sigma_bm));
  // jacobian of transform to correlation matrix

  invGamma_bm <- inverse_spd(multiply_lower_tri_self_transpose(Gamma_bm));

  llik_BM  <- - 0.5 * (trace_gen_quad_form(invC, invGamma_bm, eXbar') +
		       k * ldetC + (2 * (m - 1) * ldet_BM)); // puta que pariu 2

  increment_log_prob(llik_BM);
  
  /** pops **/

  for(i in 1:m)
    for(j in 1:(ni[i]))
      {
	eX[i,j] <- (X[i,j] - Xbar[i]) ./ sigma;
	eX[i,j] ~ multi_normal_cholesky(zero_vector, Gamma);
      }

  increment_log_prob(- m * sum (ni) * sum(sigma)); 
  // jacobian of transform to correlation matrix
}

generated quantities {
  cov_matrix[k] Sigma;
  cov_matrix[k] Sigma_bm;
  
  Sigma <- quad_form_diag(multiply_lower_tri_self_transpose(Gamma), 
			  sigma;
  Sigma_bm <- quad_form_diag(multiply_lower_tri_self_transpose(Gamma_bm), 
			     sigma_bm);
}
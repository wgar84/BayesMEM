data {
  int k; // traits
  int m; // taxa
  cov_matrix[m] C; // phylo
  int ni[m]; // amostras
  int ni_max; // máximo amostra
  vector[k] X[m,ni_max]; // dados
  //vector[k] priorX;
  //cov_matrix[k] priorS;
}

transformed data {
  vector[k] zero_vector;
    
  for (i in 1:k)
    zero_vector[i] <- 0;
 
}

parameters {
  vector[k] terminal[m]; 
  vector[k] root;
  
  cholesky_factor_corr[k] GammaW; // Sigma = sGG's
  vector<lower=0>[k] sigmaW;

  cholesky_factor_corr[k] GammaB; 
  vector<lower=0>[k] sigmaB;

}

model {

  real ldetB;
  real llikB; 
  
  matrix[k,k] invGammaB;
  vector[k] ex[m,ni_max];
  matrix[m,k] exbar;

  /** priors **/
    
  /** brownian **/

  for(i in 1:m)
    exbar[i] <- to_row_vector ((terminal[i] - root) ./ sigmaB);
  
  ldetB <- 2 * (sum (log (diagonal (GammaB))) + 
		sum(sigmaB));
  // jacobian of transform to correlation matrix
  
  invGammaB <- inverse_spd(multiply_lower_tri_self_transpose(GammaB));
  
  llikB <- - 0.5 * (trace_gen_quad_form(invGammaB, C, exbar) +
		    m * ldetB);
  
  increment_log_prob(llikB);
  
  /** pops **/

  for(i in 1:m)
    for(j in 1:(ni[i]))
      {
	ex[i,j] <- (X[i,j] - terminal[i]) ./ sigmaW;
	ex[i,j] ~ multi_normal_cholesky(zero_vector, GammaW);
      }

  increment_log_prob(- m * sum (ni) * sum(sigmaW)); 
  // jacobian of transform to correlation matrix
}

generated quantities {
  cov_matrix[k] SigmaW;
  cov_matrix[k] SigmaB;
  
  SigmaW <- quad_form_diag(multiply_lower_tri_self_transpose(GammaW), 
			   sigmaW);
  SigmaB <- quad_form_diag(multiply_lower_tri_self_transpose(GammaB), 
			   sigmaB);
}
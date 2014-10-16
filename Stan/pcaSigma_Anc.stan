data {
  int k; // traits
  int m; // taxa
  cov_matrix[2 * (m-1)] C; // phylo
  int ni[m]; // amostras
  int ni_max; // máximo amostra
  vector[k] X[m,ni_max]; // dados
  int n_pc;
  matrix[k,n_pc] eigen_Wml; //
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
  vector[k] ancestor[m-2];
  vector[k] root;
  
  cholesky_factor_corr[k] GammaW; // Sigma = sGG's
  vector<lower=0>[k] sigmaW;

  cholesky_factor_corr[n_pc] GammaB; 
  vector<lower=0>[n_pc] sigmaB;
}

model {
  real ldetB;
  real llikB; 
  
  matrix[n_pc,n_pc] invGammaB;
  vector[k] ex[m,ni_max];
  matrix[2 * (m-1),k] exbar;
  matrix[2 * (m-1),n_pc] exbar_W;

  /** priors **/
 
  // GammaB ~ lkj_corr_cholesky(2) // need it?
 
  /** pops **/
  
  for(i in 1:m)
    for(j in 1:(ni[i]))
      {
	ex[i,j] <- (X[i,j] - terminal[i]) ./ sigmaW;
	ex[i,j] ~ multi_normal_cholesky(zero_vector, GammaW);
      }
  
  increment_log_prob(- m * sum (ni) * sum(sigmaW)); 
  // jacobian of transform to correlation matrix
  
  /** means **/

  for(i in 1:m)
    exbar[i] <- to_row_vector(terminal[i] - root);
    
  for(i in 1:(m-2))
    exbar[i+m] <- to_row_vector(ancestor[i] - root);

  exbar_W <- exbar * eigen_Wml;
  
  for(i in 1:(2*(m-1)))
    exbar_W[i] <- exbar_W[i] ./ to_row_vector(sigmaB);
  
  ldetB <- 2 * (sum (log (diagonal (GammaB))) + sum(sigmaB));
  // jacobian of transform to correlation matrix
  
  invGammaB <- inverse_spd(multiply_lower_tri_self_transpose(GammaB));
  
  llikB <- - 0.5 * (trace_gen_quad_form(invGammaB, C, exbar_W));
  llikB <- llikB - ((m - 1) * ldetB); ### * 2 * 0.5
  
  increment_log_prob(llikB);
  
}

generated quantities {
  cov_matrix[n_pc] SigmaB_W;
  matrix[k,k] SigmaB;
  cov_matrix[k] SigmaW;
  
  SigmaB_W <- quad_form_diag(multiply_lower_tri_self_transpose(GammaB), sigmaB);
  SigmaB <- quad_form_sym(SigmaB_W, eigen_Wml');
  SigmaW <- quad_form_diag(multiply_lower_tri_self_transpose(GammaW), sigmaW);
}
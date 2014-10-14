data {
  int k; // traits
  int m; // taxa
  cov_matrix[2 * (m-1)] C; // phylo
  int ni[m]; // amostras
  int ni_max; // máximo amostra
  vector[k] X[m,ni_max]; // dados
  //vector[k] priorX;
  //cov_matrix[k] priorS;
}

transformed data {
  vector[k] zero_vector;
  int dim_flag;
  
  for (i in 1:k)
    zero_vector[i] <- 0;

  if (k >= m)
    {dim_flag <- m - 1;}
  else
    {dim_flag <- k;}
}


parameters {
  vector[k] terminal[m];
  vector[k] ancestor[m-2];
  vector[k] root;
  
  cholesky_factor_corr[k] GammaW; // Sigma = sGG's
  vector<lower=0>[k] sigmaW;

  cholesky_factor_corr[dim_flag] GammaB; 
  vector<lower=0>[dim_flag] sigmaB;
}


transformed parameters {

  cov_matrix[k] SigmaW;
  
  matrix[k,k] eigen_vec;
  matrix[k,dim_flag] eigen_vec_trimmed;
  
  SigmaW <- quad_form_diag(multiply_lower_tri_self_transpose(GammaW), sigmaW);

  eigen_vec <- eigenvectors_sym(SigmaW);
  eigen_vec_trimmed <- block(eigen_vec, 1, 1, k, dim_flag);
}


model {
  real ldetB;
  real llikB; 
  
  matrix[dim_flag,dim_flag] invGammaB;
  vector[k] ex[m,ni_max];
  matrix[2 * (m-1),k] exbar;

  matrix[2 * (m-1),dim_flag] exbar_W;
    
  /** priors **/
  
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

  exbar_W <- exbar * eigen_vec_trimmed;
  
  for(i in 1:(2*(m-1)))
    exbar_W[i] <- exbar_W[i] ./ to_row_vector(sigmaB);
  
  ldetB <- 2 * (sum (log (diagonal (GammaB))) + sum(sigmaB));
  // jacobian of transform to correlation matrix
  
  invGammaB <- inverse_spd(multiply_lower_tri_self_transpose(GammaB));
  
  llikB <- - 0.5 * (trace_gen_quad_form(invGammaB, C, exbar_W));
  llikB <- llikB - ((m - 1) * ldetB);
  
  increment_log_prob(llikB);
  
}

generated quantities {
  cov_matrix[dim_flag] SigmaB_W;
  matrix[k,k] SigmaB;
  
  SigmaB_W <- quad_form_diag(multiply_lower_tri_self_transpose(GammaB), sigmaB);
  SigmaB <- quad_form_sym(SigmaB_W, eigen_vec_trimmed');
}
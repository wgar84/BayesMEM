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

parameters {
  vector[k] terminal[m];
  vector[k] ancestor[m-2];
  vector[k] root;
  
  cholesky_factor_corr[k] GammaW;
  vector<lower=0>[k] sigmaW;

  real<lower=0> drift;

  cholesky_factor_corr[k] Gamma_beta; 
  vector<lower=0>[k] sigma_beta;

}

transformed parameters {
  cov_matrix[k] SigmaW;
  cov_matrix[k] Sigma_beta;
  cov_matrix[k] SigmaB;

  SigmaW <- quad_form_diag(multiply_lower_tri_self_transpose(GammaW), sigmaW);
  
  Sigma_beta <- quad_form_diag(multiply_lower_tri_self_transpose(Gamma_beta), sigma_beta);
  SigmaB <- drift * SigmaW + quad_form_sym(Sigma_beta, SigmaW); 
}

model {
  
  real ldetB;
  real llikB; 
  
  matrix[k,k] invSigmaB;
  matrix[2 * (m-1),k] eXbar;
  
  /** priors **/
  
  Gamma_beta ~ lkj_corr_cholesky(2); // só pra dar um pouco mais de densidade pra I_k
  
  /** brownian **/

  for(i in 1:m)
    eXbar[i] <- to_row_vector (terminal[i] - root);
  
  for(i in 1:(m-2))
    eXbar[i+m] <- to_row_vector (ancestor[i] - root);
  
  
  ldetB <- log_determinant(SigmaB);
  invSigmaB <- inverse_spd(SigmaB);
  
  llikB  <- - 0.5 * (trace_gen_quad_form(invSigmaB, C, eXbar) +
		     2 * (m - 1) * ldetB);

  increment_log_prob(llikB);
  
  /** pops **/

  for(i in 1:m)
    for(j in 1:(ni[i]))
      X[i,j] ~ multi_normal(terminal[i], SigmaW);
      

  // increment_log_prob(- m * sum (ni) * sum(sigma)); 
  // jacobian of transform to correlation matrix
}


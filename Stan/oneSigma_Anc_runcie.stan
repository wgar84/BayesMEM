data {
  int k; // traits
  int m; // taxa
  cov_matrix[2 * (m-1)] C; // phylo
  int ni[m]; // amostras
  int ni_max; // máximo amostra
  vector[k] X[m,ni_max]; // dados
  // runcie decomposition
  int n_fac;
  real a1W;
  real a2W;
  real b1W;
  real b2W;
  real asW;
  real bsW;
  real niW;
  real a1B;
  real a2B;
  real b1B;
  real b2B;
  real asB;
  real bsB;
  real niB;
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
  // factors
  matrix[k,n_fac] LambdaW;
  matrix[k,n_fac] LambdaB;
  // trait specific
  vector<lower=0>[k] PsiW;
  vector<lower=0>[k] PsiB;
  // precision for lambda
  matrix<lower=0>[k,n_fac] phiW;
  matrix<lower=0>[k,n_fac] phiB;
  // shrinkage
  vector<lower=0>[n_fac] deltaW;
  vector<lower=0>[n_fac] deltaB;

}

model {
  
  real ldetB;
  real llikB; 

  // precision scaling
  vector<lower=0>[n_fac] tauW;
  vector<lower=0>[n_fac] tauB;
  
  matrix[k,k] SigmaW;
  matrix[k,k] SigmaB;
  
  matrix[2 * (m-1),k] exbar;
  
  /** hyper **/

  deltaW[1] ~ gamma(a1W, b1W);
  deltaB[1] ~ gamma(a1B, b1B);
  tauW[1] <- deltaW[1];
  tauB[1] <- deltaB[1];
  
  for(i in 2:n_fac)
    {
      deltaW[i] ~ gamma(a2W, b2W);
      tauW[i] <- prod(head(deltaW, i));
      deltaB[i] ~ gamma(a2B, b2B);
      tauB[i] <- prod(head(deltaB, i));
    }

  for(i in 1:k)
    {
      PsiW[i] ~ inv_gamma(asW, bsW);
      PsiB[i] ~ inv_gamma(asB, bsB);
      for(j in 1:n_fac)
	{
	  phiW[i,j] ~ gamma(0.5 * niW, 0.5 * niW);
	  LambdaW[i,j] ~ normal(0, (phiW[i,j] * tauW[j])^(-1));
	  phiB[i,j] ~ gamma(0.5 * niB, 0.5 * niB);
	  LambdaB[i,j] ~ normal(0, (phiB[i,j] * tauB[j])^(-1));
	}
    }
  
  SigmaW <- tcrossprod(LambdaW) + diag_matrix(PsiW);
  SigmaB <- tcrossprod(LambdaB) + diag_matrix(PsiB);

  /** brownian **/
  
  for(i in 1:m)
    exbar[i] <- to_row_vector (terminal[i] - root);

  for(i in 1:(m-2))
    exbar[i+m] <- to_row_vector (ancestor[i] - root);
  
  llikB <- - 0.5 * (trace_gen_quad_form(inverse_spd(SigmaW), C, exbar) +
		    2 * (m - 1) * log_determinant(SigmaW));
  
  increment_log_prob(llikB);
  
  /** pops **/

  for(i in 1:m)
    for(j in 1:(ni[i]))
      X[i,j] ~ multi_normal_cholesky(terminal[i], SigmaW);
}

generated quantities {
  cov_matrix[k] SigmaW;
  cov_matrix[k] SigmaB;
  
  SigmaW <- tcrossprod(LambdaW) + diag_matrix(PsiW);
  SigmaB <- tcrossprod(LambdaB) + diag_matrix(PsiB);

}
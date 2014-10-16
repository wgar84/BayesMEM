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
  matrix[n_fac,k] LambdaW;
  matrix[n_fac,k] LambdaB;
  // trait specific
  vector<lower=0>[k] PsiW;
  vector<lower=0>[k] PsiB;
  // precision for lambda
  matrix<lower=0>[n_fac,k] phiW;
  matrix<lower=0>[n_fac,k] phiB;
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

  PsiW ~ inv_gamma(asW, bsW);
  PsiB ~ inv_gamma(asB, bsB);
  for(i in 1:n_fac)
    {
      if (i == 1)
	{
	  deltaW[i] ~ gamma(a1W, b1W);
	  deltaB[i] ~ gamma(a1B, b1B);
	}
      else
	{
	  deltaW[i] ~ gamma(a2W, b2W);
	  deltaB[i] ~ gamma(a2B, b2B);
	}
      tauW[i] <- prod(head(deltaW, i));
      tauB[i] <- prod(head(deltaB, i));
      phiW[i] ~ gamma(0.5 * niW, 0.5 * niW);
      phiB[i] ~ gamma(0.5 * niB, 0.5 * niB);
      LambdaW[i] ~ normal(0, exp(-log(phiW[i] * tauW[i]))); // invertendo, esquisito
      LambdaB[i] ~ normal(0, exp(-log(phiB[i] * tauB[i])));
    }
  SigmaW <- crossprod(LambdaW) + diag_matrix(PsiW);
  SigmaB <- crossprod(LambdaB) + diag_matrix(PsiB);

  /** brownian **/
  
  for(i in 1:m)
    {
      exbar[i] <- to_row_vector (terminal[i] - root);
      if(i < m-1)
	exbar[i+m] <- to_row_vector (ancestor[i] - root);
    }
 
  llikB <- - 0.5 * (trace_gen_quad_form(inverse_spd(SigmaW), C, exbar) +
		    2 * (m - 1) * log_determinant(SigmaW));
  
  increment_log_prob(llikB);
  
  /** pops **/

  for(i in 1:m)
    head(X[i], ni[i]) ~ multi_normal(terminal[i], SigmaW); //subset não-nulo de X[i]
}

generated quantities {
  cov_matrix[k] SigmaW;
  cov_matrix[k] SigmaB;
  
  SigmaW <- crossprod(LambdaW) + diag_matrix(PsiW);
  SigmaB <- crossprod(LambdaB) + diag_matrix(PsiB);

}
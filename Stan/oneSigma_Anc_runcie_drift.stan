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
  // trait specific
  vector<lower=0>[k] PsiW;
  // precision for lambda
  matrix<lower=0>[n_fac,k] phiW;
  // shrinkage
  vector<lower=0>[n_fac] deltaW;
  real<lower=0> variance_brownian;
  
}

model {
  
  real llikB; 

  // precision scaling
  vector<lower=0>[n_fac] tauW;
  
  matrix[k,k] SigmaW;
  matrix[k,k] SigmaB;
  
  matrix[2 * (m-1),k] exbar;
  
  /** hyper **/

  PsiW ~ inv_gamma(asW, bsW);
  for(i in 1:n_fac)
    {
      if (i == 1)
	deltaW[i] ~ gamma(a1W, b1W);
      else
	deltaW[i] ~ gamma(a2W, b2W);
      tauW[i] <- prod(head(deltaW, i));
      phiW[i] ~ gamma(0.5 * niW, 0.5 * niW);
      LambdaW[i] ~ normal(0, exp(-log(phiW[i] * tauW[i]))); // invertendo, esquisito
    }

  SigmaW <- crossprod(LambdaW) + diag_matrix(PsiW);
  SigmaB <- variance_brownian * SigmaW;

  /** brownian **/
  
  for(i in 1:m)
    {
      exbar[i] <- to_row_vector (terminal[i] - root);
      if(i < m-1)
	exbar[i+m] <- to_row_vector (ancestor[i] - root);
    }
 
  llikB <- - 0.5 * (trace_gen_quad_form(inverse_spd(SigmaB), C, exbar) +
		    2 * (m - 1) * log_determinant(SigmaB));
  
  increment_log_prob(llikB);
  
  /** pops **/

  for(i in 1:m)
    head(X[i], ni[i]) ~ multi_normal(terminal[i], SigmaW); //subset não-nulo de X[i]
}

generated quantities {
  cov_matrix[k] SigmaW;
  
  SigmaW <- crossprod(LambdaW) + diag_matrix(PsiW);
  
}
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
  real a1S;
  real a2S;
  real b1S;
  real b2S;
  real asS;
  real bsS;
  real niS;
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
  matrix[n_fac,k] LambdaS;
  // trait specific
  vector<lower=0>[k] PsiW;
  vector<lower=0>[k] PsiS;
  // precision for lambda
  matrix<lower=0>[n_fac,k] phiW;
  matrix<lower=0>[n_fac,k] phiS;
  // shrinkage
  vector<lower=0>[n_fac] deltaW;
  vector<lower=0>[n_fac] deltaS;

  real variance_brownian;

}

model {
  
  real llikB; 

  // precision scaling
  vector<lower=0>[n_fac] tauW;
  vector<lower=0>[n_fac] tauS;
  
  matrix[k,k] SigmaW;
  matrix[k,k] SigmaS;
  matrix[k,k] SigmaB;
  
  matrix[2 * (m-1),k] exbar;
  
  /** hyper **/

  PsiW ~ inv_gamma(asW, bsW);
  PsiS ~ inv_gamma(asS, bsS);
  for(i in 1:n_fac)
    {
      if (i == 1)
	{
	  deltaW[i] ~ gamma(a1W, b1W);
	  deltaS[i] ~ gamma(a1S, b1S);
	}
      else
	{
	  deltaW[i] ~ gamma(a2W, b2W);
	  deltaS[i] ~ gamma(a2S, b2S);
	}
      tauW[i] <- prod(head(deltaW, i));
      tauS[i] <- prod(head(deltaS, i));
      phiW[i] ~ gamma(0.5 * niW, 0.5 * niW);
      phiS[i] ~ gamma(0.5 * niS, 0.5 * niS);
      LambdaW[i] ~ normal(0, exp(-log(phiW[i] * tauW[i]))); // invertendo, esquisito
      LambdaS[i] ~ normal(0, exp(-log(phiS[i] * tauS[i])));
    }
  SigmaW <- crossprod(LambdaW) + diag_matrix(PsiW);
  SigmaS <- crossprod(LambdaS) + diag_matrix(PsiS);
  
  SigmaB <- variance_brownian * SigmaW + quad_form(SigmaS, SigmaW);

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
  cov_matrix[k] SigmaB;
  cov_matrix[k] SigmaS;
  
  SigmaW <- crossprod(LambdaW) + diag_matrix(PsiW);
  SigmaS <- crossprod(LambdaS) + diag_matrix(PsiS);
  SigmaB <- variance_brownian * SigmaW + quad_form(SigmaS, SigmaW);
}
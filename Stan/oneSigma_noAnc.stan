data {
  int k; // traits
  int m; // taxa
  cov_matrix[m] C; // phylo
  int ni[m]; // amostras
  int ni_max; // máximo amostra
  vector[k] X[m,ni_max]; // dados
  cov_matrix[k] id; 
  vector[k] zeroes;
  vector[k] priorX;
  cov_matrix[k] priorS;
}

transformed data {
  cov_matrix[m] inv_C;
  real ldet_C;
  cov_matrix[k] precS;
  
  ldet_C <- log_determinant(C);
  inv_C <- inverse_spd(C);

  precS <- inverse_spd(priorS);

}

parameters {
  vector[k] Xbar[m]; 
  vector[k] alpha;
  cov_matrix[k] invSigma_bm; // Brownian Motion Matrix
  cov_matrix[k] invSigma;
  // int<lower=0,upper=100> lambdacent; // later
}

model {
  matrix[m,k] Xbar_center;
  
  real ldet_BM;

  real ldet_current;
  real llik_cum; 

  // priors
  
  invSigma_bm ~ wishart(k, precS);
  invSigma ~ wishart(k, precS);
  
  for (i in 1:m)
    Xbar[i] ~ multi_normal_prec(priorX, precS);
  
  alpha ~ multi_normal_prec(priorX, precS);
  
  for (i in 1:m)
    Xbar_center[i] <- to_row_vector(Xbar[i] - alpha);
  

  ldet_BM <- - log_determinant(invSigma_bm);
  
  llik_cum  <- -0.5 * (trace_gen_quad_form(invSigma_bm, inv_C, Xbar_center) +
		       k * ldet_C + m * ldet_BM);
  
  ldet_current <- - log_determinant(invSigma);
  

  for (i in 1:m)
    {
      for (j in 1:(ni[i]))
	llik_cum <- llik_cum - 0.5 * (ldet_current + 
				      quad_form_sym(invSigma, X[i][j] - Xbar[i]));
    }
     
  increment_log_prob(llik_cum);
}

generated quantities{
  cov_matrix[k] Sigma;
  cov_matrix[k] Sigma_bm;
  
  Sigma_bm <- inverse_spd(invSigma_bm);
  Sigma <- inverse_spd(invSigma);
  
}


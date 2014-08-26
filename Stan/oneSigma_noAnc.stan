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
  matrix[m,m] inv_C;
  real ldet_C;

  ldet_C <- log_determinant(C);
  inv_C <- inverse_spd(C);
}

parameters {
  vector[k] Xbar[m]; 
  vector[k] alpha;
  cov_matrix[k] Sigma_bm; // Brownian Motion Matrix
  cov_matrix[k] Sigma;
  // int<lower=0,upper=100> lambdacent; // later
}

model {
  matrix[m,k] Xbar_center;
  
  real ldet_BM;

  matrix[k,k] inv_current; // reciclável
  
  real ldet_current;
  real llik_cum; 

  // priors
  
  Sigma_bm ~ inv_wishart(k, priorS);
  Sigma ~ inv_wishart(k, priorS);
  
  for (i in 1:m)
    Xbar[i] ~ multi_normal(priorX, priorS);
  
  alpha ~ multi_normal(priorX, priorS);
  
  for (i in 1:m)
    Xbar_center[i] <- to_row_vector(Xbar[i] - alpha);
  

  ldet_BM <- log_determinant(Sigma_bm);
  
  llik_cum  <- -0.5 * (trace_gen_quad_form(inverse_spd(Sigma_bm), 
					   inv_C, Xbar_center) +
		       k * ldet_C + m * ldet_BM);
  
  ldet_current <- log_determinant(Sigma);
  inv_current <- inverse_spd(Sigma);

  for (i in 1:m)
    {
      for (j in 1:(ni[i]))
	llik_cum <- llik_cum - 0.5 * (ldet_current + 
				      quad_form_sym(inv_current, 
						    X[i][j] - Xbar[i]));
    }
     
  increment_log_prob(llik_cum);
  
}

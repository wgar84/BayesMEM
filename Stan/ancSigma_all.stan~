data {
  int k; // traits
  int m; // taxa
  cov_matrix[2 * (m - 1)] C; // phylo
  int ni[m]; // amostras
  int ni_max; // máximo amostra
  vector[k] X[m,ni_max]; // dados
  cov_matrix[k] id; 
  vector[k] zeroes;
  vector[k] priorX;
  cov_matrix[k] priorS;
}

parameters {
  vector[k] Xbar[2 * (m - 1)]; 
  vector[k] alpha;
  cov_matrix[k] Sigma_bm; // Brownian Motion Matrix
  cov_matrix[k] Sigma;
  // int<lower=0,upper=100> lambdacent; // later
}

model {
  matrix[2*(m-1),k] Xbar_center;
  real ldet_C;
  matrix[2*(m-1),2*(m-1)] inv_C;
  real ldet_BM;

  matrix[k,k] inv_current; // reciclável
  real ldet_current;
  real llik_cum; 

  // priors
  
  Sigma_bm ~ inv_wishart(k, priorS);
  Sigma ~ inv_wishart(k, priorS);
  
  for (i in 1:(2*(m-1)))
    Xbar[i] ~ multi_normal(priorX, priorS);
  
  alpha ~ multi_normal(priorX, priorS);
  
  /**
  for (i in 1:m)
    {
      Sigma[i] ~ inv_wishart(k+1, 0.002*id);
      Xbar[i] ~ multi_normal_prec(zeroes, 10000 * id);
    }

  Sigma_bm ~ inv_wishart(k+1, 0.002*id);

  alpha ~ multi_normal_prec(zeroes, 1000 * id);
  **/

  // brownian
  
  for (i in 1:(2*(m-1)))
    Xbar_center[i] <- to_row_vector(Xbar[i] - alpha);
  
  ldet_C <- log_determinant(C);
  ldet_BM <- log_determinant(Sigma_bm);
  
  inv_C <- inverse_spd(C);

  llik_cum  <- -0.5 * (trace_gen_quad_form(inverse_spd(Sigma_bm), 
					   inv_C, Xbar_center) +
		       k * ldet_C + m * ldet_BM);
  
  // multinormal
  
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

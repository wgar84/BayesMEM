data {
  int k; // no. de taxa
  int p; // no. de caracteres
  // sem querer, troquei o k pelo p =P
  cov_matrix[k] C; // matriz de comprimentos de ramo
  vector[p] geoffroyi[45];
  vector[p] kuhlii[129];
  vector[p] jacchus[32];
  vector[p] penicillata[68];
  vector[p] argentata[43];
  vector[p] humeralifera[21];
  vector[p] pygmaea[74];
  cov_matrix[p] priorV;
  vector[p] priorB;
}

parameters {
  vector[p] B[k]; // medias
  vector[p] alpha; // raiz
  cov_matrix[p] V; // covariancia do processo browniano
  cov_matrix[p] sigma_geoffroyi;
  cov_matrix[p] sigma_kuhlii;
  cov_matrix[p] sigma_jacchus;
  cov_matrix[p] sigma_penicillata;
  cov_matrix[p] sigma_argentata;
  cov_matrix[p] sigma_humeralifera;
  cov_matrix[p] sigma_pygmaea;
}

model {
  real loglik_geo;
  real loglik_kuh;
  real loglik_jac;
  real loglik_pen;
  real loglik_arg;
  real loglik_hum;
  real loglik_pyg;
  
  real logdet_geo;
  real logdet_kuh;
  real logdet_jac;
  real logdet_pen;
  real logdet_arg;
  real logdet_hum;
  real logdet_pyg;
  
  matrix[p,p] inv_geoffroyi;
  matrix[p,p] inv_kuhlii;
  matrix[p,p] inv_jacchus;
  matrix[p,p] inv_penicillata;
  matrix[p,p] inv_argentata;
  matrix[p,p] inv_humeralifera;
  matrix[p,p] inv_pygmaea;
  
  matrix[k,p] B_center;
  real loglik_bm;
  real logdet_C;
  real logdet_V;
  
  real llsum;

  // priors
  
  V ~ inv_wishart(p, priorV);

  for (i in 1:k)
    B[i] ~ multi_normal(priorB, priorV);

  alpha ~ multi_normal(priorB, priorV);

  sigma_geoffroyi ~ inv_wishart(p, priorV);
  sigma_kuhlii ~ inv_wishart(p, priorV);
  sigma_jacchus ~ inv_wishart(p, priorV);
  sigma_penicillata ~ inv_wishart(p, priorV);
  sigma_argentata ~ inv_wishart(p, priorV);
  sigma_humeralifera ~ inv_wishart(p, priorV);
  sigma_pygmaea ~ inv_wishart(p, priorV);


  // loglik BM

  for (i in 1:k)
    B_center[i] <- to_row_vector(B[i] - alpha);
  
  logdet_C <- log(determinant(C));
  logdet_V <- log(determinant(V));
  
  loglik_bm  <- -0.5 * (trace_gen_quad_form(inverse_spd(V), inverse_spd(C), B_center) +
			p * logdet_C + k * logdet_V);
  // trace_gen_quad_form é o traço do produto V * Bc' * C * Bc
  // dá na mesma a ordem dessa bosta pq tr(ABCD) = tr(DABC) = tr(CDAB) = tr(BCDA)


  // loglik sps
  
  logdet_geo <- log_determinant(sigma_geoffroyi);
  inv_geoffroyi <- inverse_spd(sigma_geoffroyi);
  
  loglik_geo <- 0;
  for (i in 1:45)
    loglik_geo <- loglik_geo - 0.5 * (logdet_geo + 
				      quad_form_sym(inv_geoffroyi, geoffroyi[i] - B[1]));

  logdet_kuh <- log_determinant(sigma_kuhlii);
  inv_kuhlii <- inverse_spd(sigma_kuhlii);
  
  loglik_kuh <- 0;
  for (i in 1:129)
    loglik_kuh <- loglik_kuh - 0.5 * (logdet_kuh + 
				      quad_form_sym(inv_kuhlii, kuhlii[i] - B[2]));
  
  logdet_jac <- log_determinant(sigma_jacchus);
  inv_jacchus <- inverse_spd(sigma_jacchus);
  
  loglik_jac <- 0;
  for (i in 1:32)
    loglik_jac <- loglik_jac - 0.5 * (logdet_jac + 
				      quad_form_sym(inv_jacchus, jacchus[i] - B[3]));

  logdet_pen <- log_determinant(sigma_penicillata);
  inv_penicillata <- inverse_spd(sigma_penicillata);
  
  loglik_pen <- 0;
  for (i in 1:68)
    loglik_pen <- loglik_pen - 0.5 * (logdet_pen + 
				      quad_form_sym(inv_penicillata, 
						    penicillata[i] - B[4]));

  logdet_arg <- log_determinant(sigma_argentata);
  inv_argentata <- inverse_spd(sigma_argentata);
  
  loglik_arg <- 0;
  for (i in 1:43)
    loglik_arg <- loglik_arg - 0.5 * (logdet_arg + 
				      quad_form_sym(inv_argentata, argentata[i] - B[5]));

  logdet_hum <- log_determinant(sigma_humeralifera);
  inv_humeralifera <- inverse_spd(sigma_humeralifera);
  
  loglik_hum <- 0;
  for (i in 1:21)
    loglik_hum <- loglik_hum - 0.5 * (logdet_hum + 
				      quad_form_sym(inv_humeralifera, 
						    humeralifera[i] - B[6]));
  
  logdet_pyg <- log_determinant(sigma_pygmaea);
  inv_pygmaea <- inverse_spd(sigma_pygmaea);
  
  loglik_pyg <- 0;
  for (i in 1:74)
    loglik_pyg <- loglik_pyg - 0.5 * (logdet_pyg + 
				      quad_form_sym(inv_pygmaea, pygmaea[i] - B[7]));

  llsum <- loglik_bm;
  llsum <- llsum + loglik_geo;
  llsum <- llsum + loglik_kuh;
  llsum <- llsum + loglik_jac;
  llsum <- llsum + loglik_pen;
  llsum <- llsum + loglik_arg;
  llsum <- llsum + loglik_hum;
  llsum <- llsum + loglik_pyg;
  
  increment_log_prob(llsum);
  
}

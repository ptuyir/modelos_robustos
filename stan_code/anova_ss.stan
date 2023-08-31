data {
  int<lower=1> n; // number of observations 
  int<lower=1> n_prev; // number of observations forecast
  int<lower=1> m; // number of individuals 
  int<lower=1> p; // number of time periods
  int<lower=1, upper=m> t[n]; // row indices
  int<lower=1, upper=p> h[n]; // column indices
  int<lower=1, upper=m> t_prev[n_prev];
  int<lower=1, upper=p> h_prev[n_prev]; 
  //matrix[n,2] id_time; // matrix of individual and time identifiers
  vector[n] y; // response variable
}

parameters {
  real u; 
  vector[m-1] alpha; 
  vector[m-1] beta; 
  row_vector[m-1] beta_line1;
  real<lower=0> sigma_u; 
  real<lower=0> sigma_alpha; 
  real<lower=0> sigma_beta;

}

transformed parameters {
  matrix[m,p] u_ij;
  //vector [m] alpha_final = append_row(0,cumulative_sum(alpha));
  vector [m] alpha_final;
  matrix [m,p]  beta_final;
  beta_final[:,1] = rep_vector(0,m);
  beta_final[1,2:] = beta_line1;
  //beta_final[1, :] =
  
  for( l in 1:(p-1)){
  	if( l==1){
  	 alpha_final[l] = 0;
  	}
  	alpha_final[l+1] = alpha_final[l] + alpha[l];
  	
  }
  for( l in 2:p){
  for( r in 2:p){
  	beta_final[l,r] = beta_final[l-1,r] + beta[l-1];
  	}
  }
  
  for (i in 1:m) {
    for (j in 1:p) {
      u_ij[i,j] = u + alpha_final[i] + beta_final[i,j];
    }
  }
}

model {
  // priors
  u ~ normal(0, 100);
  alpha ~ normal(0, sqrt(sigma_alpha));
  beta_line1 ~ normal(0, 100);
  for (i in 1:(m-1)) {
      beta[i] ~ normal(0, sqrt(sigma_beta));
    }
  sigma_alpha ~ inv_gamma(0.001,0.001);
  sigma_beta ~ inv_gamma(0.001,0.001);
  sigma_u ~ inv_gamma(0.001,0.001);

  // likelihood
  for (k in 1:n) {
    //int i = t[k];
    //int j = h[k];
    y[k] ~ normal((u_ij[t[k],h[k]]), sqrt(sigma_u));
  }
}

generated quantities {
 vector[n] log_lik;
 vector[n] y_pred;
 vector[n_prev] y_new;
 real reserva;
 for (k in 1:n_prev) {
    //int i = t_prev[k];
    //int j = h_prev[k];
    y_new[k] = exp(u + alpha_final[t_prev[k]] + beta_final[t_prev[k],h_prev[k]]);
  }
  for (k in 1:n) {
    //int i = t[k];
    //int j = h[k];
    y_pred[k] = normal_rng(exp(u_ij[t[k],h[k]]), sqrt(sigma_u));
    log_lik[k] = normal_lpdf(y[k]|u_ij[t[k],h[k]], sqrt(sigma_u));
  }
  reserva = sum(y_new);
}


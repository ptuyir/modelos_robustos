functions {
}
data {
  int<lower=1> N;  // numero de valores observados do triangulo
  vector[N] Y;  // variavel respostas log claims
  int<lower=1> K_alpha;  // numero de efeitos para alpha
  int<lower=1> K_beta;  // numero de efeitos para beta
  matrix[N, K_alpha] X_alpha;  // matriz desenho para alpha
  matrix[N, K_beta] X_beta;  // matrix desenho para beta
  int<lower=1> N_prev;  // numero de observacoes faltantes do triangulo para fazer previsao
  int<lower=1> K_alpha_prev;  // numero de efeitos para alpha para previsao
  int<lower=1> K_beta_prev;  // numero de efeitos para beta para previsao
  matrix[N_prev, K_alpha_prev] X_alpha_prev;  // matriz desenho de alpha para as previsoes
  matrix[N_prev, K_beta_prev] X_beta_prev;  // matrix desenho de beta  para as previsoes
}
transformed data {
}

parameters {
  vector[K_alpha-1] a;  // efeitos de alpha
  vector[K_beta-1] b; // efeitos de beta
  real mu;  // media mu
  real<lower=0> sigma2;  // variancia
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  vector[K_alpha] alpha = append_row(a, -sum(a)); // fazendo a soma dos alphas ser 0
  vector[K_beta] beta = append_row(b,-sum(b)); // fazendo a soma dos betas ser 0
  //real sigma = sqrt(1/tau);
  //real sigma2 = sigma^2;
  vector[N] muij; // muij = mu+X*beta+X2*alpha
  muij=mu+X_alpha*alpha+X_beta*beta;
  lprior += normal_lpdf(a | 0, 10);
  lprior += normal_lpdf(b |0,10);
  lprior += normal_lpdf(mu | 0, 100);
  lprior += inv_gamma_lpdf(sigma2 | 0.001, 0.001);
}
model {
  // likelihood including constants
  target += normal_lpdf(Y | muij, sqrt(sigma2));
 
  // priors including constants
  target += lprior;
}
  

generated quantities {
  vector[N] log_lik;
  vector[N] y_pred;
  vector[N_prev] muij_new;
  vector[N_prev] y_new;
  
  real reserva;
  muij_new=exp(mu+X_alpha_prev*alpha+X_beta_prev*beta);
  for (n in 1:N) {
    y_pred[n] = normal_rng(exp(muij[n]), sqrt(sigma2));
    log_lik[n]= normal_lpdf(Y[n] | muij[n], sqrt(sigma2));
  }
  for (k in 1:N_prev) {
    y_new[k]= normal_rng(muij_new[k], sqrt(sigma2));
  }
  reserva = sum(y_new);

}

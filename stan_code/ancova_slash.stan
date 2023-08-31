functions {
}
data {
  int<lower=1> N;  // numero de valores observados do triangulo
  vector[N] Y;  // variavel respostas log claims
  int<lower=1> K_beta;  // numero de efeitos para beta
  vector[N] X_alpha;  // vetor desenho para alpha que tem efeito linear
  matrix[N, K_beta] X_beta;  // matrix desenho para beta
  int<lower=1> N_prev;  // numero de observacoes faltantes do triangulo para fazer previsao
  int<lower=1> K_beta_prev;  // numero de efeitos para beta para previsao
  vector[N_prev] X_alpha_prev;   // vetor desenho para alpha que tem efeito linear das previsoes
  matrix[N_prev, K_beta_prev] X_beta_prev;  // matrix desenho de beta  para as previsoes
}
transformed data {
}

parameters {
  real alpha;  // efeitos de alpha
  vector[K_beta-1] b; // efeitos de beta
  real mu; 
  real<lower=0> sigma2;
  vector<lower=0,upper=1>[N] lambda;
  real <lower=0> ni;
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  vector[K_beta] beta = append_row(b,-sum(b));
  vector[N] muij;
  vector<lower=0>[N] lambda2;
  muij=mu+X_alpha*alpha+X_beta*beta;
  for (j in 1:N) {
    lambda2[j] = sqrt(sigma2)/(lambda[j]);
  }
  
  lprior += normal_lpdf(alpha | 0, 10);
  lprior += normal_lpdf(b |0,10);
  lprior += normal_lpdf(mu | 0, 100);
  lprior += inv_gamma_lpdf(sigma2 | 0.001, 0.001);
  lprior += beta_lpdf(lambda| ni,1);

}
model {
  // likelihood including constants
  target += gamma_lpdf(ni|1,0.1);
  target += normal_lpdf(Y | muij, lambda2);

  target += lprior;
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_pred;
  vector[N_prev] lambda_prev;
  vector[N_prev] muij_new;
  vector[N_prev] y_new;
  real reserva;
  muij_new=exp(mu+X_alpha_prev*alpha+X_beta_prev*beta);
  for (n in 1:N) {
    y_pred[n] = normal_rng(exp(muij[n]), lambda2[n]);
    log_lik[n]= normal_lpdf(Y[n] | muij[n], lambda2[n]);
  }
  for (m in 1:N_prev){
    lambda_prev[m] = beta_rng(ni,1);
    y_new[m] = normal_rng(muij_new[m],sqrt(sigma2)/(lambda_prev[m]));
    }
  reserva = sum(y_new);
}

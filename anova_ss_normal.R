##### DADOS SIMULADOS SS MODELO ANOVA ERRO NORMAL #####
rm(list=ls())
library(dplyr)
library(lattice)
library(brms)
library(ChainLadder)
library(rstanarm)
library(rstan)
#semente
set.seed(10)

# media, variancia e graus de liberdade
mu<- 200;sigma2<-1;ni<-3
# numero de linhas e colunas do triangulo de run-off
m<-20
# numero de infos que vou gerar para preencer o triangulo
n<-400

## gerando os coeficientes dos efeitos das linhas N(0,4)
alpha <-NA
alpha[1]<- 0
erro_alpha<-rnorm(m-1,0,2)
for (i in 2:m) {
  alpha[i]<- alpha[i-1]+erro_alpha[i-1]
  
}

## gerando os coeficientes dos efeito das colunas N(0,16)
beta<- matrix(NA,nrow=20,ncol = 20)
beta[,1]<-0
erro<-rnorm(m-1,0,4)
beta[1,2:20]<-rnorm(m-1,0,4)
for (j in 2:20) {
  for(l in 2:20) {
    beta[j,l] <- beta[j-1,l]+erro[j-1]
  }
}

## erro N(0,1)
epsilon<- rnorm(400,0,1)

## gerando os muij segundo o modelo log anova
muij<-mu+matrix(rep(alpha,20),ncol=m)+beta

## gerando os valores de zij zij= muij+ erro
zij<-muij+matrix(epsilon,ncol=i)

# salvando o triangulo total
triangulo_total <- zij

## torando NA os valores abaixo da diagonal do triangulo que e desconhecido
zij[row(zij) > i - col(zij) + 1]<-NA

## colocando os dados do formato que estava no livro para fazer as estimacoes
n<-20
Claims <- data.frame(originf = factor(rep(1:n, n:1)),
                     dev=sequence(n:1),
                     inc.paid= na.omit(as.vector(t(zij))))
## triangulo de run-off
inc.triangle <- with(Claims, {
  M <- matrix(nrow=n, ncol=n,
              dimnames=list(origin=levels(originf), dev=1:n))
  M[cbind(originf, dev)] <- inc.paid
  M})

## triangulo acumulado
cum.triangle <- t(apply(inc.triangle, 1, cumsum))

## pegando os dados da diagonal principal que sao as ultimas obs do triangulo acumulado
latest.paid <- cum.triangle[row(cum.triangle) == n - col(cum.triangle) + 1]

## adicionando os valores acumulados nos dados
Claims$cum.paid <- cum.triangle[with(Claims, cbind(originf, dev))]

## usava isso em outros modelos ver se posso exluir
names(Claims)[3:4] <- c("inc.paid.k", "cum.paid.k")
ids <- with(Claims, cbind(originf, dev))
Claims <- within(Claims,{
  cum.paid.kp1 <- cbind(cum.triangle[,-1], NA)[ids]
  inc.paid.kp1 <- cbind(inc.triangle[,-1], NA)[ids]
  devf <- factor(dev)
}
)

## fatores de desenvolvimento
f <- sapply((n-1):1, function(i) {
  sum( cum.triangle[1:i, n-i+1] ) / sum( cum.triangle[1:i, n-i] )
})

## adicionando 1 para estimar o triangulo completo acumulado
tail <- 1
(f <- c(f, tail))

full.triangle <- cum.triangle
for(k in 1:(n-1)){
  full.triangle[(n-k+1):n, k+1] <- full.triangle[(n-k+1):n,k]*f[k]
}
full.triangle

## pegando os ultimos valores do triangulo acumulado completo 
(ultimate.paid <- full.triangle[,n])
## subtrai o valor acumulado final do valor acumulado do triangulo de vdd 
## tem o valor da reserva
sum(ultimate.paid - latest.paid)

##### DADOS PARA FAZER PREVISAO #####
## gerando matrix com dados para previsao
dados_prev<-matrix(1,nrow=n,ncol=n)

# zerando o que esta abaixo da diagonal principal que e conhecido
dados_prev[row(dados_prev) <= n - col(dados_prev) + 1]<-NA
dados_prev

### colocando no formato padrao do livro para fazer o modelo

allClaims <- data.frame(originf =factor(sort(rep(1:n, n))),
                        dev=rep(1:n,n),
                        inc.paid= as.vector(t(dados_prev)))
## remove na
allClaims<-allClaims %>% na.omit()


#### MODELO SS LOG ANOVA BAYESIANO ERRO NORMAL DADOS SIMULADOS ####
## lista com os dados usados no modelo
data_list <- list(n = length((Claims$inc.paid.k)),
                  n_prev = length(allClaims$inc.paid),
                  m = 20,
                  p = 20,
                  t = as.integer(Claims$originf),
                  h =as.integer(Claims$dev) ,
                  t_prev = as.integer(allClaims$originf),
                  h_prev =as.integer(allClaims$dev) ,
                  #id_time = as.matrix(id_time),
                  y = Claims$inc.paid.k)

## rodando o modelo no stan
anova_normal_ss <- stan("/home//Pessoal/Modelos_tcc_teste_validation/MODELOS USADOS/anova_ss_simulacao.stan",
                         data=data_list,
                         warmup=150000,iter=300000,thin=50,
                         chains=2, control=list(adapt_delta=0.99),cores = getOption("mc.cores", 2))  

## ANALISE BASICA PARA VER CONVERGENCIA E SE CAPTURAMOS OS PARAMETROS
summary(anova_normal_ss)
resultado<-as.matrix(anova_normal_ss)
coefs <-apply(resultado, 2, median)
coefs <-apply(resultado, 2, quantile, probs=c(0.025,0.975))
traceplot(anova_normal_ss,pars=c("sigma_beta","sigma_u","sigma_alpha"))
traceplot(anova_normal_ss,pars=c("ni"))
mu<-extract(anova_normal_ss)$u
sigma2<-extract(anova_normal_ss)$sigma_u
quantile(mu, probs = c(0.025,0.5,0.975))
quantile(sigma2, probs = c(0.025,0.5,0.975))
teste<-extract(anova_normal_ss)$log_lik
soma<-apply(teste,MARGIN=1,FUN=sum)
# lppd <-sum(log(colMeans(exp(log_lik))))
lppd <-sum(log(mean(exp(soma))))
p_waic <- sum(stats::var(soma))
waic <- -2*lppd +2*p_waic

## SALVANDO OS DADOS SIMULADOS E O RESULTADO DO MODELO
dados_simulados<- list(alpha=alpha,beta=beta,epsilon=epsilon,muij=muij,zij=zij,Claims,triangulo_total=triangulo_total)
saveRDS(anova_normal_ss,file="/home//Pessoal/Modelos_tcc_teste_validation/modelos_erro_normal//simulacao_anova_ss_normal.rds")
saveRDS(dados_simulados,file="/home//Pessoal/Modelos_tcc_teste_validation/dados_simulados_anova_ss_normal.rds")

#### MODELO ANOVA SS ERRO NORMAL DADOS CHOY ####
n<-18
Claims <- data.frame(originf = factor(rep(1978:1995, 18:1)),
                     dev=sequence(n:1),
                     inc.paid= c(3323,8332,9572,10172,7631,3855,3252,4433,2188,333,
                                 199,692,311,0.01,405,293,76,14,
                                 3785,10342,8330,7849,2839,3577,1404,1721,1065,
                                 156,35,259,250,420,6,1,0.01,4677,9989,8746,10228,
                                 8572,5787,3855,1445,1612,626,1172,589,438,473,
                                 370,31,5288,8089,12839,11829,7560,6383,4118,3016,1575,1985,2645,
                                 266,38,45,115,2294,9869,10242,13808,8775,5419,2424,
                                 1597,4149,1296,917,295,428,359,3600,7514,8247,9327,
                                 8584,4245,4096,3216,2014,593,1188,691,368,3642,7394,
                                 9838,9733,6377,4884,11920,4188,4492,1760,944,921,
                                 2463,5033,6980,7722,6702,7834,5579,3622,1300,3069,1370,
                                 2267,5959,6175,7051,8102,6339,6978,4396,3107,903,
                                 2009,3700,5298,6885,6477,7570,5855,5751,3871,1860,
                                 5282,3640,7538,5157,5766,6862,2572,2331,3517,5310,
                                 6066,10149,9265,5262,2314,4487,4112,7000,11163,10057,
                                 2607,3952,8228,7895,9317,2595,5403,6579,15546,3155,
                                 4974,7961,2626,5704,2827))

## triangulo de run-off
inc.triangle <- with(Claims, {
  M <- matrix(nrow=n, ncol=n,
              dimnames=list(origin=levels(originf), dev=1:n))
  M[cbind(originf, dev)] <- inc.paid
  M})

## REMOVEMOS 2 DIAGONAIS PARA USAR COMO CONJUNTO DE VALIDACAO DO MODELO
diag1<-c()
m<-18
for (i in 1:18) {
  diag1[i]<-inc.triangle[i,m]
  inc.triangle[i,m] <-NA
  m<-m-1
}

m<-17
diag2<-c()
for (i in 1:17) {
  diag2[i]<-inc.triangle[i,m]
  inc.triangle[i,m] <-NA
  m<-m-1
}

## NOVO TRIANGULO VAI SER 16 X 16
n <- 16
inc.triangle <- inc.triangle[1:n,1:n]

Claims <- data.frame(originf = factor(rep(1978:1993, n:1)),
                     dev=sequence(n:1),
                     inc.paid= (na.omit(as.vector(t(inc.triangle)))))


## triangulo acumulado
cum.triangle <- t(apply(inc.triangle, 1, cumsum))

## pegando os dados da diagonal principal que sao as ultimas obs do triangulo acumulado
latest.paid <- cum.triangle[row(cum.triangle) == n - col(cum.triangle) + 1]

## adicionando os valores acumulados nos dados
Claims$cum.paid <- cum.triangle[with(Claims, cbind(originf, dev))]

## usava isso em outros modelos ver se posso exluir
names(Claims)[3:4] <- c("inc.paid.k", "cum.paid.k")
ids <- with(Claims, cbind(originf, dev))
Claims <- within(Claims,{
  cum.paid.kp1 <- cbind(cum.triangle[,-1], NA)[ids]
  inc.paid.kp1 <- cbind(inc.triangle[,-1], NA)[ids]
  devf <- factor(dev)
}
)

## fatores de desenvolvimento
f <- sapply((n-1):1, function(i) {
  sum( cum.triangle[1:i, n-i+1] ) / sum( cum.triangle[1:i, n-i] )
})

## adicionando 1 para estimar o triangulo completo acumulado
tail <- 1
(f <- c(f, tail))

full.triangle <- cum.triangle
for(k in 1:(n-1)){
  full.triangle[(n-k+1):n, k+1] <- full.triangle[(n-k+1):n,k]*f[k]
}
full.triangle

## pegando os ultimos valores do triangulo acumulado completo 
(ultimate.paid <- full.triangle[,n])
## subtrai o valor acumulado final do valor acumulado do triangulo de vdd 
## tem o valor da reserva
sum(ultimate.paid - latest.paid)

##### DADOS PARA FAZER PREVISAO #####
## gerando matrix com dados para previsao
dados_prev<-matrix(1,nrow=16,ncol=16)

# zerando o que esta abaixo da diagonal principal que e conhecido
dados_prev[row(dados_prev) <= n - col(dados_prev) + 1]<-NA
dados_prev

### colocando no formato padrao do livro para fazer o modelo

allClaims <- data.frame(originf =factor(sort(rep(1978:1993, 16))),
                        dev=rep(1:16,16),
                        inc.paid= as.vector(t(dados_prev)))
## remove na
allClaims<-allClaims %>% na.omit()

#### MODELO SS LOG ANOVA BAYESIANO ERRO NORMAL DADOS CHOY ####
## PASSANDO LOG NOS PAGAMENTOS
Claims$logClaims<-log(Claims$inc.paid.k)

## LISTA COM OS DADOS USADOS NO MODELO 
data_list <- list(n = length((Claims$inc.paid.k)),
                  n_prev = length(allClaims$inc.paid),
                  m = 16,
                  p = 16,
                  t = as.integer(Claims$originf),
                  h =as.integer(Claims$dev) ,
                  t_prev = as.integer(allClaims$originf),
                  h_prev =as.integer(allClaims$dev) ,
                  #id_time = as.matrix(id_time),
                  y = Claims$logClaims)

## rodando o modelo no stan
anova_normal_ss <- stan("/home//Pessoal/Modelos_tcc_teste_validation/MODELOS USADOS/anova_ss_choy.stan",
                         data=data_list,
                         warmup=150000,iter=300000,thin=50,
                         chains=2, control=list(adapt_delta=0.99),cores = getOption("mc.cores", 2))  

## ANALISE BASICA DO RESULTADO PARA VER CONVERGENCIA E VALORES DOS PARAMETROS
summary(anova_normal_ss)
resultado<-as.matrix(anova_normal_ss)
y_new<-extract(anova_normal_ss)$reserva
quantile((y_new))
coefs <-apply(resultado, 2, median)
coefs <-apply(resultado, 2, quantile, probs=c(0.025,0.975))
traceplot(anova_normal_ss,pars=c("u","sigma_u"))
traceplot(anova_normal_ss,pars=c("sigma_beta","sigma_u","sigma_alpha"))

teste<-extract(anova_normal_ss)$log_lik
soma<-apply(teste,MARGIN=1,FUN=sum)
# lppd <-sum(log(colMeans(exp(log_lik))))
lppd <-sum(log(mean(exp(soma))))
p_waic <- sum(stats::var(soma))
waic <- -2*lppd +2*p_waic

waic(extract(anova_normal_ss)$log_lik)

## salvando resultado do modelo
saveRDS(anova_normal_ss,file="/home//Pessoal/Modelos_tcc_teste_validation/modelos_erro_normal/result_choy_ss_normal.rds")

#### MODELO ANOVA SS ERRO NORMAL DADOS COTE ####
n<-10
Claims <- data.frame(originf = factor(rep(1981:1990, 10:1)),
                     dev=sequence(n:1),
                     inc.paid=c(3488 ,11071 ,12690 ,10730, 11582 ,6396 ,2449, 2456, 2418, 584,
                                1169 ,11612 , 7769, 10997, 11261, 4577 ,2866,  727 , 294 ,
                                1478  ,9310 ,14711 , 8780 , 8778 ,6303 ,2969 , 215,
                                1186, 10666, 11061 , 9624 , 9287, 6181, 4537 ,
                                1737 ,12144 ,11640 ,12516 , 5647, 4071,
                                1571, 10582, 15176 ,14503,  9947,
                                1199 ,15878 ,12799 ,14273,
                                1263, 14810 ,12176,
                                986 , 9017 ,
                                683))


## triangulo de run-off
inc.triangle <- with(Claims, {
  M <- matrix(nrow=n, ncol=n,
              dimnames=list(origin=levels(originf), dev=1:n))
  M[cbind(originf, dev)] <- inc.paid
  M})

##  REMOVENDO DUAS DIAGONAIS PARA USAR COMO CONJUNTO DE VALIDACAO DO MODELO
diag1<-c()
m<-10
for (i in 1:10) {
  diag1[i]<-inc.triangle[i,m]
  inc.triangle[i,m] <-NA
  m<-m-1
}

m<-9
diag2<-c()
for (i in 1:9) {
  diag2[i]<-inc.triangle[i,m]
  inc.triangle[i,m] <-NA
  m<-m-1
}

## NOVO TRIANGULO SERA 8X8
n <- 8
inc.triangle <- inc.triangle[1:n,1:n]

Claims <- data.frame(originf = factor(rep(1:n, n:1)),
                     dev=sequence(n:1),
                     inc.paid= (na.omit(as.vector(t(inc.triangle)))))

## triangulo acumulado
cum.triangle <- t(apply(inc.triangle, 1, cumsum))

## pegando os dados da diagonal principal que sao as ultimas obs do triangulo acumulado
latest.paid <- cum.triangle[row(cum.triangle) == n - col(cum.triangle) + 1]

## adicionando os valores acumulados nos dados
Claims$cum.paid <- cum.triangle[with(Claims, cbind(originf, dev))]

## usava isso em outros modelos ver se posso exluir
names(Claims)[3:4] <- c("inc.paid.k", "cum.paid.k")
ids <- with(Claims, cbind(originf, dev))
Claims <- within(Claims,{
  cum.paid.kp1 <- cbind(cum.triangle[,-1], NA)[ids]
  inc.paid.kp1 <- cbind(inc.triangle[,-1], NA)[ids]
  devf <- factor(dev)
}
)

## fatores de desenvolvimento
f <- sapply((n-1):1, function(i) {
  sum( cum.triangle[1:i, n-i+1] ) / sum( cum.triangle[1:i, n-i] )
})

## adicionando 1 para estimar o triangulo completo acumulado
tail <- 1
(f <- c(f, tail))

full.triangle <- cum.triangle
for(k in 1:(n-1)){
  full.triangle[(n-k+1):n, k+1] <- full.triangle[(n-k+1):n,k]*f[k]
}
full.triangle

## pegando os ultimos valores do triangulo acumulado completo 
(ultimate.paid <- full.triangle[,n])
## subtrai o valor acumulado final do valor acumulado do triangulo de vdd 
## tem o valor da reserva
sum(ultimate.paid - latest.paid)

##### DADOS PARA FAZER PREVISAO #####
## gerando matrix com dados para previsao
dados_prev<-matrix(1,nrow=8,ncol=8)

# zerando o que esta abaixo da diagonal principal que e conhecido
dados_prev[row(dados_prev) <= n - col(dados_prev) + 1]<-NA
dados_prev

### colocando no formato padrao do livro para fazer o modelo

allClaims <- data.frame(originf =factor(sort(rep(1:8, 8))),
                        dev=rep(1:8,8),
                        inc.paid= as.vector(t(dados_prev)))
## remove na
allClaims<-allClaims %>% na.omit()


#### MODELO SS LOG ANOVA BAYESIANO ERRO NORMAL ####
## PASSANDO LOG NOS PAGAMENTOS
Claims$logClaims<-log(Claims$inc.paid.k)

#### MODELO SS LOG ANOVA BAYESIANO ERRO NORMAL DADOS COTE ####

## LISTA COM OS DADOS QUE SERAO USADOS NO MODELO
data_list <- list(n = length((Claims$inc.paid.k)),
                  n_prev = length(allClaims$inc.paid),
                  m = 8,
                  p = 8,
                  t = as.integer(Claims$originf),
                  h =as.integer(Claims$dev) ,
                  t_prev = as.integer(allClaims$originf),
                  h_prev =as.integer(allClaims$dev) ,
                  #id_time = as.matrix(id_time),
                  y = Claims$logClaims)

## rodando o modelo no stan
anova_normal_ss <- stan("/home//Pessoal/Modelos_tcc_teste_validation/MODELOS USADOS/anova_ss_choy.stan",
                        data=data_list,
                        warmup=150000,iter=300000,thin=50,
                        chains=2, control=list(adapt_delta=0.99),cores = getOption("mc.cores", 2))

## ANALISE BASICA  PARA VER CONVERGENCIA E VALORES DOS PARAMETROS
summary(anova_normal_ss)
resultado<-as.matrix(anova_normal_ss)
y_new<-rstan::extract(anova_normal_ss)$reserva
quantile((y_new))
coefs <-apply(resultado, 2, median)
coefs <-apply(resultado, 2, quantile, probs=c(0.025,0.975))
traceplot(anova_normal_ss,pars=c("u","sigma_u"))
traceplot(anova_normal_ss,pars=c("sigma_beta","sigma_u","sigma_alpha"))

teste<-extract(anova_normal_ss)$log_lik
soma<-apply(teste,MARGIN=1,FUN=sum)
# lppd <-sum(log(colMeans(exp(log_lik))))
lppd <-sum(log(mean(exp(soma))))
p_waic <- sum(stats::var(soma))
waic <- -2*lppd +2*p_waic
waic(extract(anova_normal_ss)$log_lik)

## SALVANDO RESULTADO DO MODELO 
saveRDS(anova_normal_ss,file="/home//Pessoal/Modelos_tcc_teste_validation/modelos_erro_normal/result_cote_ss_normal.rds")


#### MODELO ANOVA SS ERRO NORMAL DADOS SOA ####

n<-12
Claims <- data.frame(originf = factor(rep(1:12, 12:1)),
                     dev=sequence(n:1),
                     inc.paid= c(1750,2500,1050,700,125,70,35,25,15,30,5,0,
                                 1500,2750,1550,725,100,75,30,10,25,25,5,
                                 1600,2800,1650,675,125,50,25,40,20,10,
                                 1125,2650,1550,740,70,45,55,15,5,
                                 1900,2200,1475,650,60,60,40,15,
                                 1500,2600,1350,650,100,60,30,
                                 2000,1950,1570,750,90,60,
                                 1700,2300,1650,700,115,
                                 1400,2900,1440,725,
                                 2225,2400,1525,
                                 1700,1950,
                                 1575))

## triangulo de run-off
inc.triangle <- with(Claims, {
  M <- matrix(nrow=n, ncol=n,
              dimnames=list(origin=levels(originf), dev=1:n))
  M[cbind(originf, dev)] <- inc.paid
  M})

diag1<-c()
m<-12
for (i in 1:12) {
  diag1[i]<-inc.triangle[i,m]
  inc.triangle[i,m] <-NA
  m<-m-1
}

m<-11
diag2<-c()
for (i in 1:11) {
  diag2[i]<-inc.triangle[i,m]
  inc.triangle[i,m] <-NA
  m<-m-1
}

## NOVO TRIANGULO SERA 8X8
n <- 10
inc.triangle <- inc.triangle[1:n,1:n]

Claims <- data.frame(originf = factor(rep(1:n, n:1)),
                     dev=sequence(n:1),
                     inc.paid= (na.omit(as.vector(t(inc.triangle)))))

## triangulo acumulado
cum.triangle <- t(apply(inc.triangle, 1, cumsum))

## pegando os dados da diagonal principal que sao as ultimas obs do triangulo acumulado
latest.paid <- cum.triangle[row(cum.triangle) == n - col(cum.triangle) + 1]

## adicionando os valores acumulados nos dados
Claims$cum.paid <- cum.triangle[with(Claims, cbind(originf, dev))]

## usava isso em outros modelos ver se posso exluir
names(Claims)[3:4] <- c("inc.paid.k", "cum.paid.k")
ids <- with(Claims, cbind(originf, dev))
Claims <- within(Claims,{
  cum.paid.kp1 <- cbind(cum.triangle[,-1], NA)[ids]
  inc.paid.kp1 <- cbind(inc.triangle[,-1], NA)[ids]
  devf <- factor(dev)
}
)

## fatores de desenvolvimento
f <- sapply((n-1):1, function(i) {
  sum( cum.triangle[1:i, n-i+1] ) / sum( cum.triangle[1:i, n-i] )
})

## adicionando 1 para estimar o triangulo completo acumulado
tail <- 1
(f <- c(f, tail))

full.triangle <- cum.triangle
for(k in 1:(n-1)){
  full.triangle[(n-k+1):n, k+1] <- full.triangle[(n-k+1):n,k]*f[k]
}
full.triangle

## pegando os ultimos valores do triangulo acumulado completo 
(ultimate.paid <- full.triangle[,n])
## subtrai o valor acumulado final do valor acumulado do triangulo de vdd 
## tem o valor da reserva
sum(ultimate.paid - latest.paid)

##### DADOS PARA FAZER PREVISAO #####
## gerando matrix com dados para previsao
dados_prev<-matrix(1,nrow=10,ncol=10)

# zerando o que esta abaixo da diagonal principal que e conhecido
dados_prev[row(dados_prev) <= n - col(dados_prev) + 1]<-NA
dados_prev

### colocando no formato padrao do livro para fazer o modelo

allClaims <- data.frame(originf =factor(sort(rep(1:10, 10))),
                        dev=rep(1:10,10),
                        inc.paid= as.vector(t(dados_prev)))
## remove na
allClaims<-allClaims %>% na.omit()


#### MODELO SS LOG ANOVA BAYESIANO ERRO NORMAL ####
## PASSANDO LOG NOS PAGAMENTOS
Claims$logClaims<-log(Claims$inc.paid.k)

#### MODELO SS LOG ANOVA BAYESIANO ERRO NORMAL DADOS SOA ####

## LISTA COM OS DADOS QUE SERAO USADOS NO MODELO
data_list <- list(n = length((Claims$inc.paid.k)),
                  n_prev = length(allClaims$inc.paid),
                  m = 10,
                  p = 10,
                  t = as.integer(Claims$originf),
                  h =as.integer(Claims$dev) ,
                  t_prev = as.integer(allClaims$originf),
                  h_prev =as.integer(allClaims$dev) ,
                  #id_time = as.matrix(id_time),
                  y = Claims$logClaims)

## rodando o modelo no stan
anova_normal_ss <- stan("/home//Pessoal/Modelos_tcc_teste_validation/MODELOS USADOS/anova_ss_choy.stan",
                        data=data_list,
                        warmup=150000,iter=300000,thin=50,
                        chains=2, control=list(adapt_delta=0.99),cores = getOption("mc.cores", 2))

## ANALISE BASICA  PARA VER CONVERGENCIA E VALORES DOS PARAMETROS
summary(anova_normal_ss)
resultado<-as.matrix(anova_normal_ss)
y_new<-rstan::extract(anova_normal_ss)$reserva
quantile((y_new))
coefs <-apply(resultado, 2, median)
coefs <-apply(resultado, 2, quantile, probs=c(0.025,0.5,0.975))
traceplot(anova_normal_ss,pars=c("u","sigma_u"))
traceplot(anova_normal_ss,pars=c("sigma_beta","sigma_u","sigma_alpha"))

teste<-extract(anova_normal_ss)$log_lik
soma<-apply(teste,MARGIN=1,FUN=sum)
# lppd <-sum(log(colMeans(exp(log_lik))))
lppd <-sum(log(mean(exp(soma))))
p_waic <- sum(stats::var(soma))
waic <- -2*lppd +2*p_waic
waic(extract(anova_normal_ss)$log_lik)

## SALVANDO RESULTADO DO MODELO 
saveRDS(anova_normal_ss,file="/home//Pessoal/Modelos_tcc_teste_validation/modelos_erro_normal/result_soa_ss_normal.rds")

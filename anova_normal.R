#### MODELO ANOVA ERRO NORMAL COM DADOS SIMULADOS ####

##### DADOS SIMULADOS #####
rm(list=ls())
library(dplyr)
library(lattice)
library(brms)
library(ChainLadder)
library(rstanarm)
library(rstan)
#semente
set.seed(123)

# media, variancia
mu<- 200;sigma2<-1
# numero de linhas e colunas do triangulo de run-off
i<-20
# numero de infos que vou gerar para preencer o triangulo
n<-400

## gerando os coeficientes dos efeitos das linhas N(0,4)
alpha<-rnorm(i-1,0,2)

## gerando os coeficientes dos efeito das colunas N(0,16)
beta<-rnorm(i-1,0,4)

## restricao que a soma dos coeficientes tem que ser 0
alpha[i]<- -sum(alpha)
beta[i]<- -sum(beta)

## erro N(0,1)
epsilon<- rnorm(400,0,1)

## gerando os muij segundo o modelo log anova
muij<-mu+matrix(rep(alpha,20),ncol=i)+t(matrix(rep(beta,20),ncol=i))

## gerando os valores de zij zij= muij+ erro
zij<-muij+matrix(epsilon,ncol=i)

triangulo_completo <- zij
## tornando NA os valores abaixo da diagonal do triangulo que sao desconhecidos
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
allClaims <- data.frame(originf = sort(rep(1:n, n)),
                        dev = rep(1:n,n),
                        inc.paid= (as.vector(t(dados_prev))))
## remove na
allClaims<-allClaims %>% na.omit()
## colocando como fator
allClaims$dev<-as.factor(allClaims$dev)
allClaims$originf<-as.factor(allClaims$originf)

## matriz com 0 e 1 para fazer as previsoes para as colunas
X_b_prev <- as.data.frame(model.matrix(~ dev, data=allClaims,
                                  contrasts.arg = list(dev = contrasts(allClaims$dev, contrasts = FALSE))))[,-1]
X_b_prev<-cbind(0,X_b_prev)

## matrix com 0 e 1 para alfa das linhas
X_a_prev<-as.data.frame( model.matrix(~originf, data=allClaims,
                                 contrasts.arg = list(originf = contrasts(allClaims$originf, contrasts = FALSE))))[,-1]
X_a_prev<-cbind(0,X_a_prev)

## infos para o modelo
K_prev<- ncol(X_b_prev)
N_prev<- nrow(X_b_prev)

##### MODELO LOG ANOVA ERRO NORMAL #####

## colocando dev como fator
Claims$dev<-as.factor(Claims$dev)

## Model matrix com todos os coeficientes
## beta das colunas
X_b <- as.data.frame(model.matrix(~ dev, data=Claims,
                                  contrasts.arg = list(dev = contrasts(Claims$dev, contrasts = FALSE))))[,-1]

## alfa das linhas
X_a<-as.data.frame( model.matrix(~originf, data=Claims,
                                 contrasts.arg = list(originf = contrasts(Claims$originf, contrasts = FALSE))))[,-1]
## infos para o modelo
Y <- Claims$inc.paid.k
K<- ncol(X_b)
N<- nrow(X_b)

## rodando o modelo no stan 
anova_normal <- stan("/home//Pessoal/Modelos_tcc_teste_validation/MODELOS USADOS/anova_normal_simulacao.stan",
                     data=list(X_alpha=X_a,X_beta=X_b,Y=Y,K_alpha=K,K_beta=K,N=N,X_alpha_prev=X_a_prev,
                               X_beta_prev=X_b_prev, K_alpha_prev=K_prev,N_prev=N_prev,K_beta_prev=K_prev),
                     warmup=10000,iter=40000,thin=50,chains=2, control=list(adapt_delta=0.99),
                     cores = getOption("mc.cores", 2))  

modelo<-extract(anova_normal)

betas<- extract(anova_normal)$beta
alphas<- extract(anova_normal)$alpha
boxplot(betas)
boxplot(alphas)
teste<-extract(anova_normal)$log_lik
soma<-apply(teste,MARGIN=1,FUN=sum)
# lppd <-sum(log(colMeans(exp(log_lik))))
lppd <-sum(log(mean(exp(soma))))
p_waic <- sum(stats::var(soma))
waic <- -2*lppd +2*p_waic
quantile(modelo$reserva)
waic(extract(anova_normal)$log_lik)


View(summary(anova_normal))
resultado<-as.matrix(anova_normal)
coefs <-apply(resultado, 2, median)

dados_simulados<- list(alpha=alpha,beta=beta,epsilon=epsilon,muij=muij,zij=zij,Claims,triangulo_completo=triangulo_completo)
saveRDS(anova_normal,file="/home//Pessoal/Modelos_tcc_teste_validation/modelos_erro_normal/result_simulacao_anova_normal.rds")
saveRDS(dados_simulados,file="/home//Pessoal/Modelos_tcc_teste_validation/dados_simulados_anova_normal.rds")
y_new_pred<-data.frame(as.matrix(extract(anova_normal)$y_new))

final.triangle<-apply(y_new_pred,2,median)

#### MODELO ANOVA ERRO NORMAL DADOS CHOY ####
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
inc.triangle <- inc.triangle[1:16,1:16]

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


###### DADOS PARA FAZER PREVISAO ######
dados_prev<-matrix(1,nrow=16,ncol=16)

# zerando o que esta abaixo da diagonal principal que e conhecido
dados_prev[row(dados_prev) <= n - col(dados_prev) + 1]<-NA

allClaims <- data.frame(originf =factor(sort(rep(1978:1993, 16))),
                        dev=rep(1:16,16),
                        inc.paid= as.vector(t(dados_prev)))
allClaims<-allClaims %>% na.omit()
allClaims$dev<-as.factor(allClaims$dev)
allClaims$originf<-as.factor(allClaims$originf)
X_b_prev <- as.data.frame(model.matrix(~ dev, data=allClaims,
                                       contrasts.arg = list(dev = contrasts(allClaims$dev, contrasts = FALSE))))[,-1]
X_b_prev<-cbind(0,X_b_prev)
## alfa das linhas
X_a_prev<-as.data.frame( model.matrix(~originf, data=allClaims,
                                      contrasts.arg = list(originf = contrasts(allClaims$originf, contrasts = FALSE))))[,-1]

## infos para o modelo
K_prev<- ncol(X_b_prev)
N_prev<- nrow(X_b_prev)

##### MODELO LOG ANOVA ERRO NORMAL #####
Claims$logClaims<-log(Claims$inc.paid.k)
Claims$dev<-as.factor(Claims$dev)

## teste que eu estava fazendo apagar depois
anov<-lm(Claims$logClaims~dev+originf, data=Claims,contrasts=list(originf="contr.sum",dev="contr.sum"))
summary(anov)

## Model matrix com todos os coeficientes
## beta das colunas
X_b <- as.data.frame(model.matrix(~ dev, data=Claims,
                                  contrasts.arg = list(dev = contrasts(Claims$dev, contrasts = FALSE))))[,-1]

## alfa das linhas
X_a<-as.data.frame( model.matrix(~originf, data=Claims,
                                 contrasts.arg = list(originf = contrasts(Claims$originf, contrasts = FALSE))))[,-1]
X_a <- X_a[,-c(17,18)]
## infos para o modelo
Y <- Claims$logClaims
K<- ncol(X_b)
N<- nrow(X_b)

## rodando o modelo no stan
anova_normal <- stan("/home//Pessoal/Modelos_tcc_teste_validation/MODELOS USADOS/anova_normal_choy.stan",
                     data=list(X_alpha=X_a,X_beta=X_b,Y=Y,K_alpha=K,K_beta=K,N=N,X_alpha_prev=X_a_prev,
                               X_beta_prev=X_b_prev, K_alpha_prev=K_prev,N_prev=N_prev,K_beta_prev=K_prev),
                     warmup=10000,iter=40000,thin=50,
                     chains=2, control=list(adapt_delta=0.99),cores = getOption("mc.cores", 2))  


modelo<-extract(anova_normal)
alpha<-extract(anova_normal)$alpha
beta<-extract(anova_normal)$beta
boxplot(alpha)
boxplot(beta)


teste<-extract(anova_normal)$log_lik
soma<-apply(teste,MARGIN=1,FUN=sum)
# lppd <-sum(log(colMeans(exp(log_lik))))
lppd <-sum(log(mean(exp(soma))))
p_waic <- sum(stats::var(soma))
waic <- -2*lppd +2*p_waic
quantile((modelo$reserva))
waic(extract(anova_normal)$log_lik)


View(summary(anova_normal))
resultado<-as.matrix(anova_normal)
coefs <-apply(resultado, 2, median)
y_new<-(extract(anova_normal)$y_new)
View(extract(anova_normal)$y_new)
y_novo<-apply(y_new,MARGIN=2,FUN=median)
saveRDS(anova_normal,file="/home//Pessoal/Modelos_tcc_teste_validation/modelos_erro_normal/result_choy_anova_normal.rds")


#### MODELO ANOVA ERRO NORMAL DADOS COTE ####
n<-10
Claims <- data.frame(originf = factor(rep(1:10, 10:1)),
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


###### DADOS PARA FAZER PREVISAO ######
dados_prev<-matrix(1,nrow=8,ncol=8)

# zerando o que esta abaixo da diagonal principal que e conhecido
dados_prev[row(dados_prev) <= n - col(dados_prev) + 1]<-NA

allClaims <- data.frame(originf =factor(sort(rep(1:8, 8))),
                        dev=rep(1:8,8),
                        inc.paid= as.vector(t(dados_prev)))
allClaims<-allClaims %>% na.omit()
allClaims$dev<-as.factor(allClaims$dev)
allClaims$originf<-as.factor(allClaims$originf)
X_b_prev <- as.data.frame(model.matrix(~ dev, data=allClaims,
                                       contrasts.arg = list(dev = contrasts(allClaims$dev, contrasts = FALSE))))[,-1]
X_b_prev<-cbind(0,X_b_prev)
## alfa das linhas
X_a_prev<-as.data.frame( model.matrix(~originf, data=allClaims,
                                      contrasts.arg = list(originf = contrasts(allClaims$originf, contrasts = FALSE))))[,-1]

## infos para o modelo
K_prev<- ncol(X_b_prev)
N_prev<- nrow(X_b_prev)

##### MODELO LOG ANOVA ERRO NORMAL #####
Claims$logClaims<-log(Claims$inc.paid.k)
Claims$dev<-as.factor(Claims$dev)

## teste que eu estava fazendo apagar depois
anov<-lm(Claims$logClaims~dev+originf, data=Claims,contrasts=list(originf="contr.sum",dev="contr.sum"))
summary(anov)

## Model matrix com todos os coeficientes
## beta das colunas
X_b <- as.data.frame(model.matrix(~ dev, data=Claims,
                                  contrasts.arg = list(dev = contrasts(Claims$dev, contrasts = FALSE))))[,-1]

## alfa das linhas
X_a<-as.data.frame( model.matrix(~originf, data=Claims,
                                 contrasts.arg = list(originf = contrasts(Claims$originf, contrasts = FALSE))))[,-1]
X_a <- X_a[,-c(9,10)]
## infos para o modelo
Y <- Claims$logClaims
K<- ncol(X_b)
N<- nrow(X_b)

## rodando o modelo no stan
anova_normal <- stan("/home//Pessoal/Modelos_tcc_teste_validation/MODELOS USADOS/anova_normal_choy.stan",
                     data=list(X_alpha=X_a,X_beta=X_b,Y=Y,K_alpha=K,K_beta=K,N=N,X_alpha_prev=X_a_prev,
                               X_beta_prev=X_b_prev, K_alpha_prev=K_prev,N_prev=N_prev,K_beta_prev=K_prev),
                     warmup=10000,iter=40000,thin=50,
                     chains=2, control=list(adapt_delta=0.99),cores = getOption("mc.cores", 2))  


modelo<-extract(anova_normal)
alpha<-extract(anova_normal)$alpha
beta<-extract(anova_normal)$beta
boxplot(alpha)
boxplot(beta)


teste<-extract(anova_normal)$log_lik
soma<-apply(teste,MARGIN=1,FUN=sum)
# lppd <-sum(log(colMeans(exp(log_lik))))
lppd <-sum(log(mean(exp(soma))))
p_waic <- sum(stats::var(soma))
waic <- -2*lppd +2*p_waic
quantile((modelo$reserva))
waic(extract(anova_normal)$log_lik)


View(summary(anova_normal))
resultado<-as.matrix(anova_normal)
coefs <-apply(resultado, 2, median)
y_new<-(extract(anova_normal)$y_new)
View(extract(anova_normal)$y_new)
y_novo<-apply(y_new,MARGIN=2,FUN=median)
saveRDS(anova_normal,file="/home//Pessoal/Modelos_tcc_teste_validation/modelos_erro_normal/result_cote_anova_normal.rds")


#### MODELO ANOVA ERRO NORMAL DADOS SOA ####
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


###### DADOS PARA FAZER PREVISAO ######
dados_prev<-matrix(1,nrow=10,ncol=10)

# zerando o que esta abaixo da diagonal principal que e conhecido
dados_prev[row(dados_prev) <= n - col(dados_prev) + 1]<-NA

allClaims <- data.frame(originf =factor(sort(rep(1:10, 10))),
                        dev=rep(1:10,10),
                        inc.paid= as.vector(t(dados_prev)))
allClaims<-allClaims %>% na.omit()
allClaims$dev<-as.factor(allClaims$dev)
allClaims$originf<-as.factor(allClaims$originf)
X_b_prev <- as.data.frame(model.matrix(~ dev, data=allClaims,
                                       contrasts.arg = list(dev = contrasts(allClaims$dev, contrasts = FALSE))))[,-1]
X_b_prev<-cbind(0,X_b_prev)
## alfa das linhas
X_a_prev<-as.data.frame( model.matrix(~originf, data=allClaims,
                                      contrasts.arg = list(originf = contrasts(allClaims$originf, contrasts = FALSE))))[,-1]

## infos para o modelo
K_prev<- ncol(X_b_prev)
N_prev<- nrow(X_b_prev)

##### MODELO LOG ANOVA ERRO NORMAL #####
Claims$logClaims<-log(Claims$inc.paid.k)
Claims$dev<-as.factor(Claims$dev)

## teste que eu estava fazendo apagar depois
anov<-lm(Claims$logClaims~dev+originf, data=Claims,contrasts=list(originf="contr.sum",dev="contr.sum"))
summary(anov)

## Model matrix com todos os coeficientes
## beta das colunas
X_b <- as.data.frame(model.matrix(~ dev, data=Claims,
                                  contrasts.arg = list(dev = contrasts(Claims$dev, contrasts = FALSE))))[,-1]

## alfa das linhas
X_a<-as.data.frame( model.matrix(~originf, data=Claims,
                                 contrasts.arg = list(originf = contrasts(Claims$originf, contrasts = FALSE))))[,-1]
# X_a <- X_a[,-c(10)]
## infos para o modelo
Y <- Claims$logClaims
K<- ncol(X_b)
N<- nrow(X_b)

## rodando o modelo no stan
anova_normal <- stan("/home//Pessoal/Modelos_tcc_teste_validation/MODELOS USADOS/anova_normal_choy.stan",
                     data=list(X_alpha=X_a,X_beta=X_b,Y=Y,K_alpha=K,K_beta=K,N=N,X_alpha_prev=X_a_prev,
                               X_beta_prev=X_b_prev, K_alpha_prev=K_prev,N_prev=N_prev,K_beta_prev=K_prev),
                     warmup=10000,iter=40000,thin=50,
                     chains=2, control=list(adapt_delta=0.99),cores = getOption("mc.cores", 2))  


modelo<-extract(anova_normal)
alpha<-extract(anova_normal)$alpha
beta<-extract(anova_normal)$beta
boxplot(alpha)
boxplot(beta)


teste<-extract(anova_normal)$log_lik
soma<-apply(teste,MARGIN=1,FUN=sum)
# lppd <-sum(log(colMeans(exp(log_lik))))
lppd <-sum(log(mean(exp(soma))))
p_waic <- sum(stats::var(soma))
waic <- -2*lppd +2*p_waic
quantile((modelo$reserva))
waic(extract(anova_normal)$log_lik)


View(summary(anova_normal))
resultado<-as.matrix(anova_normal)
coefs <-apply(resultado, 2, median)
y_new<-(extract(anova_normal)$y_new)
View(extract(anova_normal)$y_new)
y_novo<-apply(y_new,MARGIN=2,FUN=median)
saveRDS(anova_normal,file="/home//Pessoal/Modelos_tcc_teste_validation/modelos_erro_normal/result_soa_anova_normal.rds")

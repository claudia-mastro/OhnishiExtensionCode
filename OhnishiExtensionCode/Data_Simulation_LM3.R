source("~/OhnishiExtensionCode/effects_LM3.R")
##############
#Seed
##############
set.seed(1222)
################
#Global Settings
################

N<-rep(Nj,
       times = J)
Ns <- list(0)
Ns[[1]] <- 1:N[1]
for (j in 2:J) {
  Ns[[j]] <- (Ns[[j-1]][N[j-1]]+1):(Ns[[j-1]][N[j-1]]+N[j])
}

###########
#Predictors
###########
x_long<-c(rnorm(n = sum(N)/6, mean=-12, sd=0.5),
          rnorm(n = sum(N)/6, mean=-8, sd=0.5),
          rnorm(n = sum(N)/6, mean=-4, sd=1),
          rnorm(n = sum(N)/6, mean=4, sd=1),
          rnorm(n = sum(N)/6, mean=8, sd=0.5),
          rnorm(n = sum(N)/6, mean=12, sd=0.5))
x_long <- sample(x_long, sum(N), replace=TRUE)
x_long <- (x_long - mean(x_long))/sd(x_long)
# x_long <- cbind(x_long, rbinom(sum(N), 1, 0.4))
# x_long <- cbind(x_long, rbinom(sum(N), 1, 0.8))
q_long<-cbind(1, x_long)

################
#True Parameters
################
beta_true<-matrix(NA,
                  nrow = 1,
                  ncol = 12)
beta_true <- c(2, 4, 6, -9, 8, -6, -4, -2, 5, 7, 10, -8)
sigma2_true <- 0.2

delta_h_true<-c(-0.8, 1.6)
tau2_h_true<-1
logit_h_true<-rnorm(n = sum(N),
                     mean = (q_long%*%delta_h_true),
                     sd = sqrt(tau2_h_true))
etah_true<-exp(logit_h_true)/(1.00 + exp(logit_h_true))

delta_l_true<-c(0.2, -0.6)
tau2_l_true<-1
logit_l_true<-rnorm(n = sum(N),
                     mean = (q_long%*%delta_l_true),
                     sd = sqrt(tau2_l_true))
etal_true<-(etah_true + exp(logit_l_true))/(1.00 + exp(logit_l_true))

ph0_true <- 0.86
ph1_true <- 0.94
pl0_true <- 0.92
pl1_true <- 0.82

lamh0_true <- rbinom(sum(N), 1, ph0_true)
lamh1_true <- rbinom(sum(N), 1, ph1_true)
laml0_true <- rbinom(sum(N), 1, pl0_true)
laml1_true <- rbinom(sum(N), 1, pl1_true)
laml0_true[lamh0_true==1 & lamh1_true==0] <- 1
laml1_true[lamh0_true==1 & lamh1_true==0] <- 0

h_true <- lamh0_true * etah_true ^ lamh1_true
l_true <- laml0_true * etal_true ^ laml1_true


#####
#Data
#####
set.seed(id)
M<-J
a<-sample(seq(0, 1, 1/(nalpha-1)),
          size = J,
          replace = TRUE)
Z<-list(0)
T<-rep(NA,
       times = J)
for(j in 1:J){
  
  Z[[j]]<-rbinom(n = N[j],
                 size = 1,
                 prob = a[j])
  T[j]<-sum(Z[[j]])
  a[j]<-T[j]/N[j]
  
}
Z_long<-Z[[1]]
T_long<-rep(a[1], 
            times = N[1])
for(j in 2:J){
  
  Z_long<-c(Z_long,
            Z[[j]])
  T_long<-c(T_long,
            rep(a[j],
                times = N[j]))
  
}

D<-list(0)
S<-rep(NA,
       times = J)
G_true <- rep(NA, sum(N))
G_true[T_long > l_true] <- 1
G_true[T_long < h_true] <- 2
G_true[T_long >= h_true & T_long <= l_true] <- 3

Ji <- sort(rep(1:J, times=Nj))
dta <- data.table(Ji, T_long, Z_long, h_true, l_true)

dta[T_long > l_true, G_long_true := 1]
dta[T_long < h_true, G_long_true := 2]
dta[T_long >= h_true & T_long <= l_true, G_long_true := 3]

dta[G_long_true==1, D_long := 1]
dta[G_long_true==2, D_long := 0]
dta[G_long_true==3, D_long := Z_long]

dta[,S_long := mean(D_long), by=Ji]

eff.a <- 0.8
eff.s <- 0.4
eff.sp <- 0.8

dta[eff.a > l_true, G_eff := 1]
dta[eff.a < h_true, G_eff := 2]
dta[eff.a >= h_true & eff.a <= l_true, G_eff := 3]
attach(dta, warn.conflicts=FALSE)

W_true<-cbind(1, S_long, T_long, Z_long, h_true, l_true, S_long*h_true, S_long*l_true,
         T_long*h_true, T_long*l_true, Z_long*h_true, Z_long*l_true)

Y0 <- Y1 <- Y0p <- Y_long <- rep(NA, sum(N))
for(j in 1:sum(N)) {
  W0p <- W0 <- W1 <- W_true[j,]
  
  W0[2] <- W1[2] <- eff.s
  W0p[2] <- eff.sp
  
  W0[7] <- W1[7] <- eff.s*h_true[j]
  W0p[7]  <- eff.sp*h_true[j]
  
  W0[8] <- W1[8] <- eff.s*l_true[j]
  W0p[8] <- eff.sp*l_true[j]
  
  W0[3] <- W1[3] <- W0p[3] <- eff.a
  W0[9] <- W1[9] <- W0p[9] <- eff.a*h_true[j]
  W0[10] <- W1[10] <- W0p[10] <- eff.a*l_true[j]
  
  W0[4] <- W0[11] <- W0p[12] <- W0p[4] <- W0p[11] <- W0p[12] <- 0
  W1[4] <- 1
  W1[11] <- h_true[j]
  W1[12] <- l_true[j]
  
  mu0 <- W0 %*% beta_true
  mu1 <- W1 %*% beta_true
  mu0p <- W0p %*% beta_true
  
  eps <- rnorm(1, 0, sqrt(sigma2_true))
  
  ## Potential Outcomes
  Y0[j]  <- mu0  + eps
  Y1[j]  <- mu1  + eps
  Y0p[j] <- mu0p  + eps
  
  ## Observed Data
  mu <- W_true[j,]%*%beta_true
  Y_long[j] <- mu + eps
}


CADE.true <- mean((Y1 - Y0)[G_eff == 3])
CASE.true <- mean((Y0 - Y0p)[G_eff == 3])


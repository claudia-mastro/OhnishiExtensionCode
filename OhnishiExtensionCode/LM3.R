###############
#Packages
###############
library(mnormt)
library(matrixStats)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
id <- as.integer(args[1])
print(id)
J <- as.integer(args[2])
print(J)
Nj <- as.integer(args[3])
print(Nj)
nalpha <- as.integer(args[4])
print(nalpha)
v <- paste0("LM3", J*Nj)
print(v)
source("~/OhnishiExtensionCode/Data_Simulation_LM3.R")

#########
#Function
#########
logSumExp<-function(x){
  m<-max(x)
  val<-m + 
    log(sum(exp(x - m)))
  return(val)
  
}

##############
#Seed
##############
set.seed(id)

################
#Global Settings
################
mcmc_samples<-50000
burnin <- 25000
thin <- 10
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

sigma2_beta<-100.00^2
a_sigma2<-3
b_sigma2<-2
sigma2_delta<-100.00^2
a_tau2<-3
b_tau2<-2
shape_tau2_update<-sum(N)/2.00 +
  a_tau2
a_pi <- 1
b_pi <- 1

###########
#Parameters
###########
beta <- rep(NA, 12)
beta <- list(beta)[rep(1,mcmc_samples)]

sigma2 <- NA
sigma2 <- list(sigma2)[rep(1,mcmc_samples)]

h <- l <- rep(NA, sum(N))
h <- list(h)[rep(1,mcmc_samples)]
l <- list(l)[rep(1,mcmc_samples)]

etah <- etal <- rep(NA, sum(N))
etah <- list(etah)[rep(1,mcmc_samples)]
etal <- list(etal)[rep(1,mcmc_samples)]

delta_h <- delta_l <- rep(NA, ncol(q_long))
delta_h <- list(delta_h)[rep(1,mcmc_samples)]
delta_l <- list(delta_l)[rep(1,mcmc_samples)]

tau2_h <- tau2_l <- NA
tau2_h <- list(tau2_h)[rep(1,mcmc_samples)]
tau2_l <- list(tau2_l)[rep(1,mcmc_samples)]

lamh0 <- lamh1 <- laml0 <- laml1 <- rep(NA, sum(N))
lamh0 <- list(lamh0)[rep(1,mcmc_samples)]
lamh1 <- list(lamh1)[rep(1,mcmc_samples)]
laml0 <- list(laml0)[rep(1,mcmc_samples)]
laml1 <- list(laml1)[rep(1,mcmc_samples)]

ph0 <- ph1 <- pl0 <- pl1 <- NA
ph0 <- list(ph0)[rep(1,mcmc_samples)]
ph1 <- list(ph1)[rep(1,mcmc_samples)]
pl0 <- list(pl0)[rep(1,mcmc_samples)]
pl1 <- list(pl1)[rep(1,mcmc_samples)]

CADE<-rep(NA, length(iters))
CASE<-rep(NA, length(iters))
pC <- rep(NA, length(iters))

###############
#Initial Values
###############
beta[[1]] <- rep(0.00, 12)

sigma2[[1]] <- 1

delta_h[[1]] <- delta_l[[1]] <- rep(0.00, ncol(q_long))
tau2_h[[1]] <- tau2_l[[1]] <- 1

logit_etah <- rnorm(n = sum(N),
                  mean = (q_long%*%delta_h[[1]]),
                  sd = sqrt(tau2_h[[1]]))
etah[[1]] <- 1.00/(1.00 + exp(-logit_etah))

logit_etal <- rnorm(n = sum(N),
                  mean = (q_long%*%delta_l[[1]]),
                  sd = sqrt(tau2_l[[1]]))
etal[[1]] <- (etah[[1]] + exp(logit_etal))/(1.00 + exp(logit_etal))

ph0[[1]] <- pl0[[1]] <- 0.9
ph1[[1]] <- pl1[[1]] <- 0.9

lamh0[[1]] <- rbinom(sum(N), 1, ph0[[1]])
lamh1[[1]] <- rbinom(sum(N), 1, ph1[[1]])
laml0[[1]] <- rbinom(sum(N), 1, pl0[[1]])
laml1[[1]] <- rbinom(sum(N), 1, pl1[[1]])

h[[1]] <- lamh0[[1]]*etah[[1]]^lamh1[[1]]
l[[1]] <- laml0[[1]]*etal[[1]]^laml1[[1]]

####################
#Metropolis Settings
####################
metrop_sd_h <- 0.05
acctot_h <- rep(1, times = sum(N))
acctot_h <- list(acctot_h)[rep(1,mcmc_samples)]

metrop_sd_l <- 0.05
acctot_l <- rep(1, times = sum(N))
acctot_l <- list(acctot_l)[rep(1,mcmc_samples)]

W <- cbind(1, S_long, T_long, Z_long, h[[1]], l[[1]], 
           S_long*h[[1]], S_long*l[[1]],
           T_long*h[[1]], T_long*l[[1]],
           Z_long*h[[1]], Z_long*l[[1]])

###################
#Main Sampling Loop
###################
for(s in 2:mcmc_samples){ 
  print(s)
  #############
  #beta, sigma2
  #############
  cov_beta <- chol2inv(chol(crossprod(W)/sigma2[[s-1]] + diag(12)/sigma2_beta))
  mu_beta <- cov_beta%*%(crossprod(W, Y_long))/sigma2[[s-1]]
  beta[[s]] <- rmnorm(n = 1, mean = mu_beta, varcov = cov_beta)
  
  rate <- crossprod(Y_long - W%*%beta[[s]])/2.00 + b_sigma2
  shape <- length(Y_long)/2.00 + a_sigma2
  sigma2[[s]] <- 1.00/rgamma(n = 1, shape = shape, rate = 1/rate)
  
  ###
  #pi
  ###
  ph0[[s]] <- rbeta(1, a_pi + sum(lamh0[[s-1]]), b_pi + sum(N) - sum(lamh0[[s-1]]))
  ph1[[s]] <- rbeta(1, a_pi + sum(lamh1[[s-1]]), b_pi + sum(N) - sum(lamh1[[s-1]]))
  pl0[[s]] <- rbeta(1, a_pi + sum(laml0[[s-1]]), b_pi + sum(N) - sum(laml0[[s-1]]))
  pl1[[s]] <- rbeta(1, a_pi + sum(laml1[[s-1]]), b_pi + sum(N) - sum(laml1[[s-1]]))
  
  #######
  #lambda
  #######
  lamh0[[s]] <- sapply(ph0[[s]]^(1-(Z_long==1 & D_long==0)), function(p) rbinom(1, 1, p))
  lamh1[[s]] <- sapply(ph1[[s]]^(1-(Z_long==0 & D_long==1)), function(p) rbinom(1, 1, p))
  laml0[[s]] <- sapply(ph0[[s]]^(1-(Z_long==1 & D_long==0)), function(p) rbinom(1, 1, p))
  laml1[[s]] <- sapply(ph1[[s]]^(1-(Z_long==0 & D_long==1)), function(p) rbinom(1, 1, p))
  
  laml0[[s]][lamh0[[s]]==1 & lamh1[[s]]==0] <- 1
  laml1[[s]][lamh0[[s]]==1 & lamh1[[s]]==0] <- 0
   
  ####
  #eta
  ####
  logit_etah_old <- logit_etal_old <- rep(NA, sum(N))
  etah_old <- etal_old <- rep(NA, sum(N))
  h_old <- l_old <- rep(NA, sum(N))
  
  ## Never-takers
  sub <- Z_long==1 & D_long==0
  
  ## eta_h
  logit_etah_old[sub] <- logit_etah[sub]
  logit_etah[sub] <- rnorm(sum(sub), logit_etah_old[sub], metrop_sd_h)
  
  etah_old[sub] <- etah[[s-1]][sub]
  etah[[s]][sub] <- (T_long[sub] + exp(logit_etah[sub]))/(1.00 + exp(logit_etah[sub]))
  
  h_old[sub] <- h[[s-1]][sub]
  h[[s]][sub] <- lamh0[[s]][sub]*etah[[s]][sub]^lamh1[[s]][sub]
  
  W_old <- W
  W[sub,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s-1]], 
             S_long*h[[s]], S_long*l[[s-1]],
             T_long*h[[s]], T_long*l[[s-1]],
             Z_long*h[[s]], Z_long*l[[s-1]])[sub,]
  
  denom <- sum(dnorm(logit_etah_old[sub], q_long[sub,]%*%delta_h[[s-1]], 
                 tau2_h[[s-1]], log=TRUE) +
    dnorm(Y_long[sub], W_old[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  numer <- sum(dnorm(logit_etah[sub], q_long[sub,]%*%delta_h[[s-1]], 
                     tau2_h[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  ratio <- exp(numer - denom)
  uni_draw <- runif(1, min = 0.00, max = 1.00)
  
  if (ratio < uni_draw) {
    logit_etah[sub] <- logit_etah_old[sub]
    etah[[s]][sub] <- etah_old[sub]
    h[[s]][sub] <- h_old[sub]
    W[sub,] <- W_old[sub,]
  }
  
  ##eta_l
  logit_etal_old[sub] <- logit_etal[sub]
  logit_etal[sub] <- rnorm(sum(sub), logit_etal_old[sub], metrop_sd_l)
  
  etal_old[sub] <- etal[[s-1]][sub]
  etal[[s]][sub] <- (etah[[s]][sub] + exp(logit_etal[sub]))/(1.00 + exp(logit_etal[sub]))
  
  l_old[sub] <- l[[s-1]][sub]
  l[[s]][sub] <- laml0[[s]][sub]*etal[[s]][sub]^laml1[[s]][sub]
  
  W_old <- W
  W[sub,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s]], 
             S_long*h[[s]], S_long*l[[s]],
             T_long*h[[s]], T_long*l[[s]],
             Z_long*h[[s]], Z_long*l[[s]])[sub,]
  
  denom <- sum(dnorm(logit_etal_old[sub], q_long[sub,]%*%delta_l[[s-1]], 
                     tau2_l[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W_old[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  numer <- sum(dnorm(logit_etal[sub], q_long[sub,]%*%delta_l[[s-1]], 
                     tau2_l[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  ratio <- exp(numer - denom)
  uni_draw <- runif(1, min = 0.00, max = 1.00)
  
  if (ratio < uni_draw) {
    logit_etal[sub] <- logit_etal_old[sub]
    etal[[s]][sub] <- etal_old[sub]
    l[[s]][sub] <- l_old[sub]
    W[sub,] <- W_old[sub,]
  }
  
  ## Always-takers
  sub <- Z_long==0 & D_long==1

  ##eta_h
  logit_etah_old[sub] <- logit_etah[sub]
  logit_etah[sub] <- rnorm(sum(sub), logit_etah_old[sub], metrop_sd_h)
  
  etah_old[sub] <- etah[[s-1]][sub]
  etah[[s]][sub] <- (T_long[sub]*exp(logit_etah[sub]))/(1.00 + exp(logit_etah[sub]))
  
  h_old[sub] <- h[[s-1]][sub]
  h[[s]][sub] <- lamh0[[s]][sub]*etah[[s]][sub]^lamh1[[s]][sub]
  
  W_old <- W
  W[sub,] <- cbind(1, S_long, T_long, Z_long, h[[s]][sub], l[[s-1]], 
             S_long*h[[s]], S_long*l[[s-1]],
             T_long*h[[s]], T_long*l[[s-1]],
             Z_long*h[[s]], Z_long*l[[s-1]])[sub,]
  
  denom <- sum(dnorm(logit_etah_old[sub], q_long[sub,]%*%delta_h[[s-1]], 
                     tau2_h[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W_old[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  numer <- sum(dnorm(logit_etah[sub], q_long[sub,]%*%delta_h[[s-1]], 
                     tau2_h[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  ratio <- exp(numer - denom)
  uni_draw <- runif(1, min = 0.00, max = 1.00)
  
  if (ratio < uni_draw) {
    logit_etah[sub] <- logit_etah_old[sub]
    etah[[s]][sub] <- etah_old[sub]
    h[[s]][sub] <- h_old[sub]    
    W[sub,] <- W_old[sub,]
  }
  
  ##eta_l
  logit_etal_old[sub] <- logit_etal[sub]
  logit_etal[sub] <- rnorm(sum(sub), logit_etal_old[sub], metrop_sd_l)
  
  etal_old[sub] <- etal[[s-1]][sub]
  etal[[s]][sub] <- (etah[[s]][sub] + T_long[sub]*exp(logit_etal[sub]))/(1.00 + exp(logit_etal[sub]))
  
  l_old[sub] <- l[[s-1]][sub]
  l[[s]][sub] <- laml0[[s]][sub]*etal[[s]][sub]^laml1[[s]][sub]
  
  W_old <- W
  W[sub,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s]], 
             S_long*h[[s]], S_long*l[[s]],
             T_long*h[[s]], T_long*l[[s]],
             Z_long*h[[s]], Z_long*l[[s]])[sub,]
  
  denom <- sum(dnorm(logit_etal_old[sub], q_long[sub,]%*%delta_l[[s-1]], 
                     tau2_l[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W_old[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  numer <- sum(dnorm(logit_etal[sub], q_long[sub,]%*%delta_l[[s-1]], 
                     tau2_l[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  ratio <- exp(numer - denom)
  uni_draw <- runif(1, min = 0.00, max = 1.00)
  
  if (ratio < uni_draw) {
    logit_etal[sub] <- logit_etal_old[sub]
    etal[[s]][sub] <- etal_old[sub]
    l[[s]][sub] <- l_old[sub]
    W[sub,] <- W_old[sub,]
  }
  
  ## Compliers
  sub <- Z_long==D_long
  
  ##eta_h
  logit_etah_old[sub] <- logit_etah[sub]
  logit_etah[sub] <- rnorm(sum(sub), logit_etah_old[sub], metrop_sd_h)
  
  etah_old[sub] <- etah[[s-1]][sub]
  etah[[s]][sub] <- (T_long[sub]*exp(logit_etah[sub]))/(1.00 + exp(logit_etah[sub]))
  
  h_old[sub] <- h[[s-1]][sub]
  h[[s]][sub] <- lamh0[[s]][sub]*etah[[s]][sub]^lamh1[[s]][sub]
  
  W_old <- W
  W[sub,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s-1]], 
             S_long*h[[s]], S_long*l[[s-1]],
             T_long*h[[s]], T_long*l[[s-1]],
             Z_long*h[[s]], Z_long*l[[s-1]])[sub,]
  
  denom <- sum(dnorm(logit_etah_old[sub], q_long[sub,]%*%delta_h[[s-1]], 
                     tau2_h[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W_old[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  numer <- sum(dnorm(logit_etah[sub], q_long[sub,]%*%delta_h[[s-1]], 
                     tau2_h[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  ratio <- exp(numer - denom)
  uni_draw <- runif(1, min = 0.00, max = 1.00)
  
  if (ratio < uni_draw) {
    logit_etah[sub] <- logit_etah_old[sub]
    etah[[s]][sub] <- etah_old[sub]
    h[[s]][sub] <- h_old[sub]
    W[sub,] <- W_old[sub,]
  }
  
  ##eta_l
  logit_etal_old[sub] <- logit_etal[sub]
  logit_etal[sub] <- rnorm(sum(sub), logit_etal_old[sub], metrop_sd_l)
  
  etal_old[sub] <- etal[[s-1]][sub]
  etal[[s]][sub] <- (T_long[sub] + exp(logit_etal[sub]))/(1.00 + exp(logit_etal[sub]))
  
  l_old[sub] <- l[[s-1]][sub]
  l[[s]][sub] <- laml0[[s]][sub]*etal[[s]][sub]^laml1[[s]][sub]
  
  W_old <- W
  W[sub,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s]], 
             S_long*h[[s]], S_long*l[[s]],
             T_long*h[[s]], T_long*l[[s]],
             Z_long*h[[s]], Z_long*l[[s]])[sub,]
  
  denom <- sum(dnorm(logit_etal_old[sub], q_long[sub,]%*%delta_l[[s-1]], 
                     tau2_l[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W_old[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  numer <- sum(dnorm(logit_etal[sub], q_long[sub,]%*%delta_l[[s-1]], 
                     tau2_l[[s-1]], log=TRUE) +
                 dnorm(Y_long[sub], W[sub,]%*%beta[[s]], sigma2[[s]], log=TRUE))
  
  ratio <- exp(numer - denom)
  uni_draw <- runif(1, min = 0.00, max = 1.00)
  
  if (ratio < uni_draw) {
    logit_etal[sub] <- logit_etal_old[sub]
    etal[[s]][sub] <- etal_old[sub]
    l[[s]][sub] <- l_old[sub]
    W[sub,] <- W_old[sub,]
  }

  ########
  #delta_h
  ########
  cov_delta_h <- chol2inv(chol(crossprod(q_long)/tau2_h[[s-1]] + 
                                 diag(ncol(q_long))/sigma2_delta))
  mu_delta_h <- cov_delta_h%*%(crossprod(q_long, logit_etah))/tau2_h[[s-1]]
  delta_h[[s]] <- rmnorm(n = 1, mean = mu_delta_h, varcov = cov_delta_h)
  
  #######
  #tau2_h
  #######
  rate <- crossprod(logit_etah - q_long%*%delta_h[[s]])/2.00 + b_tau2
  tau2_h[[s]] <- 1.00/rgamma(n = 1, shape = shape_tau2_update, rate = 1/rate)
  
  ########
  #delta_l
  ########
  cov_delta_l <- chol2inv(chol(crossprod(q_long)/tau2_l[[s-1]] + 
                                 diag(ncol(q_long))/sigma2_delta))
  mu_delta_l <- cov_delta_l%*%(crossprod(q_long, logit_etal))/tau2_l[[s-1]]
  delta_l[[s]] <- rmnorm(n = 1, mean = mu_delta_l, varcov = cov_delta_l)
  
  #######
  #tau2_l
  #######
  rate <- crossprod(logit_etal - q_long%*%delta_l[[s]])/2.00 + b_tau2
  tau2_l[[s]] <- 1.00/rgamma(n = 1, shape = shape_tau2_update, rate = 1/rate)
  
  ##########
  #Estimands
  ##########
  
  if (s %in% iters) { 
    Y1 <- rep(NA, sum(N))
    Y0 <- rep(NA, sum(N))
    Y0p <- rep(NA, sum(N))
    C <- rep(NA, sum(N))
    
    eff.a <- 0.8
    eff.s <- 0.4
    eff.sp <- 0.8
    for (ij in 1:sum(N)) {
      ## Need to figure out G(eff.a)
      if (h[[s]][ij] <= T_long[ij] & l[[s]][ij] >= T_long[ij]) {
        C[ij] <- 1
      } else {
        C[ij] <- 0
      }
      if (C[ij]==1) {
        W0p <- W0 <- W1 <- W[ij,]
        
        W0[2] <- W1[2] <- eff.s
        W0[7] <- W1[7] <- eff.s*h[[s]][ij]
        W0[8] <- W1[8] <- eff.s*l[[s]][ij]
        
        W0p[2] <- eff.sp
        W0p[7] <- eff.sp*h[[s]][ij]
        W0p[8] <- eff.sp*l[[s]][ij]
        
        W0[3] <- W0p[3] <- W1[3] <- eff.a
        W0[9] <- W0p[9] <- W1[9] <- eff.a*h[[s]][ij]
        W0[10] <- W0p[10] <- W1[10] <- eff.a*l[[s]][ij]
        
        W0[4] <- W0[11] <- W0p[12] <- W0p[4] <- W0p[11] <- W0p[12] <- 0
        W1[4] <- 1
        W1[11] <- h[[s]][ij]
        W1[12] <- l[[s]][ij]
        
        mu0<-W0%*%beta[[s]]
        var0<-sigma2[[s]]
        Y0[ij]<-rnorm(n = 1,
                      mean = mu0,
                      sd = sqrt(var0)) 
        
        mu1<-W1%*%beta[[s]]
        var1<-sigma2[[s]]
        Y1[ij]<-rnorm(n = 1,
                      mean = mu1,
                      sd = sqrt(var1))
        
        mu0p<-W0p%*%beta[[s]]
        var0p<-sigma2[[s]]
        Y0p[ij]<-rnorm(n = 1,
                       mean = mu0p,
                       sd = sqrt(var0p)) 
      }
    }
    CADE[which(iters==s)] <- sum(C*(Y1-Y0), na.rm=TRUE)/sum(C)
    CASE[which(iters==s)] <- sum(C*(Y0-Y0p), na.rm=TRUE)/sum(C)
    pC[which(iters==s)] <- mean(C)
  }
}


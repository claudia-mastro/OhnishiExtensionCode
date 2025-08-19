###############
#Packages
###############
library(mnormt)
library(matrixStats)
library(data.table)
library(invgamma)
args <- commandArgs(trailingOnly = TRUE)
id <- as.integer(args[1])
print(id)
J <- as.integer(args[2])
print(J)
Nj <- as.integer(args[3])
print(Nj)
nalpha <- as.integer(args[4])
print(nalpha)
v <- paste0("LM3_8.16")
print(v)
source("~/OhnishiExtensionCode/effects_LM3.R")
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
mcmc_samples<-100000
burnin <- 50000
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

delta_h[[1]] <- delta_l[[1]] <- rep(-1.00, ncol(q_long))
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
metrop_sd_h <- rep(0.1, sum(N))
acctot_h <- matrix(1, nrow=mcmc_samples, ncol=sum(N))

metrop_sd_l <- rep(0.1, sum(N))
acctot_l <- matrix(1, nrow=mcmc_samples, ncol=sum(N))

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
  #beta[[s]] <- beta_true
  
  rate <- crossprod(Y_long - W%*%beta[[s]])/2.00 + b_sigma2
  shape <- length(Y_long)/2.00 + a_sigma2
  sigma2[[s]] <- 1/rgamma(n = 1, shape = shape, rate = rate)
  # sigma2[[s]] <- sigma2_true
  
  ###
  #pi
  ###
  ph0[[s]] <- rbeta(1, a_pi + sum(lamh0[[s-1]]), b_pi + sum(N) - sum(lamh0[[s-1]]))
  ph1[[s]] <- rbeta(1, a_pi + sum(lamh1[[s-1]]), b_pi + sum(N) - sum(lamh1[[s-1]]))
  pl0[[s]] <- rbeta(1, a_pi + sum(laml0[[s-1]]), b_pi + sum(N) - sum(laml0[[s-1]]))
  pl1[[s]] <- rbeta(1, a_pi + sum(laml1[[s-1]]), b_pi + sum(N) - sum(laml1[[s-1]]))
  # ph0[[s]] <- ph0_true
  # ph1[[s]] <- ph1_true
  # pl0[[s]] <- pl0_true
  # pl1[[s]] <- pl1_true
  
  #######
  #lambda
  #######
  
  ph0p <- ph0[[s]]^(1-(Z_long==1 & D_long==0))
  W0_temp <- cbind(1, S_long, T_long, Z_long,
                   0, l[[s-1]],
                   0, S_long*l[[s-1]],
                   0, T_long*l[[s-1]],
                   0, Z_long*l[[s-1]])
  W1_temp <- cbind(1, S_long, T_long, Z_long,
                   etah[[s-1]]^lamh1[[s-1]], l[[s-1]],
                   S_long*etah[[s-1]]^lamh1[[s-1]], S_long*l[[s-1]],
                   T_long*etah[[s-1]]^lamh1[[s-1]], T_long*l[[s-1]],
                   Z_long*etah[[s-1]]^lamh1[[s-1]], Z_long*l[[s-1]])
  d1 <- dnorm(Y_long, W1_temp%*%beta[[s]], sqrt(sigma2[[s]]))+1e-16
  d0 <- dnorm(Y_long, W0_temp%*%beta[[s]], sqrt(sigma2[[s]]))+1e-16

  lamh0[[s]] <- sapply((d1*ph0p)/(d1*ph0p + d0*(1-ph0p)), function(p) rbinom(1, 1, p))

  ph1p <- ph1[[s]]^(1-(Z_long==0 & D_long==1))
  W0_temp <- cbind(1, S_long, T_long, Z_long,
                   lamh0[[s]], l[[s-1]],
                   S_long*lamh0[[s]], S_long*l[[s-1]],
                   T_long*lamh0[[s]], T_long*l[[s-1]],
                   Z_long*lamh0[[s]], Z_long*l[[s-1]])
  W1_temp <- cbind(1, S_long, T_long, Z_long,
                   lamh0[[s]]*etah[[s-1]], l[[s-1]],
                   S_long*lamh0[[s]]*etah[[s-1]], S_long*l[[s-1]],
                   T_long*lamh0[[s]]*etah[[s-1]], T_long*l[[s-1]],
                   Z_long*lamh0[[s]]*etah[[s-1]], Z_long*l[[s-1]])
  d1 <- dnorm(Y_long, W1_temp%*%beta[[s]], sqrt(sigma2[[s]]))+1e-16
  d0 <- dnorm(Y_long, W0_temp%*%beta[[s]], sqrt(sigma2[[s]]))+1e-16
  lamh1[[s]] <- sapply((d1*ph1p)/(d1*ph1p + d0*(1-ph1p)), function(p) rbinom(1, 1, p))

  pl0p <- (pl0[[s]]^(1-(Z_long==1 & D_long==0)))^(1-(lamh0[[s]]==1 & lamh1[[s]]==0))
  W0_temp <- cbind(1, S_long, T_long, Z_long,
                   lamh0[[s]]*etah[[s-1]]^lamh1[[s]], 0,
                   S_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], 0,
                   T_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], 0,
                   Z_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], 0)
  W1_temp <- cbind(1, S_long, T_long, Z_long,
                   lamh0[[s]]*etah[[s-1]]^lamh1[[s]], l[[s-1]],
                   S_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], S_long*etal[[s-1]]^laml1[[s-1]],
                   T_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], T_long*etal[[s-1]]^laml1[[s-1]],
                   Z_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], Z_long*etal[[s-1]]^laml1[[s-1]])
  d1 <- dnorm(Y_long, W1_temp%*%beta[[s]], sqrt(sigma2[[s]]))+1e-16
  d0 <- dnorm(Y_long, W0_temp%*%beta[[s]], sqrt(sigma2[[s]]))+1e-16
  laml0[[s]] <- sapply((d1*pl0p)/(d1*pl0p + d0*(1-pl0p)), function(p) rbinom(1, 1, p))

  pl1p <- (1-(lamh0[[s]]==1 & lamh1[[s]]==0))*(pl1[[s]]^(1-(Z_long==0 & D_long==1)))
  W0_temp <- cbind(1, S_long, T_long, Z_long,
                   lamh0[[s]]*etah[[s-1]]^lamh1[[s]], laml0[[s]],
                   S_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], laml0[[s]],
                   T_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], laml0[[s]],
                   Z_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], laml0[[s]])
  W1_temp <- cbind(1, S_long, T_long, Z_long,
                   lamh0[[s]]*etah[[s-1]]^lamh1[[s]], laml0[[s]]*etal[[s-1]],
                   S_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], S_long*laml0[[s]]*etal[[s-1]],
                   T_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], T_long*laml0[[s]]*etal[[s-1]],
                   Z_long*lamh0[[s]]*etah[[s-1]]^lamh1[[s]], Z_long*laml0[[s]]*etal[[s-1]])
  d1 <- dnorm(Y_long, W1_temp%*%beta[[s]], sqrt(sigma2[[s]]))+1e-16
  d0 <- dnorm(Y_long, W0_temp%*%beta[[s]], sqrt(sigma2[[s]]))+1e-16
  laml1[[s]] <- sapply(abs((d1*pl1p)/(d1*pl1p + d0*(1-pl1p))), function(p) rbinom(1, 1, p))

  # lamh0[[s]] <- lamh0_true
  # lamh1[[s]] <- lamh1_true
  # laml0[[s]] <- laml0_true
  # laml1[[s]] <- laml1_true
  
   
  ####
  #eta
  ####
  block_size <- 1
  logit_etah_old <- logit_etal_old <- rep(NA, sum(N))
  etah_old <- etal_old <- rep(NA, sum(N))
  h_old <- l_old <- rep(NA, sum(N))

  ## Never-takers
  sub <- which(Z_long == 1 & D_long == 0)
  if (length(sub) > 0) {
    sub <- sample(sub)
    for (start_idx in seq(1, length(sub), by = block_size)) {
      block <- sub[start_idx:min(start_idx + block_size - 1, length(sub))]

      ## eta_h
      logit_etah_old[block] <- logit_etah[block]
      logit_etah[block] <- rnorm(length(block), logit_etah_old[block], metrop_sd_h[block])
      logit_etah[logit_etah > 709.7827] <- 709.7827

      etah_old[block] <- etah[[s-1]][block]
      etah[[s]][block] <- (T_long[block] + exp(logit_etah[block]))/(1.00 + exp(logit_etah[block]))

      h_old[block] <- h[[s-1]][block]
      h[[s]][block] <- lamh0[[s]][block]*etah[[s]][block]^lamh1[[s]][block]

      W_old <- W
      W[block,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s-1]],
                         S_long*h[[s]], S_long*l[[s-1]],
                         T_long*h[[s]], T_long*l[[s-1]],
                         Z_long*h[[s]], Z_long*l[[s-1]])[block,]

      denom <- sum(dnorm(logit_etah_old[block], q_long[block,]%*%delta_h[[s-1]],
                         sqrt(tau2_h[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W_old[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      numer <- sum(dnorm(logit_etah[block], q_long[block,]%*%delta_h[[s-1]],
                         sqrt(tau2_h[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      ratio <- numer - denom
      uni_draw <- log(runif(1, min = 0.00, max = 1.00))

      if (ratio < uni_draw) {
        logit_etah[block] <- logit_etah_old[block]
        etah[[s]][block] <- etah_old[block]
        h[[s]][block] <- h_old[block]
        W[block,] <- W_old[block,]
        acctot_h[s, block] <- 0
      }
    }

    for (start_idx in seq(1, length(sub), by = block_size)) {
      block <- sub[start_idx:min(start_idx + block_size - 1, length(sub))]

      ##eta_l
      logit_etal_old[block] <- logit_etal[block]
      logit_etal[block] <- rnorm(length(block), logit_etal_old[block], metrop_sd_l)
      logit_etal[logit_etal > 709.7827] <- 709.7827

      etal_old[block] <- etal[[s-1]][block]
      etal[[s]][block] <- (etah[[s]][block] + exp(logit_etal[block]))/(1.00 + exp(logit_etal[block]))

      l_old[block] <- l[[s-1]][block]
      l[[s]][block] <- laml0[[s]][block]*etal[[s]][block]^laml1[[s]][block]

      W_old <- W
      W[block,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s]],
                       S_long*h[[s]], S_long*l[[s]],
                       T_long*h[[s]], T_long*l[[s]],
                       Z_long*h[[s]], Z_long*l[[s]])[block,]

      denom <- sum(dnorm(logit_etal_old[block], q_long[block,]%*%delta_l[[s-1]],
                         sqrt(tau2_l[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W_old[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      numer <- sum(dnorm(logit_etal[block], q_long[block,]%*%delta_l[[s-1]],
                         sqrt(tau2_l[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      ratio <- numer - denom
      uni_draw <- log(runif(1, min = 0.00, max = 1.00))

      if (ratio < uni_draw) {
        logit_etal[block] <- logit_etal_old[block]
        etal[[s]][block] <- etal_old[block]
        l[[s]][block] <- l_old[block]
        W[block,] <- W_old[block,]
        acctot_l[s, block] <- 0
      }
    }
  }

  ## Always-takers
  sub <- which(Z_long==0 & D_long==1)
  if (length(sub) > 0) {
    sub <- sample(sub)
    for (start_idx in seq(1, length(sub), by = block_size)) {
      block <- sub[start_idx:min(start_idx + block_size - 1, length(sub))]

      ##eta_h
      logit_etah_old[block] <- logit_etah[block]
      logit_etah[block] <- rnorm(length(block), logit_etah_old[block], metrop_sd_h[block])
      logit_etah[logit_etah > 709.7827] <- 709.7827

      etah_old[block] <- etah[[s-1]][block]
      etah[[s]][block] <- (T_long[block]*exp(logit_etah[block]))/(1.00 + exp(logit_etah[block]))

      h_old[block] <- h[[s-1]][block]
      h[[s]][block] <- lamh0[[s]][block]*etah[[s]][block]^lamh1[[s]][block]

      W_old <- W
      W[block,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s-1]],
                       S_long*h[[s]], S_long*l[[s-1]],
                       T_long*h[[s]], T_long*l[[s-1]],
                       Z_long*h[[s]], Z_long*l[[s-1]])[block,]

      denom <- sum(dnorm(logit_etah_old[block], q_long[block,]%*%delta_h[[s-1]],
                         sqrt(tau2_h[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W_old[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      numer <- sum(dnorm(logit_etah[block], q_long[block,]%*%delta_h[[s-1]],
                         sqrt(tau2_h[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      ratio <- numer - denom
      uni_draw <- log(runif(1, min = 0.00, max = 1.00))

      if (ratio < uni_draw) {
        logit_etah[block] <- logit_etah_old[block]
        etah[[s]][block] <- etah_old[block]
        h[[s]][block] <- h_old[block]
        W[block,] <- W_old[block,]
        acctot_h[s, block] <- 0
      }
    }

    for (start_idx in seq(1, length(sub), by = block_size)) {
      block <- sub[start_idx:min(start_idx + block_size - 1, length(sub))]

      ##eta_l
      logit_etal_old[block] <- logit_etal[block]
      logit_etal[block] <- rnorm(length(block), logit_etal_old[block], metrop_sd_l)
      logit_etal[logit_etal > 709.7827] <- 709.7827

      etal_old[block] <- etal[[s-1]][block]
      etal[[s]][block] <- (etah[[s]][block] + T_long[block]*exp(logit_etal[block]))/(1.00 + exp(logit_etal[block]))

      l_old[block] <- l[[s-1]][block]
      l[[s]][block] <- laml0[[s]][block]*etal[[s]][block]^laml1[[s]][block]

      W_old <- W
      W[block,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s]],
                       S_long*h[[s]], S_long*l[[s]],
                       T_long*h[[s]], T_long*l[[s]],
                       Z_long*h[[s]], Z_long*l[[s]])[block,]

      denom <- sum(dnorm(logit_etal_old[block], q_long[block,]%*%delta_l[[s-1]],
                         sqrt(tau2_l[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W_old[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      numer <- sum(dnorm(logit_etal[block], q_long[block,]%*%delta_l[[s-1]],
                         sqrt(tau2_l[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      ratio <- numer - denom
      uni_draw <- log(runif(1, min = 0.00, max = 1.00))

      if (ratio < uni_draw) {
        logit_etal[block] <- logit_etal_old[block]
        etal[[s]][block] <- etal_old[block]
        l[[s]][block] <- l_old[block]
        W[block,] <- W_old[block,]
        acctot_l[s, block] <- 0
      }
    }
  }


  ## Complier/NT
  sub <- which(Z_long==0 & D_long==0)

  if (length(sub) > 0) {
    sub <- sample(sub)
    for (start_idx in seq(1, length(sub), by = block_size)) {
      block <- sub[start_idx:min(start_idx + block_size - 1, length(sub))]
      ##eta_h
      logit_etah_old[block] <- logit_etah[block]
      logit_etah[block] <- rnorm(length(block), logit_etah_old[block], metrop_sd_h[block])
      logit_etah[logit_etah > 709.7827] <- 709.7827

      etah_old[block] <- etah[[s-1]][block]
      etah[[s]][block] <- (exp(logit_etah[block]))/(1.00 + exp(logit_etah[block]))


      h_old[block] <- h[[s-1]][block]
      h[[s]][block] <- lamh0[[s]][block]*etah[[s]][block]^lamh1[[s]][block]

      W_old <- W
      W[block,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s-1]],
                       S_long*h[[s]], S_long*l[[s-1]],
                       T_long*h[[s]], T_long*l[[s-1]],
                       Z_long*h[[s]], Z_long*l[[s-1]])[block,]

      denom <- sum(dnorm(logit_etah_old[block], q_long[block,]%*%delta_h[[s-1]],
                         sqrt(tau2_h[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W_old[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      numer <- sum(dnorm(logit_etah[block], q_long[block,]%*%delta_h[[s-1]],
                         sqrt(tau2_h[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      ratio <- numer - denom
      uni_draw <- log(runif(1, min = 0.00, max = 1.00))

      if (ratio < uni_draw) {
        logit_etah[block] <- logit_etah_old[block]
        etah[[s]][block] <- etah_old[block]
        h[[s]][block] <- h_old[block]
        W[block,] <- W_old[block,]
        acctot_h[s, block] <- 0
      }
    }

    for (start_idx in seq(1, length(sub), by = block_size)) {
      block <- sub[start_idx:min(start_idx + block_size - 1, length(sub))]

      ##eta_l
      logit_etal_old[block] <- logit_etal[block]
      logit_etal[block] <- rnorm(length(block), logit_etal_old[block], metrop_sd_l)
      logit_etal[logit_etal > 709.7827] <- 709.7827

      etal_old[block] <- etal[[s-1]][block]
      etal[[s]][block] <- (max(T_long[block], etah[[s]][block]) + exp(logit_etal[block]))/(1.00 + exp(logit_etal[block]))

      l_old[block] <- l[[s-1]][block]
      l[[s]][block] <- laml0[[s]][block]*etal[[s]][block]^laml1[[s]][block]

      W_old <- W
      W[block,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s]],
                       S_long*h[[s]], S_long*l[[s]],
                       T_long*h[[s]], T_long*l[[s]],
                       Z_long*h[[s]], Z_long*l[[s]])[block,]

      denom <- sum(dnorm(logit_etal_old[block], q_long[block,]%*%delta_l[[s-1]],
                         sqrt(tau2_l[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W_old[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      numer <- sum(dnorm(logit_etal[block], q_long[block,]%*%delta_l[[s-1]],
                         sqrt(tau2_l[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))

      ratio <- numer - denom
      uni_draw <- log(runif(1, min = 0.00, max = 1.00))

      if (ratio < uni_draw) {
        logit_etal[block] <- logit_etal_old[block]
        etal[[s]][block] <- etal_old[block]
        l[[s]][block] <- l_old[block]
        W[block,] <- W_old[block,]
        acctot_l[s, block] <- 0
      }
    }
  }
  
  ## Complier/AT
  sub <- which(Z_long==1 & D_long==1)
  
  if (length(sub) > 0) {
    sub <- sample(sub)
    for (start_idx in seq(1, length(sub), by = block_size)) {
      block <- sub[start_idx:min(start_idx + block_size - 1, length(sub))]
      ##eta_h
      logit_etah_old[block] <- logit_etah[block]
      logit_etah[block] <- rnorm(length(block), logit_etah_old[block], metrop_sd_h[block])
      logit_etah[logit_etah > 709.7827] <- 709.7827
      
      etah_old[block] <- etah[[s-1]][block]
      etah[[s]][block] <- (T_long[block]*exp(logit_etah[block]))/(1.00 + exp(logit_etah[block]))
      
      
      h_old[block] <- h[[s-1]][block]
      h[[s]][block] <- lamh0[[s]][block]*etah[[s]][block]^lamh1[[s]][block]
      
      W_old <- W
      W[block,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s-1]],
                         S_long*h[[s]], S_long*l[[s-1]],
                         T_long*h[[s]], T_long*l[[s-1]],
                         Z_long*h[[s]], Z_long*l[[s-1]])[block,]
      
      denom <- sum(dnorm(logit_etah_old[block], q_long[block,]%*%delta_h[[s-1]],
                         sqrt(tau2_h[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W_old[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))
      
      numer <- sum(dnorm(logit_etah[block], q_long[block,]%*%delta_h[[s-1]],
                         sqrt(tau2_h[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))
      
      ratio <- numer - denom
      uni_draw <- log(runif(1, min = 0.00, max = 1.00))
      
      if (ratio < uni_draw) {
        logit_etah[block] <- logit_etah_old[block]
        etah[[s]][block] <- etah_old[block]
        h[[s]][block] <- h_old[block]
        W[block,] <- W_old[block,]
        acctot_h[s, block] <- 0
      }
    }
    
    for (start_idx in seq(1, length(sub), by = block_size)) {
      block <- sub[start_idx:min(start_idx + block_size - 1, length(sub))]
      
      ##eta_l
      logit_etal_old[block] <- logit_etal[block]
      logit_etal[block] <- rnorm(length(block), logit_etal_old[block], metrop_sd_l)
      logit_etal[logit_etal > 709.7827] <- 709.7827
      
      etal_old[block] <- etal[[s-1]][block]
      etal[[s]][block] <- (etah[[s]][block] + exp(logit_etal[block]))/(1.00 + exp(logit_etal[block]))
      
      l_old[block] <- l[[s-1]][block]
      l[[s]][block] <- laml0[[s]][block]*etal[[s]][block]^laml1[[s]][block]
      
      W_old <- W
      W[block,] <- cbind(1, S_long, T_long, Z_long, h[[s]], l[[s]],
                         S_long*h[[s]], S_long*l[[s]],
                         T_long*h[[s]], T_long*l[[s]],
                         Z_long*h[[s]], Z_long*l[[s]])[block,]
      
      denom <- sum(dnorm(logit_etal_old[block], q_long[block,]%*%delta_l[[s-1]],
                         sqrt(tau2_l[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W_old[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))
      
      numer <- sum(dnorm(logit_etal[block], q_long[block,]%*%delta_l[[s-1]],
                         sqrt(tau2_l[[s-1]]), log=TRUE) +
                     dnorm(Y_long[block], W[block,]%*%beta[[s]], sqrt(sigma2[[s]]), log=TRUE))
      
      ratio <- numer - denom
      uni_draw <- log(runif(1, min = 0.00, max = 1.00))
      
      if (ratio < uni_draw) {
        logit_etal[block] <- logit_etal_old[block]
        etal[[s]][block] <- etal_old[block]
        l[[s]][block] <- l_old[block]
        W[block,] <- W_old[block,]
        acctot_l[s, block] <- 0
      }
    }
  }
  # logit_etah <- logit_h_true
  # logit_etal <- logit_l_true
  # etah[[s]] <- etah_true
  # etal[[s]] <- etal_true
  # h[[s]] <- h_true
  # l[[s]] <- l_true
  # W <- cbind(1, S_long, T_long, Z_long, h[[1]], l[[1]], 
  #            S_long*h[[1]], S_long*l[[1]],
  #            T_long*h[[1]], T_long*l[[1]],
  #            Z_long*h[[1]], Z_long*l[[1]])

  ########
  #delta_h
  ########
  cov_delta_h <- chol2inv(chol(crossprod(q_long)/tau2_h[[s-1]] +
                                 diag(ncol(q_long))/sigma2_delta))
  mu_delta_h <- cov_delta_h%*%(crossprod(q_long, logit_etah))/tau2_h[[s-1]]
  delta_h[[s]] <- rmnorm(n = 1, mean = mu_delta_h, varcov = cov_delta_h)
  # delta_h[[s]] <- delta_h_true
  
  #######
  #tau2_h
  #######
  rate <- crossprod(logit_etah - q_long%*%delta_h[[s]])/2.00 + b_tau2
  tau2_h[[s]] <- 1/rgamma(n = 1, shape = shape_tau2_update, rate = rate)
  # tau2_h[[s]] <- tau2_h_true
  
  ########
  #delta_l
  ########
  cov_delta_l <- chol2inv(chol(crossprod(q_long)/tau2_l[[s-1]] +
                                 diag(ncol(q_long))/sigma2_delta))
  mu_delta_l <- cov_delta_l%*%(crossprod(q_long, logit_etal))/tau2_l[[s-1]]
  delta_l[[s]] <- rmnorm(n = 1, mean = mu_delta_l, varcov = cov_delta_l)
  # delta_l[[s]] <- delta_l_true
  
  #######
  #tau2_l
  #######
  rate <- crossprod(logit_etal - q_long%*%delta_l[[s]])/2.00 + b_tau2
  tau2_l[[s]] <- 1/rgamma(n = 1, shape = shape_tau2_update, rate = rate)
  # tau2_l[[s]] <- tau2_l_true
  
  ##########
  #Estimands
  ##########
  
  if (s %in% iters) { 
    eff.a <- 0.8
    eff.s <- 0.4
    eff.sp <- 0.8
    
    effects <- CADE.CASE(eff.a, eff.s, eff.sp, h[[s]], l[[s]],
                         beta[[s]], sigma2[[s]], W)
    CADE[which(iters==s)] <- effects[[1]]
    CASE[which(iters==s)] <- effects[[2]]
  }
  if (s <= burnin & s%%100 ==0) {
    acch <- colMeans(acctot_h[(s-100):s,])
    accl <- colMeans(acctot_l[(s-100):s,])
    
    metrop_sd_h[acch < 0.15] <- metrop_sd_h[acch < 0.15]*0.9
    metrop_sd_h[acch > 0.60] <- metrop_sd_h[acch > 0.60]*1.1

    metrop_sd_l[accl < 0.15] <- metrop_sd_l[accl < 0.15]*0.9
    metrop_sd_l[accl > 0.60] <- metrop_sd_l[accl > 0.60]*1.1
  }
}
print("saving")
print(paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/CADE",
             id, ".rds"))
saveRDS(CADE, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/CADE",
                     id, ".rds"))
saveRDS(CASE, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/CASE",
                     id, ".rds"))
saveRDS(beta, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/beta",
                     id, ".rds"))
saveRDS(sigma2, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/sigma2",
                       id, ".rds"))
saveRDS(ph0, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v, "/ph0",
                  id, ".rds"))
saveRDS(ph1, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v, "/ph1",
                    id, ".rds"))
saveRDS(pl0, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v, "/pl0",
                    id, ".rds"))
saveRDS(pl1, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v, "/pl1",
                    id, ".rds"))
saveRDS(lamh0, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v, "/lamh0",
                    id, ".rds"))
saveRDS(lamh1, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v, "/lamh1",
                    id, ".rds"))
saveRDS(laml0, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v, "/laml0",
                    id, ".rds"))
saveRDS(laml1, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v, "/laml1",
                    id, ".rds"))
saveRDS(etah, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/etah",
                   id, ".rds"))
saveRDS(etal, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/etal",
                   id, ".rds"))
saveRDS(h, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/h",
                   id, ".rds"))
saveRDS(l, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/l",
                   id, ".rds"))
saveRDS(delta_h, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/deltah",
                         id, ".rds"))
saveRDS(delta_l, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/deltal",
                         id, ".rds"))
saveRDS(tau2_h, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/tau2h",
                        id, ".rds"))
saveRDS(tau2_l, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/tau2l",
                        id, ".rds"))
saveRDS(acctot_h, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/acctot_h",
                          id, ".rds"))
saveRDS(acctot_l, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/acctot_l",
                          id, ".rds"))
saveRDS(pC, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/propcomp",
                   id, ".rds"))


###############
#Packages
###############
rm(list=ls())
library(mnormt)
library(matrixStats)
source("~/project/OhnishiExtension/JWCode/Data_Simulation(2).R")

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
set.seed(5489)

################
#Global Settings
################
mcmc_samples<-100000

p_s<-ncol(v_long)

sigma2_beta<-100.00^2
sigma2_alpha<-1.00^2
sigma2_delta<-1.00^2

###########
#Parameters
###########
beta<-list(0)
for(k in 1:3){
   beta[[k]]<-matrix(NA,
                     nrow = mcmc_samples,
                     ncol = 5)
   }

sigma2<-matrix(NA,
               nrow = mcmc_samples,
               ncol = 3)

G<-list(0)
for(j in 1:J){
   G[[j]]<-matrix(NA,
                  nrow = mcmc_samples,
                  ncol = N[j])
   }

alpha<-list(0)
for(k in 1:5){
   alpha[[k]]<-matrix(NA,
                      nrow = mcmc_samples,
                      ncol = 2)
}

logith0 <- matrix(NA, nrow = mcmc_samples, ncol = sum(N))
logitl0 <- matrix(NA, nrow = mcmc_samples, ncol = sum(N))
logith1 <- matrix(NA, nrow = mcmc_samples, ncol = sum(N))
logitl1 <- matrix(NA, nrow = mcmc_samples, ncol = sum(N))

delta_h0<-matrix(NA,
                 nrow = mcmc_samples,
                 ncol = p_s)

delta_l0<-matrix(NA,
                 nrow = mcmc_samples,
                 ncol = p_s)

delta_h1<-matrix(NA,
                 nrow = mcmc_samples,
                 ncol = p_s)

delta_l1<-matrix(NA,
                 nrow = mcmc_samples,
                 ncol = p_s)

tau2_h0<-rep(NA, mcmc_samples)

tau2_l0<-rep(NA, mcmc_samples)

tau2_h1<-rep(NA, mcmc_samples)

tau2_l1<-rep(NA, mcmc_samples)


CADE<-rep(NA, mcmc_samples)

###############
#Initial Values
###############
for(k in 1:3){
   beta[[k]][1,]<-0.00
   }

sigma2[1,]<-1.00

log_pi_mat_temp<-matrix(NA,
                        nrow = sum(N),
                        ncol = 6)
log_pi_mat_temp[,6]<-0.00
for(k in 1:5){

   alpha[[k]][1,]<-0.00
   log_pi_mat_temp[,k]<-q_long%*%alpha[[k]][1,]

   }
log_pi_mat<-
pi_mat<-matrix(NA,
               nrow = sum(N),
               ncol = 6)
for(k in 1:6){

   log_pi_mat[,k]<-log_pi_mat_temp[,k] -
                   apply(log_pi_mat_temp, 1, logSumExp) 
   pi_mat[,k]<-1.00/rowSums(exp(log_pi_mat_temp - log_pi_mat_temp[,k]))

   }

for(j in 1:J){
   for(i in 1:N[j]){

      if(Z[[j]][i] == 1 & D[[j]][i] == 1){
        G[[j]][1,i]<-sample(c(1, 3, 4, 5, 6),
                            size = 1,
                            replace = TRUE)
        }

      if(Z[[j]][i] == 1 & D[[j]][i] == 0){
        G[[j]][1,i]<-sample(c(2, 4, 6),
                            size = 1,
                            replace = TRUE)
        }

      if(Z[[j]][i] == 0 & D[[j]][i] == 1){
        G[[j]][1,i]<-sample(c(1, 5, 6),
                            size = 1,
                            replace = TRUE)
        }

      if(Z[[j]][i] == 0 & D[[j]][i] == 0){
        G[[j]][1,i]<-sample(c(2, 3, 5, 6),
                            size = 1,
                            replace = TRUE)
        }

      }
   }
G_long<-G[[1]][1,]
for(j in 2:J){
   G_long<-c(G_long,
             G[[j]][1,])
   }

delta_h0[1,]<-0.00
delta_l0[1,]<-0.00
delta_h1[1,]<-0.00
delta_l1[1,]<-0.00

tau2_h0[1] <- 1
tau2_l0[1] <- 1
tau2_h1[1] <- 1
tau2_l1[1] <- 1

logith0[1,] <- rnorm(sum(N), v_long%*%delta_h0[1,], sqrt(tau2_h0[1]))
logitl0[1,] <- rnorm(sum(N), v_long%*%delta_l0[1,], sqrt(tau2_l0[1]))
logith1[1,] <- rnorm(sum(N), v_long%*%delta_h1[1,], sqrt(tau2_h1[1]))
logitl1[1,] <- rnorm(sum(N), v_long%*%delta_l1[1,], sqrt(tau2_l1[1]))

h0<-1.00/(1.00 + exp(-logith0[1,]))
l0<-1.00/(1.00 + exp(-logitl0[1,]))
h1<-1.00/(1.00 + exp(-logith1[1,]))
l1<-(h1 + 1.00)/(1.00 + exp(-logitl1[1,]))


G_a_long<-rep(NA,
              times = sum(N))
G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0)) | ((G_long == 6) & (a_long >= l1))]<-1
G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < h0)) | ((G_long == 6) & (a_long < h1))]<-2
G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= h0)) | ((G_long == 5) & (a_long < l0)) | ((G_long == 6) & (a_long >= h1) & (a_long < l1))]<-3

####################
#Metropolis Settings
####################
metrop_sd_alpha<-matrix(0.50,
                        nrow = 5,
                        ncol = 2)
acctot_alpha<-matrix(1,
                     nrow = 5,
                     ncol = 2)

metrop_sd_h0<-0.5
acctot_h0 <- 1

metrop_sd_l0<-0.5
acctot_l0 <- 1

metrop_sd_h1<-0.5
acctot_h1 <- 1

metrop_sd_l1<-0.5
acctot_l1 <- 1

###################
#Main Sampling Loop
###################
for(s in 2:mcmc_samples){

   #############
   #beta, sigma2
   #############
   for(k in 1:3){

      Y_k<-Y_long[G_a_long == k]
      W_k<-W[(G_a_long == k),]
      cov_beta<-chol2inv(chol(crossprod(W_k)/sigma2[(s-1), k] + diag(5)/sigma2_beta))
      mu_beta<-cov_beta%*%(crossprod(W_k, Y_k))/sigma2[(s-1), k]
      beta[[k]][s,]<-rmnorm(n = 1,
                            mean = mu_beta,
                            varcov = cov_beta)

      rate<-crossprod(Y_k - W_k%*%beta[[k]][s,])/2.00 +
            0.01
      shape<-length(Y_k)/2.00 +
             0.01
      sigma2[s,k]<-1.00/rgamma(n = 1,
                               shape = shape,
                               rate = rate)

      }
 
   #######
   #G, G_a
   #######
   num<-matrix(NA,
               nrow = sum(N),
               ncol = 6)
   for(k in 1:6){
      
      G_long_temp<-rep(k,
                       times = sum(N))
      G_a_long_temp<-rep(NA,
                         times = sum(N))
      G_a_long_temp[(G_long_temp == 1) | ((G_long_temp == 5) & (a_long >= l0)) | ((G_long_temp == 6) & (a_long >= l1))]<-1
      G_a_long_temp[(G_long_temp == 2) | ((G_long_temp == 4) & (a_long < h0)) | ((G_long_temp == 6) & (a_long < h1))]<-2
      G_a_long_temp[(G_long_temp == 3) | ((G_long_temp == 4) & (a_long >= h0)) | ((G_long_temp == 5) & (a_long < l0)) | ((G_long_temp == 6) & (a_long >= h1) & (a_long < l1))]<-3

      mu<-rep(NA,
              times = sum(N))
      var<-rep(NA, 
               times = sum(N))
      for(j in 1:sum(N)){

         mu[j]<-W[j,]%*%beta[[G_a_long_temp[j]]][s,]
         var[j]<-sigma2[s, G_a_long_temp[j]]
         
         }

      num[,k]<-dnorm(x = Y_long,
                     mean = mu,
                     sd = sqrt(var),
                     log = TRUE)
      if(k < 6){
        num[,k]<-num[,k] +
                 q_long%*%alpha[[k]][(s-1),]
        }

      num[((D_long == 0) & (G_a_long_temp == 1)), k]<--Inf
      num[((D_long == 1) & (G_a_long_temp == 2)), k]<--Inf
      num[((D_long != Z_long) & (G_a_long_temp == 3)), k]<--Inf
     
      }

   probs<-matrix(NA,
                 nrow = sum(N),
                 ncol = 6)
   for(k in 1:6){
      probs[,k]<-1.00/rowSums(exp(num - num[,k]))
      }
   probs[is.na(probs) == 1]<-0.00
          
   G_long<-sapply(c(1:sum(N)), 
                  function(idx){
                          sample(c(1:6), 
                                 size = 1, 
                                 prob = probs[idx,],
                                 replace = TRUE)
                          })

   G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0)) | ((G_long == 6) & (a_long >= l1))]<-1
   G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < h0)) | ((G_long == 6) & (a_long < h1))]<-2
   G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= h0)) | ((G_long == 5) & (a_long < l0)) | ((G_long == 6) & (a_long >= h1) & (a_long < l1))]<-3

   G[[1]][s,]<-G_long[1:N[1]]
   for(j in 2:J){
      G[[j]][s,]<-G_long[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
      }

   ######
   #alpha
   ######
   for(k in 1:5){

      alpha[[k]][s,]<-alpha[[k]][(s-1),]

      for(l in 1:2){

         log_pi_mat_temp_old<-log_pi_mat_temp
         log_pi_mat_old<-log_pi_mat
         pi_mat_old<-pi_mat
         denom<-sum(log_pi_mat_old[cbind(seq_along(G_long), G_long)]) +
                dnorm(x = alpha[[k]][(s-1), l],
                      mean = 0.00,
                      sd = sqrt(sigma2_alpha),
                      log = TRUE)
         
         alpha[[k]][s,l]<-rnorm(n = 1,
                                mean = alpha[[k]][(s-1), l],
                                sd = metrop_sd_alpha[k,l])

         log_pi_mat_temp[,k]<-q_long%*%alpha[[k]][s,]
         for(m in 1:6){

            log_pi_mat[,m]<-log_pi_mat_temp[,m] -
                            apply(log_pi_mat_temp, 1, logSumExp)
            pi_mat[,m]<-1.00/rowSums(exp(log_pi_mat_temp - log_pi_mat_temp[,m]))

            }
         numer<-sum(log_pi_mat[cbind(seq_along(G_long), G_long)]) +
           dnorm(x = alpha[[k]][s,l],
                 mean = 0.00,
                 sd = sqrt(sigma2_alpha),
                 log = TRUE)
         
         accept<-1
         ratio<-exp(numer - denom)
         uni_draw<-runif(n = 1,
                         min = 0.00,
                         max = 1.00)
         if(ratio < uni_draw){
           
           log_pi_mat_temp<-log_pi_mat_temp_old
           log_pi_mat<-log_pi_mat_old
           pi_mat<-pi_mat_old
           alpha[[k]][s,l]<-alpha[[k]][(s-1), l]
           accept<-0
           
         }
         
         acctot_alpha[k,l]<-acctot_alpha[k,l] + 
           accept
         
      }
      
   }
   
   ###
   #h0
   ###
   logith0[s,] <- logith0[(s-1),]
   h0<-1.00/(1.00 + exp(-logith0[s,]))
   h0_old<-h0
   
   G_a_long_old<-G_a_long
   mu_old<-sapply(1:sum(N), 
                  function(j){
                    W[j,]%*%beta[[G_a_long_old[j]]][s,]
                  })
   var_old<-sigma2[s, G_a_long_old]
   denom<-dnorm(x = Y_long,
                mean = mu_old,
                sd = sqrt(var_old),
                log = TRUE) +
     log(as.numeric((D_long == 1 & G_a_long_old == 1) | (D_long == 0 & G_a_long_old == 2) | (D_long == Z_long & G_a_long_old == 3))) +
     dnorm(x = logith0[s,],
           mean = (v_long%*%delta_h0[(s-1),]),
           sd = sqrt(tau2_h0[s-1]),
           log = TRUE)
   
   logith0[s,]<-rnorm(n = sum(N),
                     mean = logith0[(s-1),],
                     sd = metrop_sd_h0)
   h0<-1.00/(1.00 + exp(-logith0[s,]))
   G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0)) | ((G_long == 6) & (a_long >= l1))]<-1
   G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < h0)) | ((G_long == 6) & (a_long < h1))]<-2
   G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= h0)) | ((G_long == 5) & (a_long < l0)) | ((G_long == 6) & (a_long >= h1) & (a_long < l1))]<-3
   mu<-sapply(1:sum(N), 
              function(j){
                W[j,]%*%beta[[G_a_long[j]]][s,]
              })
   var<-sigma2[s, G_a_long]
   numer<-dnorm(x = Y_long,
                mean = mu_old,
                sd = sqrt(var_old),
                log = TRUE) +
     log(as.numeric((D_long == 1 & G_a_long == 1) | (D_long == 0 & G_a_long == 2) | (D_long == Z_long & G_a_long == 3))) +
     dnorm(x = logith0[s,],
           mean = (v_long%*%delta_h0[(s-1),]),
           sd = sqrt(tau2_h0[s-1]),
           log = TRUE)
   accept<-1
   ratio<-exp(numer - denom)
   uni_draw<-runif(n = 1,
                   min = 0.00,
                   max = 1.00)
   for (n in 1:sum(N)) {
     if(ratio[n] < uni_draw){
      
       G_a_long[n]<-G_a_long_old[n]
       logith0[s,n]<-logith0[(s-1),n]
       h0[n]<-h0_old[n]
       accept<-0
       
     }
     
     acctot_h0<-acctot_h0 + accept
   }
   ###
   #l0
   ###
   logitl0[s,] <- logitl0[(s-1),]
   l0<-1.00/(1.00 + exp(-logitl0[s,]))
   l0_old<-l0
   
   G_a_long_old<-G_a_long
   mu_old<-sapply(1:sum(N), 
                  function(j){
                    W[j,]%*%beta[[G_a_long_old[j]]][s,]
                  })
   var_old<-sigma2[s, G_a_long_old]
   denom<-dnorm(x = Y_long,
                mean = mu_old,
                sd = sqrt(var_old),
                log = TRUE) +
     log(as.numeric((D_long == 1 & G_a_long_old == 1) | (D_long == 0 & G_a_long_old == 2) | (D_long == Z_long & G_a_long_old == 3))) +
     dnorm(x = logitl0[s,],
           mean = (v_long%*%delta_l0[(s-1),]),
           sd = sqrt(tau2_l0[s-1]),
           log = TRUE)
   
   logitl0[s,]<-rnorm(n = sum(N),
                     mean = logitl0[(s-1),],
                     sd = metrop_sd_l0)
   l0<-1.00/(1.00 + exp(-logitl0[s,]))
   G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0)) | ((G_long == 6) & (a_long >= l1))]<-1
   G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < l0)) | ((G_long == 6) & (a_long < h1))]<-2
   G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= l0)) | ((G_long == 5) & (a_long < l0)) | ((G_long == 6) & (a_long >= h1) & (a_long < l1))]<-3
   mu<-sapply(1:sum(N), 
              function(j){
                W[j,]%*%beta[[G_a_long[j]]][s,]
              })
   var<-sigma2[s, G_a_long]
   numer<-dnorm(x = Y_long,
                mean = mu_old,
                sd = sqrt(var_old),
                log = TRUE) +
     log(as.numeric((D_long == 1 & G_a_long == 1) | (D_long == 0 & G_a_long == 2) | (D_long == Z_long & G_a_long == 3))) +
     dnorm(x = logitl0[s,],
           mean = (v_long%*%delta_l0[(s-1),]),
           sd = sqrt(tau2_l0[s-1]),
           log = TRUE)
   accept<-1
   ratio<-exp(numer - denom)
   uni_draw<-runif(n = 1,
                   min = 0.00,
                   max = 1.00)
   for (n in 1:sum(N)) {
     if(ratio[n] < uni_draw){
       
       G_a_long[n]<-G_a_long_old[n]
       logitl0[s,n]<-logitl0[(s-1),n]
       l0[n]<-l0_old[n]
       accept<-0
       
     }
     acctot_l0<-acctot_l0 + accept
   }
   ###
   #h1
   ###
   logith1[s,] <- logith1[(s-1),]
   h1<-1.00/(1.00 + exp(-logith1[s,]))
   h1_old<-h1
   
   G_a_long_old<-G_a_long
   mu_old<-sapply(1:sum(N), 
                  function(j){
                    W[j,]%*%beta[[G_a_long_old[j]]][s,]
                  })
   var_old<-sigma2[s, G_a_long_old]
   denom<-dnorm(x = Y_long,
                mean = mu_old,
                sd = sqrt(var_old),
                log = TRUE) +
     log(as.numeric((D_long == 1 & G_a_long_old == 1) | (D_long == 0 & G_a_long_old == 2) | (D_long == Z_long & G_a_long_old == 3))) +
     dnorm(x = logith1[s,],
           mean = (v_long%*%delta_h1[(s-1),]),
           sd = sqrt(tau2_h1[s-1]),
           log = TRUE)
   
   logith1[s,]<-rnorm(n = sum(N),
                     mean = logith1[(s-1),],
                     sd = metrop_sd_h1)
   h1<-1.00/(1.00 + exp(-logith1[s,]))
   G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0)) | ((G_long == 6) & (a_long >= l1))]<-1
   G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < h1)) | ((G_long == 6) & (a_long < h1))]<-2
   G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= h1)) | ((G_long == 5) & (a_long < l0)) | ((G_long == 6) & (a_long >= h1) & (a_long < l1))]<-3
   mu<-sapply(1:sum(N), 
              function(j){
                W[j,]%*%beta[[G_a_long[j]]][s,]
              })
   var<-sigma2[s, G_a_long]
   numer<-dnorm(x = Y_long,
                mean = mu_old,
                sd = sqrt(var_old),
                log = TRUE) +
     log(as.numeric((D_long == 1 & G_a_long == 1) | (D_long == 0 & G_a_long == 2) | (D_long == Z_long & G_a_long == 3))) +
     dnorm(x = logith1[s,],
           mean = (v_long%*%delta_h1[(s-1),]),
           sd = sqrt(tau2_h1[s-1]),
           log = TRUE)
   accept<-1
   ratio<-exp(numer - denom)
   uni_draw<-runif(n = 1,
                   min = 0.00,
                   max = 1.00)
   for (n in 1:sum(N)) {
     if(ratio[n] < uni_draw){
       
       G_a_long[n]<-G_a_long_old[n]
       logith1[s,n]<-logith1[(s-1),n]
       h1[n]<-h1_old[n]
       accept<-0
       
     }
     
     acctot_h1<-acctot_h1 + accept
   }
   
   ###
   #l1
   ###
   logitl1[s,] <- logitl1[(s-1),]
   l1<-(h1 + 1.00)/(1.00 + exp(-logitl1[s,]))
   l1_old<-l1
   
   G_a_long_old<-G_a_long
   mu_old<-sapply(1:sum(N), 
                  function(j){
                    W[j,]%*%beta[[G_a_long_old[j]]][s,]
                  })
   var_old<-sigma2[s, G_a_long_old]
   denom<-dnorm(x = Y_long,
                mean = mu_old,
                sd = sqrt(var_old),
                log = TRUE) +
     log(as.numeric((D_long == 1 & G_a_long_old == 1) | (D_long == 0 & G_a_long_old == 2) | (D_long == Z_long & G_a_long_old == 3))) +
     dnorm(x = logitl1[s,],
           mean = (v_long%*%delta_l1[(s-1),]),
           sd = sqrt(tau2_l1[s-1]),
           log = TRUE)
   
   logitl1[s,]<-rnorm(n = sum(N),
                     mean = logitl1[(s-1),],
                     sd = metrop_sd_l1)
   l1<-(h1 + 1.00)/(1.00 + exp(-logitl1[s,]))
   G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0)) | ((G_long == 6) & (a_long >= l1))]<-1
   G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < l1)) | ((G_long == 6) & (a_long < h1))]<-2
   G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= l1)) | ((G_long == 5) & (a_long < l0)) | ((G_long == 6) & (a_long >= h1) & (a_long < l1))]<-3
   mu<-sapply(1:sum(N), 
              function(j){
                W[j,]%*%beta[[G_a_long[j]]][s,]
              })
   var<-sigma2[s, G_a_long]
   numer<-dnorm(x = Y_long,
                mean = mu_old,
                sd = sqrt(var_old),
                log = TRUE) +
     log(as.numeric((D_long == 1 & G_a_long == 1) | (D_long == 0 & G_a_long == 2) | (D_long == Z_long & G_a_long == 3))) +
     dnorm(x = logitl1[s,],
           mean = (v_long%*%delta_l1[(s-1),]),
           sd = sqrt(tau2_l1[s-1]),
           log = TRUE)
   accept<-1
   ratio<-exp(numer - denom)
   uni_draw<-runif(n = 1,
                   min = 0.00,
                   max = 1.00)
   for (n in 1:sum(N)) {
     if(ratio[n] < uni_draw){
       
       G_a_long[n]<-G_a_long_old[n]
       logitl1[s,n]<-logitl1[(s-1),n]
       l1[n]<-l1_old[n]
       accept<-0
       
     }
     
     acctot_l1<-acctot_l1 + accept
   }
   
   #########
   #delta_h0
   #########
   ## MVN update
   
   cov_dh0 <- chol2inv(chol(crossprod(v_long)/tau2_h0[(s-1)] + diag(2)/sigma2_delta))
   mu_dh0 <- cov_dh0%*%(crossprod(v_long, logith0[s,]))/tau2_h0[(s-1)]
   delta_h0[s,] <- rmnorm(n = 1,
                         mean = mu_dh0,
                         varcov = cov_dh0)   
   #########
   #delta_l0
   #########
   ## MVN update
   
   cov_dl0 <- chol2inv(chol(crossprod(v_long)/tau2_l0[(s-1)] + diag(2)/sigma2_delta))
   mu_dl0 <- cov_dl0%*%(crossprod(v_long, logitl0[s,]))/tau2_l0[(s-1)]   
   delta_l0[s,] <- rmnorm(n = 1,
                          mean = mu_dl0,
                          varcov = cov_dl0)  
   
   #########
   #delta_h1
   #########
   ## MVN update
   
   cov_dh1 <- chol2inv(chol(crossprod(v_long)/tau2_h1[(s-1)] + diag(2)/sigma2_delta))
   mu_dh1 <- cov_dh1%*%(crossprod(v_long, logith1[s,]))/tau2_h1[(s-1)]  
   delta_h1[s,] <- rmnorm(n = 1,
                          mean = mu_dh1,
                          varcov = cov_dh1) 
   
   #########
   #delta_l1
   #########
   ## MVN update
   cov_dl1 <- chol2inv(chol(crossprod(v_long)/tau2_h1[(s-1)] + diag(2)/sigma2_delta))
   mu_dl1 <- cov_dl1%*%(crossprod(v_long, logith1[s,]))/tau2_h1[(s-1)] 
   delta_l1[s,] <- rmnorm(n = 1,
                          mean = mu_dl1,
                          varcov = cov_dl1)
   
   ########
   #tau2_h0
   ######## 
   rate<-crossprod(logith0[s,] - v_long%*%delta_h0[s,])/2.00 +
     0.01
   shape<-sum(N)/2.00 +
     0.01
   tau2_h0[s]<-1.00/rgamma(n = 1,
                            shape = shape,
                            rate = rate)

   ########
   #tau2_l0
   ######## 
   rate<-crossprod(logitl0[s,] - v_long%*%delta_l0[s,])/2.00 +
     0.01
   shape<-sum(N)/2.00 +
     0.01
   tau2_l0[s]<-1.00/rgamma(n = 1,
                            shape = shape,
                            rate = rate)
   
   ########
   #tau2_h1
   ######## 
   rate<-crossprod(logith1[s,] - v_long%*%delta_h1[s,])/2.00 +
     0.01
   shape<-sum(N)/2.00 +
     0.01
   tau2_h1[s]<-1.00/rgamma(n = 1,
                            shape = shape,
                            rate = rate)
   
   ########
   #tau2_l1
   ######## 
   rate<-crossprod(logitl1[s,] - v_long%*%delta_l1[s,])/2.00 +
     0.01
   shape<-sum(N)/2.00 +
     0.01
   tau2_l1[s]<-1.00/rgamma(n = 1,
                            shape = shape,
                            rate = rate)
   
   print("#####")
   print(round(100*(s/mcmc_samples), 2))
   # print(round(100*(min(acctot_alpha)/s), 2))
   # print(round(100*(max(acctot_alpha)/s), 2))
   # print(round(100*(min(acctot_delta_h0)/s), 2))
   # print(round(100*(max(acctot_delta_h0)/s), 2))
   # print(round(100*(min(acctot_delta_l0)/s), 2))
   # print(round(100*(max(acctot_delta_l0)/s), 2))
   # print(round(100*(min(acctot_delta_h1)/s), 2))
   # print(round(100*(max(acctot_delta_h1)/s), 2))
   # print(round(100*(min(acctot_delta_l1)/s), 2))
   # print(round(100*(max(acctot_delta_l1)/s), 2))
   
   ##########
   #Estimands
   ##########
   
   Y1 <- rep(NA, sum(N))
   Y0 <- rep(NA, sum(N))
   C <- rep(NA, sum(N))
   
   eff.a <- a[1]
   eff.s <- S[1]
   ij <- 0
   for (j in 1:J) {
     for (i in 1:N[J]) {
       ij <- ij + 1
       ## Need to figure out G(eff.a)
       if (G[[j]][s,i] %in% 1:3) {
         G.eff.a <- G[[j]][s,i]
       }  else if (G[[j]][s,i]==4) {
         if (eff.a < h0[ij]) {
           G.eff.a <- 2
         } else {
           G.eff.a <- 3
         }
       } else if (G[[j]][s,i]==5) {
         if (eff.a < l0[ij]) {
           G.eff.a <- 3
         } else {
           G.eff.a <- 1
         }
       } else if (G[[j]][s,i]==6) {
         if (eff.a < h1[ij]) {
           G.eff.a <- 2
         } else if (eff.a < l1[ij]) {
           G.eff.a <- 3
         } else {
           G.eff.a <- 1
         }
       }
       
       if (G.eff.a == 3) {
         C[ij] <- 1
       } else {
         C[ij] <- 0
       }
       
       W0 <- W1 <- W[j,]
       W0[3] <- W1[3] <- eff.s
       W0[4] <- W1[4] <- eff.a
       W0[5] <- 0
       W1[5] <- 1
       mu0<-W0%*%beta[[G.eff.a]][s,]
       var0<-sigma2[G.eff.a]
       Y0[ij]<-rnorm(n = 1,
                     mean = mu0,
                     sd = sqrt(var0)) 
       mu1<-W1%*%beta[[G.eff.a]][s,]
       var1<-sigma2[G.eff.a]
       Y1[ij]<-rnorm(n = 1,
                     mean = mu1,
                     sd = sqrt(var1))
     }
   }
   CADE[s] <- sum(C*(Y1-Y0))/sum(C)
}


#write.csv(CADE, "/home/cim24/project/OhnishiExtension/Results/JWCode_CADE_100k.csv")
saveRDS(list(beta, sigma2, alpha, delta_h0, delta_l0, delta_h1, delta_l1,
             tau2_h0, tau2_l0, tau2_h1, tau2_l1), 
        "/home/cim24/project/OhnishiExtension/Results/NEWCode_params_100k.RDS")

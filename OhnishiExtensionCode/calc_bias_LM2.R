library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
J <- as.integer(args[2])
Nj <- as.integer(args[3])
nalpha <- as.integer(args[4])

v <- paste0("LM2_6.7N", J*Nj)
id <- i
source("~/project/OhnishiExtension/JWCode/Data_Simulation_LM2.R")

beta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/beta", id, ".rds"))
beta_bias <- matrix(NA, nrow=6, ncol=5)
for (r in 1:6) {
  for (c in 1:5) {
    beta_bias[r,c] <- median(sapply(unlist(lapply(beta, function(x) x[r,c]))[5000:10000], 
    function(x) (x-beta_true[r,c])/abs(beta_true[r,c])*100))
  }
}


sig2 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/sigma2", id, ".rds"))
sig2_bias <- rep(NA, 6)
for (k in 1:6) {
    sig2_bias[k] <- median(sapply(unlist(lapply(sig2, function(x) x[k]))[5000:10000],
                                  function(x) (x-sigma2_true[k])/abs(sigma2_true[k])*100))
}

alpha <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/alpha", id, ".rds"))
alpha_bias <- matrix(NA, nrow=2, ncol=5)
for (r in 1:2) {
  for (c in 1:5) {
    alpha_bias[r,c] <- median(sapply(unlist(lapply(alpha, function(x) x[r,c]))[5000:10000], 
                                    function(x) (x-alpha_true[c,r])/abs(alpha_true[c,r])*100))
  }
}

deltah0 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/deltah0", id, ".rds"))
deltal0 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/deltal0", id, ".rds"))
deltah1 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/deltah1", id, ".rds"))
deltal1 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/deltal1", id, ".rds"))
delta_bias <- matrix(NA, nrow=4, ncol=2)

for (c in 1:2) {
  delta_bias[1,c] <- median(sapply(unlist(lapply(deltah0, function(x) x[c]))[5000:10000], 
                                  function(x) (x-delta_h0_true[c])/abs(delta_h0_true[c])*100))
  delta_bias[2,c] <- median(sapply(unlist(lapply(deltal0, function(x) x[c]))[5000:10000], 
                                   function(x) (x-delta_l0_true[c])/abs(delta_l0_true[c])*100))
  delta_bias[3,c] <- median(sapply(unlist(lapply(deltah1, function(x) x[c]))[5000:10000], 
                                  function(x) (x-delta_h1_true[c])/abs(delta_h1_true[c])*100))
  delta_bias[4,c] <- median(sapply(unlist(lapply(deltal1, function(x) x[c]))[5000:10000], 
                                   function(x) (x-delta_l1_true[c])/abs(delta_l1_true[c])*100))
}

tau2h0 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/tau2h0", id, ".rds"))
tau2l0 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/tau2l0", id, ".rds"))
tau2h1 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/tau2h1", id, ".rds"))
tau2l1 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/tau2l1", id, ".rds"))
tau2_bias <- matrix(NA, nrow=4, ncol=1)

  tau2_bias[1] <- median(sapply(unlist(lapply(tau2h0, function(x) x))[5000:10000], 
                                   function(x) (x-tau2_h0_true)/abs(tau2_h0_true)*100))
  tau2_bias[2] <- median(sapply(unlist(lapply(tau2l0, function(x) x))[5000:10000], 
                                   function(x) (x-tau2_l0_true)/abs(tau2_l0_true)*100))
  tau2_bias[3] <- median(sapply(unlist(lapply(tau2h1, function(x) x))[5000:10000], 
                                   function(x) (x-tau2_h1_true)/abs(tau2_h1_true)*100))
  tau2_bias[4] <- median(sapply(unlist(lapply(tau2l1, function(x) x))[5000:10000], 
                                   function(x) (x-tau2_l1_true)/abs(tau2_l1_true)*100))




saveRDS(c(c(t(beta_bias)), sig2_bias, c(alpha_bias), c(t(delta_bias)), tau2_bias),
paste0("~/project/OhnishiExtension/Results/", v, "/param_bias", id, ".rds"))


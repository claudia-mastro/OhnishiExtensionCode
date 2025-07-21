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

mcmc_samples<-100000
burnin <- 50000
thin <- 10
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

v <- paste0("LM2_new", J*Nj, "_2")
id <- i
source("~/OhnishiExtensionCode/Data_Simulation_LM2.R")

beta <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/beta", id, ".rds"))
beta_bias <- matrix(NA, nrow=6, ncol=5)
for (r in 1:6) {
  for (c in 1:5) {
    beta_bias[r,c] <- median(sapply(unlist(lapply(beta, function(x) x[r,c]))[burnin:mcmc_samples], 
    function(x) (x-beta_true[r,c])))
  }
}


sig2 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/sigma2", id, ".rds"))
sig2_bias <- rep(NA, 6)
for (k in 1:6) {
    sig2_bias[k] <- median(sapply(unlist(lapply(sig2, function(x) x[k]))[burnin:mcmc_samples],
                                  function(x) (x-sigma2_true[k])))
}

alpha <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/alpha", id, ".rds"))
alpha_bias <- matrix(NA, nrow=2, ncol=5)
for (r in 1:2) {
  for (c in 1:5) {
    alpha_bias[r,c] <- median(sapply(unlist(lapply(alpha, function(x) x[r,c]))[burnin:mcmc_samples], 
                                    function(x) (x-alpha_true[c,r])))
  }
}

deltah0 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/deltah0", id, ".rds"))
deltal0 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/deltal0", id, ".rds"))
deltah1 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/deltah1", id, ".rds"))
deltal1 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/deltal1", id, ".rds"))
delta_bias <- matrix(NA, nrow=4, ncol=2)

for (c in 1:2) {
  delta_bias[1,c] <- median(sapply(unlist(lapply(deltah0, function(x) x[c]))[burnin:mcmc_samples], 
                                  function(x) (x-delta_h0_true[c])))
  delta_bias[2,c] <- median(sapply(unlist(lapply(deltal0, function(x) x[c]))[burnin:mcmc_samples], 
                                   function(x) (x-delta_l0_true[c])))
  delta_bias[3,c] <- median(sapply(unlist(lapply(deltah1, function(x) x[c]))[burnin:mcmc_samples], 
                                  function(x) (x-delta_h1_true[c])))
  delta_bias[4,c] <- median(sapply(unlist(lapply(deltal1, function(x) x[c]))[burnin:mcmc_samples], 
                                   function(x) (x-delta_l1_true[c])))
}

tau2h0 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/tau2h0", id, ".rds"))
tau2l0 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/tau2l0", id, ".rds"))
tau2h1 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/tau2h1", id, ".rds"))
tau2l1 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/tau2l1", id, ".rds"))
tau2_bias <- matrix(NA, nrow=4, ncol=1)

  tau2_bias[1] <- median(sapply(unlist(lapply(tau2h0, function(x) x))[burnin:mcmc_samples], 
                                   function(x) (x-tau2_h0_true)))
  tau2_bias[2] <- median(sapply(unlist(lapply(tau2l0, function(x) x))[burnin:mcmc_samples], 
                                   function(x) (x-tau2_l0_true)))
  tau2_bias[3] <- median(sapply(unlist(lapply(tau2h1, function(x) x))[burnin:mcmc_samples], 
                                   function(x) (x-tau2_h1_true)))
  tau2_bias[4] <- median(sapply(unlist(lapply(tau2l1, function(x) x))[burnin:mcmc_samples], 
                                   function(x) (x-tau2_l1_true)))




saveRDS(c(c(t(beta_bias)), sig2_bias, c(alpha_bias), c(t(delta_bias)), tau2_bias),
paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/param_bias", id, ".rds"))


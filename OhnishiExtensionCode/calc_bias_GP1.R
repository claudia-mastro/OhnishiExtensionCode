library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

v <- "GP1"
id <- i
source("~/OhnishiExtensionCode/Data_Simulation_GP.R")

mu <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/mu", id, ".rds"))
mu_bias <- matrix(NA, nrow=3, ncol=5)
for (r in 1:3) {
  for (c in 1:5) {
    mu_bias[r,c] <- median(sapply(unlist(lapply(mu, function(x) x[r,c])), 
    function(x) (x-mu456[r,c])/abs(mu456[r,c])*100))
  }
}


sig2 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/sig2", id, ".rds"))
sig2_bias <- rep(NA, 6)
for (k in 1:6) {
    sig2_bias[k] <- median(sapply(unlist(lapply(sig2, function(x) x[k])), 
                               function(x) (x-sigma2_true[k])/abs(sigma2_true[k])*100))
}

beta <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/beta", id, ".rds"))
beta_bias <- matrix(NA, nrow=3, ncol=5)
for (r in 1:3) {
  for (c in 1:5) {
    beta_bias[r,c] <- median(sapply(unlist(lapply(beta[[r]], function(x) x[c])), 
                                  function(x) (x-beta_true[[r]][c])/abs(beta_true[[r]][c])*100))
  }
}

saveRDS(c(beta_bias, mu_bias, sig2_bias), paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/param_bias", id, ".rds"))

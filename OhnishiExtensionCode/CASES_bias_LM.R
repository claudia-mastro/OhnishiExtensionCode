library(mnormt)
library(matrixStats)
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
J <- as.integer(args[2])
Nj <- as.integer(args[3])
nalpha <- as.integer(args[4])

mcmc_samples<-10000
burnin <- 5000
thin <- 10
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

v <- paste0("LM_5.2_J", J, "Nj", Nj, "Nalpha", nalpha)
#v <- "LM_5.1_nospace"

CASE.S.bias <- rep(NA, length(iters))

source("~/project/OhnishiExtension/JWCode/effects_LM.R")
R <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/G", i, ".rds"))
h0 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/h0", i, ".rds"))
h1 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/h1", i, ".rds"))
l0 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l0", i, ".rds"))
l1 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l1", i, ".rds"))
beta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/beta", i, ".rds"))
sig2 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/sigma2", i, ".rds"))
id <- i

source("~/project/OhnishiExtension/JWCode/Data_Simulation_LM.R")

CASE.S.true <- CASE.S(0.4, 0.8, 0, 3, G_long_true, h0_true, l0_true, h1_true, l1_true, beta_true, sigma2_true)

for (j in 1:length(iters)) {
  s <- iters[j]
  CASE.S.est <- CASE.S(0.4, 0.8, 0, 3, R[[s]], h0[[s]], l0[[s]], h1[[s]], l1[[s]], beta[[s]], sig2[[s]])
  CASE.S.bias[j] <- (CASE.S.est-CASE.S.true)
}

saveRDS(CASE.S.bias, paste0("~/project/OhnishiExtension/Results/", v, "/CASE_S_bias", i, ".rds"))

v <- paste0("LM_5.2_J", J, "Nj", Nj, "Nalpha", nalpha)
source("~/project/OhnishiExtension/JWCode/Data_Simulation_LM.R")
CASE.S.true <- CASE.S(0.4, 0.8, 0, 3, G_long_true, h0_true, l0_true, h1_true, l1_true, beta_true, sigma2_true)
CASE.S.bias <- matrix(NA, nrow=length(iters), ncol=500)
for (i in 1:500) {
  if (file.exists(paste0("~/project/OhnishiExtension/Results/", v, "/CASE_S_bias", i, ".rds"))) {
    CASE.S.bias[,i] <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/CASE_S_bias", i, ".rds"))
  }
}
mean(apply(CASE.S.bias, 2, median))
mean(apply(CASE.S.bias/CASE.S.true*100, 2, median))
mean(apply(CASE.S.bias, 2, var))
lower <- apply(CASE.S.bias+CASE.S.true, 2, quantile, 0.025)
upper <- apply(CASE.S.bias+CASE.S.true, 2, quantile, 0.975)
mean((CASE.S.true >= lower) & (CASE.S.true <= upper))
mean(lower)
mean(upper)
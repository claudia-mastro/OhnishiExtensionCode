library(DescTools)

J <- 100
Nj <- 25
nalpha <- 11
v <- paste0("LM2_N", J*Nj)

mcmc_samples<-50000
burnin <- 25000
thin <- 10
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

Y_pred <- matrix(NA, nrow=sum(N), ncol=200)
for (i in 1:200) {
  if (file.exists(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CADE", i, ".rds"))) {
  id <- i
  source("~/OhnishiExtensionCode/Data_Simulation_LM2.R")
  beta <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/beta", i, ".rds"))
  G <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/G", i, ".rds"))
  Ys <- matrix(NA, nrow=sum(N), ncol=length(iters))
  si <- 1
  for (s in iters) {
    Ysi <- rep(NA, sum(N))
    for (k in 1:6) {
     Ysi[c(G[[s]])==k] <- W[c(G[[s]])==k,] %*% beta[[s]][k,]
    }
    Ys[,si] <- Ysi
    si <- si + 1
  }
  Y_pred[,i] <- (Y_long - rowMeans(Ys))
  }
}
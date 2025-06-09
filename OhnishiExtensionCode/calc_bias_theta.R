library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
id <- 1
source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2.R")

theta_med <- matrix(NA, nrow=500, ncol=200)
for (j in 1:200) {
  theta <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/theta", j, ".rds"))
  for (i in 1:500) {
    Ri <- R_long_true[i]
    theta_med[i,j] <- (median(unlist(lapply(theta, function(x) x[i, Ri]))) - theta_true[i, Ri])/
      theta_true[i, Ri]*100
  }
  rm(theta)
}
saveRDS(rowMeans(theta_med, na.rm=TRUE), "~/project/OhnishiExtension/Results/GP2_4.15_ERT/theta_bias_trueR_byi.rds")

theta_bias <- matrix(NA, nrow=500, ncol=200)
for (j in 1:200) {
  theta <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/theta", j, ".rds"))
  R <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/R", j, ".rds"))
  for (i in 1:500) {
    Ri <- unlist(lapply(R, function(x) x[i]))
    thetai_bias <- rep(NA, 1000)
    for (s in 1:1000) {
      thetai_bias[s] <- (theta[[s]][i, Ri[s]] - theta_true[i, Ri[s]])/theta_true[i, Ri[s]]*100
    }
    theta_bias[i, j] <- mean(thetai_bias)
  }
  rm(theta)
  rm(R)
}
saveRDS(theta_bias, "~/project/OhnishiExtension/Results/GP2_4.15_ERT/theta_bias_predR.rds")
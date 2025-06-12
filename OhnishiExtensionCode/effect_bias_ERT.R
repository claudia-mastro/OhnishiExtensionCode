library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)

v <- "GP2_4.15_ERT"
source("~/project/OhnishiExtension/JWCode/effects.R")

burnin <- 5000
thin <- 10
iters <- burnin:10000
iters <- iters[seq(1, 10000-burnin + thin, thin)]

bias_CADE <- matrix(NA, length(iters), 200)
bias_CASE <- matrix(NA, length(iters), 200)
for (i in 1:200) {
  id <- i
  print(i)
  if (file.exists(paste0("~/project/OhnishiExtension/Results/", v, "/CADECASE", id, ".rds"))) {
    source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2_ERT.R")
    phi_true_mat <- matrix(0, 5, 6)
    phi_true_mat[1:3,1] <- phi_true[[1]]
    phi_true_mat[1:3,2] <- phi_true[[2]]
    phi_true_mat[1:3,3] <- phi_true[[3]]
    phi_true_mat[1:4,4] <- phi_true[[4]]
    phi_true_mat[1:4,5] <- phi_true[[5]]
    phi_true_mat[1:5,6] <- phi_true[[6]]
    
    eff <- CADE.CASE.ERT(0.8, 0.4, 0.8, 0, R_long_true, h4_true, l5_true, h6_true, l6_true, 
              phi_true_mat, theta_true, mu_true, sigma2_true, psi2_true)
    est <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/CADECASE", id, ".rds"))
    bias_CADE[,i] <- (est[,1]-eff[[1]])/eff[[1]]*100
    bias_CASE[,i] <- (est[,2]-eff[[2]])/eff[[2]]*100
  }
}
saveRDS(bias_CADE, paste0("~/project/OhnishiExtension/Results/", v, "/CADE_bias.rds"))
saveRDS(bias_CASE, paste0("~/project/OhnishiExtension/Results/", v, "/CASE_bias.rds"))
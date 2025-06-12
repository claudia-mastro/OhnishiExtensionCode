library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
id <- as.integer(args[1])
source("~/project/OhnishiExtension/JWCode/effects.R")
source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2.R")
phi_true_mat <- matrix(0, 6, 6)
phi_true_mat[1:4,1] <- phi_true[[1]]
phi_true_mat[1:4,2] <- phi_true[[2]]
phi_true_mat[1:4,3] <- phi_true[[3]]
phi_true_mat[1:5,4] <- phi_true[[4]]
phi_true_mat[1:5,5] <- phi_true[[5]]
phi_true_mat[1:6,6] <- phi_true[[6]]

effd <- effs <- rep(NA, 1000)
for (i in 1:1000) {
  print(i)
  # eff <- CADE.CASE(0.8, 0.4, 0.8, 0, G_long_true, h0_true, l0_true, h1_true, l1_true,
  #                  beta_true, sigma2_true)
  eff <- CADE.CASE(0.8, 0.4, 0.8, 0, R_long_true, h4_true, l5_true, h6_true, l6_true,
            phi_true_mat, theta_true, mu_true, sigma2_true, psi2_true)
  effd[i] <- eff[[1]]
  effs[i] <- eff[[2]]
}
write.csv(mean((effd-CADE.true)/abs(CADE.true)),
          paste0("~/project/OhnishiExtension/Results/fixallbias", id, ".csv"))
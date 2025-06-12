library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

v <- "GP2_fixall"
source("~/project/OhnishiExtension/JWCode/effects.R")

burnin <- 5000
thin <- 10
iters <- burnin:10000
iters <- iters[seq(1, 10000-burnin + thin, thin)]

bias_CADE <- rep(NA, length(iters))
bias_CASE <- rep(NA, length(iters))
id <- i
print(i)
source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2.R")
phi_true_mat <- matrix(0, 6, 6)
phi_true_mat[1:4,1] <- phi_true[[1]]
phi_true_mat[1:4,2] <- phi_true[[2]]
phi_true_mat[1:4,3] <- phi_true[[3]]
phi_true_mat[1:5,4] <- phi_true[[4]]
phi_true_mat[1:5,5] <- phi_true[[5]]
phi_true_mat[1:6,6] <- phi_true[[6]]

eff <- CADE.CASE(0.8, 0.4, 0.8, 0, R_long_true, h4_true, l5_true, h6_true, l6_true, 
                 phi_true_mat, theta_true, mu_true, sigma2_true, psi2_true)

# R <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/R", id, ".rds"))
# h4 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/h4", id, ".rds"))
# l5 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l5", id, ".rds"))
# h6 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/h6", id, ".rds"))
# l6 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l6", id, ".rds"))
# phi <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/phi", id, ".rds"))
# theta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/theta", id, ".rds"))
# mu <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/mu", id, ".rds"))
# sigma2 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/sig2", id, ".rds"))
# psi2 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/psi2", id, ".rds"))

j <- 0
for (s in iters) {
  j <- j+1
  est <- CADE.CASE(0.8, 0.4, 0.8, 0, R_long_true, h4_true, l5_true, h6_true, l6_true, 
                   phi_true_mat, theta_true, mu_true, sigma2_true, psi2_true)
  bias_CADE[j] <- (est[[1]]-eff[[1]])/abs(eff[[1]])*100
  bias_CASE[j] <- (est[[2]]-eff[[2]])/abs(eff[[2]])*100
}
saveRDS(bias_CADE, paste0("~/project/OhnishiExtension/Results/", v, "/CADE_bias", id, ".rds"))
saveRDS(bias_CASE, paste0("~/project/OhnishiExtension/Results/", v, "/CASE_bias", id, ".rds"))
library(coda)
library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)

mcmc_samples<-10000
burnin <- 5000
thin <- 10
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

J <- 100
Nj <- 25
nalpha <- 11
v <- paste0("GP2")
#v <- "LM_5.5_nospace"

CADE.bias <- rep(NA, 500)
CASE.bias <- rep(NA, 500)
CADE.abs.bias <- rep(NA, 500)
CASE.abs.bias <- rep(NA, 500)
CADE.var <- rep(NA, 500)
CASE.var <- rep(NA, 500)
CADE.true.vec <- rep(NA, 500)
CASE.true.vec <- rep(NA, 500)
CADE.lower <- rep(NA, 500)
CADE.upper <- rep(NA, 500)
CASE.lower <- rep(NA, 500)
CASE.upper <- rep(NA, 500)
CADE.coverage <- rep(NA, 500)
CASE.coverage <- rep(NA, 500)

# CADE.convergence <- rep(NA, 500)
# CASE.convergence <- rep(NA, 500)
# for (i in 1:500) {
#   if (file.exists(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/beta", i, ".rds"))) {
#     CADE <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CADE", i, ".rds"))
#     CASE <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CASE", i, ".rds"))
#     CADE.convergence[i] <- pnorm(abs(geweke.diag(CADE[-1], frac1=0.5)$z), lower.tail=FALSE)*2 > 0.05
#     CASE.convergence[i] <- pnorm(abs(geweke.diag(CASE[-1], frac1=0.5)$z), lower.tail=FALSE)*2 > 0.05
#   }
# }
# mean(CADE.convergence, na.rm=TRUE)
# mean(CASE.convergence, na.rm=TRUE)

for (i in 310:500) {
if (file.exists(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CADE", i, ".rds"))) {
  id <- i
  source("~/OhnishiExtensionCode/Data_Simulation_GP2.R")
  phi_true_mat <- matrix(0, 6, 6)
  phi_true_mat[1:4,1] <- phi_true[[1]]
  phi_true_mat[1:4,2] <- phi_true[[2]]
  phi_true_mat[1:4,3] <- phi_true[[3]]
  phi_true_mat[1:5,4] <- phi_true[[4]]
  phi_true_mat[1:5,5] <- phi_true[[5]]
  phi_true_mat[1:6,6] <- phi_true[[6]]
  
  eff <- CADE.CASE(0.8, 0.4, 0.8, 0, R_long_true, h4_true, l5_true, h6_true, l6_true, 
                   phi_true_mat, theta_true, mu_true, sigma2_true, psi2_true)
  CADE.true.vec[i] <- eff[[1]]
  CASE.true.vec[i] <- eff[[2]]
  
  CADE <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CADE", i, ".rds"))
  CASE <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CASE", i, ".rds"))
  
  CADE.lower[i] <- quantile(CADE[iters], 0.025, na.rm=TRUE)
  CADE.upper[i] <- quantile(CADE[iters], 0.975, na.rm=TRUE)
  CASE.lower[i] <- quantile(CASE[iters], 0.025, na.rm=TRUE)
  CASE.upper[i] <- quantile(CASE[iters], 0.975, na.rm=TRUE)
  CADE.coverage[i] <- (eff[[1]] >= quantile(CADE[iters], 0.025, na.rm=TRUE)) & (eff[[1]] <= quantile(CADE[iters], 0.975, na.rm=TRUE))
  CASE.coverage[i] <- (eff[[2]] >= quantile(CASE[iters], 0.025, na.rm=TRUE)) & (eff[[2]] <= quantile(CASE[iters], 0.975, na.rm=TRUE))
  
  CADE.abs.bias[i] <- median(CADE[iters]-eff[[1]])
  CASE.abs.bias[i] <- median(CASE[iters]-eff[[2]])
  CADE.bias[i] <- median(CADE[iters]-eff[[1]])/eff[[1]] * 100
  CASE.bias[i] <- median(CASE[iters]-eff[[2]])/eff[[2]] * 100
  CADE.var[i] <- var(CADE[iters])
  CASE.var[i] <- var(CASE[iters])
}
}
paste(round(eff[[1]], 2),
round(mean(CADE.abs.bias, na.rm=TRUE), 2),
round(mean(CADE.bias, na.rm=TRUE),2),
round(mean(CADE.var, na.rm=TRUE), 4),
round(mean(CADE.coverage, na.rm=TRUE),2),
round(mean(CADE.lower, na.rm=TRUE),2),
round(mean(CADE.upper, na.rm=TRUE),2), sep=" & ")

paste(round(eff[[2]], 2),
round(mean(CASE.abs.bias, na.rm=TRUE), 2),
round(mean(CASE.bias, na.rm=TRUE),2),
round(mean(CASE.var, na.rm=TRUE), 4),
round(mean(CASE.coverage, na.rm=TRUE),2),
round(mean(CASE.lower, na.rm=TRUE),2),
round(mean(CASE.upper, na.rm=TRUE),2), sep=" & ")

print(sum(!is.na(CADE.bias)))

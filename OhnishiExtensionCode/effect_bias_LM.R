library(LaplacesDemon)
library(coda)

mcmc_samples<-100000
burnin <- 50000
thin <- 10
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

J <- 100
Nj <- 25
nalpha <- 11
v <- paste0("LM2_N", J*Nj)
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

CADE.convergence <- rep(NA, 500)
CASE.convergence <- rep(NA, 500)
for (i in 1:500) {
  if (file.exists(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/beta", i, ".rds"))) {
    CADE <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CADE", i, ".rds"))
    CASE <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CASE", i, ".rds"))
    CADE.convergence[i] <- pnorm(abs(geweke.diag(CADE[-1], frac1=0.5)$z), lower.tail=FALSE)*2 > 0.05
    CASE.convergence[i] <- pnorm(abs(geweke.diag(CASE[-1], frac1=0.5)$z), lower.tail=FALSE)*2 > 0.05
  }
}
mean(CADE.convergence, na.rm=TRUE)
mean(CASE.convergence, na.rm=TRUE)

for (i in 1:500) {
if (file.exists(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CADE", i, ".rds"))) {
  id <- i
  source("~/OhnishiExtensionCode/Data_Simulation_LM2.R")
  CADE.true.vec[i] <- CADE.true
  CASE.true.vec[i] <- CASE.true
  
  CADE <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CADE", i, ".rds"))
  CASE <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/CASE", i, ".rds"))
  
  CADE.lower[i] <- quantile(CADE[iters], 0.025)
  CADE.upper[i] <- quantile(CADE[iters], 0.975)
  CASE.lower[i] <- quantile(CASE[iters], 0.025)
  CASE.upper[i] <- quantile(CASE[iters], 0.975)
  CADE.coverage[i] <- (CADE.true >= quantile(CADE[iters], 0.025)) & (CADE.true <= quantile(CADE[iters], 0.975))
  CASE.coverage[i] <- (CASE.true >= quantile(CASE[iters], 0.025)) & (CASE.true <= quantile(CASE[iters], 0.975))
  
  CADE.abs.bias[i] <- median(CADE[iters]-CADE.true)
  CASE.abs.bias[i] <- median(CASE[iters]-CASE.true)
  CADE.bias[i] <- median(CADE[iters]-CADE.true)/CADE.true * 100
  CASE.bias[i] <- median(CASE[iters]-CASE.true)/CASE.true * 100
  CADE.var[i] <- var(CADE[iters])
  CASE.var[i] <- var(CASE[iters])
  
}
}
paste(round(CADE.true, 2),
round(mean(CADE.abs.bias, na.rm=TRUE), 2),
round(mean(CADE.bias, na.rm=TRUE),2),
round(mean(CADE.var, na.rm=TRUE), 4),
round(mean(CADE.coverage, na.rm=TRUE),2),
round(mean(CADE.lower, na.rm=TRUE),2),
round(mean(CADE.upper, na.rm=TRUE),2), sep=" & ")

paste(round(CASE.true, 2),
round(mean(CASE.abs.bias, na.rm=TRUE), 2),
round(mean(CASE.bias, na.rm=TRUE),2),
round(mean(CASE.var, na.rm=TRUE), 4),
round(mean(CASE.coverage, na.rm=TRUE),2),
round(mean(CASE.lower, na.rm=TRUE),2),
round(mean(CASE.upper, na.rm=TRUE),2), sep=" & ")

print(sum(!is.na(CADE.bias)))

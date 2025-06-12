library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
id <- as.integer(args[1])
v <- "GP1_ig32"

source("~/project/OhnishiExtension/JWCode/effects_GP1.R")
source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP.R")

# burnin <- 0
# thin <- 1
# iters <- burnin:4000
# iters <- iters[seq(1, 4000-burnin + thin, thin)+1]
iters <- 1:4000

eff <- CADE.CASE(0.8, 0.4, 0.8, 0, G_long_true, h0_true, l0_true, h1_true, l1_true,
                          beta_true, sigma2_true)
Gid <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/G", id, ".rds"))
h0id <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/h4", id, ".rds"))
l0id <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l5", id, ".rds"))
h1id <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/h6", id, ".rds"))
l1id <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l6", id, ".rds"))
betaid <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/beta", id, ".rds"))
sigma2id <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/sig2", id, ".rds"))


bias_CADE <- rep(NA, length(iters))
bias_CASE <- rep(NA, length(iters))
j <- 0
ests <- matrix(NA, nrow=2, ncol=4000)
for (i in iters) {
  j <- j+1
  beta_list <- list(betaid[[1]][[i]], betaid[[2]][[i]], betaid[[3]][[i]],
                    betaid[[4]][[i]], betaid[[5]][[i]], betaid[[6]][[i]])
  est <- CADE.CASE(0.8, 0.4, 0.8, 0, Gid[[i]], h0id[[i]], l0id[[i]], h1id[[i]], l1id[[i]],
                   beta_list, sigma2id[[i]])
  ests[1,j] <- est[[1]]
  ests[2,j] <- est[[2]]
  bias_CADE[j] <- (est[[1]]-eff[[1]])/abs(eff[[1]])*100
  bias_CASE[j] <- (est[[2]]-eff[[2]])/abs(eff[[2]])*100
}
saveRDS(bias_CADE, paste0("~/project/OhnishiExtension/Results/", v, "/CADE_bias", id, ".rds"))
saveRDS(bias_CASE, paste0("~/project/OhnishiExtension/Results/", v, "/CASE_bias", id, ".rds"))
library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
id <- as.integer(args[1])

v <- "LM_4.24"
source("~/project/OhnishiExtension/JWCode/effects_LM.R")
burnin <- 50000
thin <- 50
iters <- burnin:100000
iters <- iters[seq(1, 100000-burnin + thin, thin)]
CADEa <- CADEs <- CASEs <- rep(NA, length(iters))

source("~/project/OhnishiExtension/JWCode/Data_Simulation(2).R")
CADEa_true <- CADE.A(0.8, G_long_true, h0_true, l0_true, h1_true, l1_true, beta_true, sigma2_true, 200)

R <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/G", id, ".rds"))
h4 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/h0", id, ".rds"))
l5 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l0", id, ".rds"))
h6 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/h1", id, ".rds"))
l6 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l1", id, ".rds"))
beta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/beta", id, ".rds"))
sig2 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/sigma2", id, ".rds"))

j <- 1
for (i in iters) {
  print(j)
  
  CADEa[j] <- CADE.A(0.8, c(t(R[[i]])), h4[[i]], l5[[i]], h6[[i]], l6[[i]], beta[[i]], sig2[[i]], 200)
  CADEs[j] <- CADE.S(0.4, c(t(R[[i]])), h4i, l5i, h6i, l6i, beta[[i]], sig2[[i]])
  CASEs[j] <- CASE.S(0.4, 0.8, 0, 3, c(t(R[[i]])), h4i, l5i, h6i, l6i, beta[[i]], sig2[[i]])
  
  j <- j +1
}

saveRDS(cbind(CADE, CASE), paste0("~/project/OhnishiExtension/Results/", v, "/CADECASE", id, ".rds"))

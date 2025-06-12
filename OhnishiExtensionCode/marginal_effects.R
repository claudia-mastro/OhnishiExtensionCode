library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
id <- as.integer(args[1])

v <- "GP2_4.10"
source("~/project/OhnishiExtension/JWCode/effects.R")
burnin <- 5000
thin <- 10
iters <- burnin:10000
iters <- iters[seq(1, 10000-burnin + thin, thin)]
CADEa <- CADEs <- CASEs <- rep(NA, length(iters))

source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2.R")
R <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/R", id, ".rds"))
h4 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/h4", id, ".rds"))
l5 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l5", id, ".rds"))
h6 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/h6", id, ".rds"))
l6 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l6", id, ".rds"))
l6 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/l6", id, ".rds"))
phi <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/phi", id, ".rds"))
theta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/theta", id, ".rds"))
mu <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/mu", id, ".rds"))
sig2 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/sig2", id, ".rds"))
psi2 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/psi2", id, ".rds"))

j <- 1
for (i in iters) {
  print(j)
  h4i <- invlogit(h4[[i]])
  l5i <- invlogit(l5[[i]])
  h6i <- invlogit(h6[[i]])
  l6i <- (h6i + exp(l6[[i]]))/(1+exp(l6[[i]]))
  
  CADEa[j] <- CADE.A(0.8, R[[i]], h4i, l5i, h6i, l6i, phi[[i]], theta[[i]], mu[[i]], sig2[[i]], psi2[[i]], 200)
  CADEs[j] <- CADE.S(0.4, R[[i]], h4i, l5i, h6i, l6i, phi[[i]], theta[[i]], mu[[i]], sig2[[i]], psi2[[i]])
  CASEs[j] <- CASE.S(0.4, 0.8, 0, 3, R[[i]], h4i, l5i, h6i, l6i, phi[[i]], theta[[i]], mu[[i]], sig2[[i]], psi2[[i]])
  
  j <- j +1
}

saveRDS(cbind(CADE, CASE), paste0("~/project/OhnishiExtension/Results/", v, "/CADECASE", id, ".rds"))

#install.packages("label.switching", repos = "http://cran.r-project.org")
library(mnormt)
library(matrixStats)
library(label.switching)
args <- commandArgs(trailingOnly = TRUE)
id <- as.integer(args[1])
print(id)
J <- as.integer(args[2])
print(J)
Nj <- as.integer(args[3])
print(Nj)
nalpha <- as.integer(args[4])
print(nalpha)
v <- paste0("LM2_sepN", J*Nj, "_IG23")

mcmc_samples<-50000
burnin <- 25000
thin <- 10
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

source("~/OhnishiExtensionCode/Data_Simulation_LM2.R")

beta <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/beta", id, ".rds"))
sig <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/sigma2", id, ".rds"))
alpha <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/alpha", id, ".rds"))

beta1.df <- do.call("rbind", lapply(beta, function(x) x[,1]))
beta2.df <- do.call("rbind", lapply(beta, function(x) x[,2]))
beta3.df <- do.call("rbind", lapply(beta, function(x) x[,3]))
beta4.df <- do.call("rbind", lapply(beta, function(x) x[,4]))
beta5.df <- do.call("rbind", lapply(beta, function(x) x[,5]))
sig.df <- do.call("rbind", sig)
alpha1.df <- cbind(do.call("rbind", lapply(alpha, function(x) x[1,])), 0)
alpha2.df <- cbind(do.call("rbind", lapply(alpha, function(x) x[2,])), 0)
G <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/G", id, ".rds"))
G.df <- do.call("rbind", (lapply(G, function(x) c(x))))
mcmc.pars <- array(data = NA, dim = c(mcmc_samples, 6, 8))
mcmc.pars[,,1] <- beta1.df
mcmc.pars[,,2] <- beta2.df
mcmc.pars[,,3] <- beta3.df
mcmc.pars[,,4] <- beta4.df
mcmc.pars[,,5] <- beta5.df
mcmc.pars[,,6] <- sig.df
mcmc.pars[,,7] <- alpha1.df
mcmc.pars[,,8] <- alpha2.df

#p <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/pi", id, ".rds"))
labs <- label.switching(c("AIC"), z=G.df, K=6, mcmc=mcmc.pars, data=Y_long)
pm <- permute.mcmc(mcmc.pars, labs$permutations$`AIC`)

beta <- pm$output[,,1:5]
sigma2 <- pm$output[,,6]

h0 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/h0", id, ".rds"))
l0 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/l0", id, ".rds"))
h1 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/h1", id, ".rds"))
l1 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/l1", id, ".rds"))


CADE <- rep(NA, mcmc_samples)
CASE <- rep(NA, mcmc_samples)
Y0_mat <- Y1_mat <- matrix(NA, nrow=length(iters), ncol=sum(N))
si <- 0
for (s in iters) {
  Y1 <- rep(NA, sum(N))
  Y0 <- rep(NA, sum(N))
  Y0p <- rep(NA, sum(N))
  C <- rep(NA, sum(N))

  si <- si + 1

  betas <- beta[si,,]
  G_long <- labs$clusters
  eff.a <- 0.8
  eff.s <- 0.4
  eff.sp <- 0.8
  for (ij in 1:sum(N)) {
    ## Need to figure out G(eff.a)
    if (G_long[ij] %in% 1:3) {
      G.eff.a <- G_long[ij]
    }  else if (G_long[ij]==4) {
      if (eff.a < h0[[s]][ij]) {
        G.eff.a <- 2
      } else {
        G.eff.a <- 3
      }
    } else if (G_long[ij]==5) {
      if (eff.a < l0[[s]][ij]) {
        G.eff.a <- 3
      } else {
        G.eff.a <- 1
      }
    } else if (G_long[ij]==6) {
      if (eff.a < h1[[s]][ij]) {
        G.eff.a <- 2
      } else if (eff.a < l1[[s]][ij]) {
        G.eff.a <- 3
      } else {
        G.eff.a <- 1
      }
    }

    if (G.eff.a == 3) {
      C[ij] <- 1
    } else {
      C[ij] <- 0
    }
    if (a_long[ij] == eff.a & S_long[ij] == eff.s) {
      if (Z_long[ij] == 0) {
        Y0[ij] <- Y_long[ij]
        W0p <- W1 <- W[ij,]
        W1[3] <- eff.s
        W0p[3] <- eff.sp
        W0p[4] <- W1[4] <- eff.a
        W0p[5] <- 0
        W1[5] <- 1

        mu1<-W1%*%betas[G_long[ij],]
        var1<-sigma2[si,G_long[ij]]
        Y1[ij]<-rnorm(n = 1,
                      mean = mu1,
                      sd = sqrt(var1))

        mu0p<-W0p%*%betas[G_long[ij],]
        var0p<-sigma2[si,G_long[ij]]
        Y0p[ij]<-rnorm(n = 1,
                       mean = mu0p,
                       sd = sqrt(var0p))
      } else if (Z_long[ij]==1) {
        Y1[ij] <- Y_long[ij]
        W0p <- W0 <- W[ij,]
        W0[3] <- eff.s
        W0p[3] <- eff.sp
        W0p[4] <- W0[4] <- eff.a
        W0p[5] <- W0[5] <- 0

        mu0<-W0%*%betas[G_long[ij],]
        var0<-sigma2[si,G_long[ij]]
        Y0[ij]<-rnorm(n = 1,
                      mean = mu0,
                      sd = sqrt(var0))

        mu0p<-W0p%*%betas[G_long[ij],]
        var0p<-sigma2[si,G_long[ij]]
        Y0p[ij]<-rnorm(n = 1,
                       mean = mu0p,
                       sd = sqrt(var0p))
      }
    } else {
      W0p <- W0 <- W1 <- W[ij,]
      W0[3] <- W1[3] <- eff.s
      W0p[3] <- eff.sp
      W0p[4] <- W0[4] <- W1[4] <- eff.a
      W0p[5] <- W0[5] <- 0
      W1[5] <- 1

      mu0<-W0%*%betas[G_long[ij],]
      var0<-sigma2[si,G_long[ij]]
      Y0[ij]<-rnorm(n = 1,
                    mean = mu0,
                    sd = sqrt(var0))

      mu1<-W1%*%betas[G_long[ij],]
      var1<-sigma2[si,G_long[ij]]
      Y1[ij]<-rnorm(n = 1,
                    mean = mu1,
                    sd = sqrt(var1))

      mu0p<-W0p%*%betas[G_long[ij],]
      var0p<-sigma2[si,G_long[ij]]
      Y0p[ij]<-rnorm(n = 1,
                     mean = mu0p,
                     sd = sqrt(var0p))
    }
  }
  CADE[s] <- sum(C*(Y1-Y0))/sum(C)
  CASE[s] <- sum(C*(Y0-Y0p))/sum(C)
  Y0_mat[si,] <- Y0
  Y1_mat[si,] <- Y1
}

# sigma2 <- sig
# for (s in iters) {
#   Y1 <- rep(NA, sum(N))
#   Y0 <- rep(NA, sum(N))
#   Y0p <- rep(NA, sum(N))
#   C <- rep(NA, sum(N))
#   si <- si + 1
#   eff.a <- 0.8
#   eff.s <- 0.4
#   eff.sp <- 0.8
#   G_long <- c(G[[s]])
#   for (ij in 1:sum(N)) {
#     ## Need to figure out G(eff.a)
#     if (G_long[ij] %in% 1:3) {
#       G.eff.a <- G_long[ij]
#     }  else if (G_long[ij]==4) {
#       if (eff.a < h0[[s]][ij]) {
#         G.eff.a <- 2
#       } else {
#         G.eff.a <- 3
#       }
#     } else if (G_long[ij]==5) {
#       if (eff.a < l0[[s]][ij]) {
#         G.eff.a <- 3
#       } else {
#         G.eff.a <- 1
#       }
#     } else if (G_long[ij]==6) {
#       if (eff.a < h1[[s]][ij]) {
#         G.eff.a <- 2
#       } else if (eff.a < l1[[s]][ij]) {
#         G.eff.a <- 3
#       } else {
#         G.eff.a <- 1
#       }
#     }
# 
#     if (G.eff.a == 3) {
#       C[ij] <- 1
#     } else {
#       C[ij] <- 0
#     }
#     if (a_long[ij] == eff.a & S_long[ij] == eff.s) {
#       if (Z_long[ij] == 0) {
#         Y0[ij] <- Y_long[ij]
#         W0p <- W1 <- W[ij,]
#         W1[3] <- eff.s
#         W0p[3] <- eff.sp
#         W0p[4] <- W1[4] <- eff.a
#         W0p[5] <- 0
#         W1[5] <- 1
# 
#         mu1<-W1%*%beta[[s]][G_long[ij],]
#         var1<-sigma2[[s]][G_long[ij]]
#         Y1[ij]<-rnorm(n = 1,
#                       mean = mu1,
#                       sd = sqrt(var1))
# 
#         mu0p<-W0p%*%beta[[s]][G_long[ij],]
#         var0p<-sigma2[[s]][G_long[ij]]
#         Y0p[ij]<-rnorm(n = 1,
#                        mean = mu0p,
#                        sd = sqrt(var0p))
#       } else if (Z_long[ij]==1) {
#         Y1[ij] <- Y_long[ij]
#         W0p <- W0 <- W[ij,]
#         W0[3] <- eff.s
#         W0p[3] <- eff.sp
#         W0p[4] <- W0[4] <- eff.a
#         W0p[5] <- W0[5] <- 0
# 
#         mu0<-W0%*%beta[[s]][G_long[ij],]
#         var0<-sigma2[[s]][G_long[ij]]
#         Y0[ij]<-rnorm(n = 1,
#                       mean = mu0,
#                       sd = sqrt(var0))
# 
#         mu0p<-W0p%*%beta[[s]][G_long[ij],]
#         var0p<-sigma2[[s]][G_long[ij]]
#         Y0p[ij]<-rnorm(n = 1,
#                        mean = mu0p,
#                        sd = sqrt(var0p))
#       }
#     } else {
#       W0p <- W0 <- W1 <- W[ij,]
#       W0[3] <- W1[3] <- eff.s
#       W0p[3] <- eff.sp
#       W0p[4] <- W0[4] <- W1[4] <- eff.a
#       W0p[5] <- W0[5] <- 0
#       W1[5] <- 1
# 
#       mu0<-W0%*%beta[[s]][G_long[ij],]
#       var0<-sigma2[[s]][G_long[ij]]
#       Y0[ij]<-rnorm(n = 1,
#                     mean = mu0,
#                     sd = sqrt(var0))
# 
#       mu1<-W1%*%beta[[s]][G_long[ij],]
#       var1<-sigma2[[s]][G_long[ij]]
#       Y1[ij]<-rnorm(n = 1,
#                     mean = mu1,
#                     sd = sqrt(var1))
# 
#       mu0p<-W0p%*%beta[[s]][G_long[ij],]
#       var0p<-sigma2[[s]][G_long[ij]]
#       Y0p[ij]<-rnorm(n = 1,
#                      mean = mu0p,
#                      sd = sqrt(var0p))
#     }
#   }
#   CADE[s] <- sum(C*(Y1-Y0))/sum(C)
#   CASE[s] <- sum(C*(Y0-Y0p))/sum(C)
#   Y0_mat[si,] <- Y0
#   Y1_mat[si,] <- Y1
# }

saveRDS(CADE, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/CADEaic",
                     id, ".rds"))
saveRDS(CASE, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/CASEaic",
                     id, ".rds"))
saveRDS(Y0_mat, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/Y0aic",
                     id, ".rds"))
saveRDS(Y1_mat, paste0("/home/cim24/palmer_scratch/OhnishiExtension/Results/", v,"/Y1aic",
                     id, ".rds"))
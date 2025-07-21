library(matrixStats)
library(xtable)

mcmc_samples<-100000
burnin <- 50000
thin <- 100
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

J <- 40
Nj <- 25
nalpha <- 11

v <- paste0("LM2_new", J*Nj, "_2")

acctot <- matrix(NA, nrow=500, ncol=10)
for (i in 1:500) {
  if (file.exists(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/acctot_alpha", i, ".rds"))) {
    alpha <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/acctot_alpha", i, ".rds"))
    rc <- 0
    for (r in 1:5) {
      for (c in 1:2) {
        rc <- rc + 1
        acctot[i, rc] <- mean(unlist(lapply(alpha, function(x) x[r,c]))[-1])
      }
    }
  }
}

alpha.table <- data.frame(round(colMeans(acctot, na.rm=TRUE), 4),
                          round(colMins(acctot, na.rm=TRUE), 4),
                          round(colMaxs(acctot, na.rm=TRUE), 4))
rownames(alpha.table) <- c("$\\alpha_{AT0}$", "$\\alpha_{AT1}$", "$\\alpha_{NT0}$",  
           "$\\alpha_{NT1}$", "$\\alpha_{C0}$","$\\alpha_{C1}$",
           "$\\alpha_{NT-C0}$", "$\\alpha_{NT-C1}$","$\\alpha_{C-AT0}$", 
           "$\\alpha_{C-AT1}$")
print(xtable(alpha.table), sanitize.text.function = function(x){x})

acctothl <- matrix(NA, nrow=500, ncol=4)
for (i in 25:500) {
  if (file.exists(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/h0", i, ".rds"))) {
    h0 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/h0", i, ".rds"))
    l0 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/l0", i, ".rds"))
    h1 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/h1", i, ".rds"))
    l1 <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/l1", i, ".rds"))
    
    acctothl[i, 1] <- mean(unlist(lapply(h0, function(x) mean(x, na.rm=TRUE))), na.rm=TRUE)
    acctothl[i, 2] <- mean(unlist(lapply(l0, function(x) mean(x, na.rm=TRUE))), na.rm=TRUE)
    acctothl[i, 3] <- mean(unlist(lapply(h1, function(x) mean(x, na.rm=TRUE))), na.rm=TRUE)
    acctothl[i, 4] <- mean(unlist(lapply(l1, function(x) mean(x, na.rm=TRUE))), na.rm=TRUE)
  }
}

hl.table <- data.frame(round(colMeans(acctothl, na.rm=TRUE), 4),
                          round(colMins(acctothl, na.rm=TRUE), 4),
                          round(colMaxs(acctothl, na.rm=TRUE), 4))
rownames(hl.table) <- c("$h_0$", "$\\ell_0$", "$h_1$",  
                           "$\\ell_1$")
print(xtable(hl.table), sanitize.text.function = function(x){x})

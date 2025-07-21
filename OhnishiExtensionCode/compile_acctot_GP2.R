library(matrixStats)
library(xtable)

mcmc_samples<-10000
burnin <- 5000
thin <- 10
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

J <- 40
Nj <- 25
nalpha <- 11

v <- "GP2_N1000"

acctot_phi <- matrix(NA, nrow=500)
acctot_gamma <- matrix(NA, nrow=500)
acctot_h4 <- matrix(NA, nrow=500)
acctot_l5 <- matrix(NA, nrow=500)
acctot_h6 <- matrix(NA, nrow=500)
acctot_l6 <- matrix(NA, nrow=500)
for (i in 1:500) {
  if (file.exists(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/acctots", i, ".rds"))) {
    acctots <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/acctots", i, ".rds"))
    acctot_phi[i] <- mean(unlist(lapply(acctots[[1]], function(x) mean(x))))
    acctot_gamma[i] <- mean(unlist(lapply(acctots[[2]], function(x) mean(x))))
    acctot_h4[i] <- mean(unlist(lapply(acctots[[3]], function(x) mean(x))))
    acctot_l5[i] <- mean(unlist(lapply(acctots[[4]], function(x) mean(x))))
    acctot_h6[i] <- mean(unlist(lapply(acctots[[5]], function(x) mean(x))))
    acctot_l6[i] <- mean(unlist(lapply(acctots[[6]], function(x) mean(x))))
  }
}

acctot.table <- cbind(acctot_phi, acctot_gamma, acctot_h4, acctot_l5, acctot_h6, acctot_l6)
acc.table <- data.frame(round(colMeans(acctot.table, na.rm=TRUE), 4),
                          round(colMins(acctot.table, na.rm=TRUE), 4),
                          round(colMaxs(acctot.table, na.rm=TRUE), 4))
rownames(acc.table) <- c("$\\phi$", "$\\gamma$", "$h_4$",  
           "$\\ell_5$", "$h_6$","$\\ell_6$")
print(xtable(acc.table), sanitize.text.function = function(x){x})
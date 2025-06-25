i <- 3
v <- "GP2"
mu <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/mu", i, ".rds"))
id <- i
# J <- 100
# Nj <- 25
# nalpha <- 11
source("~/OhnishiExtensionCode/Data_Simulation_GP2.R")
par(mfrow=c(1,1))

plot(0, xlim=c(1, 10000), ylim=c(-20,20), type='n', ylab="mu")
lines(unlist(lapply(mu, function(x) x[5])), type='l', col='blue3')
lines(unlist(lapply(mu, function(x) x[6])), type='l', col='purple3')
lines(unlist(lapply(mu, function(x) x[4])), type='l', col='green3')
lines(unlist(lapply(mu, function(x) x[3])), type='l', col='yellow3')
lines(unlist(lapply(mu, function(x) x[2])), type='l', col='orange3')
lines(unlist(lapply(mu, function(x) x[1])), type='l', col='red3')
abline(h=mu_true[6], col='purple', lwd=2, lty=2)
abline(h=mu_true[5], col='blue', lwd=2, lty=2)
abline(h=mu_true[4], col='green', lwd=2, lty=2)
abline(h=mu_true[3], col='yellow', lwd=2, lty=2)
abline(h=mu_true[2], col='orange', lwd=2, lty=2)
abline(h=mu_true[1], col='red', lwd=2, lty=2)


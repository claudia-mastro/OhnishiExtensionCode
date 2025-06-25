i <- 4
v <- "LM_500"
beta <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/beta", i, ".rds"))
id <- i
J <- 20
Nj <- 25
nalpha <- 11
source("~/OhnishiExtensionCode/Data_Simulation_LM.R")
par(mfrow=c(3,2))

plot(0, xlim=c(1, 10000), ylim=c(-2,2), type='n', ylab="Beta1")
lines(unlist(lapply(beta, function(x) x[3,1])), type='l', col='yellow3')
lines(unlist(lapply(beta, function(x) x[2,1])), type='l', col='blue3')
lines(unlist(lapply(beta, function(x) x[1,1])), type='l', col='red3')
abline(h=beta_true[3,1], col='yellow', lwd=2, lty=2)
abline(h=beta_true[2,1], col='blue', lwd=2, lty=2)
abline(h=beta_true[1,1], col='red', lwd=2, lty=2)

plot(0, xlim=c(1, 10000), ylim=c(-2,2), type='n', ylab="Beta2")
lines(unlist(lapply(beta, function(x) x[3,2])), type='l', col='yellow3')
lines(unlist(lapply(beta, function(x) x[2,2])), type='l', col='blue3')
lines(unlist(lapply(beta, function(x) x[1,2])), type='l', col='red3')
abline(h=beta_true[3,2], col='yellow', lwd=2, lty=2)
abline(h=beta_true[2,2], col='blue', lwd=2, lty=2)
abline(h=beta_true[1,2], col='red', lwd=2, lty=2)

plot(0, xlim=c(1, 10000), ylim=c(-1,3), type='n', ylab="Beta3")
lines(unlist(lapply(beta, function(x) x[1,3])), type='l', col=('red3'))
lines(unlist(lapply(beta, function(x) x[2,3])), type='l', col=('blue3'))
lines(unlist(lapply(beta, function(x) x[3,3])), type='l', col=('yellow3'))
abline(h=beta_true[3,3], col='yellow', lwd=2, lty=2)
abline(h=beta_true[2,3], col='blue', lwd=2, lty=2)
abline(h=beta_true[1,3], col='red', lwd=2, lty=2)

plot(0, xlim=c(1, 10000), ylim=c(-2,2), type='n', ylab="Beta4")
lines(unlist(lapply(beta, function(x) x[1,4])), type='l', col=('red3'))
lines(unlist(lapply(beta, function(x) x[2,4])), type='l', col=('blue3'))
lines(unlist(lapply(beta, function(x) x[3,4])), type='l', col=('yellow3'))
abline(h=beta_true[3,4], col='yellow', lwd=2, lty=2)
abline(h=beta_true[2,4], col='blue', lwd=2, lty=2)
abline(h=beta_true[1,4], col='red', lwd=2, lty=2)

plot(0, xlim=c(1, 10000), ylim=c(-1,3), type='n', ylab="Beta5")
lines(unlist(lapply(beta, function(x) x[1,5])), type='l', col=('red3'))
lines(unlist(lapply(beta, function(x) x[2,5])), type='l', col=('blue3'))
lines(unlist(lapply(beta, function(x) x[3,5])), type='l', col=('yellow3'))
abline(h=beta_true[3,5], col='yellow', lwd=2, lty=2)
abline(h=beta_true[2,5], col='blue', lwd=2, lty=2)
abline(h=beta_true[1,5], col='red', lwd=2, lty=2)



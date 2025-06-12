library(DescTools)
library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)

v <- "GP1_fixc"

mcmc_samples<-10000
burnin <- 5000
thin <- 10
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

# Y_pred <- matrix(NA, nrow=500, ncol=10)
# Y_bias <- matrix(NA, nrow=500, ncol=10)
# for (id in 1:10) {
#   source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP.R")
#   R <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/G", id, ".rds"))
#   beta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/beta", id, ".rds"))
#   for (i in 1:500) {
#     Ri <- Mode(unlist(lapply(R, '[[', i))[iters])
#     if (length(Ri)>1) {
#       Ri <- Ri[1]
#     }
#     if (is.na(Ri)) {
#       Ri <- mean(unlist(lapply(R, '[[', i))[iters])
#     }
#           if (Ri %in% 1:3) {
#             betai <- sapply(1:5, function(i) mean(unlist(lapply(beta[[Ri]], function(x) x[[i]]))[iters]))
#           } else if (k %in% 4:6) {
#             betai <- sapply(1:5, function(i) mean(unlist(lapply(beta[[Ri]], function(x) x[,i]))[iters]))
#           }
#     
#     Y_pred[i,id] <- W[i,]%*%thetai
#     Y_bias[,i] <- (Y_pred[,i] - Y_long)/Y_long * 100
#     
#   }
# }
# 
# Y_bias <- matrix(NA, nrow=500, ncol=200)
# for (i in 1:200) {
#   id <- i
#   source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2.R")
#   Y_bias[,i] <- (Y_pred[,i] - Y_long)/Y_long * 100
# }
# round(mean(colMeans(Y_bias)), 2)

Y_data <- matrix(NA, nrow=500, ncol=200)
for (ds in 2:69) {
  id <- ds
  source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP.R")
  Y_data[,ds] <- Y_long
}
saveRDS(Y_data, paste0("~/project/OhnishiExtension/Results/", v, "/Y_true.RDs"))

Y_data <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/Y_true.RDs"))

Y_bias <- matrix(NA, nrow=500, ncol=200)
for (j in 2:200) {
  print(j)
  R <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/G", j, ".rds"))
  beta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/beta", j, ".rds"))
  for (n in 1:500) {
    preds <- rep(NA, length(iters))
    for (s in 1:length(iters)) {
      if (R[[iters[s]]][n] %in% 1:3) {
        betai <- sapply(1:5, function(i) mean(unlist(lapply(beta[[R[[iters[s]]][n]]], function(x) x[i]))[iters]))
      } else if (R[[iters[s]]][n] %in% 4:6) {
        betai <- sapply(1:5, function(i) mean(unlist(lapply(beta[[R[[iters[s]]][n]]], function(x) x[n,i]))[iters]))
      }
      preds[s] <- W[i,]%*%betai
    }
    Y_bias[i,j] <- median((preds - Y_data[i,j])/Y_data[i,j])
  }
}

saveRDS(Y_bias, paste0("~/project/OhnishiExtension/Results/", v, "/Y_bias.RDs"))
mean(rowMeans(Y_bias, na.rm=TRUE))*100

# par(mfrow = c(4,4))
# for (i in 1:16) {
#   id <- i
#   source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2_ERT.R")
#   plot(Y_long,Y_pred[,i])
#   abline(a=0, b=1, col='red')
# }
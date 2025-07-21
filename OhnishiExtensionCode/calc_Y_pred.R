library(DescTools)
library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)

v <- "GP2"

Y_data <- matrix(NA, nrow=500, ncol=500)
for (i in 1:500) {
  id <- i
  source("~/OhnishiExtensionCode/Data_Simulation_GP2.R")
  Y_data[,i] <- Y_long
}
saveRDS(Y_data, paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/Y_true.RDs"))

Y_bias <- matrix(NA, nrow=500, ncol=500)
for (j in 285:500) {
  print(j)
  R <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/R", j, ".rds"))
  mu <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/mu", j, ".rds"))
  theta <- readRDS(paste0("~/palmer_scratch/OhnishiExtension/Results/", v, "/theta", j, ".rds"))
  for (i in 1:500) {
    Ri <- unlist(lapply(R, '[[', i))[5000:10000]
    mui <- matrix(unlist(mu), ncol=6, byrow = TRUE)[5000:10000,]
    mui <- mui[cbind(seq_along(Ri), Ri)]    
    thetai <- matrix(unlist(lapply(theta, function(x) x[i,])), ncol=6, byrow=TRUE)[5000:10000,]
    thetai <- thetai[cbind(seq_along(Ri), Ri)]    
    
    preds <- mui + thetai
    Y_bias[i,j] <- median((preds - Y_data[i,j])/Y_data[i,j])
  }
}

mean(rowMeans(Y_bias, na.rm=TRUE))*100

par(mfrow = c(4,4))
for (i in 1:16) {
  id <- i
  source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2_ERT.R")
  plot(Y_long,Y_pred[,i])
  abline(a=0, b=1, col='red')
}
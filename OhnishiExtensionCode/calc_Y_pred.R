library(DescTools)

v <- "GP2_4.15_ERT"

Y_pred <- matrix(NA, nrow=500, ncol=200)
for (j in 1:200) {
  R <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/R", j, ".rds"))
  mu <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/mu", j, ".rds"))
  theta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/theta", j, ".rds"))
  for (i in 1:500) {
    Ri <- Mode(unlist(lapply(R, '[[', i))[5000:10000])
    if (length(Ri)>1) {
      Ri <- Ri[1]
    }
    if (is.na(Ri)) {
      Ri <- mean(unlist(lapply(R, '[[', i))[5000:10000])
    }
    mui <- mean(unlist(lapply(mu, '[[', Ri))[5000:10000])
    thetai <- mean(unlist(lapply(theta, function(x) x[i, Ri]))[5000:10000])
    
    Y_pred[i,j] <- mui + thetai
  }
}

Y_bias <- matrix(NA, nrow=500, ncol=200)
for (i in 1:200) {
  id <- i
  source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2.R")
  Y_bias[,i] <- (Y_pred[,i] - Y_long)/Y_long * 100
}
round(mean(colMeans(Y_bias)), 2)

Y_data <- matrix(NA, nrow=500, ncol=200)
for (i in 163:200) {
  id <- i
  source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2.R")
  Y_data[,i] <- Y_long
}
saveRDS(Y_data, paste0("~/project/OhnishiExtension/Results/", v, "/Y_true.RDs"))

Y_bias <- matrix(NA, nrow=500, ncol=200)
for (j in 1:200) {
  R <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/R", j, ".rds"))
  mu <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/mu", j, ".rds"))
  theta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/theta", j, ".rds"))
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
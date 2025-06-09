library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
id <- 1
source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2.R")

h4_bias <- matrix(NA, nrow=500, ncol=200)
l5_bias <- matrix(NA, nrow=500, ncol=200)
h6_bias <- matrix(NA, nrow=500, ncol=200)
l6_bias <- matrix(NA, nrow=500, ncol=200)
for (j in 1:200) {
  h4 <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/h4", j, ".rds"))
  for (i in which(R_long_true==4)) {
    h4_bias[i,j] <- mean((unlist(lapply(h4, function(x) x[i])) - h4_true[i])/h4_true[i]*100)
  }
  l5 <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/l5", j, ".rds"))
  for (i in which(R_long_true==5)) {
    l5_bias[i,j] <- mean((unlist(lapply(l5, function(x) x[i])) - l5_true[i])/l5_true[i]*100)
  }
  h6 <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/h6", j, ".rds"))
  for (i in which(R_long_true==6)) {
    h6_bias[i,j] <- mean((unlist(lapply(h6, function(x) x[i])) - h6_true[i])/h6_true[i]*100)
  }
  l6 <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/l6", j, ".rds"))
  for (i in which(R_long_true==6)) {
    l6_bias[i,j] <- mean((unlist(lapply(l6, function(x) x[i])) - l6_true[i])/l6_true[i]*100)
  }
}
saveRDS(h4_bias, "~/project/OhnishiExtension/Results/GP2_4.15_ERT/h4_bias_trueR_byi.rds")
saveRDS(l5_bias, "~/project/OhnishiExtension/Results/GP2_4.15_ERT/l5_bias_trueR_byi.rds")
saveRDS(h6_bias, "~/project/OhnishiExtension/Results/GP2_4.15_ERT/h6_bias_trueR_byi.rds")
saveRDS(l6_bias, "~/project/OhnishiExtension/Results/GP2_4.15_ERT/l6_bias_trueR_byi.rds")

round(mean(rowMeans(h4_bias_trueR_byi[R_long_true==4,], na.rm=TRUE)), 2)
round(median(rowMeans(h4_bias_trueR_byi[R_long_true==4,], na.rm=TRUE)), 2)
round(sd(rowMeans(h4_bias_trueR_byi[R_long_true==4,], na.rm=TRUE)), 2)

round(mean(rowMeans(l5_bias_trueR_byi[R_long_true==5,], na.rm=TRUE)),2)
round(median(rowMeans(l5_bias_trueR_byi[R_long_true==5,], na.rm=TRUE)),2)
round(sd(rowMeans(l5_bias_trueR_byi[R_long_true==5,], na.rm=TRUE)),2)

round(mean(rowMeans(h6_bias_trueR_byi[R_long_true==6,], na.rm=TRUE)),2)
round(median(rowMeans(h6_bias_trueR_byi[R_long_true==6,], na.rm=TRUE)),2)
round(sd(rowMeans(h6_bias_trueR_byi[R_long_true==6,], na.rm=TRUE)),2)

round(mean(rowMeans(l6_bias_trueR_byi[R_long_true==6,], na.rm=TRUE)),2)
round(median(rowMeans(l6_bias_trueR_byi[R_long_true==6,], na.rm=TRUE)),2)
round(sd(rowMeans(l6_bias_trueR_byi[R_long_true==6,], na.rm=TRUE)),2)

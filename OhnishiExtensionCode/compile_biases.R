v <- "GP1_ig32"
CADE_bias <- c()
for (i in 1:100) {
  if (file.exists(paste0("~/project/OhnishiExtension/Results/", v,"/CADE_bias", i, ".rds"))) {
  eff <- readRDS(paste0(paste0("~/project/OhnishiExtension/Results/", v, "/CADE_bias", i, ".rds")))
  CADE_bias <- cbind(CADE_bias, eff)
  }
}
mean(colMeans(CADE_bias))
mean(colVars(CADE_bias))

CASE_bias <- c()
for (i in 1:100) {
  if (file.exists(paste0("~/project/OhnishiExtension/Results/", v,"/CASE_bias", i, ".rds"))) {
    eff <- readRDS(paste0(paste0("~/project/OhnishiExtension/Results/", v, "/CASE_bias", i, ".rds")))
    CASE_bias <- cbind(CASE_bias, eff)
  }
}
mean(colMeans(CASE_bias))
mean(colVars(CASE_bias))
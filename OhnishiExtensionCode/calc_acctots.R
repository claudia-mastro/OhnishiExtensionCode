res <- matrix(NA, nrow=200, ncol=36+12+6+6+6+6)
for (s in 1:200) {
  acctots <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.8/acctots", s, ".rds"))
  acctot_phi <- acctots[[1]]
  acctot_gamma <- acctots[[2]]
  acctot_h4 <- acctots[[3]]
  acctot_l5 <- acctots[[4]]
  acctot_h6 <- acctots[[5]]
  acctot_l6 <- acctots[[6]]
  k <- 1
  for (i in 1:6) {
    for (j in 1:6) {
      res[s, k] <- mean(unlist(lapply(acctot_phi, function(x) x[i,j])))
      k <- k + 1
    }
  }
  
  for (i in 1:6) {
    for (j in 1:2) {
      res[s, k] <- mean(unlist(lapply(acctot_gamma, function(x) x[i,j])))
      k <- k + 1
    }
  }
  
  for (i in 1:6) {
    res[s, k] <- mean(unlist(lapply(acctot_h4, function(x) x[i])))
    k <- k + 1
  }
  for (i in 1:6) {
    res[s, k] <- mean(unlist(lapply(acctot_l5, function(x) x[i])))
    k <- k + 1
  }
  for (i in 1:6) {
    res[s, k] <- mean(unlist(lapply(acctot_h6, function(x) x[i])))
    k <- k + 1
  }
  for (i in 1:6) {
    res[s, k] <- mean(unlist(lapply(acctot_l6, function(x) x[i])))
    k <- k + 1
  }
}
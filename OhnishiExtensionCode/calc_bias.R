mu1_med <- mu2_med <- mu3_med <- mu4_med <- mu5_med <- mu6_med <- rep(NA, 200)
for (i in 1:200) {
  mu <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/mu", i, ".rds"))
  mu1_med[i] <- median(unlist(lapply(mu, function(x) x[1]))[5000:10000])
  mu2_med[i] <- median(unlist(lapply(mu, function(x) x[2]))[5000:10000])
  mu3_med[i] <- median(unlist(lapply(mu, function(x) x[3]))[5000:10000])
  mu4_med[i] <- median(unlist(lapply(mu, function(x) x[4]))[5000:10000])
  mu5_med[i] <- median(unlist(lapply(mu, function(x) x[5]))[5000:10000])
  mu6_med[i] <- median(unlist(lapply(mu, function(x) x[6]))[5000:10000])
}
bias_mu1 <- (mu1_med[!is.na(mu1_med)] - mu_true[1])/mu_true[1] * 100
bias_mu2 <- (mu2_med[!is.na(mu2_med)] - mu_true[2])/mu_true[2] * 100
bias_mu3 <- (mu3_med[!is.na(mu3_med)] - mu_true[3])/mu_true[3] * 100
bias_mu4 <- (mu4_med[!is.na(mu4_med)] - mu_true[4])/mu_true[4] * 100
bias_mu5 <- (mu5_med[!is.na(mu5_med)] - mu_true[5])/mu_true[5] * 100
bias_mu6 <- (mu6_med[!is.na(mu6_med)] - mu_true[6])/mu_true[6] * 100

mean_mu <- c(round(mean(bias_mu1), 2), round(mean(bias_mu2), 2),
              round(mean(bias_mu3), 2), round(mean(bias_mu4), 2), 
              round(mean(bias_mu5), 2), round(mean(bias_mu6), 2))

med_mu <- c(round(median(bias_mu1), 2), round(median(bias_mu2), 2),
             round(median(bias_mu3), 2), round(median(bias_mu4), 2),
             round(median(bias_mu5), 2), round(median(bias_mu6), 2))

sd_mu <- c(round(sd(bias_mu1), 2), round(sd(bias_mu2), 2), 
            round(sd(bias_mu3), 2), round(sd(bias_mu4), 2), 
            round(sd(bias_mu5), 2), round(sd(bias_mu6), 2))

sig21_med <- sig22_med <- sig23_med <- sig24_med <- sig25_med <- sig26_med <- rep(NA, 200)
for (i in c(1:195, 210:214)) {
  sig2 <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/sig2", i, ".rds"))
  sig21_med[i] <- median(unlist(lapply(sig2, function(x) x[1]))[5000:10000])
  sig22_med[i] <- median(unlist(lapply(sig2, function(x) x[2]))[5000:10000])
  sig23_med[i] <- median(unlist(lapply(sig2, function(x) x[3]))[5000:10000])
  sig24_med[i] <- median(unlist(lapply(sig2, function(x) x[4]))[5000:10000])
  sig25_med[i] <- median(unlist(lapply(sig2, function(x) x[5]))[5000:10000])
  sig26_med[i] <- median(unlist(lapply(sig2, function(x) x[6]))[5000:10000])
}
bias_sig21 <- (sig21_med[!is.na(sig21_med)] - sigma2_true[1])/sigma2_true[1] * 100
bias_sig22 <- (sig22_med[!is.na(sig22_med)] - sigma2_true[2])/sigma2_true[2] * 100
bias_sig23 <- (sig23_med[!is.na(sig23_med)] - sigma2_true[3])/sigma2_true[3] * 100
bias_sig24 <- (sig24_med[!is.na(sig24_med)] - sigma2_true[4])/sigma2_true[4] * 100
bias_sig25 <- (sig25_med[!is.na(sig25_med)] - sigma2_true[5])/sigma2_true[5] * 100
bias_sig26 <- (sig26_med[!is.na(sig26_med)] - sigma2_true[6])/sigma2_true[6] * 100

mean_sig <- c(round(mean(bias_sig21), 2), round(mean(bias_sig22), 2), 
               round(mean(bias_sig23), 2), round(mean(bias_sig24), 2), 
               round(mean(bias_sig25), 2), round(mean(bias_sig26), 2))

med_sig <- c(round(median(bias_sig21), 2), round(median(bias_sig22), 2), 
              round(median(bias_sig23), 2), round(median(bias_sig24), 2), 
              round(median(bias_sig25), 2), round(median(bias_sig26), 2))

sd_sig <- c(round(sd(bias_sig21), 2), round(sd(bias_sig22), 2),
             round(sd(bias_sig23), 2), round(sd(bias_sig24), 2), 
             round(sd(bias_sig25), 2), round(sd(bias_sig26), 2))

psi21_med <- psi22_med <- psi23_med <- psi24_med <- psi25_med <- psi26_med <- rep(NA, 200)
for (i in c(1:195, 210:214)) {
  psi2 <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/psi2", i, ".rds"))
  psi21_med[i] <- median(unlist(lapply(psi2, function(x) x[1]))[5000:10000])
  psi22_med[i] <- median(unlist(lapply(psi2, function(x) x[2]))[5000:10000])
  psi23_med[i] <- median(unlist(lapply(psi2, function(x) x[3]))[5000:10000])
  psi24_med[i] <- median(unlist(lapply(psi2, function(x) x[4]))[5000:10000])
  psi25_med[i] <- median(unlist(lapply(psi2, function(x) x[5]))[5000:10000])
  psi26_med[i] <- median(unlist(lapply(psi2, function(x) x[6]))[5000:10000])
}
bias_psi21 <- (psi21_med[!is.na(psi21_med)] - psi2_true[1])/psi2_true[1] * 100
bias_psi22 <- (psi22_med[!is.na(psi22_med)] - psi2_true[2])/psi2_true[2] * 100
bias_psi23 <- (psi23_med[!is.na(psi23_med)] - psi2_true[3])/psi2_true[3] * 100
bias_psi24 <- (psi24_med[!is.na(psi24_med)] - psi2_true[4])/psi2_true[4] * 100
bias_psi25 <- (psi25_med[!is.na(psi25_med)] - psi2_true[5])/psi2_true[5] * 100
bias_psi26 <- (psi26_med[!is.na(psi26_med)] - psi2_true[6])/psi2_true[6] * 100

mean_psi2 <- c(round(mean(bias_psi21), 2), round(mean(bias_psi22), 2), 
               round(mean(bias_psi23), 2), round(mean(bias_psi24), 2), 
               round(mean(bias_psi25), 2), round(mean(bias_psi26), 2))

med_psi2 <- c(round(median(bias_psi21), 2), round(median(bias_psi22), 2), 
              round(median(bias_psi23), 2), round(median(bias_psi24), 2), 
              round(median(bias_psi25), 2), round(median(bias_psi26), 2))

sd_psi2 <- c(round(sd(bias_psi21), 2), round(sd(bias_psi22), 2),
             round(sd(bias_psi23), 2), round(sd(bias_psi24), 2),
             round(sd(bias_psi25), 2), round(sd(bias_psi26), 2))

delta_h4_0 <- delta_h4_1 <- delta_l5_0 <- delta_l5_1 <- delta_h6_0 <- delta_h6_1 <-
  delta_l6_0 <- delta_l6_1 <- rep(NA, 200)
for (i in c(1:195, 210:214)) {
  delta <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/delta", i, ".rds"))
  delta_h4_0[i] <- median(unlist(lapply(delta, function(x) x[1,1]))[5000:10000])
  delta_h4_1[i] <- median(unlist(lapply(delta, function(x) x[1,2]))[5000:10000])
  delta_l5_0[i] <- median(unlist(lapply(delta, function(x) x[2,1]))[5000:10000])
  delta_l5_1[i] <- median(unlist(lapply(delta, function(x) x[2,2]))[5000:10000])
  delta_h6_0[i] <- median(unlist(lapply(delta, function(x) x[3,1]))[5000:10000])
  delta_h6_1[i] <- median(unlist(lapply(delta, function(x) x[3,2]))[5000:10000])
  delta_l6_0[i] <- median(unlist(lapply(delta, function(x) x[4,1]))[5000:10000])
  delta_l6_1[i] <- median(unlist(lapply(delta, function(x) x[4,2]))[5000:10000])  
}

bias_delta_h40 <- (delta_h4_0[!is.na(delta_h4_0)] - delta_h4_true[1])/delta_h4_true[1] * 100
bias_delta_h41 <- (delta_h4_1[!is.na(delta_h4_1)] - delta_h4_true[2])/delta_h4_true[2] * 100
bias_delta_l50 <- (delta_l5_0[!is.na(delta_l5_0)] - delta_l5_true[1])/delta_l5_true[1] * 100
bias_delta_l51 <- (delta_l5_1[!is.na(delta_l5_1)] - delta_l5_true[2])/delta_l5_true[2] * 100
bias_delta_h60 <- (delta_h6_0[!is.na(delta_h6_0)] - delta_h6_true[1])/delta_h6_true[1] * 100
bias_delta_h61 <- (delta_h6_1[!is.na(delta_h6_1)] - delta_h6_true[2])/delta_h6_true[2] * 100
bias_delta_l60 <- (delta_l6_0[!is.na(delta_l6_0)] - delta_l6_true[1])/delta_l6_true[1] * 100
bias_delta_l61 <- (delta_l6_1[!is.na(delta_l6_1)] - delta_l6_true[2])/delta_l6_true[2] * 100

mean_delta <- c(round(mean(bias_delta_h40), 2), round(mean(bias_delta_h41), 2), 
                   round(mean(bias_delta_l50), 2), round(mean(bias_delta_l51), 2), 
                   round(mean(bias_delta_h60), 2), round(mean(bias_delta_h61), 2),
                   round(mean(bias_delta_l60), 2), round(mean(bias_delta_l61), 2))

med_delta <- c(round(median(bias_delta_h40), 2), round(median(bias_delta_h41), 2), 
               round(median(bias_delta_l50), 2), round(median(bias_delta_l51), 2),
               round(median(bias_delta_h60), 2), round(median(bias_delta_h61), 2),
               round(median(bias_delta_l60), 2), round(median(bias_delta_l61), 2))

sd_delta <- c(round(sd(bias_delta_h40), 2), round(sd(bias_delta_h41), 2),
              round(sd(bias_delta_l50), 2), round(sd(bias_delta_l51), 2), 
              round(sd(bias_delta_h60), 2), round(sd(bias_delta_h61), 2),
              round(sd(bias_delta_l60), 2), round(sd(bias_delta_l61), 2))

tau2_h4 <- tau2_l5 <- tau2_h6 <- tau2_l6 <- rep(NA, 200)
for (i in c(1:195, 210:214)) {
  tau2 <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/tau2", i, ".rds"))
  tau2_h4[i] <- median(unlist(lapply(tau2, function(x) x[1]))[5000:10000])
  tau2_l5[i] <- median(unlist(lapply(tau2, function(x) x[2]))[5000:10000])
  tau2_h6[i] <- median(unlist(lapply(tau2, function(x) x[3]))[5000:10000])
  tau2_l6[i] <- median(unlist(lapply(tau2, function(x) x[4]))[5000:10000])
}

bias_tau2_h4 <- (tau2_h4[!is.na(tau2_h4)] - tau2_h4_true)/tau2_h4_true * 100
bias_tau2_l5 <- (tau2_l5[!is.na(tau2_l5)] - tau2_l5_true)/tau2_l5_true * 100
bias_tau2_h6 <- (tau2_h6[!is.na(tau2_h6)] - tau2_h6_true)/tau2_h6_true * 100
bias_tau2_l6 <- (tau2_l6[!is.na(tau2_l6)] - tau2_l6_true)/tau2_l6_true * 100

mean_tau2 <- c(round(mean(bias_tau2_h4), 2),  
                   round(mean(bias_tau2_l5), 2),
                   round(mean(bias_tau2_h6), 2), 
                   round(mean(bias_tau2_l6), 2))

med_tau2 <- c(round(median(bias_tau2_h4), 2), 
               round(median(bias_tau2_l5), 2), 
               round(median(bias_tau2_h6), 2), 
               round(median(bias_tau2_l6), 2))

sd_tau2 <- c(round(sd(bias_tau2_h4), 2),
              round(sd(bias_tau2_l5), 2), 
              round(sd(bias_tau2_h6), 2), 
              round(sd(bias_tau2_l6), 2))

gamma_10 <- gamma_11 <- gamma_20 <- gamma_21 <- gamma_30 <- gamma_31 <-
  gamma_40 <- gamma_41 <- gamma_50 <- gamma_51 <- gamma_60 <- gamma_61 <- rep(NA, 200)
for (i in c(1:195, 210:214)) {
  gamma <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/gamma", i, ".rds"))
  gamma_10[i] <- median(unlist(lapply(gamma, function(x) x[1,1]))[5000:10000])
  gamma_11[i] <- median(unlist(lapply(gamma, function(x) x[1,2]))[5000:10000])
  gamma_20[i] <- median(unlist(lapply(gamma, function(x) x[2,1]))[5000:10000])
  gamma_21[i] <- median(unlist(lapply(gamma, function(x) x[2,2]))[5000:10000])
  gamma_30[i] <- median(unlist(lapply(gamma, function(x) x[3,1]))[5000:10000])
  gamma_31[i] <- median(unlist(lapply(gamma, function(x) x[3,2]))[5000:10000])
  gamma_40[i] <- median(unlist(lapply(gamma, function(x) x[4,1]))[5000:10000])
  gamma_41[i] <- median(unlist(lapply(gamma, function(x) x[4,2]))[5000:10000])  
  gamma_50[i] <- median(unlist(lapply(gamma, function(x) x[5,1]))[5000:10000])
  gamma_51[i] <- median(unlist(lapply(gamma, function(x) x[5,2]))[5000:10000])  
}

bias_gamma10 <- (gamma_10[!is.na(gamma_10)] - gamma_true[1,1])/gamma_true[1,1] * 100
bias_gamma11 <- (gamma_11[!is.na(gamma_11)] - gamma_true[1,2])/gamma_true[1,2] * 100
bias_gamma20 <- (gamma_20[!is.na(gamma_20)] - gamma_true[2,1])/gamma_true[2,1] * 100
bias_gamma21 <- (gamma_21[!is.na(gamma_21)] - gamma_true[2,2])/gamma_true[2,2] * 100
bias_gamma30 <- (gamma_30[!is.na(gamma_30)] - gamma_true[3,1])/gamma_true[3,1] * 100
bias_gamma31 <- (gamma_31[!is.na(gamma_31)] - gamma_true[3,2])/gamma_true[3,2] * 100
bias_gamma40 <- (gamma_40[!is.na(gamma_40)] - gamma_true[4,1])/gamma_true[4,1] * 100
bias_gamma41 <- (gamma_41[!is.na(gamma_41)] - gamma_true[4,2])/gamma_true[4,2] * 100
bias_gamma50 <- (gamma_50[!is.na(gamma_50)] - gamma_true[5,1])/gamma_true[5,1] * 100
bias_gamma51 <- (gamma_51[!is.na(gamma_51)] - gamma_true[5,2])/gamma_true[5,2] * 100

mean_gam <- c(round(mean(bias_gamma10), 2), round(mean(bias_gamma11), 2),
              round(mean(bias_gamma20), 2), round(mean(bias_gamma21), 2),
              round(mean(bias_gamma30), 2), round(mean(bias_gamma31), 2),
              round(mean(bias_gamma40), 2), round(mean(bias_gamma41), 2),
              round(mean(bias_gamma50), 2), round(mean(bias_gamma51), 2))

med_gam <- c(round(median(bias_gamma10), 2), round(median(bias_gamma11), 2), 
             round(median(bias_gamma20), 2), round(median(bias_gamma21), 2),
             round(median(bias_gamma30), 2), round(median(bias_gamma31), 2),
             round(median(bias_gamma40), 2), round(median(bias_gamma41), 2),
             round(median(bias_gamma50), 2), round(median(bias_gamma51), 2))

sd_gam <- c(round(sd(bias_gamma10), 2), round(sd(bias_gamma11), 2),
            round(sd(bias_gamma20), 2), round(sd(bias_gamma21), 2),
            round(sd(bias_gamma30), 2), round(sd(bias_gamma31), 2),
            round(sd(bias_gamma40), 2), round(sd(bias_gamma41), 2), 
            round(sd(bias_gamma50), 2), round(sd(bias_gamma51), 2))

phi <- matrix(NA, nrow=200, ncol=3*3 + 4 + 4 + 5)
for (j in c(1:195)) {
  phii <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/phi", j, ".rds"))
  phi[j,1:3] <- (colMeans(do.call(rbind, lapply(phii, function(x) x[,1])))[1:3] - phi_true[[1]])/
    phi_true[[1]]*100
  phi[j,4:6] <- (colMeans(do.call(rbind, lapply(phii, function(x) x[,2])))[1:3] - phi_true[[2]])/
    phi_true[[2]]*100
  phi[j,7:9] <- (colMeans(do.call(rbind, lapply(phii, function(x) x[,3])))[1:3] - phi_true[[3]])/
    phi_true[[3]]*100
  phi[j,10:13] <- (colMeans(do.call(rbind, lapply(phii, function(x) x[,4])))[1:4] - phi_true[[4]])/
    phi_true[[4]]*100
  phi[j,14:17] <- (colMeans(do.call(rbind, lapply(phii, function(x) x[,5])))[1:4] - phi_true[[5]])/
    phi_true[[5]]*100
  phi[j,18:22] <- (colMeans(do.call(rbind, lapply(phii, function(x) x[,6])))[1:5] - phi_true[[6]])/
    phi_true[[6]]*100
}  
phi_table <- cbind(colMeans(phi, na.rm=TRUE), colMedians(phi, na.rm=TRUE), colSds(phi, na.rm=TRUE))
names <- c("$\\phi_{1_{X}}$", "$\\phi_{1_{S}}$", "$\\phi_{1_{Z}}$",
                   "$\\phi_{2_{X}}$", "$\\phi_{2_{S}}$", "$\\phi_{2_{Z}}$",
                   "$\\phi_{3_{X}}$", "$\\phi_{3_{S}}$", "$\\phi_{3_{Z}}$",
                   "$\\phi_{4_{X}}$", "$\\phi_{4_{S}}$", "$\\phi_{4_{Z}}$", "$\\phi_{4_{h_4}}$",
                   "$\\phi_{4_{X}}$", "$\\phi_{4_{S}}$", "$\\phi_{4_{Z}}$", "$\\phi_{4_{l_5}}$",
                   "$\\phi_{4_{X}}$", "$\\phi_{4_{S}}$", "$\\phi_{4_{Z}}$", "$\\phi_{4_{h_6}}$", 
                   "$\\phi_{4_{l_6}}$")
row.names(phi_table) <- names
print(xtable(phi_table), sanitize.text.function = function(x){x})

  
R_prop <- matrix(NA, nrow=200, ncol=6)
for (i in 1:195) {
  R <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/R", i, ".rds"))
  R_prop[i,] <- colMeans(do.call(rbind, lapply(R, function(x) table(x))))/500
}
R_prop_true <- table(R_long_true)/500
R_prop[1,] - R_prop_true
round(rowMeans(apply(R_prop, 1, function(x) (x*100 - R_prop_true*100)), na.rm=TRUE), 2)

bias_R_prop <- apply(R_prop, 1, function(x) (x - R_prop_true)/R_prop_true *100)
round(rowMeans(bias_R_prop, na.rm=TRUE), 2)

names <- c("$\\mu_{AT}$", "$\\mu_{NT}$", "$\\mu_{C}$", "$\\mu_{NT-C}$", "$\\mu_{C-AT}$", 
           "$\\mu_{NT-C-AT}$", "$\\sigma^2_{AT}$", "$\\sigma^2_{NT}$", "$\\sigma^2_{C}$", 
           "$\\sigma^2_{NT-C}$", "$\\sigma^2_{C-AT}$", "$\\sigma^2_{NT-C-AT}$",
           "$\\psi^2_{AT}$", "$\\psi^2_{NT}$", "$\\psi^2_{C}$", "$\\psi^2_{NT-C}$", 
           "$\\psi^2_{C-AT}$", "$\\psi^2_{NT-C-AT}$",
           "$\\delta_{{h4}_0}$", "$\\delta_{{h4}_1}$", "$\\delta_{{l5}_0}$", "$\\delta_{{l5}_1}$",
           "$\\delta_{{h6}_0}$", "$\\delta_{{h6}_1}$", "$\\delta_{{l6}_0}$", "$\\delta_{{l6}_1}$",
           "$\\tau^2_{h4}$", "$\\tau^2_{l5}$", "$\\tau^2_{h6}$", "$\\tau^2_{l6}$",
           "$\\gamma_{1_0}$", "$\\gamma_{1_1}$", "$\\gamma_{2_0}$", "$\\gamma_{2_1}$", 
           "$\\gamma_{3_0}$", "$\\gamma_{3_1}$", "$\\gamma_{4_0}$", "$\\gamma_{4_1}$",
           "$\\gamma_{5_0}$", "$\\gamma_{5_1}$")
library(xtable)
out_table <- cbind(c(mean_mu, mean_sig, mean_psi2, mean_delta, mean_tau2, mean_gam),
                   c(med_mu, med_sig, med_psi2, med_delta, med_tau2, med_gam),
                   c(sd_mu, sd_sig, sd_psi2, sd_delta, sd_tau2, sd_gam))
rownames(out_table) <- names
print(xtable(out_table), sanitize.text.function = function(x){x})

theta_med <- matrix(NA, nrow=500, ncol=200)
for (j in c(1:195, 210:214)) {
  theta <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/theta", j, ".rds"))
  R <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/R", j, ".rds"))
  for (i in 1:500) {
    Ri <- R_long_true[i]
    for (s in 1:1000) {
      theta_med[i,j] <- median(unlist(lapply(theta, function(x) x[i, Ri])))
    }
  }
  rm(theta)
}

theta_bias <- matrix(NA, nrow=500, ncol=200)
for (j in c(1:195, 210:214)) {
  theta <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/theta", j, ".rds"))
  R <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_4.15_ERT/R", j, ".rds"))
  for (i in 1:500) {
    Ri <- unlist(lapply(R, function(x) x[i]))
    thetai_bias <- rep(NA, 1000)
    for (s in 1:1000) {
      thetai_bias <- theta[[s]][i, Ri[s]] - theta_true[i, Ri[s]]
    }
  }
  rm(theta)
}

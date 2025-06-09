library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

v <- "GP2_fixphi"
id <- i
source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2.R")

mu <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/mu", id, ".rds"))
sig2 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/sig2", id, ".rds"))
psi2 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/psi2", id, ".rds"))

bias_mu <- bias_sig2 <- bias_psi2 <- rep(NA, 6)
for(k in 1:6) {
  bias_mu[k] <- mean(sapply((unlist(lapply(mu, function(x) x[k]))[5000:10000]),
                          function(m) (m-mu_true[k])/abs(mu_true[k])*100))
  
  bias_sig2[k] <- mean(sapply((unlist(lapply(sig2, function(x) x[k]))[5000:10000]),
                           function(s) (s-sigma2_true[k])/abs(sigma2_true[k])*100))
  
  bias_psi2[k] <- mean(sapply((unlist(lapply(psi2, function(x) x[k]))[5000:10000]),
                              function(p) (p-psi2_true[k])/abs(psi2_true[k])*100))
}

delta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/delta", id, ".rds"))
bias_delta_h40 <- mean(sapply(unlist(lapply(delta, function(x) x[1,1]))[5000:10000],
                          function(d) (d - delta_h4_true[1])/abs(delta_h4_true[1])*100))
bias_delta_h41 <- mean(sapply(unlist(lapply(delta, function(x) x[1,2]))[5000:10000],
                              function(d) (d - delta_h4_true[2])/abs(delta_h4_true[2]) * 100))
bias_delta_l50 <- mean(sapply(unlist(lapply(delta, function(x) x[2,1]))[5000:10000],
                              function(d) (d - delta_l5_true[1])/abs(delta_l5_true[1])*100))
bias_delta_l51 <- mean(sapply(unlist(lapply(delta, function(x) x[2,2]))[5000:10000],
                              function(d) (d - delta_l5_true[2])/abs(delta_l5_true[2])*100))
bias_delta_h60 <- mean(sapply(unlist(lapply(delta, function(x) x[3,1]))[5000:10000],
                              function(d) (d - delta_h6_true[1])/abs(delta_h6_true[1])*100))
bias_delta_h61 <- mean(sapply(unlist(lapply(delta, function(x) x[3,2]))[5000:10000],
                              function(d) (d - delta_h6_true[2])/abs(delta_h6_true[2])*100))
bias_delta_l60 <- mean(sapply(unlist(lapply(delta, function(x) x[4,1]))[5000:10000],
                             function(d) (d - delta_l6_true[1])/abs(delta_l6_true[1])*100))
bias_delta_l61 <- mean(sapply(unlist(lapply(delta, function(x) x[4,2]))[5000:10000],
                              function(d) (d - delta_l6_true[2])/abs(delta_l6_true[2])*100))

tau2 <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/tau2", id, ".rds"))
bias_tau2_h4 <- mean(sapply(unlist(lapply(delta, function(x) x[1]))[5000:10000],
                            function(t) (t - tau2_h4_true)/abs(tau2_h4_true)*100))
bias_tau2_l5 <- mean(sapply(unlist(lapply(delta, function(x) x[2]))[5000:10000],
                            function(t) (t - tau2_l5_true)/abs(tau2_l5_true)*100))
bias_tau2_h6 <- mean(sapply(unlist(lapply(delta, function(x) x[3]))[5000:10000],
                            function(t) (t - tau2_h6_true)/abs(tau2_h6_true)*100))
bias_tau2_l6 <- mean(sapply(unlist(lapply(delta, function(x) x[4]))[5000:10000],
                            function(t) (t - tau2_l6_true)/abs(tau2_l6_true)*100))

gamma <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/gamma", id, ".rds"))
bias_gamma <- matrix(NA, nrow=5, ncol=2)
for (r in 1:5) {
  for (c in 1:2) {
    bias_gamma[r,c] <- mean(sapply(unlist(lapply(gamma, function(x) x[r,c]))[5000:10000],
                                   function(g) (g - gamma_true[r,c])/abs(gamma_true[r,c])*100))
  }
}
bias_gamma <- c(t(bias_gamma))

phi <- rep(NA, 4*3 + 5 + 5 + 6)
phi_sd <- rep(NA, 4*3 + 5 + 5 + 6)
phii <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/phi", id, ".rds"))
phi[1:4] <- rowMeans(apply(do.call(rbind, lapply(phii, function(x) x[,1]))[5000:10000,1:4], 1,
                  function(p) (p - phi_true[[1]])/abs(phi_true[[1]])*100))
phi[5:8] <- rowMeans(apply(do.call(rbind, lapply(phii, function(x) x[,2]))[5000:10000,1:4], 1,
                     function(p) (p - phi_true[[2]])/abs(phi_true[[2]])*100))
phi[9:12] <- rowMeans(apply(do.call(rbind, lapply(phii, function(x) x[,3]))[5000:10000,1:4], 1,
                             function(p) (p - phi_true[[3]])/abs(phi_true[[3]])*100))
phi[13:17] <- rowMeans(apply(do.call(rbind, lapply(phii, function(x) x[,4]))[5000:10000,1:5], 1,
                             function(p) (p - phi_true[[4]])/abs(phi_true[[4]])*100))
phi[18:22] <- rowMeans(apply(do.call(rbind, lapply(phii, function(x) x[,5]))[5000:10000,1:5], 1,
                       function(p) (p - phi_true[[5]])/abs(phi_true[[5]])*100))
phi[23:28] <- rowMeans(apply(do.call(rbind, lapply(phii, function(x) x[,6]))[5000:10000,1:6], 1,
                             function(p) (p - phi_true[[6]])/abs(phi_true[[6]])*100))


theta_bias <- rep(NA, sum(N))
theta <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/theta", j, ".rds"))
R <- readRDS(paste0("~/project/OhnishiExtension/Results/", v, "/R", j, ".rds"))
for (j in 1:500) {
  Ri <- unlist(lapply(R, function(x) x[j]))
  thetai_bias <- rep(NA, 1000)
  for (s in 1:1000) {
    thetai_bias[s] <- (theta[[s]][j, Ri[s]] - theta_true[j, Ri[s]])/abs(theta_true[j, Ri[s]])*100
  }
  theta_bias[j] <- mean(thetai_bias)
}

theta_table <- c(mean(theta_bias), sd(theta_bias))

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
           "$\\gamma_{5_0}$", "$\\gamma_{5_1}$",
           "$\\phi_{1_{X}}$", "$\\phi_{1_{S}}$", "$\\phi_{1_{T}}$", "$\\phi_{1_{Z}}$",
           "$\\phi_{2_{X}}$", "$\\phi_{2_{S}}$", "$\\phi_{2_{T}}$", "$\\phi_{2_{Z}}$",
           "$\\phi_{3_{X}}$", "$\\phi_{3_{S}}$", "$\\phi_{3_{T}}$", "$\\phi_{3_{Z}}$",
           "$\\phi_{4_{X}}$", "$\\phi_{4_{S}}$", "$\\phi_{4_{T}}$", "$\\phi_{4_{Z}}$", "$\\phi_{4_{h_4}}$",
           "$\\phi_{5_{X}}$", "$\\phi_{5_{S}}$", "$\\phi_{5_{T}}$", "$\\phi_{5_{Z}}$", "$\\phi_{5_{l_5}}$",
           "$\\phi_{6_{X}}$", "$\\phi_{6_{S}}$", "$\\phi_{6_{T}}$", "$\\phi_{6_{Z}}$", "$\\phi_{6_{h_6}}$", 
           "$\\phi_{6_{l_6}}$", "$\\theta$")

out_table <- c(bias_mu, bias_sig2, bias_psi2, bias_delta_h40, bias_delta_h41,
               bias_delta_l50, bias_delta_l51, bias_delta_h60, bias_delta_h61,
               bias_delta_l60, bias_delta_l61, bias_tau2_h4, bias_tau2_l5, bias_tau2_h6,
               bias_tau2_l6, bias_gamma, phi, mean(theta_bias))
out_table <- data.frame(out_table)
rownames(out_table) <- names

saveRDS(out_table, paste0("~/project/OhnishiExtension/Results/", v, "/param_bias", id, ".rds"))

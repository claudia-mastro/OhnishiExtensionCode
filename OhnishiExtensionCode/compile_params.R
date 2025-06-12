means <- c()
for (i in 1:500) {
  if (file.exists(paste0("~/project/OhnishiExtension/Results/GP2_fixphi/param_bias", i, ".rds"))) {
    f <- readRDS(paste0("~/project/OhnishiExtension/Results/GP2_fixphi/param_bias", i, ".rds"))
    means <- cbind(means, f[,1])
  }
}
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
out_table <- data.frame(rowMeans(means, na.rm=TRUE), row.names = names)
print(xtable(out_table), sanitize.text.function = function(x){x})


means <- c()
for (i in 1:500) {
  if (file.exists(paste0("~/project/OhnishiExtension/Results/GP1_fixc/param_bias", i, ".rds"))) {
    f <- readRDS(paste0("~/project/OhnishiExtension/Results/GP1_fixc/param_bias", i, ".rds"))
    means <- cbind(means, f)
  }
}
names <- c("$\\mu_{AT1}$","$\\mu_{NT1}$", "$\\mu_{C1}$", 
           "$\\mu_{ATX}$","$\\mu_{NTX}$", "$\\mu_{CX}$", 
           "$\\mu_{ATS}$","$\\mu_{NTS}$", "$\\mu_{CS}$", 
           "$\\mu_{ATT}$","$\\mu_{NTT}$", "$\\mu_{CT}$", 
           "$\\mu_{ATZ}$","$\\mu_{NTZ}$", "$\\mu_{CZ}$", 
           "$\\sigma^2_{AT}$", "$\\sigma^2_{NT}$", "$\\sigma^2_{C}$", 
           "$\\sigma^2_{NT-C}$", "$\\sigma^2_{C-AT}$", "$\\sigma^2_{NT-C-AT}$",
           "$\\psi^2_{AT}$", "$\\psi^2_{NT}$", "$\\psi^2_{C}$", "$\\psi^2_{NT-C}$", 
           "$\\psi^2_{C-AT}$", "$\\psi^2_{NT-C-AT}$",)

means <- c()
for (i in 1:500) {
  if (file.exists(paste0("~/project/OhnishiExtension/Results/LM2_5.30/param_bias", i, ".rds"))) {
    f <- readRDS(paste0("~/project/OhnishiExtension/Results/LM2_5.30/param_bias", i, ".rds"))
    means <- cbind(means, f)
  }
}
names <- c("$\\beta_{AT0}$", "$\\beta_{AT1}$", "$\\beta_{AT2}$", "$\\beta_{AT3}$", "$\\beta_{AT4}$",
           "$\\beta_{NT0}$", "$\\beta_{NT1}$", "$\\beta_{NT2}$", "$\\beta_{NT3}$", "$\\beta_{NT4}$",
           "$\\beta_{C0}$", "$\\beta_{C1}$", "$\\beta_{C2}$", "$\\beta_{C3}$", "$\\beta_{C4}$",
           "$\\beta_{NT-C0}$", "$\\beta_{NT-C1}$", "$\\beta_{NT-C2}$", "$\\beta_{NT-C3}$", "$\\beta_{NT-C4}$",
           "$\\beta_{C-AT0}$", "$\\beta_{C-AT1}$", "$\\beta_{C-AT2}$", "$\\beta_{C-AT3}$", "$\\beta_{C-AT4}$",
           "$\\beta_{NT-C-AT0}$", "$\\beta_{NT-C-AT1}$", "$\\beta_{NT-C-AT2}$", "$\\beta_{NT-C-AT3}$", "$\\beta_{NT-C-AT4}$",
           "$\\sigma^2_{AT}$", "$\\sigma^2_{NT}$", "$\\sigma^2_{C}$", 
           "$\\sigma^2_{NT-C}$", "$\\sigma^2_{C-AT}$", "$\\sigma^2_{NT-C-AT}$",
           "$\\alpha_{AT0}$", "$\\alpha_{AT1}$", "$\\alpha_{NT0}$", "$\\alpha_{NT1}$",
           "$\\alpha_{C0}$", "$\\alpha_{C1}$", "$\\alpha_{NT-C0}$", "$\\alpha_{NT-C1}$",
           "$\\alpha_{C-AT0}$", "$\\alpha_{C-AT1}$", "$\\delta_{h0_0}$", "$\\delta_{h0_1}$",
           "$\\delta_{l0_0}$", "$\\delta_{l0_1}$", "$\\delta_{h1_0}$", "$\\delta_{h1_1}$",
           "$\\delta_{l1_0}$", "$\\delta_{l1_1}$", "$\\tau2_{h0}$", "$\\tau2_{l0}$",
           "$\\tau2_{h1}$", "$\\tau2_{l1}$")
out_table <- data.frame(rowMeans(means, na.rm=TRUE), row.names = names)
print(xtable(out_table), sanitize.text.function = function(x){x})
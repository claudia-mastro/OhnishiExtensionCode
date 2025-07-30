##############
#Seed
##############
set.seed(1222)
################
#Global Settings
################

N<-rep(Nj,
       times = J)
Ns <- list(0)
Ns[[1]] <- 1:N[1]
for (j in 2:J) {
  Ns[[j]] <- (Ns[[j-1]][N[j-1]]+1):(Ns[[j-1]][N[j-1]]+N[j])
}

M<-J
a<-sample(seq(0, 1, 1/(nalpha-1)),
          size = J,
          replace = TRUE)
Z<-list(0)
T<-rep(NA,
       times = J)
for(j in 1:J){

   Z[[j]]<-rbinom(n = N[j],
                  size = 1,
                  prob = a[j])
   T[j]<-sum(Z[[j]])
   a[j]<-T[j]/N[j]

   }
Z_long<-Z[[1]]
T_long<-rep(a[1], 
            times = N[1])
for(j in 2:J){

   Z_long<-c(Z_long,
             Z[[j]])
   T_long<-c(T_long,
             rep(a[j],
                 times = N[j]))

   }

###########
#Predictors
###########
x_long<-c(rnorm(n = sum(N)/6, mean=-12, sd=0.5),
          rnorm(n = sum(N)/6, mean=-8, sd=0.5),
          rnorm(n = sum(N)/6, mean=-4, sd=1),
          rnorm(n = sum(N)/6, mean=4, sd=1),
          rnorm(n = sum(N)/6, mean=8, sd=0.5),
          rnorm(n = sum(N)/6, mean=12, sd=0.5))
x_long <- sample(x_long, sum(N), replace=TRUE)
x_long <- (x_long - mean(x_long))/sd(x_long)
# x_long <- cbind(x_long, rbinom(sum(N), 1, 0.4))
# x_long <- cbind(x_long, rbinom(sum(N), 1, 0.8))
q_long<-cbind(1, x_long)

################
#True Parameters
################
beta_true<-matrix(NA,
                  nrow = 1,
                  ncol = 12)
beta_true <- c(2, 4, 6, -9, 8, -6, -4, -2, 5, 7, 10, -8)
sigma2_true <- 0.2

delta_h_true<-c(-0.8, 1.6)
tau2_h_true<-0.1
logit_h_true<-rnorm(n = sum(N),
                     mean = (q_long%*%delta_h_true),
                     sd = sqrt(tau2_h_true))
etah_true<-1.00/(1.00 + exp(-logit_h_true))

delta_l_true<-c(0.8, 1.6)
tau2_l_true<-0.1
logit_l_true<-rnorm(n = sum(N),
                     mean = (q_long%*%delta_l_true),
                     sd = sqrt(tau2_l_true))
etal_true<-1.00/(1.00 + exp(-logit_l_true))

ph0_true <- 0.86
ph1_true <- 0.94
pl0_true <- 0.92
pl1_true <- 0.82

lamh0_true <- rbinom(sum(N), 1, ph0_true)
lamh1_true <- rbinom(sum(N), 1, ph1_true)
laml0_true <- rbinom(sum(N), 1, pl0_true)
laml1_true <- rbinom(sum(N), 1, pl1_true)
laml0_true[lamh0_true==1 & lamh1_true==0] <- 1
laml1_true[lamh0_true==1 & lamh1_true==0] <- 0

h_true <- lamh0_true * etah_true ^ lamh1_true
l_true <- laml0_true * etal_true ^ laml1_true


#####
#Data
#####
set.seed(id)
D<-list(0)
S<-rep(NA,
       times = J)
G_true <- rep(NA, sum(N))
G_true[T_long > l_true] <- 1
G_true[T_long < h_true] <- 2
G_true[T_long >= h_true & T_long <= l_true] <- 3

Ji <- sort(rep(1:J, times=Nj))
dta <- data.table(Ji, T_long, Z_long, h_true, l_true)

dta[T_long > l_true, G_long_true := 1]
dta[T_long < h_true, G_long_true := 2]
dta[T_long >= h_true & T_long <= l_true, G_long_true := 3]

dta[G_long_true==1, D_long := 1]
dta[G_long_true==2, D_long := 0]
dta[G_long_true==3, D_long := Z_long]

dta[,S_long := mean(D_long), by=Ji]

attach(dta)

W<-cbind(1, S_long, T_long, Z_long, h_true, l_true, S_long*h_true, S_long*l_true,
         T_long*h_true, T_long*l_true, Z_long*h_true, Z_long*l_true)
   
Y0_long <- rep(NA, sum(N))
Y1_long <- rep(NA, sum(N))
Y0p_long <- rep(NA, sum(N))
C <- rep(NA, sum(N))
CADE.G <- rep(NA, sum(N))
eff.a <- 0.8
eff.s <- 0.4
eff.sp <- 0.8

for(j in 1:sum(N)) {
  if (G_long_true[j] == 3) {
    C[j] <- 1
  } else {
    C[j] <- 0
  }
  W0p <- W0 <- W1 <- W[j,]
  
  W0[2] <- W1[2] <- eff.s
  W0[7] <- W1[7] <- eff.s*h_true[j]
  W0[8] <- W1[8] <- eff.s*l_true[j]
  
  W0p[2] <- eff.sp
  W0p[7] <- eff.sp*h_true[j]
  W0p[8] <- eff.sp*l_true[j]
  
  W0[3] <- W0p[3] <- W1[3] <- eff.a
  W0[9] <- W0p[9] <- W1[9] <- eff.a*h_true[j]
  W0[10] <- W0p[10] <- W1[10] <- eff.a*l_true[j]
  
  W0[4] <- W0[11] <- W0p[12] <- W0p[4] <- W0p[11] <- W0p[12] <- 0
  W1[4] <- 1
  W1[11] <- h_true[j]
  W1[12] <- l_true[j]
  
  mu0 <- W0 %*% beta_true
  var0 <- sigma2_true
  Y0_long[j] <- rnorm(n = 1,
                      mean = mu0,
                      sd = sqrt(var0))
  mu1 <- W1 %*% beta_true
  var1 <- sigma2_true
  Y1_long[j] <- rnorm(n = 1,
                      mean = mu1,
                      sd = sqrt(var1))
  mu0p <- W0p %*% beta_true
  var0p <- sigma2_true
  Y0p_long[j] <- rnorm(n = 1,
                       mean = mu0p,
                       sd = sqrt(var0p))
}

Y_long<-rep(NA,
            time = sum(N))

for(j in 1:sum(N)){
  if (T_long[j] == eff.a & S_long[j] == eff.s & G_long_true[j] == 3) {
    if (Z_long[j]==0) {
      Y_long[j] <- Y0_long[j]
    } else if (Z_long[j]==1) {
      Y_long[j] <- Y1_long[j]
    }
  } else if (T_long[j] == eff.a & S_long[j] == eff.sp & G_long_true[j] == 3 & Z_long[j]==0) {
    Y_long[j] <- Y0p_long[j]
  } else {
    mu<-W[j,]%*%beta_true
    var<-sigma2_true
    Y_long[j]<-rnorm(n = 1,
                     mean = mu,
                     sd = sqrt(var))
  }
}
Y<-list(0)
Y[[1]]<-Y_long[1:N[1]]
for(j in 2:J){
  Y[[j]]<-Y_long[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
}

CADE.true <- sum((G_long_true==3)*(Y1_long - Y0_long))/(sum(G_long_true==3))
CASE.true <- sum((G_long_true==3)*(Y0_long - Y0p_long))/(sum(G_long_true==3))










                    

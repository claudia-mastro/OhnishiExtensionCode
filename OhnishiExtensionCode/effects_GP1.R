CADE.CASE <- function(alpha, S, Sp, Z, R, h4, l5, h6, l6, beta, sig2) {
  Y1 <- rep(NA, sum(N))
  Y0 <- rep(NA, sum(N))
  Y0p <- rep(NA, sum(N))
  C <- rep(NA, sum(N))

  eff.a <- alpha
  eff.s <- S
  eff.sp <- Sp
  ij <- 0
  for (j in 1:J) {
    for (i in 1:N[j]) {
      ij <- ij + 1
      if (R[ij] %in% 1:3) {
        beta_ij<-beta[[R[ij]]]
      } else if (R[ij] %in% 4:6) {
        beta_ij<-beta[[R[ij]]][ij,]
      }
      ## Need to figure out G(eff.a)
      if (R[ij] %in% 1:3) {
        G.eff.a <- R[ij]
      }  else if (R[ij]==4) {
        if (eff.a < h4[ij]) {
          G.eff.a <- 2
        } else {
          G.eff.a <- 3
        }
      } else if (R[ij]==5) {
        if (eff.a < l5[ij]) {
          G.eff.a <- 3
        } else {
          G.eff.a <- 1
        }
      } else if (R[ij]==6) {
        if (eff.a < h6[ij]) {
          G.eff.a <- 2
        } else if (eff.a < l6[ij]) {
          G.eff.a <- 3
        } else {
          G.eff.a <- 1
        }
      }
      
      if (G.eff.a == 3) {
        C[ij] <- 1
      } else {
        C[ij] <- 0
      }
      if (C[ij]==1) {
        W0p <- W0 <- W1 <- W[ij,]
        W0[3] <- W1[3] <- eff.s
        W0p[3] <- eff.sp
        W0p[4] <- W0[4] <- W1[4] <- eff.a
        W0p[5] <- W0[5] <- 0
        W1[5] <- 1
        mu0 <- W0 %*% beta_ij
        var0 <- sig2[R[j]]
        Y0[ij] <- rnorm(n = 1,
                        mean = mu0,
                        sd = 0)
        mu1 <- W1 %*% beta_ij
        var1 <- sigma2_true[R[ij]]
        Y1[ij] <- rnorm(n = 1,
                        mean = mu1,
                        sd = 0)
        
        W0p[5] <- 0
        mu0p <- W0p %*% beta_ij
        var0p <- sig2[R[ij]]
        Y0p[ij] <- rnorm(n = 1,
                         mean = mu0p,
                         sd = 0)
        
        # if (a_long[ij] == eff.a & S_long[ij]/N[j] == eff.s) {
        #   if (Z_long[j] == 0) {
        #     Y0[ij] <- Y_long[ij]
        #   } else if (Z_long[j] == 1) {
        #     Y1[ij] <- Y_long[ij]
        #   }
        # }
        # if (a_long[ij] == eff.a & S_long[ij]/N[j] == eff.sp) {
        #   Y0p[ij] <- Y_long[ij] 
        # }
      }
    }
  }
  CADE <- sum(C*(Y1-Y0), na.rm=TRUE)/sum(C)
  CASE <- sum(C*(Y0-Y0p), na.rm=TRUE)/sum(C)
  
  return(list(CADE, CASE))
}
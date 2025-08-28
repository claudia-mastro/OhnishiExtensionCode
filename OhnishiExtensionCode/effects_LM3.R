CADE.CASE <- function(alpha, s, Sp, h, l, beta, W) {
  mu <- rep(NA, sum(N))
  mu1 <- rep(NA, sum(N))
  mu0p <- rep(NA, sum(N))
  C <- rep(NA, sum(N))

  for(ij in 1:sum(N)) {
    if (h[ij] <= alpha & l[ij] >= alpha) {
      C[ij] <- 1
    } else {
      C[ij] <- 0
    }
    W0p <- W0 <- W1 <- W[ij,]
    
    W0[2] <- W1[2] <- s
    W0[7] <- W1[7] <- s*h[ij]
    W0[8] <- W1[8] <- s*l[ij]
    
    W0p[2] <- Sp
    W0p[7] <- Sp*h[ij]
    W0p[8] <- Sp*l[ij]
    
    W0[3] <- W0p[3] <- W1[3] <- alpha
    W0[9] <- W0p[9] <- W1[9] <- alpha*h[ij]
    W0[10] <- W0p[10] <- W1[10] <- alpha*l[ij]
    
    W0[4] <- W0[11] <- W0p[12] <- W0p[4] <- W0p[11] <- W0p[12] <- 0
    W1[4] <- 1
    W1[11] <- h[ij]
    W1[12] <- l[ij]
    
    mu0[ij] <- W0 %*% beta

    mu1[ij] <- W1 %*% beta

    mu0p[ij] <- W0p %*% beta
  }
  CADE <- sum(C*(mu1-mu0), na.rm=TRUE)/sum(C)
  CASE <- sum(C*(mu0-mu0p), na.rm=TRUE)/sum(C)

  return(list(CADE, CASE))
}


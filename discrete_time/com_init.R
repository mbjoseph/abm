# com_init() initializes community and environment parameters
com_init <- function(A, Emin, Emax, N, 
                     pEmin, pEmax, P,
                     ERmin, ERmax, 
                     pERmin, pERmax,
                     sig, sig.p){
  E <- runif(A, Emin, Emax) # habitat patch environment type
  pE <- runif(N, pEmin, pEmax) # within-host environmental conditions
  mu.i <- runif(N, ERmin, ERmax) # optimum environment for hosts
  mu.p <- runif(P, pERmin, pERmax) # parasite optima
  sigma <- rep(sig, N) # runif(N, 1, 10) # # niche width
  sigma.p <- rep(sig.p, P)
  
  Z <- array(dim=c(N)) # normalization constant
  
  # now calculating all elements of Z coefficient 
  # (ensures all species have equal Pr(establishment) in regional pool)
  for (i in 1:N){
    integrand <- function(E) {
      exp(-((E - mu.i[i]) ^ 2) / (2 * sigma[i] ^ 2))
    }
    res <- integrate(integrand, lower=ERmin, upper=ERmax)
    Z[i] <- 1 / res$value
  }
  
  # probability of establishment
  Pcol <- array(dim=c(A, N))
  for (i in 1:A){
    for (j in 1:N){
      Pcol[i, j] <- Z[j] * exp(-((E[i] - mu.i[j]) ^ 2) / (2*sigma[j] ^ 2))
    }
  }
  
  # calculate parasite Z
  pZ <- rep(NA, P)
  for (i in 1:P){
    integrand <- function(pE) {
      exp(-((pE - mu.p[i]) ^ 2) / (2 * sigma.p[i] ^ 2))
    }
    res <- integrate(integrand, lower=pERmin, upper=pERmax)
    pZ[i] <- 1 / res$value
  }
  
  # parasite probability of establishment
  pPcol <- array(dim=c(N, P))
  for (i in 1:N){
    for (j in 1:P){
      pPcol[i, j] <- pZ[j] * exp(-((pE[i] - mu.p[j]) ^ 2) / (2*sigma.p[j] ^ 2))
    }
  }
  
  # store host niche data
  species <- rep(1:N, each=A)
  E <- rep(E, N)
  Pr.estab <- c(Pcol)
  niche.d <- data.frame(species, E, Pr.estab)
  niche.d <- niche.d[with(niche.d, order(species, E)),]
  
  # store parasite niche data
  parasite.species <- rep(1:P, each=N)
  pE <- rep(pE, P)
  pPr.estab <- c(pPcol)
  pniche.d <- data.frame(parasite.species, pE, pPr.estab)
  pniche.d <- pniche.d[with(pniche.d, order(parasite.species, pE)),]
  return(list(Pcol=Pcol, pPcol=pPcol, niche.d=niche.d, pniche.d=pniche.d))
}

static_init <- function(Y=5, modalO=2, z=0.1, 
                        pEmin=pEmin, pEmax=pEmax, P=P, 
                        pERmin=pERmin, pERmax=pERmax,
                        sig.p=sig.p){

  ## Define host community
  octaves <- 1:10
  s <- Y*exp(-(z*(octaves-modalO)^2)) 
  s <- round(s)
  plog2 <- 2^octaves
  abund <- rep(plog2, s)
  abund <- sort(abund, decreasing=T)
  rank <- rep(octaves, s)
  rank <- sort(rank, decreasing=T)
  N <- length(abund)
  A <- sum(abund)
  preston <- data.frame(s, octaves, plog2)
  
  # put into array
  hosts <- array(0, dim=c(A, N))
  for (i in 1:N){
    if (i == 1){
      end <- abund[i]
      hosts[1:end, i] <- 1
    } else {
      start <- end + 1
      end <- start + abund[i] - 1
      hosts[start:end, i] <- 1
    }
  }
  
  # plot host community
  #require(reshape2)
  #require(ggplot2)
  #x1 <- melt(hosts)
  #names(x1) = c("x", "y", "Occupancy")
  #x1$Occupancy <- factor(x1$Occupancy, labels=c("unoccupied", "occupied"))
  #ggplot(x1, aes(x=x, y=y, fill=Occupancy)) + 
  #  geom_tile(color="black") + 
  #  theme_classic() + 
  #  scale_fill_manual(values=c("lightgrey", "black")) + 
  #  xlab("Individual number") + 
  #  ylab("Species identity")

  ## Define parasite community
  pE <- runif(N, pEmin, pEmax) # within-host environmental conditions
  mu.p <- runif(P, pERmin, pERmax) # parasite optima
  sigma.p <- rep(sig.p, P) # runif(P, 1, 10) # parasite niche width
  
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
  
  # store parasite niche data
  parasite.species <- rep(1:P, each=N)
  pE <- rep(pE, P)
  pPr.estab <- c(pPcol)
  pniche.d <- data.frame(parasite.species, pE, pPr.estab)
  pniche.d <- pniche.d[with(pniche.d, order(parasite.species, pE)),]
  return(list(preston=preston, pPcol=pPcol, pniche.d=pniche.d, 
              abund=abund, rank=rank, N=N, A=A, state=hosts))
}

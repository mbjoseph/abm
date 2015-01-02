# symb_init() initializes the symbiont species traits w.r.t. environmental conditions

# Arguments:
## H = number of host species in regional pool
## sEmin = environmental condition minimum
## sEmax = environmental condition maximum
## S = number of symbionts
## sig.s = symbiont niche breadth

# Returns:
## sniche.d = data frame of symbiont niche data
## P.col = host X symbiont species array of colonization probabilities

symb_init <- function(H, S, sEmin, sEmax, sERmin, sERmax, sig.s){
  # setup ordered vector of environmental conditions
  Xe <- sort(runif(H, sEmin, sEmax))
  
  # symbiont optima
  mu.s <- runif(S, sERmin, sERmax)
  
  # symbiont niche width
  sigma.s <- rep(sig.s, S) 
  
  # calculate symbiont Z
  sZ <- rep(NA, S)
  for (i in 1:S){
    integrand <- function(Xe) {
      exp(-((Xe - mu.s[i]) ^ 2) / (2 * sigma.s[i] ^ 2))
    }
    res <- integrate(integrand, lower=sERmin, upper=sERmax)
    sZ[i] <- 1 / res$value
  }
  
  # symbiont probability of establishment
  Pcol <- array(dim=c(H, S))
  for (i in 1:S){
    for (j in 1:H){
      Pcol[j, i] <- sZ[i] * exp(-((Xe[j] - mu.s[i]) ^ 2) / (2 * sigma.s[i] ^ 2))
    }
  }
  
  # store symbiont niche data
  symbiont.species <- rep(1:S, each=H)
  host.condition <- rep(Xe, S)
  Pr.estab <- c(Pcol)
  sniche.d <- data.frame(symbiont.species, host.condition, Pr.estab)
  return(list(sniche.d = sniche.d, Pcol = Pcol))
}

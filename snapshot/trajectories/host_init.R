# host_init() initializes the host species traits w.r.t. environmental conditions

# Arguments:
## cells = number of cells in environment (max number of hosts)
## Emin = environmental condition minimum
## Emax = environmental condition maximum
## H = number of hosts
## sig.h = host niche breadth

# Returns:
## hniche.d = data frame of host niche data
## P.col = cell X host species array of colonization probabilities

host_init <- function(cells, Emin, Emax, H, sig.h){
  # setup ordered vector of environmental conditions
  Xe <- sort(runif(cells, Emin, Emax))
  
  # assume regional min & max = local min & max
  ERmin <- Emin 
  ERmax <- Emax
  
  # host optima
  mu.h <- runif(H, ERmin, ERmax)
  
  # host niche width
  sigma.h <- rep(sig.h, H) 
  
  # calculate host Z
  hZ <- rep(NA, H)
  for (i in 1:H){
    integrand <- function(Xe) {
      exp(-((Xe - mu.h[i]) ^ 2) / (2 * sigma.h[i] ^ 2))
    }
    res <- integrate(integrand, lower=ERmin, upper=ERmax)
    hZ[i] <- 1 / res$value
  }
  
  # host probability of establishment
  Pcol <- array(dim=c(cells, H))
  for (i in 1:H){
    for (j in 1:cells){
      Pcol[j, i] <- hZ[i] * exp(-((Xe[j] - mu.h[i]) ^ 2) / (2 * sigma.h[i] ^ 2))
    }
  }
  
  # store host niche data
  host.species <- rep(1:H, each=cells)
  environmental.condition <- rep(Xe, H)
  Pr.estab <- c(Pcol)
  hniche.d <- data.frame(host.species, environmental.condition, Pr.estab)
  return(list(hniche.d = hniche.d, Pcol = Pcol))
}

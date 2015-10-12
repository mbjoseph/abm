# Helper functions for gillespie implementation of ABM

# Rate function generates rates of all events given current state & params
# state: 2 X ncells array containing 0's if cell is unoccupied, 
#        and a number indicating which host or symbiont is present
# params: a list of the following parameters
# r = host birth rate
# d = background host death rate
# c = host colonization rate
# beta_d = effect of infection on death rate
#  (if beta_d = 0, infection has no effect, 
#   beta_d < 0 reduction in death rate, 
#   beta_d > 0 increase in death rate)
# phi = contact rate
# a_pen = among-species contact rate penalty (probability)
# rs = rate of symbiont rain
# cells = number of lanscape cells
# H = number of hosts present
# nS = number of symbionts present
# gamma = host recovery rate
# hosts = result from host_init
# symbionts = result from symb_init
ratefun <- function(state, params){
  stopifnot(params$beta_d > -1)
  res <- matrix(0, nrow=4, ncol=ncol(state))
  x_occ <- state[1, ] > 0 # host slots occupied
  s_occ <- state[2, ] > 0 # symbiont slots occupied
  res[1, x_occ] <- params$r                                   # birth
  res[2, x_occ & !s_occ] <- params$d                          # death uninfect
  res[2, x_occ & s_occ] <- params$d * (1 + params$beta_d)     # death infect
  res[3, !x_occ] <- params$c * params$H                       # colon
  res[4, s_occ & s_occ] <- params$gamma                       # recov
  
  # generate table of counts for each host species that is uninfected
  n_susc <- table(state[1, x_occ & !s_occ])
  nsv <- rep(0, params$H)
  nsv[as.numeric(names(n_susc))] <- c(n_susc)
  # put into matrix form
  N_H <- matrix(nsv, nrow=params$H, ncol=params$nS)
  
  # define symbiont colonization rates
  if (!any(x_occ & !s_occ)){
    sCmat <- matrix(0, nrow=params$H, ncol=params$nS)
  } else {
    sCmat <- params$rs * N_H * params$symbionts$Pcol
  }
  
  # define contact rates
  if (!any(x_occ) | !any(s_occ)){ # need susceptibles and infectives for transm.
    Tmat <- matrix(0, nrow=params$H, ncol=params$nS)
  } else {
    # generate table of counts for number of hosts infected with each symbiont
    n_infct <- table(state[2, ])[-1] # -1 because it returns number of 0s
    niv <- rep(0, params$nS)
    niv[as.numeric(names(n_infct))] <- c(n_infct)
    
    # put the counts into matrices
    N_I <- matrix(niv, nrow=params$H, ncol=params$nS, byrow=TRUE)
    
    # generate transmision matrix
    Tmat <- params$phi * N_H * N_I * params$symbionts$Pcol
  }
  
  # return unrolled transmission matrix
  c(t(res), c(sCmat), c(Tmat))
}

# agent-based model step function
ABMstep <- function(state, action, cell, params, event){
  stopifnot(length(cell) == 1)
  stopifnot(!is.na(cell))
  counter <- 0
  same_species <- NA
  if (action == "birth"){
    # find cell for host to disperse to
    cells <- dim(state)[2]
    disp.cell <- sample.int(cells, size=1)
    if (state[1, disp.cell] == 0){
      success <- rbinom(1, 1, prob = params$hosts$Pcol[disp.cell, state[1, cell]])
      if (success == 1){
        state[1, disp.cell] <- state[1, cell]
        counter <- counter + 1
      }
    }
  }
  if (action == "cntct"){
    # determine which symbiont species we have
    infecting_symbiont_sp <- params$symb_index[event - params$non_cntct_indices]
    state[2, cell] <- infecting_symbiont_sp
    counter <- counter + 1
  }
  if (action == "rains"){
    stopifnot(state[2, cell] == 0)
    # which symbiont species
    colonizing_symbiont_sp <- params$symb_index[event - 4*params$cells]
    state[2, cell] <- colonizing_symbiont_sp
    counter <- counter + 1
  }
  if (action == "death"){
    state[, cell] <- 0 # host and symbiont dies
    counter <- counter + 1
  }
  if (action == "recov"){
    state[2, cell] <- 0 # symbiont dies
    counter <- counter + 1
  }
  if (action == "colon"){
    if (state[1, cell] == 0){
      spnum <- sample.int(params$H, 1)
      success <- rbinom(1, 1, prob = params$hosts$Pcol[cell, spnum])
      if (success) state[1, cell] <- spnum
      #state[1, cell] <- ifelse(success == 1, spnum, 0)
      counter <- counter + success
    }
  }
  return(list(state=state, counter = counter))
}

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

host_init <- function(cells, H, pcol){
  host.species <- rep(1:H, each=cells)
  environmental.condition <- rep(0, H)
  Pcol <- array(pcol, dim=c(cells, H))  
 return(list(Pcol = Pcol))
}

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
  if (length(sig.s) == 1){
    sigma.s <- rep(sig.s, S)    
  } else {
    stopifnot(length(sig.s) == S)
    sigma.s <- sig.s
  }
  
  # calculate symbiont Z
  sZ <- rep(NA, S)
  for (i in 1:S){
    integrand <- function(Xe) {
      exp(-((Xe - mu.s[i]) ^ 2) / (2 * sigma.s[i] ^ 2))
    }
    res <- integrate(integrand, lower=sERmin, upper=sERmax, subdivisions=1000L)
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

# test 
check <- FALSE
if (check) { 
  out <- symb_init(H=100, S=100, sEmin=-100, sEmax=100, sERmin=-100, sERmax=100, sig.s=1)
#  out$Pcol
  out$sniche.d$Symbiont <- as.factor(out$sniche.d$symbiont.species)
  
  ggplot(out$sniche.d, aes(x=host.condition, y=Pr.estab)) + 
    geom_line(aes(group=Symbiont), size=1) + 
    theme_classic() + 
    geom_point(aes(x=host.condition, y=-.01), pch="|", size=4) + 
    xlab("Within-host condition") + 
    ylab("Probability of establishment") + 
    theme(legend.position="bottom") + 
    theme(legend.background = element_rect(color="black", size=.5)) 
}

# removes duplicate sequences of repeated integers from a column
sub_unique <- function(d, col){
  d <- as.data.frame(d)
  numcol <- which(names(d) == col)
  keep <- cumsum(rle(d[, col])$length)
  d[keep, ]
}

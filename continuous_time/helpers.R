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
  X <- state[1, ] # host state
  S <- state[2, ] # symbiont state
  stopifnot(params$beta_d > -1)
  with(c(as.list(params)),
       list(
         birth = ifelse(X > 0, r, 0),
         death = ifelse(X > 0, d * (1 + beta_d * (S > 0)), 0),
         colon = ifelse(X == 0, c * H, 0), 
         rains = ifelse((X > 0 & S == 0), rs * nS, 0), 
         recov = ifelse((X > 0 & S > 0), gamma, 0),
         cntct = ifelse(sum(X) == 0, 0, 
                        ifelse(mode == "density", 
                          phi * sum(X > 0 & S == 0) * sum(X>0 & S > 0), # dd
                          phi * sum(X > 0 & S == 0) * sum(X>0 & S > 0) / sum(X>0) # fd
                          )
                        )
         # could incorporate a_pen here, to avoid unrealized contacts
       )
       )
}

# helper function to update symbiont state through transmission
transmit <- function(state, cell1, cell2, params){
  X <- state[1, ]
  S <- state[2, ]
  pr.estab <- params$symbionts$Pcol[X[cell2], S[cell1]] # is a host x symb array
  success <- rbinom(1, 1, pr.estab)
  if (success){
    S[cell2] <- S[cell1]
  }
  return(list(S=S, success=success))
}


# agent-based model step function
ABMstep <- function(state, action, cell, params){
  counter <- 0
  same_species <- NA
  if (action == "birth"){
    # find cell for host to disperse to
    cells <- dim(state)[2]
    disp.cell <- sample(1:cells, size=1)
    if (state[1, disp.cell] == 0){
      success <- rbinom(1, 1, prob = params$hosts$Pcol[disp.cell, state[1, cell]])
      if (success == 1){
        state[1, disp.cell] <- state[1, cell]
        counter <- counter + 1
      }
    }
  }
  if (action == "cntct"){
    # find S and I hosts at random (proportional to abundance)
    Iindivs <- which(state[1, ] > 0 & state[2, ] > 0)
    if (length(Iindivs) == 1){
      cell <- Iindivs
    } else {
      cell <- sample(Iindivs, size=1)
    }
    
    Sindivs <- which(state[1, ] > 0 & state[2, ] == 0)
    if (length(Sindivs) == 1){
      ind_contacted <- Sindivs
    } else {
      ind_contacted <- sample(Sindivs, size=1)
    }
    sp1 <- state[1, cell]
    sp2 <- state[1, ind_contacted]
    same_species <- as.numeric(sp1 == sp2)
    contact_realized <- rbinom(1, 1, 
                               same_species + 
                                 (1 - same_species) * params$a_pen)
    # if contact realized, check for transmission
    if (contact_realized){
      counter <- counter + 1
      # call to transmission function
      t_out <- transmit(state, cell, ind_contacted, params)
      if (t_out$success){
        state[2, ] <- t_out$S
      }
    }
  }
  if (action == "rains"){
    if(state[2, cell] == 0){
      # eventually call to symbiont rain function
      attempting <- sample(1:params$nS, size=1)
      success <- rbinom(1, 1, params$symbionts$Pcol[state[1, cell], attempting])
      if (success){
        state[2, cell] <- attempting
        counter <- counter + success
      }
    }
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
      spnum <- sample(1:params$H, 1)
      success <- rbinom(1, 1, prob = params$hosts$Pcol[cell, spnum])
      state[1, cell] <- ifelse(success == 1, spnum, 0)
      counter <- counter + success
    }
  }
  return(list(state=state, counter = counter, same_species = same_species))
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
  sigma.s <- rep(sig.s, S) 
  
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


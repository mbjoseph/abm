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
  with(c(as.list(params)),
       list(
         birth = ifelse(X > 0, r, 0),
         death = ifelse(X > 0, d * (1 + beta_d * (S > 0)), 0),
         colon = ifelse(X == 0, c * H, 0), 
         cntct = ifelse(S > 0, phi * (sum(X > 0 & S == 0) - 1), 0), 
         # ^ minus one, can't contact self
         # should incorporate a_pen here, to avoid unrealized contacts
         rains = ifelse((X > 0 & S == 0), rs * nS, 0), 
         recov = ifelse((X > 0 & S > 0), gamma, 0)
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
ABMstep <- function(state, action, cell, params, transmitted){
  counter <- 0
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
    # find individual for host to contact
    Sindivs <- which(state[1, ] > 0 & state[2, ] == 0)
    if (length(Sindivs) == 1){
      ind_contacted <- Sindivs
    } else {
      ind_contacted <- sample(Sindivs, size=1)
    }
    sp1 <- state[1, cell]
    sp2 <- state[1, ind_contacted]
    same_species <- sp1 == sp2
    contact_realized <- rbinom(1, 1, ifelse(same_species, 1, params$a_pen))
    # if contact realized, check for transmission
    if (contact_realized){
      counter <- counter + 1
      # call to transmission function
      t_out <- transmit(state, cell, ind_contacted, params)
      if (t_out$success){
        state[2, ] <- t_out$S
        transmitted[sp1, sp2, state[2, cell]] <- transmitted[sp1, sp2, state[2, cell]] + 1
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
  return(list(state=state, counter = counter, 
              transmitted=transmitted))
}

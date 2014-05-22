# Helper functions for gillespie implementation of ABM

# Rate function generates rates of all events given current state & params
ratefun <- function(X, S, params){
  with(c(as.list(params)),
       list(
         birth = ifelse(X > 0, r, 0),
         death = ifelse(X > 0, d, 0),
         colon = ifelse(X == 0, c * H, 0), 
         cntct = ifelse(S > 0, phi * (sum(X > 0 & S == 0)), 0), # minus one, can't contact self
         rains = ifelse((X > 0 & S == 0), rs * nS, 0)
       )
  )
}

# helper function to update symbiont state through transmission
transmit <- function(X, S, cell1, cell2, params){
  pr.estab <- params$symbionts$Pcol[X[cell2], S[cell1]] # is a host x symb array
  success <- rbinom(1, 1, pr.estab)
  if (success){
    S[cell2] <- S[cell1]
  }
  return(S)
}

# agent-based model step function
ABMstep <- function(state, action, cell, H, S, params){
  counter <- 0
  if (action == "birth"){
    # find cell to disperse to
    disp.cell <- sample(1:cells, size=1)
    if (state[disp.cell] == 0){
      success <- rbinom(1, 1, prob = params$hosts$Pcol[disp.cell, state[cell]])
      state[disp.cell] <- state[disp.cell] + ifelse(success == 1, state[cell], 0)
      counter <- counter + success
    }
  }
  if (action == "cntct"){
    # find individual to contact
    Sindivs <- which(X > 0 & S == 0)
    if (length(Sindivs) == 1){
      ind_contacted <- Sindivs
    } else {
      ind_contacted <- sample(Sindivs, size=1)
    }
    same_species <- ifelse(state[cell] == state[ind_contacted], 1, 0)
    contact_realized <- rbinom(1, 1, 
                               prob = ifelse(same_species, phi, phi * a_pen))
    # if contact realized, check for transmission
    if (contact_realized){
      # call to transmission function
      S <- transmit(state, S, cell, ind_contacted, params)
    }
  }
  if (action == "rains"){
    if(S[cell] == 0){
      # eventually call to symbiont rain function
      attempting <- sample(1:params$nS, size=1)
      success <- rbinom(1, 1, prob=params$symbionts$Pcol[state[cell], attempting])
      if (success){
        S[cell] <- attempting
      }
    }
  }
  if (action == "death"){
    state[cell] <- 0
    S[cell] <- 0
    counter <- counter + 1
  }
  if (action == "colon"){
    if (state[cell] == 0){
      spnum <- sample(1:H, 1)
      success <- rbinom(1, 1, prob = params$hosts$Pcol[cell, spnum])
      state[cell] <- ifelse(success == 1, spnum, 0)
      counter <- counter + success
    }
  }
  return(list(state=state, S = S, counter = counter))
}

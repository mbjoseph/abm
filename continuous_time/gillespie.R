# plants on grid, reproduce with rate r, die rate d
rm(list=ls())
r <- 2.0 # reproductive rate
d <- 1 # death rate
c <- .1 # colonization rate for empty cells
cells <- 100 # number of habitat patches

maxt <- 50
t <- 0

X <- rep(0, cells)

# set initial condition
# X[sample(1:cells, size=1)] <- 1

t.int <- 1
n.ind <- sum(X)
tau <- NULL

pars <- c(r = r, d = d, c = c, cells = cells)

# Rate function
ratefun <- function(X, params){
  with(c(as.list(params)),
       list(
         birth = ifelse(X == 1, r, 0),
         death = ifelse(X == 1, d, 0),
         colon = rep(c, length(X)))
       )
}

# agent-based model step function
ABMstep <- function(state, action, cell){
  if (action == "birth"){
    # find cell to disperse to
    disp.cell <- sample(1:cells, size=1)
    if (state[disp.cell] == 0){
      state[disp.cell] <- state[disp.cell] + 1
    }
  }
  if (action == "death"){
    state[cell] <- 0
  }
  if (action == "colon"){
    # ensure that this function is not called if cell occupied
    state[cell] <- 1
  }
  return(state)
}


ntrans <- 3 * cells
transctr <- rep(0, ntrans)

while(t[t.int] < maxt){
  
  ratelist <- ratefun(X, pars)
  rates <- unlist(ratelist)
  rnames <- rep(names(ratelist), each=cells)
  ntrans <- length(rates)
  tot.rates <- sum(rates)
  tau[t.int] <- rexp(1, tot.rates)
  event <- sample((1:ntrans)[rates > 0], 
                  size=1, 
                  prob = rates[rates > 0])
  transctr[event] <- transctr[event] + 1
  f_name <- names(rates[event])
  event_type <- rnames[event]
  cell_num <- event - cells * (which(names(ratelist) == event_type) - 1)
  # update time
  t.int <- t.int + 1
  t[t.int] <- t[t.int - 1] + tau[t.int-1]
  
  # carry out action on chosen cell
  if (event_type == "colon" & X[cell_num] == 1){
    transctr[event] <- transctr[event] - 1 # uncount that event
  } else {
    X <- ABMstep(X, event_type, cell_num)
  }
  n.ind[t.int] <- sum(X)
}

par(mfrow=c(1, 2))
plot(t, n.ind, type="l")
plot(density(log(tau)))

relist(transctr, ratelist)

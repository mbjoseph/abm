# plants on grid, reproduce with rate r, die rate d
rm(list=ls())
r <- 1 # reproductive rate
d <- .01 # death rate
c <- .01 # colonization rate for empty cells
cells <- 50 # number of habitat patches

# establish environmental conditions for each cell
Emin <- -10
Emax <- 10
Xe <- sort(runif(cells, Emin, Emax))

H <- 3 # number of host species
sig.h <- runif(H, 1, 10) # host niche breadth
ERmin <- Emin
ERmax <- Emax
mu.h <- runif(H, ERmin, ERmax) # host optima
sigma.h <- rep(sig.h, H) # runif(P, 1, 10) # parasite niche width

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

par(mfrow=c(1, 3))
plot(Xe, Pcol[, 1], xlim=c(ERmin, ERmax), 
     ylim=c(0, max(Pcol)), 
     type="l")
for (i in 2:H){
  lines(Xe, Pcol[, i], col=i)
}
points(x=Xe, y=rep(0, length(Xe)))
# store parasite niche data
host.species <- rep(1:H, each=cells)
environmental.condition <- rep(Xe, H)
Pr.estab <- c(Pcol)
hniche.d <- data.frame(host.species, environmental.condition, Pr.estab)

maxt <- 2000
t <- 0

X <- rep(0, cells)

# set initial condition
# X[sample(1:cells, size=1)] <- 1

t.int <- 1
n.ind <- array(0, dim=c(1, H))
tau <- NULL

pars <- c(r = r, d = d, c = c, cells = cells, H=H)

# Rate function
ratefun <- function(X, params){
  with(c(as.list(params)),
       list(
         birth = ifelse(X > 0, r, 0),
         death = ifelse(X > 0, d, 0),
         colon = rep(c * H, cells))
       )
}

# agent-based model step function
ABMstep <- function(state, action, cell, H, Pcol){
  counter <- 0
  if (action == "birth"){
    # find cell to disperse to
    disp.cell <- sample(1:cells, size=1)
    if (state[disp.cell] == 0){
      success <- rbinom(1, 1, prob = Pcol[disp.cell, state[cell]])
      state[disp.cell] <- state[disp.cell] + ifelse(success == 1, state[cell], 0)
      counter <- counter + success
    }
  }
  if (action == "death"){
    state[cell] <- 0
    counter <- counter + 1
  }
  if (action == "colon"){
    if (state[cell] == 0){
      spnum <- sample(1:H, 1)
      success <- rbinom(1, 1, prob = Pcol[cell, spnum])
      state[cell] <- ifelse(success == 1, spnum, 0)
      counter <- counter + success
    }
  }
  return(list(state=state, counter = counter))
}

ntrans <- 3 * cells
transctr <- rep(0, ntrans) # transition counter

while(t[t.int] < maxt){
  ratelist <- ratefun(X, pars)
  rates <- unlist(ratelist)
  rnames <- rep(names(ratelist), times=c(cells, cells, cells * H))
  ntrans <- length(rates)
  tot.rates <- sum(rates)
  tau[t.int] <- rexp(1, tot.rates)
  event <- sample((1:ntrans)[rates > 0], 
                  size=1, 
                  prob = rates[rates > 0])
  f_name <- names(rates[event])
  event_type <- rnames[event]
  cell_num <- event - cells * (which(names(ratelist) == event_type) - 1)
  
  # carry out action on chosen cell
  nextstep <- ABMstep(X, event_type, cell_num, H, Pcol)
  X <- nextstep$state
  transctr[event] <- transctr[event] + nextstep$counter

  # update time
  t.int <- t.int + 1
  t[t.int] <- t[t.int - 1] + tau[t.int-1]
  n.ind <- rbind(n.ind, tabulate(X, nbins=H))
}

plot(t, n.ind[, 1], type="l", ylim=c(0, max(n.ind)))
for (i in 2:H){
  lines(t, n.ind[, i], col=i)
}
plot(density(log(tau)))

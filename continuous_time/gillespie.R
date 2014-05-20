# plants on grid, reproduce with rate r, die rate d
rm(list=ls())
r <- 1 # reproductive rate
d <- .05 # death rate
c <- .001 # colonization rate for empty cells
phi <- .01 # interspecific contact rates
a_pen <- .1 # among-species adjustment for contact rates
rs <- .001 # rate at which symbionts rain
cells <- 100 # number of habitat patches

# establish environmental conditions for each cell
Emin <- -5
Emax <- 5
Xe <- sort(runif(cells, Emin, Emax))

H <- 2 # number of host species
sig.h <- 1.5 #runif(H, .5, 2) # host niche breadth
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

# store host niche data
host.species <- rep(1:H, each=cells)
environmental.condition <- rep(Xe, H)
Pr.estab <- c(Pcol)
hniche.d <- data.frame(host.species, environmental.condition, Pr.estab)

X <- rep(0, cells) # hosts
S <- rep(0, cells) # symbionts

# set initial condition
# X[sample(1:cells, size=1)] <- 1

t.int <- 1
n.ind <- array(0, dim=c(1, H))
s.ind <- array(0, dim=c(1, 1))
tau <- NULL

pars <- c(r = r, d = d, c = c,
          phi = phi, a_pen = a_pen, rs = rs,
          cells = cells, H=H)

# Rate function
ratefun <- function(X, S, params){
  with(c(as.list(params)),
       list(
         birth = ifelse(X > 0, r, 0),
         death = ifelse(X > 0, d, 0),
         colon = rep(c * H, cells), 
         cntct = ifelse(S > 0, phi * (sum(X > 0 & S == 0)), 0), # minus one, can't contact self
         rains = ifelse(X > 0, rs, 0)
         )
       )
}

transmit <- function(X, S, cell1, cell2){
  pr.estab <- 1
  success <- rbinom(1, 1, pr.estab)
  if (success){
    S[cell2] <- S[cell1]
  }
  return(S)
}

# agent-based model step function
ABMstep <- function(state, action, cell, H, Pcol, S){
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
      S <- transmit(state, S, cell, ind_contacted)
    }
  }
  if (action == "rains"){
    if(S[cell] == 0){
      # eventually call to symbiont rain function
      S[cell] <- rbinom(1, 1, prob = 1)
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
      success <- rbinom(1, 1, prob = Pcol[cell, spnum])
      state[cell] <- ifelse(success == 1, spnum, 0)
      counter <- counter + success
    }
  }
  return(list(state=state, S = S, counter = counter))
}

ntrans <- 3 * cells
transctr <- rep(0, ntrans) # transition counter

maxt <- 5000
t <- 0
&/test&
pb <- txtProgressBar(min=0, max=maxt, style=3)
while(t[t.int] < maxt){
  ratelist <- ratefun(X, S, pars)
  rates <- unlist(ratelist)
  rnames <- rep(names(ratelist), times=c(cells, cells, 
                                         cells, 
                                         cells, cells))
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
  nextstep <- ABMstep(X, event_type, cell_num, H, Pcol, S)
  X <- nextstep$state
  S <- nextstep$S
  transctr[event] <- transctr[event] + nextstep$counter

  # update time
  t.int <- t.int + 1
  t[t.int] <- t[t.int - 1] + tau[t.int-1]
  setTxtProgressBar(pb, t[t.int])
  n.ind <- rbind(n.ind, tabulate(X, nbins=H))
  s.ind <- rbind(s.ind, tabulate(S, nbins=1))
}

plot(t, n.ind[, 1], type="l", ylim=c(0, max(c(n.ind, s.ind))))
for (i in 2:H){
  lines(t, n.ind[, i], col=i)
}
lines(t, s.ind, col="grey")
plot(log(tau), type="l")

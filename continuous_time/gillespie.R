# plants on grid, reproduce with rate r, die rate d
rm(list=ls())
setwd("/home/max/Documents/abm/continuous_time")
library(ggplot2)
library(gridExtra)
source("host_init.R")
source("symb_init.R")
source("helpers.R")

## Define most important parameters
r <- 1 # reproductive rate
d <- .05 # death rate
c <- .001 # colonization rate for empty cells
phi <- .2 # interspecific contact rates
a_pen <- .1 # among-species adjustment for contact rates
rs <- .001 # rate at which symbionts rain
cells <- 100 # number of habitat patches

# host community parameters
Emin <- -5
Emax <- 5
H <- 50
sig.h <- 50

# symbiont community parameters
sEmin <- -5
sEmax <- 5
nS <- 50
sig.s <- 50

# initialize state arrays, time objects, and abundance counters
X <- rep(0, cells) # hosts
S <- rep(0, cells) # symbionts
n.ind <- array(0, dim=c(1, H))
s.ind <- array(0, dim=c(1, nS))
rnames <- rep(c("birth", "death", "colon", "cntct", "rains"), 
              each=cells)
ntrans <- length(rnames)
transctr <- rep(0, ntrans) # transition counter

# initialize host community and visualize host niches
hosts <- host_init(cells, Emin, Emax, H, sig.h)

hniche <- ggplot(hosts$hniche.d) + 
  geom_line(aes(x=environmental.condition, y=Pr.estab, 
                col=factor(host.species), group=host.species), 
            size=2) + 
  geom_rug(aes(x=environmental.condition)) +
  theme_bw()

# initialize symbionts and view niches
symbionts <- symb_init(H, S=nS, sEmin, sEmax, sig.s)

sniche <- ggplot(symbionts$sniche.d) +
  geom_line(aes(x=host.condition, y=Pr.estab, 
                col=factor(symbiont.species), group=symbiont.species), 
            size=2) + 
  geom_rug(aes(x=host.condition)) +
  theme_bw()

grid.arrange(hniche, sniche, ncol=2)

# store parameters for later extraction
pars <- list(r = r, d = d, c = c,
             phi = phi, a_pen = a_pen, rs = rs,
             cells = cells, H=H, nS=nS, 
             symbionts=symbionts, hosts=hosts)

maxt <- 500
t <- 0
t.int <- 1
tau <- NULL


pb <- txtProgressBar(min=0, max=maxt, style=3)
while(t[t.int] < maxt){
  
  ratelist <- ratefun(X, S, pars)
  rates <- unlist(ratelist)
  tot.rates <- sum(rates)
  tau[t.int] <- rexp(1, tot.rates)
  event <- sample((1:ntrans)[rates > 0], 
                  size=1, 
                  prob = rates[rates > 0])
  f_name <- names(rates[event])
  event_type <- rnames[event]
  cell_num <- event - cells * (which(names(ratelist) == event_type) - 1)
  
  # carry out action on chosen cell
  nextstep <- ABMstep(X, event_type, cell_num, H, S, params=pars)
  X <- nextstep$state
  S <- nextstep$S
  transctr[event] <- transctr[event] + nextstep$counter

  # update time
  t.int <- t.int + 1
  t[t.int] <- t[t.int - 1] + tau[t.int-1]
  setTxtProgressBar(pb, t[t.int])
  n.ind <- rbind(n.ind, tabulate(X, nbins=H))
  s.ind <- rbind(s.ind, tabulate(S, nbins=nS))
}

par(mfrow=c(1, 3))
plot(t, n.ind[, 1], type="l", ylim=c(0, max(c(n.ind, s.ind))))
for (i in 2:H){
  lines(t, n.ind[, i], col=i)
}

plot(t, s.ind[, 1], type="l", ylim=c(0, max(c(n.ind, s.ind))))
for (i in 2:nS){
  lines(t, s.ind[, i], col=i)
}

plot(log(tau), type="l")

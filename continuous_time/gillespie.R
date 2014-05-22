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
d <- 0.01 # death rate
c <- .001 # colonization rate for empty cells
phi <- .75 # interspecific contact rates
a_pen <- 1 # among-species adjustment for contact rates
rs <- 0.01 # rate at which symbionts rain
cells <- 150 # number of habitat patches

# host community parameters
Emin <- -20 # min environmental condition
Emax <- 20 # max environmental conditionb
H <- 20 # regional host richness
sig.h <- 10 # host niche width

# symbiont community parameters
sEmin <- -20 # min w/in host condition
sEmax <- 20 # max w/in host condition
nS <- 20 # regional symbiont richness
sig.s <- 10 # symbiont niche width

# initialize state arrays, time objects, and abundance counters
state <- array(0, dim=c(2, cells))
n.ind <- array(0, dim=c(1, H))
s.ind <- array(0, dim=c(1, nS))
rnames <- rep(c("birth", "death", "colon", "cntct", "rains"), 
              each=cells)
ntrans <- length(rnames)
transctr <- rep(0, ntrans) # transition counter

transmission_events <- array(0, dim=c(H, H, nS))

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

maxt <- 300
t <- 0
t.int <- 1
tau <- NULL

pb <- txtProgressBar(min=0, max=maxt, style=3)
while(t[t.int] < maxt){
  ratelist <- ratefun(state, pars)
  rates <- unlist(ratelist)
  tot.rates <- sum(rates)
  tau[t.int] <- rexp(1, tot.rates)
  event <- sample(1:ntrans, size=1, prob=rates)
  f_name <- names(rates[event])
  event_type <- rnames[event]
  cell_num <- event - cells * (which(names(ratelist) == event_type) - 1)
  
  # carry out action on chosen cell
  nextstep <- ABMstep(state, event_type, cell_num, 
                      params=pars, transmission_events)
  state <- nextstep$state
  transctr[event] <- transctr[event] + nextstep$counter
  transmission_events <- nextstep$transmission_events

  # update time
  t.int <- t.int + 1
  t[t.int] <- t[t.int - 1] + tau[t.int-1]
  setTxtProgressBar(pb, t[t.int])
  n.ind <- rbind(n.ind, tabulate(state[1, ], nbins=H))
  s.ind <- rbind(s.ind, tabulate(state[2, ], nbins=nS))
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

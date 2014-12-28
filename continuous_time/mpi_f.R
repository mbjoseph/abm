mpi_f <- function(iter=1, nER=1, maxt=10, H=10, nS=10, 
                  a_pen=0.1, sig.s=10, rs=.1, gamma=0.1, cells=100){
  source("symb_init.R")
  source("host_init.R")
  source("helpers.R")
  # generate environmental ranges
  ER <- sample(seq(1, 100, by=.01), size = nER, replace=FALSE)
  
  # regional symbiont niche minimum and maximum
  sERmin <- -50
  sERmax <- 50
  
  # generate values for Emin and Emax, nrow=nER
  Earray <- array(dim=c(nER, 2))
  for (i in 1:nER){
    # find range of potential means
    rmin <- sERmin + ER[i] / 2
    rmax <- sERmax - ER[i] / 2
    ERmidpoint <- runif(1, rmin, rmax)
    Earray[i, ] <- c(ERmidpoint - ER[i] / 2, ERmidpoint + ER[i] / 2)
  }
  
  # initialize parameters to stay constant across runs
  # (dynamic host community)
  r <- 1
  d <- 0.1
  beta_d <- -1
  c <- 0.1
  # H <- 10 # host richness (max)
  
  # symbiont params
  phi <- 5 # contact rate
  # a_pen <- .1 # interspecific penalty
  # rs <- .01 # rain
  # nS <- 50 # regional symbiont richness
  # sig.s <- 10 # niche width
  # gamma <- .01 # recovery rate
  
  trans_bar <- array(dim=c(H, H, nS, nER, iter))
  rich_bar <- array(dim=c(nER, iter))
  
  for (i in 1:nER){
    for (j in 1:iter){
      # initialize state arrays, time objects, and abundance counters
      state <- array(0, dim=c(2, cells))
      n.ind <- array(0, dim=c(1, H))
      s.ind <- array(0, dim=c(1, nS))
      rnames <- rep(c("birth", "death", "colon", "cntct", "rains", "recov"), 
                    each=cells)
      ntrans <- length(rnames)
      transctr <- rep(0, ntrans) # transition counter
      transmission_events <- array(0, dim=c(H, H, nS))
      
      # initialize host community
      # state[1, ] <- sample(1:H, size=cells, replace=T)
      hosts <- host_init(cells, H, pcol=.1)
      
      # initialize symbionts and view niches
      symbionts <- symb_init(H, S=nS, 
                             sEmin = Earray[i, 1], sEmax = Earray[i, 2], 
                             sERmin, sERmax, sig.s)
      
      # store parameters for later extraction
      pars <- list(r = r, d = d, c = c, beta_d = beta_d,
                   phi = phi, a_pen = a_pen, rs = rs,
                   cells = cells, H=H, nS=nS,
                   gamma= gamma, 
                   hosts=hosts,
                   symbionts=symbionts)
      
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
      
      rich <- rep(NA, dim(s.ind)[1])
      for (k in 1:dim(s.ind)[1]){
        rich[k] <- sum(s.ind[k, ] > 0)
      }
      
      trans_bar[, , , i, j] <- transmission_events
      rich_bar[i, j] <- mean(rich)
    }  
  }
  return(list(trans_bar = trans_bar, 
              rich_bar = rich_bar, 
              ER = ER, 
              Earray = Earray, 
              n.ind = n.ind, 
              s.ind = s.ind,
              t=t)
  )
}

testout <- mpi_f(nER=1, iter=1, maxt=5, nS=30, sig.s=3, H=10, 
                 a_pen = .5, cells=400, gamma=0, rs=0.1)
nsteps <- dim(testout$t)

par(mfrow=c(1, 2))
plot(x=testout$t, y=testout$s.ind[, 1], type="l", ylim=c(0, max(testout$s.ind)), 
     xlab="Time", ylab="Number of infected hosts")
for (k in 1:dim(testout$s.ind)[2]){
  lines(x=testout$t, y=testout$s.ind[, k], col=k + 2)
}

plot(x=testout$t, y=testout$n.ind[, 1], type="l", 
     ylim=c(0, max(testout$n.ind)), 
     xlab="Time", 
     ylab="Number of hosts")
for (k in 1:dim(testout$n.ind)[2]){
  lines(x=testout$t, y=testout$n.ind[, k], col=k+2)
}
par(mfrow=c(1, 1))

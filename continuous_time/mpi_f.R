mpi_f <- function(iter=1, nER=1, maxt=1, H=10, nS=10, 
                  a_pen=1, sig.s=10, rs=.1, gamma=0.1, 
                  cells=100, r=2, d=.1, 
                  beta_d_min=0, beta_d_max=0,
                  c=.1, phi=5, mode="dens"){
  stopifnot(mode == "dens" | mode == "freq")

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
  
  trans_bar <- array(dim=c(H, H, nS, nER, iter))
  rich_bar <- array(dim=c(nER, iter))
  stored_pars <- c()
  timesteps <- c()
  wasted_contacts <- c()
  win_transmission <- c()
  among_transmission <- c()
  
  for (i in 1:nER){
    for (j in 1:iter){
      # initialize state arrays, time objects, and abundance counters
      state <- array(0, dim=c(2, cells))
      n.ind <- array(0, dim=c(1, H))
      s.ind <- array(0, dim=c(1, nS))
      rnames <- c(rep(c("birth", "death", "colon", "rains", "recov"), 
                    each=cells), "cntct")
      ntrans <- length(rnames)
      transctr <- rep(0, ntrans) # transition counter
      n_contacts <- 0
      successful_cntcts <- 0
      win_t <- 0
      amng_t <- 0

      # initialize host community
      hosts <- host_init(cells, H, pcol=1)
      
      # initialize symbionts and view niches
      symbionts <- symb_init(H, S=nS, 
                             sEmin = Earray[i, 1], sEmax = Earray[i, 2], 
                             sERmin, sERmax, sig.s)
      
      # generate effects of infection on hosts
      beta_d <- runif(1, beta_d_min, beta_d_max)
      
      # store parameters for later extraction
      pars <- list(r = r, d = d, c = c, beta_d = beta_d,
                   phi = phi, a_pen = a_pen, rs = rs,
                   sig.s = sig.s,
                   cells = cells, H=H, nS=nS,
                   gamma= gamma, 
                   hosts=hosts,
                   symbionts=symbionts, 
                   mode=mode)
      
      t <- 0
      t.out <- 0
      t.int <- 1
      tau <- NULL
      
      pb <- txtProgressBar(min=0, max=maxt, style=3)
      #while(t[t.int] < maxt){
      while(t.int < maxt){
        ratelist <- ratefun(state, pars)
        rates <- unlist(ratelist)
        tot.rates <- sum(rates)
        stopifnot(!is.na(tot.rates))
        stopifnot(all(rates >= 0 ))
        tau[t.int] <- rexp(1, tot.rates)
        event <- sample(1:ntrans, size=1, prob=rates)
        f_name <- names(rates[event])
        event_type <- rnames[event]
        if (event_type=="cntct"){
          cell_num <- NA
          n_contacts <- n_contacts + 1
        } else {
          cell_num <- event - cells * (which(names(ratelist) == event_type) - 1)
        }
        
        # carry out action on chosen cell
        nextstep <- ABMstep(state, event_type, cell_num, 
                            params=pars)
        
        # update time
        t.int <- t.int + 1
        t[t.int] <- t[t.int - 1] + tau[t.int-1]
        
        state_change <- !all(state == nextstep$state)
        
        if (state_change){
          if (event_type == "cntct") {
            # count number of successful contacts
            successful_cntcts <- successful_cntcts + 1
            # count among vs. within species transmission
            if (nextstep$same_species){
              win_t <- win_t + 1
            } else {
              amng_t <- amng_t + 1
            }
          }
          state <- nextstep$state
          n.ind <- rbind(n.ind, tabulate(state[1, ], nbins=H))
          s.ind <- rbind(s.ind, tabulate(state[2, ], nbins=nS))
          t.out <- c(t.out, t[t.int])
        }
        setTxtProgressBar(pb, t.int)
      }
      
      rich <- rep(NA, dim(s.ind)[1])
      for (k in 1:dim(s.ind)[1]){
        rich[k] <- sum(s.ind[k, ] > 0)
      }

      rich_bar[i, j] <- mean(rich)
      stored_pars <- c(stored_pars, pars)
      timesteps <- c(timesteps, t.int)
      wasted_contacts <- c(wasted_contacts, 1 - successful_cntcts / n_contacts)
      win_transmission <- c(win_transmission, win_t / t[t.int])
      among_transmission <- c(among_transmission, amng_t / t[t.int])
    }  
  }
  
  res <- list(trans_bar = trans_bar, 
              rich_bar = rich_bar, 
              stored_pars = stored_pars,
              timesteps=timesteps,
              ER = ER, 
              Earray = Earray, 
              n.ind = n.ind, 
              s.ind = s.ind,
              t=t.out,
              hosts=hosts, 
              symbionts=symbionts, 
              wasted_contacts = wasted_contacts, 
              win_transmission = win_transmission, 
              among_transmission = among_transmission)
  class(res) <- "symb"
  return(res)
}

# declare plotting functions
plot.symb <- function(res, ...){  
  nsteps <- dim(res$t)
  par(mfrow=c(1, 2))
  plot(x=res$t, y=res$s.ind[, 1], type="l", 
       ylim=c(0, max(res$s.ind)), 
       xlab="Time", ylab="Number of infected hosts")
  for (k in 1:dim(res$s.ind)[2]){
    lines(x=res$t, y=res$s.ind[, k], col=k + 2)
  }
  plot(x=res$t, y=res$n.ind[, 1], type="l", 
       ylim=c(0, max(res$n.ind)), 
       xlab="Time", 
       ylab="Number of hosts")
  
  for (k in 1:dim(res$n.ind)[2]){
    lines(x=res$t, y=res$n.ind[, k], col=k+2)
  }
  par(mfrow=c(1, 1))
}


check <- FALSE

if (check){
  system.time(testout <- mpi_f(maxt=10000, nS=100, H=100, sig.s=1, c=.001,
                               beta_d_min=0, beta_d_max=0, phi=100, r=1,
                               mode="dens", cells=100))
  
  # view timeseries
  plot(testout)
  str(testout)
  
  library(ggplot2)
  # view niches
  ggplot(testout$symbionts$sniche.d, 
         aes(x=host.condition, y=Pr.estab)) + 
    geom_line(aes(col=factor(symbiont.species), group=factor(symbiont.species)), size=2) + 
    theme_classic() + 
    geom_rug(sides="b") +
    xlab("Within-host condition") + 
    ylab("Probability of establishment") + 
    theme(legend.position="bottom") + 
    theme(legend.background = element_rect(color="black", size=.5))
}


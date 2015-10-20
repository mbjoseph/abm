mpi_f <- function(maxt=1, H=10, nS=10, 
                  a_pen=1, sig.s=10, rs=.1, gamma=0.1, 
                  cells=100, r=2, d=.1, sERmin=-50, sERmax=50,
                  beta_d_min=0, beta_d_max=0, vary_beta = FALSE, 
                  c=.1, phi=5, mode="dens"){
  stopifnot(mode == "dens" | mode == "freq")

  # generate environmental ranges
  #sERmin <- -50 # regional symbiont niche minimum and maximum
  #sERmax <- 50
  Earray <- sort(runif(2, sERmin, sERmax))
  
  event_names <- c("birth", "death", "colon", "rains", "recov", "cntct")
  event_index <- rep(1:cells, times=4)
  # generate host-symbiont species indices for use with transmission matrix
  host_index <- rep(1:H, nS)
  symb_index <- rep(1:nS, each=H)
  non_cntct_indices <- 4*cells + H*nS # for ratefun extraction
  
  # initialize state arrays, time objects, and abundance counters
  state <- array(0, dim=c(2, cells))
  state_unchanged <- 0
  n.ind <- array(dim=c(maxt, H))
  n.ind[1, ] <- 0
  s.ind <- array(dim=c(maxt, nS))
  s.ind[1, ] <- 0
  rnames <- c(rep(c("birth", "death", "colon", "recov"), 
                each=cells), rep("rains", H*nS), rep("cntct", H*nS))
  ev <- rep(NA, maxt) # event counter
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
                         sEmin = Earray[1], sEmax = Earray[2], 
                         sERmin, sERmax, sig.s)
  
  # generate effects of infection on hosts
  if (vary_beta){
    beta_d <- runif(nS, beta_d_min, beta_d_max)
  } else {
    beta_d <- runif(1, beta_d_min, beta_d_max)    
  }
  
  # store parameters for later extraction
  pars <- list(r = r, d = d, c = c, beta_d = beta_d,
               phi = phi, a_pen = a_pen, rs = rs,
               sig.s = sig.s,
               cells = cells, H=H, nS=nS,
               gamma= gamma, 
               hosts=hosts,
               symbionts=symbionts, 
               mode=mode, 
               host_index=host_index, symb_index=symb_index,
               non_cntct_indices = non_cntct_indices, 
               vary_beta = vary_beta)
  
  t <- 0
  t.out <- 0
  t.int <- 1
  tau <- NULL
  
  pb <- txtProgressBar(min=0, max=maxt, style=3)
  while(t.int < maxt){
    ratelist <- ratefun(state, pars)
    tot.rates <- sum(ratelist)
    stopifnot(!is.na(tot.rates))
    stopifnot(all(ratelist >= 0))
    tau[t.int] <- rexp(1, tot.rates)
    event <- sample.int(ntrans, size=1, prob=ratelist)
    event_type <- rnames[event]
    ev[t.int] <- event_type
    if (event_type=="cntct"){
      # determine which host will be infected
      h_sp <- host_index[event - non_cntct_indices]
      cell_num <- which(state[1, ] == h_sp & state[2, ] == 0)[1]
      # take first individual, they're exchangeable
      n_contacts <- n_contacts + 1
    } else {
      if (event_type == "rains"){
        # determine which host species and individual
        h_sp <- host_index[event - 4*cells]
        cell_num <- which(state[1, ] == h_sp & state[2, ] == 0)[1]
        stopifnot(!is.na(cell_num))
      } else {
        stopifnot(event <= length(event_index))
        cell_num <- event_index[event]
      }
    }
    
    # carry out action on chosen cell
    nextstep <- ABMstep(state, event_type, 
                        cell=cell_num, params=pars, event)
    
    # update time
    t.int <- t.int + 1
    t[t.int] <- t[t.int - 1] + tau[t.int-1]
    
    state_change <- !all(state == nextstep$state)
    
    if (state_change){
      state <- nextstep$state
      n.ind[t.int, ] <- tabulate(state[1, ], nbins=H)
      s.ind[t.int, ] <- tabulate(state[2, ], nbins=nS)
      t.out <- c(t.out, t[t.int])
    } else {
      state_unchanged <- state_unchanged + 1
    }
    setTxtProgressBar(pb, t.int)
  }
  
  rich <- rep(NA, dim(s.ind)[1])
  for (k in 1:dim(s.ind)[1]){
    rich[k] <- sum(s.ind[k, ] > 0)
  }

  rich_bar <- mean(rich)
  timesteps <- t.int
  
  na_steps <- is.na(n.ind[, 1])
  n.ind <- n.ind[!na_steps, ]
  s.ind <- s.ind[!na_steps, ]
  
  res <- list(rich_bar = rich_bar, 
              pars = pars,
              timesteps=timesteps,
              ER = Earray[2] - Earray[1], 
              Earray = Earray, 
              n.ind = n.ind, 
              s.ind = s.ind,
              t=t.out,
              hosts=hosts, 
              state_unchanged=state_unchanged, 
              ev=ev)
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
  system.time(testout <- mpi_f(maxt=10000, nS=100, H=100, sig.s=runif(100, .5, 50), 
                               c=.0001, 
                               beta_d_min=-1, beta_d_max=100, vary_beta=TRUE,
                               phi=10, r=.1, rs=1,
                               mode="dens", cells=100, a_pen=1, gamma=0, d=.08))
  # view timeseries
  plot(testout)
  str(testout)
  table(testout$ev)
  
  dump('mpi_f', file='mpi_f.R')
  source('mpi_f.R')
  
  library(aprof)
  tmp <- tempfile()
  Rprof(tmp, line.profiling=TRUE)
  mpi_f(maxt=10000, nS=100, H=100, sig.s=1, c=.001,
        beta_d_min=0, beta_d_max=0, phi=100, r=1,
        mode="dens", cells=100, a_pen=1)
  Rprof(append=FALSE)
  fooaprof <- aprof('mpi_f.R', tmp)
  summary(fooaprof)
  plot(fooaprof)
  profileplot(fooaprof)
  
  library(lineprof)
  p <- lineprof(mpi_f(maxt=10000, nS=100, H=100, sig.s=1, c=.001,
                                         beta_d_min=0, beta_d_max=0, phi=100, r=1,
                                         mode="dens", cells=100, a_pen=1))
  shine(p)
  
  
  library(ggplot2)
  # view niches
  ggplot(testout$pars$symbionts$sniche.d, 
         aes(x=host.condition, y=Pr.estab)) + 
    geom_line(aes(col=factor(symbiont.species), group=factor(symbiont.species)), size=2) + 
    theme_classic() + 
    geom_rug(sides="b") +
    xlab("Within-host condition") + 
    ylab("Probability of establishment") + 
    theme(legend.position="bottom") + 
    theme(legend.background = element_rect(color="black", size=.5))
}


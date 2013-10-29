# master function for static symbiont IBM
staticIBM <- function(timesteps=500,
                      P=100, pEmin=-30, pEmax=30, pERmin=-30, pERmax=30,
                      sig.p = 5, gamma = .3, pI=.01,
                      network.pow=1, network.za=1){
  require(igraph)
  
  # initialize host and symbiont communities
  inits <- static_init(pEmin=pEmin, pEmax=pEmax, P=P, 
                       pERmin=pERmin, pERmax=pERmax,
                       sig.p=sig.p)
  
  # initialize output objects
  state <- inits$state
  A <- inits$A
  pstate <- array(0, dim=c(timesteps, A, P))
  host.richness <- inits$N
  parasite.richness <- rep(NA, timesteps)
  parasite.richness[1] <- 0
  parasite.pr.occ <- rep(NA, timesteps)
  parasite.pr.occ[1] <- 0
  
  pb <- txtProgressBar(min = 0, max = timesteps, style = 3)
  for (t in 2:timesteps){
    pstate[t,,] <- pstate[t-1,,]
    transmitted <- transmit_symbionts(t, state, pstate, network.pow, network.za, static=T)
    symbiont.immigration <- rain_symbionts(state, pstate, A, P, 
                                           pI, t, pPcol = inits$pPcol, static=T)
    symbiont.attempts <- transmitted + symbiont.immigration
    resolution <- resolve(symbiont.attempts, P)
    pstate[t, , ] <- pstate[t, , ] + resolution
    
    # HOST RECOVERY / PARASITE DEATH #
    recovery <- array(rbinom(A*N, 1, state*gamma), dim=c(A, N))
    # for the hosts that recover, remove parasites
    rec.sites <- apply(recovery, 1, sum) # recovered sites
    pstate[t, rec.sites,] <- 0 # host recovers, is again susceptible
    parasite.richness[t] <- sum(apply(pstate[t,,], 2, max))
    setTxtProgressBar(pb, t)
  }
  # check to confirm there are not multiple indivduals in any sites
  stopifnot(all(pstate[t,,] < 2))
  return(list(host.richness=host.richness, 
              state=state, pstate=pstate,
              pniche.d=inits$pniche.d, 
              parasite.richness=parasite.richness))
}
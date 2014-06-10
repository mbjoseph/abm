# master function for static symbiont IBM
staticIBM <- function(timesteps=500,
                      Y=5, modalO=2, z=0.1,
                      P=100, pEmin=-30, pEmax=30, pERmin=-30, pERmax=30,
                      sig.p = 5, gamma = .3, pI=.01,
                      network.pow=1, network.za=1,
                      A=NULL, N=NULL, prest=FALSE, cij=1){
  require(igraph)
  
  # initialize host and symbiont communities
  inits <- static_init(Y=Y, modalO=modalO, z=z,
    pEmin=pEmin, pEmax=pEmax, P=P, prest=prest, 
                       pERmin=pERmin, pERmax=pERmax,
                       sig.p=sig.p, 
                       A=A, N=N)
  
  # initialize output objects
  state <- inits$state
  A <- inits$A
  N <- inits$N
  pPcol <- inits$pPcol
  pstate <- array(0, dim=c(timesteps, A, P))
  host.richness <- inits$N
  parasite.richness <- rep(NA, timesteps)
  parasite.richness[1] <- 0
  parasite.pr.occ <- rep(NA, timesteps)
  parasite.pr.occ[1] <- 0
  beta <- rep(NA, timesteps)
  all.prevalence <- rep(NA, timesteps)
  
  pb <- txtProgressBar(min = 0, max = timesteps, style = 3)
  for (t in 2:timesteps){
    pstate[t,,] <- pstate[t-1,,]
    transmitted <- transmit(t, state, pstate, network.pow, network.za, static=T, pPcol, cij=cij)
    symbiont.immigration <- rain_symbionts(state, pstate, A, P, 
                                           pI, t, pPcol = pPcol, static=T)
    symbiont.attempts <- symbiont.immigration + transmitted 
    resolution <- resolve(symbiont.attempts, P)
    pstate[t, , ] <- pstate[t, , ] + resolution
    
    # HOST RECOVERY / PARASITE DEATH #
    recovery <- rbinom(A, 1, gamma)
    # for the hosts that recover, remove parasites
    rec.sites <- which(recovery == 1) # recovered sites
    pstate[t, rec.sites,] <- 0 # host recovers, is again susceptible
    
    # calculate the number of new infections
    diff <- pstate[t,,] - pstate[t-1,,]
    beta[t] <- length(which(diff == 1))
    
    parasite.richness[t] <- sum(apply(pstate[t,,], 2, max))
    all.prevalence[t] <- sum(apply(pstate[t,,], 1, max)) / A # number of hosts occupied / A
    setTxtProgressBar(pb, t)
  }
  # check to confirm there are not multiple indivduals in any sites
  stopifnot(all(pstate[t,,] < 2))
  return(list(host.richness=host.richness, 
              state=state, pstate=pstate,
              pniche.d=inits$pniche.d, 
              parasite.richness=parasite.richness,
              beta=beta, pE = inits$pE, 
              all.prevalence = all.prevalence))
}
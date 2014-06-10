fullIBM <- function(A=100, N=100, ERmin=-30, ERmax=30, Emin=-30, Emax=30, 
                    sig=5, timesteps=500, pM=.1, pR=1, R=2, I=0.1,
                        P=100, pEmin=-30, pEmax=30, pERmin=-30, pERmax=30,
                        sig.p = 5, gamma = .3, pI=.01,
                        network.pow=1, network.za=1){
  require(igraph)

  # initialize community and environmental traits
  inits <- com_init(A=A, Emin=Emin, Emax=Emax, N=N, 
                    pEmin=pEmin, pEmax=pEmax, P=P,
                    ERmin=ERmin, ERmax=ERmax, 
                    pERmin=pERmin, pERmax=pERmax,
                    sig=sig, sig.p=sig.p)
  
  # initialize output objects
  state <- array(0, dim=c(timesteps, A, N))
  pstate <- array(0, dim=c(timesteps, A, P))
  host.richness <- rep(NA, timesteps)
  host.richness[1] <- 0
  host.pr.occ <- rep(NA, timesteps)
  host.pr.occ[1] <- 0
  parasite.richness <- rep(NA, timesteps)
  parasite.richness[1] <- 0
  parasite.pr.occ <- rep(NA, timesteps)
  parasite.pr.occ[1] <- 0
  
  pb <- txtProgressBar(min = 0, max = timesteps, style = 3)
  for (t in 2:timesteps){
    state[t,,] <- state[t-1,,]
    pstate[t,,] <- pstate[t-1,,]
    
    transmitted <- transmit_symbionts(t, state, pstate, network.pow, network.za)
    symbiont.immigration <- rain_symbionts(state, pstate, A, P, 
                                           pI, t, pPcol = inits$pPcol)
    symbiont.attempts <- transmitted + symbiont.immigration
    resolution <- resolve(symbiont.attempts, P)
    pstate[t, , ] <- pstate[t, , ] + resolution
    
    # HOST RECOVERY / PARASITE DEATH #
    recovery <- array(rbinom(A*N, 1, c(state[t,,])*gamma), dim=c(A, N))
    # for the hosts that recover, remove parasites
    rec.sites <- apply(recovery, 1, sum) # recovered sites
    pstate[t, rec.sites,] <- 0 # host recovers, is again susceptible
    
    ## HOST DEATHS ##
    deaths <- array(rbinom(A*N, 1, c(state[t,,])*pM), dim=c(A, N))
    state[t,,] <- state[t,,] - deaths # remove hosts
    clean.sites <- apply(deaths, 1, max) # what sites are now parasite-free?
    pstate[t, clean.sites,] <- 0 # if a host in some site dies, so do the parasites
    
    ## HOST BIRTHS ##
    pot.fecundity <- array(rpois(A * N, lambda = c(state[t,,] * R)), dim=c(A, N)) # potential number of offspring
    repro <- array(rbinom(A*N, 1, pR), dim=c(A, N)) # whether reproduction actually occurs
    fecundity <- repro * pot.fecundity
    sum.fec <- apply(fecundity, 2, sum) # number of offspring per species
    
    ## HOST OFFSPRING COLONIZE EMPTY SITES ##
    occupancy <- apply(state[t,,], 1, max)
    if (sum(occupancy) < A & sum(sum.fec) > 0){ # if empty sites & new offspring
      empty.sites <- which(occupancy == 0)
      occ.sites <- which(occupancy == 1)
      # randomly assign sites (empty & filled) to each offspring individual
      col.sites <- sample(1:A, sum(sum.fec), replace=T)
      col.spec <- rep(1:N, times=sum.fec) # how many of each species colonizing
      colonizing.offspring <- array(0, dim=c(A, N)) # how many of each species colonizing each site
      for(i in 1:length(col.sites)){
        colonizing.offspring[col.sites[i], col.spec[i]] <- colonizing.offspring[col.sites[i], col.spec[i]] + 1 
      }
      # offspring attempting to colonize occupied sites fail to displace
      colonizing.offspring[occ.sites,] <- 0
      
      # which colonizing offspring can actually establish?
      binom.mat <- ifelse(colonizing.offspring > 0, 1, 0)
      colonists <- array(rbinom(n = A * N, 
                                size = c(colonizing.offspring),
                                prob = c(binom.mat * inits$Pcol)),
                         dim=c(A, N))
      
      resolution <- resolve(colonists, N)
      
      # add successful colonists
      state[t,,] <- state[t,,] + resolution
    } 
    # end of offspring colonization
    
    ## HOST IMMIGRANTS COLONIZE EMPTY SITES ##
    occupancy <- apply(state[t,,], 1, max)
    if(sum(occupancy) < A){
      empty.sites <- which(occupancy == 0)
      # which species immigrate to each site?
      immigration <- array(rbinom(length(empty.sites)*N,
                                  1, I), dim=c(length(empty.sites), N))
      # which immigrants establish?
      Pest <- immigration * inits$Pcol[empty.sites, ]
      establishment <- array(rbinom(length(Pest), 1, c(Pest)),
                             dim=c(length(empty.sites), N))
      
      resolution <- resolve(establishment, N)
      
      # add establishing immigrants
      state[t, empty.sites, ] <- state[t, empty.sites,] + resolution
    }

    host.richness[t] <- sum(apply(state[t,,], 2, max))
    host.pr.occ[t] <- length(which(state[t,,] == 1)) / A
    parasite.richness[t] <- sum(apply(pstate[t,,], 2, max))
    setTxtProgressBar(pb, t)
  }
  # check to confirm there are not multiple indivduals in any sites
  stopifnot(all(apply(state[t,,], 1, sum) < 2))
  stopifnot(all(pstate[t,,] < 2))
  return(list(host.richness=host.richness, host.pr.occ = host.pr.occ, 
              state=state, pstate=pstate,
              niche.d=inits$niche.d, pniche.d=inits$pniche.d, 
              parasite.richness=parasite.richness))
}
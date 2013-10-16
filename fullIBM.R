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

    # PARASITE TRANSMISSION WITHIN COMMUNITY #
    occupancy <- apply(state[t,,], 1, max)
    n.ind <- sum(occupancy)
    
    # build scale-free undirected contact network with Barabasi-Albert algorithm
    g <- barabasi.game(n.ind, directed=F, power=network.pow, zero.appeal=network.za)
    # turn into matrix
    m <- as.matrix(get.adjacency(g))
    # reorder to account for bias towards connected early nodes
    new.order <- sample(dim(m)[1], replace=F)
    m <- m[new.order,new.order]
    # plot(graph.adjacency(m)) # visualize the contact network
    
    # Now, identify whether susceptibles have contacted infectious
    # pstate is timestep X site X parasite
    Sindivs <- ifelse(class(pstate[t, which(occupancy==1),]) == "array", 
                      which(apply(pstate[t, which(occupancy==1),], 1,sum) < 1),
                      ifelse(sum(pstate[t, which(occupancy==1),] == 0), 1, 0))
    Iindivs <- ifelse(class(pstate[t, which(occupancy==1),]) == "array", 
                      which(apply(pstate[t, which(occupancy==1),], 1,sum) > 0),
                      ifelse(sum(pstate[t, which(occupancy==1),]) == 0, 0, 1))
    if (sum(Iindivs) > 0){ # if there are infectious individuals
      # go through transmission cycle
      Ssites <- which(occupancy == 1)[Sindivs]
      Isites <- which(occupancy == 1)[Iindivs]
      # are there any S rows in the m matrix that contact I columns?
      transmitted <- array(0, dim=c(A, P))
      if(any(m[Sindivs,Iindivs] == 1)){
        # test whether parasite establishes for each S-I contact
        # for each susceptible host
        for (s in Sindivs){
          contacts <- sum(m[s, Iindivs]) # how many contacts with infectious individuals?
          if(contacts == 0){
            next
          }else{
            Ssite <- which(occupancy == 1)[s] # which site is that host in?
            Ssp <- which(state[t, Ssite,] == 1) # which host species is susceptible
            minds <- which(m[s, ] == 1) # which individuals in the m matrix are encountered?
            i.minds <- minds[minds %in% Iindivs] # which of these are infectious?
            Isites <- which(occupancy == 1)[i.minds] # what sites do they occur in?
            Isp <- which(state[t, Isites, ] == 1) # what host species are they?
            # what parasites do they have?
            Ipars <- ifelse(class(pstate[t, Isites,]) == "matrix", 
                            apply(pstate[t, Isites,], 2, sum),
                            pstate[t, Isites,])
            # above ifelse statement to account for only one I contact
            # Ipars is a vector with the sum of parsite contacts for each parasite species
            ppestab <- inits$pPcol[Ssp,] # probability of parasite establishment in that host species
            # adjust this probability to account for interspecific transmission? built-in
            ptrial <- rbinom(sum(Ipars), Ipars, ppestab) # is establishment possible
            # store potential transmission data in a host X symbiont matrix, where number
            # indicates the number of potential establishment events
            transmitted[Ssite, ] <- ptrial
          }
        }
        resolution <- resolve(transmitted, P)
        
        # add establishing parasites
        pstate[t,,] <- pstate[t,,] + resolution
      } 
    }

    # SYMBIONT INVASION FROM REGIONAL POOL #
    # call to rain_symbionts()
    # args: sites occupied by hosts, P, pI, t, pPcol (inits$pPcol), 
    
    
    poccupancy <- apply(pstate[t,,], 1, max)
    occupancy <- apply(state[t,,], 1, max)
    open.sites <- which(occupancy - poccupancy == 1)
    if(sum(poccupancy) < sum(occupancy)){ # if not every host is infected
      empty.hosts <- which(poccupancy == 0 & occupancy == 1)
      # which parasite species immigrate to each site?
      immigration <- array(rbinom(length(empty.hosts)*P,
                                  1, pI), dim=c(length(empty.hosts), P))
      # which parasite immigrants establish?
      # first determine which species occur in the empty sites
      Ssp <- apply(state[t, open.sites,], 1, which.max)
      Pest <- immigration * inits$pPcol[Ssp, ] # P(symbiont.colonization) = I(attempting to colonize)*Pr(estab)
      symbiont.colonization <- array(rbinom(length(Pest), 1, c(Pest)),
                             dim=c(length(empty.hosts), P))
      
      # end rain_symbionts()
      
      resolution <- resolve(symbiont.colonization, P)
      
      # add establishing immigrants
      pstate[t, empty.hosts, ] <- pstate[t, empty.hosts,] + resolution
    }
    
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
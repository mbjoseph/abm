fullIBM <- function(A=100, N=100, ERmin=-30, ERmax=30, Emin=-30, Emax=30, 
                    sig=5, timesteps=500, pM=.1, pR=1, R=2, I=0.1,
                        P=100, pEmin=-30, pEmax=30, pERmin=-30, pERmax=30,
                        sig.p = 5, gamma = .3, pI=.01,
                        network.pow=1, network.za=1){
  require(igraph)
  E <- runif(A, Emin, Emax) # habitat patch environment type
  pE <- runif(N, pEmin, pEmax) # within-host environmental conditions
  mu.i <- runif(N, ERmin, ERmax) # optimum environment for hosts
  mu.p <- runif(P, pERmin, pERmax) # parasite optima
  sigma <- rep(sig, N) # runif(N, 1, 10) # # niche width
  sigma.p <- rep(sig.p, P)
  
  Z <- array(dim=c(N)) # normalization constant
  
  # now calculating all elements of Z coefficient 
  # (ensures all species have equal Pr(establishment) in regional pool)
  for (i in 1:N){
    integrand <- function(E) {
      exp(-((E - mu.i[i]) ^ 2) / (2 * sigma[i] ^ 2))
    }
    res <- integrate(integrand, lower=ERmin, upper=ERmax)
    Z[i] <- 1 / res$value
  }
  
  # probability of establishment
  Pcol <- array(dim=c(A, N))
  for (i in 1:A){
    for (j in 1:N){
      Pcol[i, j] <- Z[j] * exp(-((E[i] - mu.i[j]) ^ 2) / (2*sigma[j] ^ 2))
    }
  }
  
  # calculate parasite Z
  pZ <- rep(NA, P)
  for (i in 1:P){
    integrand <- function(pE) {
      exp(-((pE - mu.p[i]) ^ 2) / (2 * sigma.p[i] ^ 2))
    }
    res <- integrate(integrand, lower=pERmin, upper=pERmax)
    pZ[i] <- 1 / res$value
  }
  
  # parasite probability of establishment
  pPcol <- array(dim=c(N, P))
  for (i in 1:N){
    for (j in 1:P){
      pPcol[i, j] <- pZ[j] * exp(-((pE[i] - mu.p[j]) ^ 2) / (2*sigma.p[j] ^ 2))
    }
  }
  
  # store host niche data
  species <- rep(1:N, each=A)
  E <- rep(E, N)
  Pr.estab <- c(Pcol)
  niche.d <- data.frame(species, E, Pr.estab)
  niche.d <- niche.d[with(niche.d, order(species, E)),]
  
  # store parasite niche data
  parasite.species <- rep(1:P, each=N)
  pE <- rep(pE, P)
  pPr.estab <- c(pPcol)
  pniche.d <- data.frame(parasite.species, pE, pPr.estab)
  pniche.d <- pniche.d[with(pniche.d, order(parasite.species, pE)),]
  
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
    
    
    # PARASITES DO THEIR THING # 
    # PARASITE TRANSMISSION WITHIN COMMUNITY #
    # for simplicity, assume contact rates are poisson distributed and density independent
    # determine number of contacts
    # for each host, determine how many other hosts were contacted
    #contacts <- array(rpois(A*N, c(state[t,,])*exp.contact), dim=c(A, N))
    # randomly determine which individuals are contacted (homogeneous mixing)
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
    #Iindivs <- which(apply(pstate[t, which(occupancy==1),], 1, sum) > 0)
    if (sum(Iindivs) > 0){ # if there are infectious individuals
      # go through transmission cycle
      Ssites <- which(occupancy == 1)[Sindivs]
      Isites <- which(occupancy == 1)[Iindivs]
      # are there any S rows in the m matrix that contact I columns?
      establishing.parasites <- array(0, dim=c(A, P))
      if(any(m[Sindivs,Iindivs] == 1)){
        # test whether parasite establishes for each S-I contact
        # for each susceptible host
        for (s in Sindivs){
          # what is the relevant matrix? m[s, Iindivs]
          contacts <- sum(m[s, Iindivs]) # how many contacts with infectious individuals?
          if(contacts == 0){
            next
          }else{
            # which site is that host in?
            Ssite <- which(occupancy == 1)[s]
            # which host species is susceptible
            Ssp <- which(state[t, Ssite,] == 1)
            # How many contacts occur with individuals infectious with different parasite spp?
            # which individuals in the m matrix are encountered?
            minds <- which(m[s, ] == 1)
            # which of these are infectious?
            i.minds <- minds[minds %in% Iindivs]
            # what sites do they occur in?
            Isites <- which(occupancy == 1)[i.minds]
            # what host species are they?
            Isp <- which(state[t, Isites, ] == 1)
            # what parasites do they have?
            Ipars <- ifelse(class(pstate[t, Isites,]) == "matrix", 
                            apply(pstate[t, Isites,], 2, sum),
                            pstate[t, Isites,])
            # above ifelse statement to account for only one I contact
            # Ipars is a vector with the sum of parsite contacts for each parasite species
            ppestab <- pPcol[Ssp,]
            # what is the probability of parasite establishment in that host species
            # adjust this probability to account for interspecific transmission?
            # bernoulli trial to see whether establishment possible
            ptrial <- rbinom(sum(Ipars), Ipars, ppestab)
            # store potential establishment data in a host X parasite matrix, where the number
            # indicates the number of potential establishment events
            establishing.parasites[Ssite, ] <- ptrial
          }
        }
        # if multiple parasites can establish in one host individual, resolve conflict
        col.attempts <- apply(establishing.parasites, 1, sum)
        if (any(col.attempts > 1)){ # if individuals are trying to simultaneously colonize
          conflicts <- which(col.attempts > 1) # which empty sites have conflicts
          for (k in conflicts){ # for each conflict
            # how many of individuals of each species are attempting to simultaneously colonize?
            attempting <- rep(1:P, times = establishing.parasites[k,])
            # successful individual randomly selected from those attempting
            successful <- sample(attempting, size=1)
            new.row <- rep(0, length.out=P)
            new.row[successful] <- 1
            establishing.parasites[k,] <- new.row
          }
        }
        # add establishing parasites
        pstate[t,,] <- pstate[t,,] + establishing.parasites
      } 
    }
    
    # otherwise, move on to parasite colonization from regional pool
    
    # PARASITE INVASION FROM REGIONAL POOL #
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
      Pest <- immigration * pPcol[Ssp, ] # P(establishment) = I(attempting to colonize)*Pr(estab)
      establishment <- array(rbinom(length(Pest), 1, c(Pest)),
                             dim=c(length(empty.hosts), P))
      # resolve conflicts arising from simultaneous colonization
      col.attempts <- apply(establishment, 1, sum)
      if (any(col.attempts > 1)){ # if individuals are trying to simultaneously colonize
        conflicts <- which(col.attempts > 1) # which empty sites have conflicts
        for (k in conflicts){ # for each conflict
          # how many of individuals of each species are attempting to simultaneously colonize?
          attempting <- rep(1:P, times = establishment[k,])
          # successful individual randomly selected from those attempting
          successful <- sample(attempting, size=1)
          new.row <- rep(0, length.out=P)
          new.row[successful] <- 1
          establishment[k,] <- new.row
        }
      }
      # add establishing immigrants
      pstate[t, empty.hosts, ] <- pstate[t, empty.hosts,] + establishment
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
                                prob = c(binom.mat * Pcol)),
                         dim=c(A, N))
      
      # are there colonization conflicts (> 1 individual trying to colonize each site?)
      attempting <- apply(colonists, 1, sum) # number indiv. attempting to colonize each site
      if (any(attempting > 1)){
        # resolve colonization conflicts
        conflicts <- which(attempting > 1) # which sites have conflicts
        for (k in conflicts){ # for each conflict
          # how many of individuals of each species are attempting to simultaneously colonize?
          n.attempting <- rep(1:N, times = colonists[k,])
          # randomly select one successful from those attempting
          successful <- sample(n.attempting, size=1)
          new.row <- rep(0, length.out=N)
          new.row[successful] <- 1
          colonists[k,] <- new.row # individual becomes the only colonist
        }
      }
      # add successful colonists
      state[t,,] <- state[t,,] + colonists
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
      Pest <- immigration * Pcol[empty.sites, ]
      establishment <- array(rbinom(length(Pest), 1, c(Pest)),
                             dim=c(length(empty.sites), N))
      
      # resolve conflicts arising from simultaneous colonization
      col.attempts <- apply(establishment, 1, sum)
      if (any(col.attempts > 1)){ # if individuals are trying to simultaneously colonize
        conflicts <- which(col.attempts > 1) # which empty sites have conflicts
        for (k in conflicts){ # for each conflict
          # how many of individuals of each species are attempting to simultaneously colonize?
          attempting <- rep(1:N, times = establishment[k,])
          # successful individual randomly selected from those attempting
          successful <- sample(attempting, size=1)
          new.row <- rep(0, length.out=N)
          new.row[successful] <- 1
          establishment[k,] <- new.row
        }
      }
      # add establishing immigrants
      state[t, empty.sites, ] <- state[t, empty.sites,] + establishment
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
              niche.d=niche.d, pniche.d=pniche.d, parasite.richness=parasite.richness))
}
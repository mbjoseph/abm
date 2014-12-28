# transmit_symbionts allows local transmission of symbionts
transmit <- function(t, state, pstate, network.pow, 
                     network.za, static=FALSE, pPcol, cij){
  A <- dim(pstate)[2]
  P <- dim(pstate)[3]
  if (static == TRUE){
    occupancy <- apply(state, 1, max)
  } else {
    occupancy <- apply(state[t,,], 1, max)
  }
  n.ind <- sum(occupancy)
  
  # build scale-free undirected contact network with Barabasi-Albert algorithm
  g <- barabasi.game(n.ind, directed=F, power=network.pow, zero.appeal=network.za)
  m <- as.matrix(get.adjacency(g))
  new.order <- sample(dim(m)[1], replace=F) # randomizes order
  m <- m[new.order,new.order]
  #plot(graph.adjacency(m)) # visualize the contact network
  
  # determine which of the contacts are between conspecifics & heterospecifics
  m.species <- apply(state, 1, which.max)
  inter.matrix <- array(NA, dim=c(A, A))
  for (i in 1:A){
    row.spec <- m.species[i]
    col.species <- m.species[1:A]
    inter.matrix[i, ] <- ifelse(col.species != row.spec, cij, 1)
  }
  
  # use bernoulli trial to determine which connections truly happened
  bern.matrix <- array(rbinom(A ^ 2, 1, inter.matrix), dim=c(A, A))
  
  # alter m accordingly
  m <- m * bern.matrix

  # identify whether susceptibles have contacted infectious
  # number of susceptibles
  infected <- apply(pstate[t , ,], 1, sum)
  Sindivs <- which(infected == 0)
  Iindivs <- which(infected == 1)
  transmitted <- array(0, dim=c(A, P))
  if (sum(Iindivs) > 0){ # if there are infectious individuals
    # go through transmission cycle
    Ssites <- which(occupancy == 1)[Sindivs]
    Isites <- which(occupancy == 1)[Iindivs]
    # are there any S rows in the m matrix that contact I columns?
    if(any(m[Sindivs,Iindivs] == 1)){
      # test whether parasite establishes for each S-I contact
      # for each susceptible host
      for (s in Sindivs){
        contacts <- sum(m[s, Iindivs]) # how many contacts with infectious individuals?
        if(contacts == 0){
          next
        }else{
          Ssite <- which(occupancy == 1)[s] # which site is that host in?
          Ssp <- which(state[Ssite,] == 1) # which host species is susceptible
          minds <- which(m[s, ] == 1) # which individuals in the m matrix are encountered?
          i.minds <- minds[minds %in% Iindivs] # which of these are infectious?
          Isites <- which(occupancy == 1)[i.minds] # what sites do they occur in?
          Isp <- which(state[Isites, ] == 1) # what host species are they?
          interspecies.penalty <- ifelse(Ssp == Isp, 1, cij)
        #  true.contact <- rbinom(length())
          # what parasites do they have?
          if (class(pstate[t, Isites, ]) == "matrix"){
            n.each <- apply(pstate[t, Isites, ], 2, sum)
            Ipars <- rep(1:P, n.each)
          } else {
            Ipars <- which(pstate[t, Isites,] == 1)
          }
          # above ifelse statement to account for only one I contact
          # Ipars is a vector with the sum of parsite contacts for each parasite species
          ppestab <- pPcol[Ssp, Ipars] # probability of parasite establishment in that host species
          # adjust this probability to account for interspecific transmission? built-in
          ptrial <- rbinom(length(Ipars), 1, ppestab * interspecies.penalty) # is establishment possible
          # store potential transmission data in a host X symbiont matrix, where number
          # indicates the number of potential establishment events
          z <- rep(NA, P)
          for (i in 1:P){
            indx <- which(Ipars == i)
            z[i] <- sum(ptrial[indx])
          }
          transmitted[Ssite, ] <- z
        }
      }
    } 
  }
  return(transmitted)
}
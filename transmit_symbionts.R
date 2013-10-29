# transmit_symbionts allows local transmission of symbionts
transmit_symbionts <- function(t, state, pstate, network.pow, network.za, static=FALSE){
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
  
  # identify whether susceptibles have contacted infectious
  Sindivs <- ifelse(class(pstate[t, which(occupancy==1),]) == "array", 
                    which(apply(pstate[t, which(occupancy==1),], 1,sum) < 1),
                    ifelse(sum(pstate[t, which(occupancy==1),] == 0), 1, 0))
  Iindivs <- ifelse(class(pstate[t, which(occupancy==1),]) == "array", 
                    which(apply(pstate[t, which(occupancy==1),], 1,sum) > 0),
                    ifelse(sum(pstate[t, which(occupancy==1),]) == 0, 0, 1))
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
    } 
  }
  return(transmitted)
}
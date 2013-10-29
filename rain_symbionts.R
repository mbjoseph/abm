# rain_symbionts() allows colonization of hosts by symbionts from regional pool
rain_symbionts <- function(state, pstate, A, P, pI, t, pPcol, static=FALSE){
  poccupancy <- apply(pstate[t,,], 1, max)
  if (static == TRUE){
    occupancy <- apply(state, 1, max)
  } else {
    occupancy <- apply(state[t,,], 1, max)
  }
  open.sites <- which(occupancy - poccupancy == 1)
  symbiont.colonization <- array(0, dim = c(A, P))
  if(sum(poccupancy) < sum(occupancy)){ # if not every host is infected
    empty.hosts <- which(poccupancy == 0 & occupancy == 1)
    # which parasite species immigrate to each site?
    immigration <- array(rbinom(length(empty.hosts)*P, 1, pI), dim=c(length(empty.hosts), P))
    # which parasite immigrants establish?
    # first determine which species occur in the empty sites
    if (static == TRUE){
      Ssp <- apply(matrix(state[open.sites,], nrow=length(open.sites)), 1, which.max)
    } else {
      Ssp <- apply(matrix(state[t, open.sites,], nrow=length(open.sites)), 1, which.max)
    }
    Pest <- immigration * pPcol[Ssp, ] # P(symbiont.colonization) = I(attempting to colonize)*Pr(estab)
    symbiont.colonization[empty.hosts, ] <- array(rbinom(length(Pest), 1, c(Pest)),
                                   dim=c(length(empty.hosts), P))
  } 
  return(symbiont.colonization)
}
  
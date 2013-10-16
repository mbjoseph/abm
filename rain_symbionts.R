# rain_symbionts() allows colonization of hosts by symbionts from regional pool
rain_symbionts() <- function(state, pstate, P, pI, t, pPcol
  #sites occupied by hosts, P, pI, t, pPcol (inits$pPcol)
  ){
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
  } else {
    symbiont.colonization <- array(0, dim = c(A, P))
  }
  
# resolve() settles simultaneous colonization conflicts for parasites & hosts
resolve <- function(trying, n.types){
  col.attempts <- apply(trying, 1, sum)
  if (any(col.attempts > 1)){ # if individuals are trying to simultaneously colonize
    conflicts <- which(col.attempts > 1) # which empty sites have conflicts
    for (k in conflicts){ # for each conflict
      # how many of individuals of each species are attempting to simultaneously colonize?
      attempting <- rep(1:n.types, times = trying[k,])
      # successful individual randomly selected from those attempting
      successful <- sample(attempting, size=1)
      new.row <- rep(0, length.out=n.types)
      new.row[successful] <- 1
      trying[k,] <- new.row
    }
  }
  return(trying)
}
# end resolve function

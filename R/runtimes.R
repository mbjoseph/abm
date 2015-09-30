# file to evaluate runtimes for various parameter combinations
source("continuous_time/mpi_f.R")
source("continuous_time/helpers.R")
n.iter <- 5

H <- 100
nS <- 100
maxt <- seq(5000, 500000, length.out=n.iter)

runtimes <- rep(NA, n.iter)
timesteps <- rep(NA, n.iter)
tmax <- rep(NA, n.iter)

for (i in 1:n.iter){
  runtimes[i] <- system.time(testout <- mpi_f(maxt=maxt[i], nS=100, H=100, sig.s=1, 
                                              phi=3)
                             )["elapsed"]
  timesteps[i] <- testout$timesteps
  tmax[i] <- max(testout$t)
  plot(testout)
}

pairs(~ maxt + runtimes + timesteps + tmax)

# file to evaluate runtimes for various parameter combinations
source("continuous_time/janus/mpi_f.R")
source("continuous_time/janus/symb_init.R")
source("continuous_time/janus/host_init.R")
source("continuous_time/janus/helpers.R")
n.iter <- 20

H <- 100
nS <- 100
maxt <- 50000

runtimes <- rep(NA, n.iter)
timesteps <- rep(NA, n.iter)
tmax <- rep(NA, n.iter)

for (i in 1:n.iter){
  runtimes[i] <- system.time(testout <- mpi_f(iter=1, nER=1, maxt=maxt, H=H, 
                                              nS=nS, a_pen=1, sig.s=3, rs=.01, 
                                              gamma=0, cells=100, r=.4, d=.3, 
                                              beta_d_min = -1, beta_d_max = 1, 
                                              c=.001, phi=5)
                             )["elapsed"]
  timesteps[i] <- testout$timesteps
  tmax[i] <- max(testout$t)
  plot(testout)
}

pairs(~ runtimes + timesteps + tmax)

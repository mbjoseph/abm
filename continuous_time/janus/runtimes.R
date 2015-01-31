# file to evaluate runtimes for various parameter combinations
n.iter <- 50

H <- 20
nS <- 20
maxt <- 2000

runtimes <- rep(NA, n.iter)
timesteps <- rep(NA, n.iter)
tmax <- rep(NA, n.iter)

for (i in 1:n.iter){
  runtimes[i] <- system.time(testout <- mpi_f(iter=1, nER=1, maxt=maxt, H=H, 
                                              nS=nS, a_pen=1, sig.s=50, rs=.1, 
                                              gamma=0, cells=100, r=.4, d=.3, 
                                              beta_d_min = -1, beta_d_max = 1, 
                                              c=.001, phi=5)
                             )["elapsed"]
  timesteps[i] <- testout$timesteps
  tmax[i] <- max(testout$t)
}

pairs(~ runtimes + timesteps + tmax)

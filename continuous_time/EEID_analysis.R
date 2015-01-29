## Two main components for EEID poster
# 1. How does host richness & functional diversity affect symbiont transmission & richness?
#    - fixed host communities of varying functional diversity & among species transmission
#    - record symbiont transmission rates & richness
# 2. How does habitat heterogeneity affect symbiont transmission & richness?

# 1. Host richness & functional diversity's effects on symbionts
# Need to create host communities that are static that vary in heterogeneity
# and for each of these, explore consequences of varying interspecies transmiss.
# for each run, record mean symbiont transmission & richness
library(parallel)
library(Rmpi)
library(snow, lib.loc = "/home/majo3748/Rpackages")

#Make the cluster
cl <- makeCluster( mpi.universe.size(), type="MPI" )

#Check the cluster - we should get a response from each worker
clusterCall( cl, function() Sys.info()[c("nodename","machine")])

# load necessary functions
source("helpers.R")
source("symb_init.R")
source("mpi_f.R")

xx <- clusterCall(cl, mpi_f, iter=1, nER=1, maxt=1, H=100, nS=100, 
                  a_pen=1, sig.s=3, rs=.01, gamma=0, cells=100, 
                  r=.4, d=.3, beta_d = 0, c=.001, phi=2)

saveRDS(xx, file=paste(format(Sys.time(), "%b%d_%H:%M:%S"), 
                        "_res.RData", sep=""))

#Shutdown
stopCluster(cl)
mpi.quit()

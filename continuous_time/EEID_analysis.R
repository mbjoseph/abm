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
clusterSetRNGStream(cl, 123123)

#Check the cluster - we should get a response from each worker
clusterCall( cl, function() Sys.info()[c("nodename","machine")])
#mpi.remote.exec(paste(Sys.info()[c("nodename")],"checking in as",mpi.comm.rank(),"of",mpi.comm.size()))

# load necessary functions
source("helpers.R")
source("symb_init.R")

xx <- clusterCall(cl, mpi_f, nER=1, iter=1, maxt=50)
save(xx,file='test.rdata')

#Shutdown
stopCluster(cl)
mpi.quit()
library(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)

mpi.bcast.cmd( id <- mpi.comm.rank() )
mpi.bcast.cmd( np <- mpi.comm.size() )
mpi.bcast.cmd( host <- mpi.get.processor.name() )
result <- mpi.remote.exec(paste("I am", id, "of", np, "running on", host)) 

print(unlist(result))

# load functions on master 
source("helpers.R")
source("mpi_f.R")

# send functions to slaves
mpi.bcast.Robj2slave(ratefun)
mpi.bcast.Robj2slave(transmit)
mpi.bcast.Robj2slave(ABMstep)
mpi.bcast.Robj2slave(symb_init)
mpi.bcast.Robj2slave(host_init)
mpi.bcast.Robj2slave(mpi_f)

# execute
res <- mpi.remote.exec(mpi_f(maxt=50000, nS=100, H=100, sig.s=1, 
                             beta_d_min=0, beta_d_max=0, phi=3))

saveRDS(res, file=paste(format(Sys.time(), "%b%d_%H:%M:%S"), 
                       "_res.RData", sep=""))

mpi.close.Rslaves(dellog = FALSE)
mpi.exit()


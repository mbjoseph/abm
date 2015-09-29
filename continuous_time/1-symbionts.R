# Parallel execution of host-symbiont analysis
source('continuous_time/mpi_f.R')
source('continuous_time/helpers.R')
library(reshape2)
library(dplyr)
library(doMC)
registerDoMC(2)

# Section 1: host diversity, symbiont richness, and transmission
iter <- 10

system.time(
  r <- foreach(icount(iter)) %dopar% {
    mpi_f(maxt=30000, nS=50, H=50, sig.s=10, c=.0001, 
            beta_d_min=0, beta_d_max=0, phi=10, r=.1, rs=1,
            mode="dens", cells=500, a_pen=1, gamma=0, d=.08)
  }
)

# calculate symbiont richness timeseries for each iteration
rich_ts <- list()
srich_ts <- list()
ts <- list()
for (i in 1:iter){
  rich_ts[[i]] <- apply(r[[i]]$n.ind, 1, FUN=function(x) sum(x > 0))
  srich_ts[[i]] <- apply(r[[i]]$s.ind, 1, FUN=function(x) sum(x > 0))
  ts[[i]] <- r[[i]]$t
}

mts <- melt(rich_ts)
names(mts) <- c('host_richness', 'iteration')
mts$t <- melt(ts) %>%
  select(value) %>%
  unlist()
mts$symbiont_richness <- melt(srich_ts) %>%
  select(value) %>%
  unlist()

sub_unique <- function(d, col){
  numcol <- which(names(d) == col)
  keep <- cumsum(rle(d[, col])$length)
  d[keep, ]
}

ggplot(sub_unique(mts, 'host_richness'), 
       aes(x=t, y=host_richness, group=iteration)) + 
  geom_line(alpha=.5)

ggplot(sub_unique(mts, 'symbiont_richness'), 
       aes(x=t, y=symbiont_richness, group=iteration)) + 
  geom_line(alpha=.5)

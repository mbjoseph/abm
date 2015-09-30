# Parallel execution of host-symbiont analysis
source('continuous_time/mpi_f.R')
source('continuous_time/helpers.R')
library(reshape2)
library(ggplot2)
library(dplyr)
library(doMC)

registerDoMC(2)

# Section 1: host diversity, symbiont richness, and transmission
iter <- 200

system.time(
  r <- foreach(icount(iter)) %dopar% {
    mpi_f(maxt=30000, nS=50, H=50, sig.s=1, c=.0001, 
            beta_d_min=0, beta_d_max=0, phi=10, r=.1, rs=1,
            mode="dens", cells=500, a_pen=1, gamma=0, d=.08)
  }
)

# calculate symbiont richness timeseries for each iteration
rich_ts <- list()
srich_ts <- list()
ts <- list()
h_condition <- list()
eff_diversity <- list()

for (i in 1:iter){
  rich_ts[[i]] <- apply(r[[i]]$n.ind, 1, FUN=function(x) sum(x > 0))
  srich_ts[[i]] <- apply(r[[i]]$s.ind, 1, FUN=function(x) sum(x > 0))
  ts[[i]] <- r[[i]]$t
  h_condition[[i]] <- r[[i]]$pars$symbionts$sniche.d %>%
    filter(symbiont.species == 1) %>%
    select(host.condition) %>%
    unlist()
  eff_diversity[[i]] <- apply(r[[i]]$n.ind, 1, FUN = function(x){
    present <- which(x > 0)
    if (!any(present)) {
      diversity <- 0
    } else {
      h_vals <- h_condition[[i]][present]
      diversity <- max(h_vals) - min(h_vals)
    }
    diversity
    }) 
}

mdiv <- melt(eff_diversity)

mts <- melt(rich_ts)
names(mts) <- c('host_richness', 'iteration')
mts$t <- melt(ts) %>%
  select(value) %>%
  unlist()
mts$symbiont_richness <- melt(srich_ts) %>%
  select(value) %>%
  unlist()
mts$functional_diversity <- mdiv$value

mts <- tbl_df(mts)

ggplot(sub_unique(mts, 'host_richness'), 
       aes(x=t, y=host_richness, group=iteration)) + 
  geom_line(alpha=.5)

ggplot(sub_unique(mts, 'symbiont_richness'), 
       aes(x=t, y=symbiont_richness, group=iteration)) + 
  geom_line(alpha=.5)

up_lim <- 600
low_lim <- 400

mts$roundtime <- round(mts$t, 0)

sum_d <- mts %>%
  filter(roundtime < up_lim + 1, roundtime > low_lim - 1) %>%
  sub_unique('roundtime') %>%
  group_by(iteration) %>%
  summarize(hmean = mean(host_richness),
            hsd = sd(host_richness), 
            smean = mean(symbiont_richness), 
            ssd = sd(symbiont_richness), 
            dmean = mean(functional_diversity),
            dsd = sd(functional_diversity),
            n=n())
alph <- .5
ggplot(sum_d, aes(x=dmean, y=smean)) + 
  geom_point(alpha=alph) + 
  geom_segment(aes(x=dmean, xend=dmean, y=smean - ssd, yend = smean + ssd), 
               alpha=alph) +
  geom_segment(aes(x=dmean - dsd, xend=dmean + dsd, y=smean, yend=smean), 
               alpha=alph) +
  xlab('Functional diversity') + 
  ylab('Symbiont richness')

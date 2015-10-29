# Parallel execution of host-symbiont analysis
source('R/mpi_f.R')
source('R/helpers.R')
library(reshape2)
library(ggplot2)
library(dplyr)
library(doMC)

registerDoMC(4)

# Section 1: host diversity, symbiont richness, and transmission
iter <- 990
dir <- paste(getwd(), "/R/1-symbionts/sim_results", sep="")

op <- options(digits.secs =7)
system.time(
foreach(icount(iter)) %dopar% {
    mpi_f(maxt=30000, nS=50, H=50, sig.s=1, c=.0001, 
            beta_d_min=0, beta_d_max=0, phi=10, r=.1, rs=1,
            mode="dens", cells=500, a_pen=1, gamma=0, d=.08) %>%
    saveRDS(file=paste0(dir, "/res-", 
                        Sys.time(), ".RData"))
  }
)
options(op)

data <- list.files(dir, pattern="res*")
iter <- length(data)

# calculate symbiont richness timeseries for each iteration
rich_ts <- list()
srich_ts <- list()
ts <- list()
h_condition <- list()
eff_diversity <- list()
trans <- rep(NA, iter)

pb <- txtProgressBar(max=iter, style=3)
for (i in 1:iter){
  d <- readRDS(paste(dir, data[i], sep="/"))
  rich_ts[[i]] <- apply(d$n.ind, 1, FUN=function(x) sum(x > 0))
  srich_ts[[i]] <- apply(d$s.ind, 1, FUN=function(x) sum(x > 0))
  ts[[i]] <- d$t
  h_condition[[i]] <- d$pars$symbionts$sniche.d %>%
    filter(symbiont.species == 1) %>%
    select(host.condition) %>%
    unlist()
  eff_diversity[[i]] <- apply(d$n.ind, 1, FUN = function(x){
    present <- which(x > 0)
    if (!any(present)) {
      diversity <- 0
    } else {
      h_vals <- h_condition[[i]][present]
      diversity <- max(h_vals) - min(h_vals)
    }
    diversity
    }) 
  trans[i] <- sum(d$ev == 'cntct', na.rm=T) / max(d$t)
  rm(d)
  setTxtProgressBar(pb, i)
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
  geom_line(alpha=.2)

ggplot(sub_unique(mts, 'symbiont_richness'), 
       aes(x=t, y=symbiont_richness, group=iteration)) + 
  geom_line(alpha=.07)

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
            cor_div = cor(host_richness, symbiont_richness),
            n=n())
sum_d$trans <- trans[match(sum_d$iteration, 1:length(trans))]
alph <- .5
ggplot(sum_d, aes(x=dmean, y=smean)) + 
  geom_point(alpha=alph) + 
  #geom_segment(aes(x=dmean, xend=dmean, y=smean - ssd, yend = smean + ssd), 
  #             alpha=alph) +
  #geom_segment(aes(x=dmean - dsd, xend=dmean + dsd, y=smean, yend=smean), 
  #             alpha=alph) +
  xlab('Functional diversity') + 
  ylab('Symbiont richness')

ggplot(sum_d, aes(x=dmean, y=trans)) + 
  geom_point(alpha=alph) + 
  xlab('Functional diversity') + 
  ylab('Transmission rate')

ggplot(sum_d, aes(x=smean, y=trans)) + 
  geom_point(alpha=alph) + 
  xlab('Symbiont richness') + 
  ylab('Transmission rate')

ggplot(sum_d, aes(x=dmean, y=cor_div)) + 
  geom_point(alpha=alph) + 
  xlab('Functional diversity') + 
  ylab('Correlation: host and symbiont richness') + 
  geom_hline(yintercept=0, linetype='dashed')

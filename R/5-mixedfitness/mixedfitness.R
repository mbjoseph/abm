# Parallel execution of host-symbiont analysis
source('R/mpi_f.R')
source('R/helpers.R')
library(reshape2)
library(ggplot2)
library(dplyr)
library(doMC)

registerDoMC(2)

# Section 1: host diversity, symbiont nichewidth, and transmission
iter <- 10
dir <- paste(getwd(), "/R/5-mixedfitness/sim_results", sep="")
beta_max <- 10
beta_min <- -.99
nS <- 50

foreach(icount(iter)) %dopar% {
  r <- mpi_f(maxt=30000, nS=nS, H=50, sig.s=1, 
             c=.0001, vary_beta=TRUE,
             beta_d_min=beta_min, beta_d_max=beta_max, phi=10, r=.1, rs=1,
             mode="dens", cells=500, a_pen=1, gamma=0, d=.08)
  saveRDS(r, file=paste0(dir, "/res-", 
                         format(Sys.time(), "%b%d_%H:%M:%S"), ".RData"))
}


data <- list.files(dir, pattern="res*")
iter <- length(data)

# calculate symbiont richness timeseries for each iteration
rich_ts <- list()
srich_ts <- list()
ts <- list()
h_condition <- list()
eff_diversity <- list()
trans <- rep(NA, iter)
trans_each <- array(dim=c(nS, iter))
beta_s <- array(dim=c(nS, iter))

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
  # calculate colonization & transmission rates for each symbiont species
  trans_each[, i] <- apply(d$s.ind, 2, FUN=function(x) sum(diff(x) == 1) / max(d$t))
  beta_s[, i] <- d$pars$beta_d
  rm(d)
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
            cor_div = cor(host_richness, symbiont_richness),
            n=n())
sum_d$trans <- trans[match(sum_d$iteration, 1:length(trans))]
ss <- melt(sigma_s)
names(ss) <- c('symbiont', 'iteration', 'niche_width')
mtrans <- melt(trans_each)
names(mtrans) <- c('symbiont', 'iteration', 'trans')
jt <- full_join(ss, mtrans)
jt$dmean <- sum_d$dmean[match(jt$iteration, sum_d$iteration)]
jt$smean <- sum_d$smean[match(jt$iteration, sum_d$iteration)]

ggplot(jt, aes(x=dmean, y=trans, color=niche_width)) + 
  geom_point(alpha=.5) +
  xlab('Functional diversity') + 
  ylab('Symbiont transmission & colonization') + 
  scale_color_gradientn(colours=rainbow(3))


library(gtools)
jt$sigma_bin <- quantcut(jt$niche_width, q=6)
alph <- .3
ggplot(jt, aes(x=dmean, y=smean)) + 
  geom_point(alpha=alph) + 
  xlab('Functional diversity') + 
  ylab('Symbiont richness')

ggplot(jt, aes(x=dmean, y=trans, color=sigma_bin)) + 
  geom_point(alpha=.2) +
  xlab('Functional diversity') + 
  ylab('Symbiont transmission & colonization') + 
  facet_wrap(~sigma_bin)

ggplot(sum_d, aes(x=dmean, y=cor_div)) + 
  geom_point(alpha=alph) + 
  xlab('Functional diversity') + 
  ylab('Correlation: host and symbiont richness') + 
  geom_abline(yintercept=0, slope=0, linetype='dashed')
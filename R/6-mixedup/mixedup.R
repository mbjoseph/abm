# Parallel execution of host-symbiont analysis
source('R/mpi_f.R')
source('R/helpers.R')
library(reshape2)
library(ggplot2)
library(dplyr)
library(doMC)

registerDoMC(2)

iter <- 1000
dir <- paste(getwd(), "/R/6-mixedup/sim_results", sep="")
beta_max <- .99
beta_min <- -.99
sigma_max <- 50
sigma_min <- .5
nS <- 50

foreach(icount(iter)) %dopar% {
  r <- mpi_f(maxt=60000, nS=nS, H=50, sig.s=runif(nS, sigma_min, sigma_max), 
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
sigma_s <- array(dim=c(nS, iter))

pb <- txtProgressBar(max=iter)
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
  sigma_s[, i] <- d$pars$sig.s
  beta_s[, i] <- d$pars$beta_d
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
  geom_line(alpha=.01)

ggplot(sub_unique(mts, 'symbiont_richness'), 
       aes(x=t, y=symbiont_richness, group=iteration)) + 
  geom_line(alpha=.01) + 
  geom_vline(xintercept=200) + 
  stat_smooth(alpha=.1, se=FALSE, fill='blue', size=.1) + 
  ylim(0, 22)

ggplot(sub_unique(mts, 'functional_diversity'), 
       aes(x=t, y=functional_diversity, group=iteration)) + 
  geom_line(alpha=.1)

up_lim <- 1000
low_lim <- 500#400

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
bs <- melt(beta_s)
names(bs) <- c('symbiont', 'iteration', 'beta_d')
mtrans <- melt(trans_each)
names(mtrans) <- c('symbiont', 'iteration', 'trans')

ss <- melt(sigma_s)
names(ss) <- c('symbiont', 'iteration', 'niche_width')
jt <- full_join(ss, mtrans)
jt <- full_join(jt, bs)
jt$dmean <- sum_d$dmean[match(jt$iteration, sum_d$iteration)]
jt$smean <- sum_d$smean[match(jt$iteration, sum_d$iteration)]

ggplot(jt, aes(x=dmean, y=trans, color=beta_d)) + 
  geom_point(alpha=.5) +
  xlab('Functional diversity') + 
  ylab('Symbiont transmission & colonization') + 
  scale_color_gradientn(colours=rainbow(3))

library(gtools)
n_bin <- 5
jt$beta_bin <- quantcut(jt$beta_d, q=n_bin)
jt$sigma_bin <- quantcut(jt$niche_width, q=n_bin)
alph <- .3

ggplot(jt, aes(x=dmean, y=smean)) + 
  geom_point(alpha=alph) + 
  xlab('Functional diversity') + 
  ylab('Symbiont richness')

ggplot(sum_d, aes(x=dmean, y=smean)) + 
  geom_point() + 
  xlab('Functional diversity') + 
  ylab('Symbiont richness')

ggplot(jt, aes(x=dmean, y=smean)) + 
  geom_point(alpha=alph) + 
  xlab('Functional diversity') + 
  ylab('Symbiont richness') + 
  facet_grid(sigma_bin~beta_bin)

ggplot(jt, aes(x=dmean, y=trans)) + 
  geom_point(alpha=.2) +
  xlab('Functional diversity') + 
  ylab('Symbiont transmission & colonization') + 
  facet_grid(sigma_bin~beta_bin)

ggplot(sum_d, aes(x=dmean, y=cor_div)) + 
  geom_point(alpha=alph) + 
  xlab('Functional diversity') + 
  ylab('Correlation: host and symbiont richness') + 
  geom_abline(yintercept=0, slope=0, linetype='dashed')

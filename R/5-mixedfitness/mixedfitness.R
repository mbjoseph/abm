# Parallel execution of host-symbiont analysis
source('R/mpi_f.R')
source('R/helpers.R')
library(reshape2)
library(ggplot2)
library(dplyr)
library(doMC)

registerDoMC(2)

dir <- paste(getwd(), "/R/5-mixedfitness/sim_results", sep="")
data <- list.files(dir, pattern="res*")
iter <- length(data)
beta_max <- 10
beta_min <- -.99
nS <- 50

if (iter < 1000) {
  iter <- 1000 - iter
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
}

if (!("res.Rdata" %in% list.files(paste(getwd(), "/R/5-mixedfitness/", sep="")))) {
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
  save(list=ls(), file=paste(getwd(), "/R/5-mixedfitness/res.Rdata", sep=""))
} else {
  load(paste(getwd(), "/R/5-mixedfitness/res.Rdata", sep=""))
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
  geom_line(alpha=.01)

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
            beta_range = range()
            n=n())
sum_d$trans <- trans[match(sum_d$iteration, 1:length(trans))]
bs <- melt(beta_s)
names(bs) <- c('symbiont', 'iteration', 'beta_d')

sum_beta <- bs %>%
  group_by(iteration) %>% 
  summarize(beta_range = max(beta_d) - min(beta_d))

sum_d <- full_join(sum_d, sum_beta)


mtrans <- melt(trans_each)
names(mtrans) <- c('symbiont', 'iteration', 'trans')
jt <- full_join(bs, mtrans)
jt$dmean <- sum_d$dmean[match(jt$iteration, sum_d$iteration)]
jt$smean <- sum_d$smean[match(jt$iteration, sum_d$iteration)]

ggplot(jt, aes(x=dmean, y=trans, color=beta_d)) + 
  geom_point(alpha=.5) +
  xlab('Functional diversity') + 
  ylab('Symbiont transmission & colonization') + 
  scale_color_gradientn(colours=rainbow(3))


library(gtools)
jt$beta_bin <- quantcut(jt$beta_d, q=6)
alph <- .5
p1 <- ggplot(sum_d, aes(x=dmean, y=smean)) + 
  geom_point(alpha=alph) + 
  xlab('Host functional diversity') + 
  ylab('Symbiont richness')
p1

p2 <- ggplot(jt, aes(x=dmean, y=trans)) + 
  geom_point(alpha=.2) +
  xlab('Host functional diversity') + 
  ylab('Symbiont transmission & colonization') + 
  facet_wrap(~beta_bin)
p2

ggplot(sum_d, aes(x=dmean, y=cor_div)) + 
  geom_point(alpha=alph) + 
  xlab('Functional diversity') + 
  ylab('Correlation: host and symbiont richness') + 
  geom_abline(yintercept=0, slope=0, linetype='dashed')

library(gridExtra)
library(grid)
xj <- .1
yj <- .6
myplot1 <- arrangeGrob(p1, top = textGrob("A", 
                                          x = unit(xj, "npc"), 
                                          y = unit(yj, "npc"), 
                                          just = c("left","top")))
myplot2 <- arrangeGrob(p2, top = textGrob("B", 
                                          x = unit(xj, "npc"), 
                                          y = unit(yj, "npc"), 
                                          just = c("left","top")))
grid.arrange(myplot1, myplot2, ncol=2, widths=c(1, 1.5))
dev.copy2pdf(file="paper/fig/fig5.pdf", width = 8, height = 4)


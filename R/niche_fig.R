# Figure to illustrate the symbiont niches
source('R/helpers.R')
library(ggplot2)

S <- 2
mu <- runif(S, -50, 50)
sigma <- 2
xvals <- seq(-50, 50, .1)
g <- expand.grid(x=xvals, Symbiont = 1:S)


out <- symb_init(H=1000, S=S, sEmin=-50, sEmax=50, 
                       sERmin=-50, sERmax=50, sig.s=10)
out$sniche.d$Symbiont <- as.factor(out$sniche.d$symbiont.species)
host_d <- data.frame(x=runif(20, -30, 30))
mu_d <- data.frame(mu=out$mu.s, Symbiont=factor(1:S), y=0)
sd_d <- data.frame(sd = out$sig.s, Symbiont=factor(1:S), 
                   low = out$mu.s - out$sig.s, 
                   high = out$mu.s + out$sig.s)

# find niche optima and standard deviations
ggplot(out$sniche.d, aes(y=Pr.estab)) + 
  geom_line(aes(x=host.condition, color=Symbiont, group=Symbiont), 
            size=3, alpha=.7) + 
  theme_classic() + 
  geom_point(aes(x=x, y=-.01), data=host_d, pch="|", size=5) + 
  xlab("Within-host condition niche axis") + 
  ylab("Probability of establishment") +
  theme(legend.position="none") + 
  theme(legend.background = element_rect(color="black", size=.5)) + 
  scale_color_brewer(palette="Set1") + 
  geom_text(aes(x=mu, y=.005, 
                label = paste("mu[", Symbiont, "]", sep = ""), 
                color=Symbiont),
            parse = TRUE, data=mu_d, size=6) + 
  geom_text(aes(x=low, y=.005, 
                label = paste0("mu[", Symbiont, "] - sigma[", Symbiont, "]"), 
                color= Symbiont),
            parse = TRUE, data=sd_d, size=4) + 
  geom_text(aes(x=high, y=.005, 
                label = paste0("mu[", Symbiont, "] + sigma[", Symbiont, "]"), 
                color= Symbiont),
            parse = TRUE, data=sd_d, size=4) + 
  scale_y_continuous(breaks=seq(0, max(out$sniche.d$Pr.estab), .01)) + 
  geom_segment(aes(x = min(x), 
                   xend = max(x), 
                   y = -.01, yend = -.01), 
               size=3, alpha=.01, data=host_d) + 
  geom_text(aes(x=min(x) + .5 * (max(x) - min(x)), y=-.005, 
                label = 'Host functional diversity'), 
            data=host_d, size=4, alpha=.03) + 
  geom_text(aes(x=mu, y=.01, 
                label = paste("Symbiont", Symbiont, sep = " "), 
                color=Symbiont), data=mu_d, size=5)
dev.copy2pdf(file="paper/fig/niche.pdf", width = 7, height = 4)

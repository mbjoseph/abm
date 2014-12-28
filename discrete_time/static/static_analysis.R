# static model analysis
# IBM analysis
setwd("/home/max/Documents/manuscripts/abm")
source("static/static_init.R")
source("resolve.R")
source("rain_symbionts.R")
source("static/transmit.R")
source("static/staticIBM.R")
output <- staticIBM(Y = 5, P=100, sig.p=1, gamma=0.1, network.pow = 1, 
                    pEmin = -5, pEmax = 5, pERmin=-5, pERmax=5, pI=.01, 
                    timesteps=2000, A=100, N=10)
output2 <- staticIBM(Y = 5, P=100, sig.p=1, gamma=0.1, network.pow = 1, 
                    pEmin = -5, pEmax = 5, pERmin=-5, pERmax=5, pI=.01, 
                    timesteps=2000, A=100, N=10)
output3 <- staticIBM(Y = 5, P=100, sig.p=1, gamma=0.1, network.pow = 1, 
                     pEmin = -5, pEmax = 5, pERmin=-5, pERmax=5, pI=.01, 
                     timesteps=2000, A=100, N=10)
output4 <- staticIBM(Y = 5, P=100, sig.p=1, gamma=0.1, network.pow = 1, 
                     pEmin = -5, pEmax = 5, pERmin=-5, pERmax=5, pI=.01, 
                     timesteps=2000, A=100, N=10)
str(output)

par(mfrow=c(1, 2))
plot(output$parasite.richness, type="l", lty=1, col="black",
     xlab="Timestep", ylab="Symbiont richness", 
     ylim=c(0, max(c(output$beta, output$parasite.richness), na.rm=T)))
abline(v=50, col="blue", lty=2, lwd=2)
abline(h=mean(output$parasite.richness[50:1000]), col="red", lty=2, lwd=2)
hist(output$parasite.richness[50:1000], xlab="Symbiont richness", breaks=seq(20, 41, 1), main="")
abline(v=mean(output$parasite.richness[50:1000]), col="red", lty=2, lwd=2)

par(mfrow=c(1,1))
plot(output$parasite.richness, type="l", lty=1, col="black",
     xlab="Timestep", ylab="Symbiont richness", 
     ylim=c(0, max(c(output$beta, output$parasite.richness), na.rm=T)))
lines(x=1:length(output4$parasite.richness), y=output4$parasite.richness, col="orange")
lines(x=1:length(output2$parasite.richness), y=output2$parasite.richness, col="green")
lines(x=1:length(output3$parasite.richness), y=output3$parasite.richness, col="purple")


lines(output$beta)
hist(output$beta)
plot(output$all.prevalence)

# visualize parasite niche data
ggplot(output$pniche.d) + 
  geom_line(aes(x=pEr, y=pPr.estab, color=factor(parasite.species))) +
  geom_point(aes(x=output$pE, y=-.01))

# visualizing transmission process
setwd("cartoon")
png(file = "cartoon%03d.png", width=700, height=300)
for (i in 1:dim(output$pstate)[1]){
  par(mfrow=c(1,2))
  head.title <- paste("t=", i, sep="")
  image(output$pstate[i,,], xlab="Host individual", ylab="Parasite species", 
        main=head.title, col=c("white", "red"))
  plot(output$parasite.richness[1:i], type="l", lty=1, col="red",
       xlab="Timestep", ylab="Parasite richness", 
       ylim=c(0, max(output$parasite.richness, na.rm=T)))
}
dev.off()
system("convert -delay 20 *.png cartoon.gif")
system("eog cartoon.gif")
file.remove(list.files(pattern=".png"))
setwd("/home/max/Documents/manuscripts/abm")

#########################################
## Host richness vs. symbiont richness ##
#########################################

A <- 81
N <- seq(1, 81, by=20)
P <- 200
timesteps <- 100
host.rich <- rep(NA, length(N))
symb.rich <- rep(NA, length(N))
for (i in N){
  output <- staticIBM(Y = 5, P=P, sig.p=10, gamma=0.1, network.pow = .3, 
                      pEmin = -500, pEmax = 500, pERmin=-500, pERmax=500, pI=.2, 
                      timesteps=100, A=A, N=i, cij=0)
  host.rich[i] <- i
  symb.rich[i] <- mean(output$parasite.richness[(.5*timesteps):timesteps])
  if (i == 1){
    plot(output$parasite.richness, type="l", col="red",
         xlab="Timestep", ylab="Parasite richness", 
         ylim=c(0, 30))
  } else {
    lines(output$parasite.richness, type="l", col="red")
  }
}

ggplot(output$pniche.d) + 
  geom_line(aes(x=pEr, y=pPr.estab, color=factor(parasite.species))) +
  geom_point(aes(x=output$pE, y=-.01))



# changing richness in parallel
require(doMC)
registerDoMC(cores=8)

parrich <- foreach(i=1:getDoParWorkers(), .combine=rbind) %dopar% {
  A <- 121
  N <- seq(1, 121, by=1)
  P <- 200
  timesteps <- 200
  host.rich <- rep(NA, length(N))
  symb.rich <- rep(NA, length(N))
  symb.rsd <- rep(NA, length(N))
  for (i in 1:length(N)){
    output <- staticIBM(Y = 5, P=P, sig.p=10, gamma=0.1, network.pow = .3, 
                        pEmin = -500, pEmax = 500, pERmin=-500, pERmax=500, pI=.2, 
                        timesteps=timesteps, A=A, N=N[i], cij=.01)
    host.rich[i] <- N[i]
    symb.rich[i] <- mean(output$parasite.richness[(.5*timesteps):timesteps])
    symb.rsd[i] <- sd(output$parasite.richness[(.5*timesteps):timesteps])
  }
  return(data.frame(host.rich, symb.rich, symb.rsd))
}

r1 <- ggplot(parrich, aes(x=host.rich, y=symb.rich)) + 
  theme_classic() + 
  geom_point(shape=1) +
  xlab("Host richness") + 
  ylab("Symbiont richness")

r2 <- ggplot(parrich, aes(x=host.rich, y=symb.rsd)) + 
  theme_classic() + 
  geom_point(shape=1) +
  xlab("Host richness") + 
  ylab(expression(paste(sigma, "(symbiont richness)", sep="")))

require(gridExtra)
grid.arrange(r1, r2, ncol=2)

# Changing how different species are, keeping richness constant
Emin <- -50
Emax <- 50
global.median <- median(c(Emin, Emax))
n.intervals <- 500
lower.limits <- seq(Emin, global.median - .5, length.out = n.intervals)
upper.limits <- seq(Emax, global.median, length.out=n.intervals)
host.het <- upper.limits - lower.limits
plot(x=c(lower.limits[1], upper.limits[1]), y=c(1,1), type="l", 
     ylim=c(0, n.intervals))
for (i in 2:n.intervals){
  lines(x=c(lower.limits[i], upper.limits[i]), y=c(i, i))
}

require(doMC)
registerDoMC(cores=8)
parsims <- foreach(i=1:getDoParWorkers(), .combine=rbind) %dopar% {
  A <- 101
  N <- 101
  P <- 100
  timesteps <- 100
  
  host.rich <- rep(NA, n.intervals)
  symb.rich <- rep(NA, n.intervals)
  heterogeneity <- rep(NA, n.intervals)
  for (i in 1:n.intervals){
    output <- staticIBM(Y = 5, P=P, sig.p=3, gamma=0.1, network.pow = .3, 
                        pEmin = lower.limits[i], pEmax = upper.limits[i], 
                        pERmin=-30, pERmax=30, pI=.01, 
                        timesteps=timesteps, A=A, N=N)
    host.rich[i] <- N
    heterogeneity[i] <- host.het[i]
    symb.rich[i] <- mean(output$parasite.richness[(.5*timesteps):timesteps])
  }
  return(data.frame(host.rich, symb.rich, heterogeneity))
}

ggplot(parsims, aes(x=heterogeneity, y=symb.rich)) + 
  theme_classic() + 
  geom_jitter(shape=1) + 
  xlab("Host heterogeneity") + 
  ylab("Equilibrium symbiont richness")

#write.csv(parsims, file="500intervals.csv")
#write.table(parsims, file = "500intervals.csv", sep = ",", col.names = T,
#            qmethod = "double")
df <- read.csv("500intervals.csv")

ggplot(df, aes(x=heterogeneity, y=symb.rich)) + 
  theme_classic() + 
  geom_point(shape=1) + 
  xlab("Host heterogeneity") + 
  ylab("E(steady state symbiont richness)") + 
  stat_smooth(method = "lm", formula = y ~ poly(x,2), color="red")


# just changing host richness, keeping range of host conditions equal
output <- staticIBM(Y = 5, P=3, sig.p=1, gamma=0.1, network.pow = 1, 
                    pEmin = -5, pEmax = 5, pERmin=-5, pERmax=5, pI=.2, 
                    timesteps=100, A=10, N=4, cij=0)

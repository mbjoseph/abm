# IBM analysis
source("com_init.R")
source("resolve.R")
source("rain_symbionts.R")
source("transmit_symbionts.R")
source("fullIBM.R")
output <- fullIBM(A=20, sig=10, sig.p=10, pI=1, N=20, P=10, timesteps=100, I=1)
str(output)

plot(output$host.richness, type="l", xlab="Timestep", ylab="Host richness")
lines(output$parasite.richness, type="l", lty=2, col="red",
      xlab="Timestep", ylab="Parasite richness")
legend("bottomright", legend=c("Hosts", "Symbionts"), col=c("black", "red"), 
       lty=c(1, 2))

# visualize niche data
ggplot(output$niche.d) + 
  geom_line(aes(x=E, y=Pr.estab, color=factor(species)))

# visualize parasite niche data
ggplot(output$pniche.d) + 
  geom_line(aes(x=pE, y=pPr.estab, color=factor(parasite.species)))

x <- -100:100
y <- -x^2
y <- scale(y)
par(mgp=c(1,1,0))
plot(x, y, "l", yaxt="n", xaxt="n",
     ylab="Species richness", xlab="Habitat heterogeneity", 
     main="Habitat area-heterogeneity trade-off")
z <- - y^2
z <- scale(z)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
z <- range01(z) ^ 2
require(scatterplot3d)
scatterplot3d(x, y, z, highlight.3d=T, type="h", 
              #main="Hypothetical richness relationship", 
              xlab="Habitat heterogeneity", 
              ylab="Host diveristy", 
              zlab="Symbiont diversity", 
              tick.marks=F, 
              label.tick.marks=F)
library(lattice)
cloud(z~x*y,
      distance=.1,
      type=c("h", "p"),
      main="Dual-trade-off diversity relationships",
      xlab="Habitat heterogeneity", 
      ylab="Hosts", 
      zlab="Symbionts", 
      shade=T, 
      pch=1,
      pch.col="black",
      aspect = c(1.2, 1))


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
par(mgp=c(1,1,0))
plot(x, y, "l", yaxt="n", xaxt="n", ylab="Species richness", xlab="Habitat heterogeneity", main="Habitat area-heterogeneity trade-off")

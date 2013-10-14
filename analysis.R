# IBM analysis
source("com_init.R")
source("fullIBM.R")
output <- fullIBM(sig.p=2, pI=0.01)
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

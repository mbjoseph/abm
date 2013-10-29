# static model analysis
# IBM analysis
source("static/static_init.R")
source("resolve.R")
source("rain_symbionts.R")
source("transmit_symbionts.R")
source("static/staticIBM.R")
output <- staticIBM()
str(output)

plot(output$parasite.richness, type="l", lty=2, col="red",
      xlab="Timestep", ylab="Parasite richness")

# visualize parasite niche data
ggplot(output$pniche.d) + 
  geom_line(aes(x=pE, y=pPr.estab, color=factor(parasite.species)))

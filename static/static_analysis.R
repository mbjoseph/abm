# static model analysis
# IBM analysis
setwd("/home/max/Documents/manuscripts/abm")
source("static/static_init.R")
source("resolve.R")
source("rain_symbionts.R")
source("static/transmit.R")
source("static/staticIBM.R")
output <- staticIBM(Y = 5, P=50, sig.p=1, gamma=0.1, network.pow = 1, 
                    pEmin = -5, pEmax = 5, pERmin=-5, pERmax=5, pI=.01, 
                    timesteps=100)
str(output)

plot(output$parasite.richness, type="l", lty=2, col="red",
      xlab="Timestep", ylab="Parasite richness", 
     ylim=c(0, max(c(output$beta, output$parasite.richness), na.rm=T)))
lines(output$beta)
hist(output$beta)

# visualize parasite niche data
ggplot(output$pniche.d) + 
  geom_line(aes(x=pE, y=pPr.estab, color=factor(parasite.species)))

# visualizing transmission process
setwd("cartoon")
png(file = "cartoon%03d.png", width=500, height=250)
for (i in 1:dim(output$pstate)[1]){
  head.title <- paste("t=", i, sep="")
  image(output$pstate[i,,], xlab="Host individual", ylab="Parasite species", 
        main=head.title)
}
dev.off()
system("convert -delay 20 *.png cartoon.gif")
system("eog cartoon.gif")
file.remove(list.files(pattern=".png"))
setwd("/home/max/Documents/manuscripts/abm")

# showing estimated transmission/contact rates between susceptibles and infected


# showing estimates of R0

# Post-JANUS analysis
setwd("continuous_time")
load("test.rdata")
str(xx)
ERvec <- rep(NA, length(xx))
rich <- rep(NA, length(xx))

for (i in 1:length(xx)){
  ERvec[i] <- xx[[i]]$ER
  rich[i] <- xx[[i]]$rich_bar
}

plot(ERvec, rich)

# Post-JANUS analysis
setwd("continuous_time")
system("rsync --update -raz --progress majo3748@login.rc.colorado.edu:/projects/majo3748/abm/ ~/Documents/manuscripts/abm/continuous_time/")
load("test.rdata")
#load("sig_s2.rdata")
str(xx[[1]])
ERvec <- rep(NA, length(xx))
rich <- rep(NA, length(xx))

for (i in 1:length(xx)){
  ERvec[i] <- xx[[i]]$ER
  rich[i] <- xx[[i]]$rich_bar
}

svg("figs/sigs2.svg", width=9, height=4)
op <- par(mfrow=c(1, 2), mai=c(1.25, 1.25, .75, .75))

# show intervals used
ERlow <- rep(NA, length(xx))
ERhigh <- rep(NA, length(xx))

for (i in 1:length(xx)){
  ERlow[i] <- xx[[i]]$Earray[, 1]
  ERhigh[i] <- xx[[i]]$Earray[, 2]
}

ERd <- data.frame(ERvec, ERlow, ERhigh)
ERd <- ERd[order(ERd$ERvec), ]

plot(x=rep(ERd$ERvec[1], 2), y=c(ERd$ERlow[1], ERd$ERhigh[1]), 
     xlim=range(ERvec), ylim=c(min(ERlow), max(ERhigh)), 
     type="l", 
     xlab="Host functional diversity", 
     ylab="Within-host environment range", 
     #xaxt="n"
     )
for (i in 2:nrow(ERd)){
  lines(x=rep(ERd$ERvec[i], 2), y=c(ERd$ERlow[i], ERd$ERhigh[i]))
}
par(mai=c(1.25, 1, .75, .75))
plot(ERvec, rich, 
     xlab="Host functional diversity", 
     ylab="Mean symbiont richness")

par(op)
dev.off()

# Calculate mean transmission rates within and among species
win_bar <- rep(NA, length(xx))
among_bar <- rep(NA, length(xx))
all_bar <- rep(NA, length(xx))

nsymb <- dim(xx[[1]]$trans_bar)[3]
nhosts <- dim(xx[[1]]$trans_bar)[1]

for (i in 1:length(xx)){
  # gather within species data
  withind <- rep(NA, nsymb)
  for (j in 1:nsymb){
    withind[j] <- sum(diag(xx[[i]]$trans_bar[, , j, ,])) / nhosts
  }
  win_bar[i] <- mean(withind)
  
  # do among species data
  amongd <- rep(NA, nsymb)
  for (j in 1:nsymb){
    submat <- xx[[i]]$trans_bar[, , j, ,]
    diag(submat) <- NA
    amongd[j] <- sum(submat, na.rm=T) / length(!is.na(submat))
  }
  among_bar[i] <- mean(amongd)
  
  # do all transmissions
  for (j in 1:nsymb){
    alld <- rep(NA, nsymb)
    for (j in 1:nsymb){
      submat <- xx[[i]]$trans_bar[, , j, ,]
      alld[j] <- sum(submat) / nhosts ^ 2
    }
    all_bar[i] <- mean(alld)
  }
}

# Characterize and plot mean prevalence for each symbiont X each iteration
prev <- array(dim=c(nsymb, length(xx)))
for (i in 1:nsymb){
  for (j in 1:length(xx)){
    prev.series <- xx[[j]]$s.ind[, i] / 100 # 100 cells
    prev[i, j] <- mean(prev.series)
  }
}

svg("figs/trans.svg", width=9, height=4)
op <- par(mfrow=c(1, 2))
plot(ERvec, all_bar, 
     xlab="Host functional diversity", 
     ylab="Mean transmission rate")

jitt <- .5
plot(x = ERvec + runif(length(xx), -jitt, jitt), 
     y = prev[1, ], 
     ylim=c(0, max(prev)), 
     col=ifelse(prev[1, ] == 0, "lightgray", "black"), 
     xlab="Host functional diversity", 
     ylab="Symbiont prevalence")
for (i in 2:nsymb){
  points(x = ERvec + runif(length(xx), -jitt, jitt), 
         y = prev[i, ], 
         col=ifelse(prev[i, ] == 0, "red", "blue"))
}
legend("topright", legend=c("Never invaded", "Successfully invaded"), 
       col=c("red", "blue"), pch=c(1, 1), bty="n", 
       cex=.8)
par(op)
dev.off()


str(xx[[1]])
plot(xx[[1]]$t, xx[[1]]$s.ind[, 1], type="l", ylim=c(0, max(xx[[1]]$s.ind)))
for (i in 1:length(xx)){
  for (j in 1:dim(xx[[1]]$s.ind)[2]){
    lines(xx[[i]]$t, xx[[i]]$s.ind[, j], col=j)
  }
}
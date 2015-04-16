# Post-JANUS analysis
library(ggplot2)
library(gridExtra)

system("rsync --update -raz --progress majo3748@login.rc.colorado.edu:/projects/majo3748/abm/ ~/Documents/manuscripts/abm/continuous_time/")

system("mv continuous_time/Mar* continuous_time/results")

# gather results
dir <- paste(getwd(), "/continuous_time/results", sep="")
data <- list.files(dir, pattern="Mar*")
source("continuous_time/mpi_f.R")
source("continuous_time/helpers.R")

# initialize vectors
ERvec <- c()
ERlow <- c()
ERhigh <- c()
rich <- c()
beta_d <- c()
sig_s <- c()
phi <- c()
mode <- c()
ncells <- c()

# read data
pindex <- 1
for (i in 1:length(data)){ # each node result
  d <- readRDS(paste("continuous_time/results/", data[i], sep=""))
  for (j in 1:length(d)){
    # catch errors (only occur at small niche widths...?)
    if (class(d[[j]]) == "try-error"){
      print(paste("try-error at node ", i, ", cpu ", j))
      next 
    } else { # otherwise, extract values
      ERvec <- c(ERvec, d[[j]]$ER)
      ERlow <- c(ERlow, d[[j]]$Earray[, 1])
      ERhigh <- c(ERhigh, d[[j]]$Earray[, 2])
      rich <- c(rich, d[[j]]$rich_bar)
      beta_d <- c(beta_d, d[[j]]$stored_pars$beta_d)
      if (is.null(d[[j]]$stored_pars$sig.s)) { # early iterations need fix
        d[[j]]$stored_pars$sig.s <- 1
      }
      sig_s <- c(sig_s, d[[j]]$stored_pars$sig.s)
      phi <- c(phi, d[[j]]$stored_pars$phi)
      mode <- c(mode, d[[j]]$stored_pars$mode)
      ncells <- c(ncells, d[[j]]$stored_pars$cells)
      pindex <- pindex + 1
    }
  }
  plot(d[[1]], main=paste0("Functional diversity = ", d[[1]]$ER))
  rm(d) # remove to conserve RAM
  print(paste("Completed", i, "of", length(data)))
}

# combine in data.frame
ERd <- data.frame(ERvec, ERlow, ERhigh, rich, beta_d, sig_s, phi, mode, ncells)
ERd <- ERd[order(ERd$ERvec), ]
library(plyr)
ERd$transmission <- mapvalues(ERd$mode, from=c("dens", "freq"), 
                              to=c("Density-dependent transmission", 
                                   "Frequency-dependent transmission"))

# plot effect of functional diversity for commensal symbionts
commensals <- subset(ERd, beta_d == 0)
p1 <- ggplot(commensals) + 
  geom_segment(aes(x=ERvec, xend=ERvec, y=ERlow, yend=ERhigh)) + 
  facet_grid(sig_s~mode) + 
  xlab("Host functional diversity") + 
  ylab("Local range of within-host environments") + 
  theme_bw()

p2 <- ggplot(commensals, aes(x=ERvec, y=rich, group=ncells)) + 
  facet_grid(.~transmission) +
  geom_point(shape=1)+ 
  xlab("Host functional diversity") + 
  ylab("Mean symbiont richness") + 
  theme_classic() + 
  stat_smooth(method="lm", formula = y ~ 0 + x + I(x^2))

grid.arrange(p2, p1, ncol=1)
#dev.off()

 
# make surface plot showing effect of functional div & infection eff
eff_d <- subset(ERd, ncells==1000)
library(scales)
ggplot(eff_d, aes(x=ERvec, y=beta_d)) + 
  geom_point(aes(color=rich), size=4, alpha=.5) + 
  scale_colour_gradient2(mid="blue", high="red") +
  theme_bw()


mod <- lm(rich ~ 0 + ERvec * beta_d + I(ERvec^2))
coefs <- coef(mod)

x1 <- 0:100
x2 <- seq(-1, 1, length.out=100)
X <- expand.grid(x1=x1, x2=x2)

mu <- with(X, coefs[1] * x1 + coefs[2] * x2 + 
             coefs[3] * x1^2 + coefs[4] * x1 * x2)

require(ggplot2)
d <- data.frame(mu=mu, x1=X$x1, x2=X$x2)
p1 <- ggplot(d, aes(x1, x2, z=mu)) + 
  theme_classic() +
  geom_tile(aes(fill=mu, color=mu)) +
  stat_contour(binwidth=1.5) +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       name="Average\nsymbiont\nrichness", midpoint=20) +
  scale_color_gradient2(low="blue", mid="white", high="red", 
                       name="Average\nsymbiont\nrichness", midpoint=20) +
  xlab("Functional diversity") + ylab("Effect on fitness") #+
#  ggtitle("Parasitism reduces the scaling of\n symbiont richness with host diversity")
p1

scatter3d(rich ~ ERvec + beta_d)


## plot diversity-richness relationship for symbionts with varying transmission rates
td <- subset(ERd, beta_d == 0)

ggplot(td, aes(x=ERvec, y=rich, color=factor(phi))) + 
  facet_wrap(~sig_s) +
  geom_point(shape=1)+ 
  xlab("Host functional diversity") + 
  ylab("Mean symbiont richness") + 
  theme_bw() + 
  stat_smooth(method="lm", formula = y ~  x + I(x^2))






# end applicable code (so far) for dynamic host communities
# Calculate mean transmission rates within and among species
win_bar <- rep(NA, length(res[[1]]))
among_bar <- rep(NA, length(res[[1]]))
all_bar <- rep(NA, length(res[[1]]))

nsymb <- dim(res[[1]][[1]]$trans_bar)[3]
nhosts <- dim(res[[1]][[1]]$trans_bar)[1]

for (i in 1:length(res[[1]])){
  # gather within species data
  withind <- rep(NA, nsymb)
  for (j in 1:nsymb){
    withind[j] <- sum(diag(res[[1]][[i]]$trans_bar[, , j, ,])) / nhosts
  }
  win_bar[i] <- mean(withind)
  
  # do among species data
  amongd <- rep(NA, nsymb)
  for (j in 1:nsymb){
    submat <- res[[1]][[i]]$trans_bar[, , j, ,]
    diag(submat) <- NA
    amongd[j] <- sum(submat, na.rm=T) / length(!is.na(submat))
  }
  among_bar[i] <- mean(amongd)
  
  # do all transmissions
  for (j in 1:nsymb){
    alld <- rep(NA, nsymb)
    for (j in 1:nsymb){
      submat <- res[[1]][[i]]$trans_bar[, , j, ,]
      alld[j] <- sum(submat) / nhosts ^ 2
    }
    all_bar[i] <- mean(alld)
  }
}

# Characterize and plot mean prevalence for each symbiont X each iteration
prev <- array(dim=c(nsymb, length(res[[1]])))
for (i in 1:nsymb){
  for (j in 1:length(res[[1]])){
    prev.series <- res[[1]][[j]]$s.ind[, i] / 100 # 100 cells
    prev[i, j] <- mean(prev.series)
  }
}

svg("figs/trans.svg", width=9, height=4)
op <- par(mfrow=c(1, 2))
plot(ERvec, all_bar, 
     xlab="Host functional diversity", 
     ylab="Mean transmission rate")

jitt <- .5
plot(x = ERvec + runif(length(res[[1]]), -jitt, jitt), 
     y = prev[1, ], 
     ylim=c(0, max(prev)), 
     col=ifelse(prev[1, ] == 0, "lightgray", "black"), 
     xlab="Host functional diversity", 
     ylab="Symbiont prevalence")
for (i in 2:nsymb){
  points(x = ERvec + runif(length(res[[1]]), -jitt, jitt), 
         y = prev[i, ], 
         col=ifelse(prev[i, ] == 0, "red", "blue"))
}
legend("topright", legend=c("Never invaded", "Successfully invaded"), 
       col=c("red", "blue"), pch=c(1, 1), bty="n", 
       cex=.8)
par(op)
dev.off()


str(res[[1]][[1]])
plot(res[[1]][[1]]$t, res[[1]][[1]]$s.ind[, 1], type="l", ylim=c(0, max(res[[1]][[1]]$s.ind)))
for (i in 1:length(res[[1]])){
  for (j in 1:dim(res[[1]][[1]]$s.ind)[2]){
    lines(res[[1]][[i]]$t, res[[1]][[i]]$s.ind[, j], col=j)
  }
}
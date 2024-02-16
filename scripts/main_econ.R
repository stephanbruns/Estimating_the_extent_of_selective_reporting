
## Estimating the extent of selective reporting: An application to economics

rm(list = ls(all = TRUE))                     #clear environment

## Required packages
library(doParallel)
library(foreach)
library(metafor)
library(openxlsx)
library(readxl)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(writexl)

## Paths for results and functions
path_results <- "G:/My Drive/InflatedSigInEco/results/Main/"
path_fcts <- "G:/My Drive/InflatedSigInEco/scripts/"

## Functions
source(paste(path_fcts,"econ_functions.r", sep=""))

## Load data
dat_final <- read_excel("G:/My Drive/InflatedSigInEco/data/data_final.xlsx")
dat_final$vi <- dat_final$se^2  #calculating variance
dim(dat_final)

## Exclusion of large t-values
dat_final$tvs <- dat_final$eff / dat_final$se
dat_final <- dat_final %>%
  dplyr::filter(abs(tvs) <= 20)
  #dplyr::filter(abs(tvs) <= 50)
  #dplyr::filter(abs(tvs) <= 100)
dim(dat_final)  #70399x13

## Check
zdat <- list()
jj <- 1
for (i in unique(dat_final$filecode)){
  zdat[[jj]] <- dat_final[which(dat_final$filecode==i), ]
  jj <- jj + 1
}

mss <- NULL
for (i in 1:dim(summary(zdat))[1]) { 
  mss[i] <- length(zdat[[i]]$se)
}

length(mss) #192
sum(mss)    #70399

## Estimating genuine effects for each research field/meta-analysis
## This produces the estimates of the genuine effects & the number of adequately powered studies

GE_wls  <- NULL;  GE_pp  <- NULL
GE_waap <- NULL;  n.waap <- NULL

for (i in 1:length(summary(zdat)[,1])) {	
  
  #wls
  wls <- sum(zdat[[i]]$eff / zdat[[i]]$se^2) / sum(1/zdat[[i]]$se^2)
  GE_wls <- c(GE_wls, rep(wls, length(zdat[[i]]$eff)))
  
  #waap
  se.crit <- abs(wls) / 2.8
  pos <- which(zdat[[i]]$se <= se.crit)
  n.waap[i] <- length(pos) / length(zdat[[i]]$se)
  
  #at least two adequately powered studies are needed
  if (length(pos) > 1) {
    GE_waap <- c(GE_waap, rep(sum(zdat[[i]]$eff[pos] / zdat[[i]]$se[pos]^2) / 
                sum(1/zdat[[i]]$se[pos]^2),	length(zdat[[i]]$eff)))
  } else {
    GE_waap <- c(GE_waap, rep(NA, length(zdat[[i]]$eff)))
  }
  
  #PET-PEESE (PET if alpha is not significant, PEESE if it is significant.)
  if (coefficients(summary(lm(zdat[[i]]$eff ~ zdat[[i]]$se, weights=(1/zdat[[i]]$se)^2)))[1,4] > 0.1) {
    GE_pp <- c(GE_pp, rep(coefficients(summary(lm(zdat[[i]]$eff ~ zdat[[i]]$se, 
              weights=(1/zdat[[i]]$se)^2)))[1,1], length(zdat[[i]]$eff)))
  } else {
    se2 <- zdat[[i]]$se^2
    GE_pp <-c(GE_pp, rep(coefficients(summary(lm(zdat[[i]]$eff ~ se2, 
             weights=(1/zdat[[i]]$se)^2)))[1,1],length(zdat[[i]]$eff)))
  }
}	

#head(data.frame(GE_wls, GE_pp, GE_waap)) #to check

## GE_wapp inserts GE_pp whenever there are not at least two adequately powered estimates.
GE_wapp <- GE_waap
for (i in 1:sum(mss)) {
  if (is.na(GE_wapp[i]) == TRUE) {
    GE_wapp[i] <- GE_pp[i]
  }
}

# Main case
case <- "GE_wapp"
GE_case <- GE_wapp

## Add estimated genuine effect to zdat[[i]]
dat.ideas <- do.call(rbind.data.frame, zdat)
ii <- 1
for (i in unique(dat.ideas$filecode)) {	
  zdat[[ii]]$GE <- GE_case[which(dat.ideas$filecode == i)]
  zdat[[ii]]$n.waap <- n.waap[[ii]]
  ii <- ii + 1
}
head(dat.ideas)
dim(dat.ideas)  #Run the above code to add GE & n.waap to the dataframe dat.ideas
#View(dat.ideas)

## Descriptive table

## --------------
## micro vs macro
## --------------

p <- 1
o <- 1
zdat.micro <- list()
zdat.macro <- list()
i.micro <- NULL
i.macro <- NULL
for (i in 1:length(summary(zdat)[,1])) {
  if(zdat[[i]]$field[1] == "Macro" | zdat[[i]]$field[1] == "macro") {
    zdat.macro[[p]] <- zdat[[i]]
    i.macro[p] <- i
    p <- p+1
  }
  if(zdat[[i]]$field[1] == "micro") {
    zdat.micro[[o]] <- zdat[[i]]
    i.micro[o] <- i
    o <- o+1
  }
}


## ----------
## exp vs obs
## ----------

p <- 1
o <- 1
zdat.obs <- list()
zdat.exp <- list()
i.obs <- NULL
i.exp <- NULL
for (i in 1:length(summary(zdat)[,1])) {
  if(zdat[[i]]$experiment[1] == 1) { 
    zdat.exp[[p]] <- zdat[[i]]
    i.exp[p] <- i
    p <- p+1
  }
  if(zdat[[i]]$experiment[1] == 0) {
    zdat.obs[[o]] <- zdat[[i]]
    i.obs[o] <- i
    o <- o+1
  }
}


## ---------------------
## power: 0 vs something
## ---------------------

p <- 1
o <- 1
zdat.pow  <- list()
zdat.npow <- list()
i.pow  <- NULL
i.npow <- NULL
for (i in 1:length(summary(zdat)[,1])) {
  if(zdat[[i]]$n.waap[1] > 0) { 
    zdat.pow[[p]] <- zdat[[i]]
    i.pow[p] <- i
    p <- p+1
  }
  if(zdat[[i]]$n.waap[1] == 0) {
    zdat.npow[[o]] <- zdat[[i]]
    i.npow[o] <- i
    o <- o+1
  }
}

length(zdat.pow)  
length(zdat.npow) 

d.tab <- matrix(ncol=8, nrow=7)

#Full
d.tab[1,1] <- length(zdat)
d.tab[1,2] <- sum(mss)
d.tab[1,3] <- round(mean(mss))
d.tab[1,4] <- min(mss)
d.tab[1,5] <- round(quantile(mss, probs=0.25))
d.tab[1,6] <- round(median(mss))
d.tab[1,7] <- round(quantile(mss, probs=0.75))
d.tab[1,8] <- max(mss)

#Micro
d.tab[2,1] <- length(zdat.micro)
d.tab[2,2] <- sum(mss[i.micro])
d.tab[2,3] <- round(mean(mss[i.micro]))
d.tab[2,4] <- min(mss[i.micro])
d.tab[2,5] <- round(quantile(mss[i.micro], probs=0.25))
d.tab[2,6] <- round(median(mss[i.micro]))
d.tab[2,7] <- round(quantile(mss[i.micro], probs=0.75))
d.tab[2,8] <- max(mss[i.micro])

#Macro
d.tab[3,1] <- length(zdat.macro)
d.tab[3,2] <- sum(mss[i.macro])
d.tab[3,3] <- round(mean(mss[i.macro]))
d.tab[3,4] <- min(mss[i.macro])
d.tab[3,5] <- round(quantile(mss[i.macro], probs=0.25))
d.tab[3,6] <- round(median(mss[i.macro]))
d.tab[3,7] <- round(quantile(mss[i.macro], probs=0.75))
d.tab[3,8] <- max(mss[i.macro])

#Exp
d.tab[4,1] <- length(zdat.exp)
d.tab[4,2] <- sum(mss[i.exp])
d.tab[4,3] <- round(mean(mss[i.exp]))
d.tab[4,4] <- min(mss[i.exp])
d.tab[4,5] <- round(quantile(mss[i.exp], probs=0.25))
d.tab[4,6] <- round(median(mss[i.exp]))
d.tab[4,7] <- round(quantile(mss[i.exp], probs=0.75))
d.tab[4,8] <- max(mss[i.exp])

#Obs
d.tab[5,1] <- length(zdat.obs)
d.tab[5,2] <- sum(mss[i.obs])
d.tab[5,3] <- round(mean(mss[i.obs]))
d.tab[5,4] <- min(mss[i.obs])
d.tab[5,5] <- round(quantile(mss[i.obs], probs=0.25))
d.tab[5,6] <- round(median(mss[i.obs]))
d.tab[5,7] <- round(quantile(mss[i.obs], probs=0.75))
d.tab[5,8] <- max(mss[i.obs])

#Pow
d.tab[6,1] <- length(zdat.pow)
d.tab[6,2] <- sum(mss[i.pow])
d.tab[6,3] <- round(mean(mss[i.pow]))
d.tab[6,4] <- min(mss[i.pow])
d.tab[6,5] <- round(quantile(mss[i.pow], probs=0.25))
d.tab[6,6] <- round(median(mss[i.pow]))
d.tab[6,7] <- round(quantile(mss[i.pow], probs=0.75))
d.tab[6,8] <- max(mss[i.pow])

#Npow
d.tab[7,1] <- length(zdat.npow)
d.tab[7,2] <- sum(mss[i.npow])
d.tab[7,3] <- round(mean(mss[i.npow]))
d.tab[7,4] <- min(mss[i.npow])
d.tab[7,5] <- round(quantile(mss[i.npow], probs=0.25))
d.tab[7,6] <- round(median(mss[i.npow]))
d.tab[7,7] <- round(quantile(mss[i.npow], probs=0.75))
d.tab[7,8] <- max(mss[i.npow])

rownames(d.tab) <- c("Total", "Microeconomics", "Macroeconomics","Experimental", 
                     "Observational",	"Share of APS > 0", "Share of APS = 0")
colnames(d.tab) <- c("N (Meta)", "N (Estimates)", "Mean", "Min", "Q25", "Q50", "Q75", "Max")#,

print(d.tab)
write.csv(d.tab,paste(path_results,case,"_desc.table.csv",sep=""))


## -----------------------------------------------------------------------------

## Grids of z-values & p-values
# z.grid is used for visualization and goes from 1 to 10 while
# p.grid is used for exact calculations and goes from z-values of 0 to Inf (absolute values).

# p.grid for main results table
p.grid.tab <- c(-Inf,
                qnorm(0.001/2, lower.tail = TRUE),
                qnorm(0.01/2,  lower.tail = TRUE),
                qnorm(0.05/2,  lower.tail = TRUE),
                qnorm(0.1/2,   lower.tail = TRUE),
                qnorm(0.2/2,   lower.tail = TRUE),
                qnorm(0.3/2,   lower.tail = TRUE),
                qnorm(0.4/2,   lower.tail = TRUE),
                qnorm(0.5/2,   lower.tail = TRUE),
                qnorm(0.6/2,   lower.tail = TRUE),
                qnorm(0.7/2,   lower.tail = TRUE),
                qnorm(0.8/2,   lower.tail = TRUE),
                qnorm(0.9/2,   lower.tail = TRUE),
                0,
                qnorm(0.9/2,   lower.tail = FALSE),
                qnorm(0.8/2,   lower.tail = FALSE),
                qnorm(0.7/2,   lower.tail = FALSE),
                qnorm(0.6/2,   lower.tail = FALSE),
                qnorm(0.5/2,   lower.tail = FALSE),
                qnorm(0.4/2,   lower.tail = FALSE),
                qnorm(0.3/2,   lower.tail = FALSE),
                qnorm(0.2/2,   lower.tail = FALSE),
                qnorm(0.1/2,   lower.tail = FALSE),
                qnorm(0.05/2,  lower.tail = FALSE),
                qnorm(0.01/2,  lower.tail = FALSE),
                qnorm(0.001/2, lower.tail = FALSE),
                Inf)

# We need to calculate for the orginal z-values also the frequencies
# and we can do this directly for the two-sided test grid
p.grid.tab2 <- c(0, 
                 qnorm(0.9/2,   lower.tail = FALSE),
                 qnorm(0.8/2,   lower.tail = FALSE),
                 qnorm(0.7/2,   lower.tail = FALSE),
                 qnorm(0.6/2,   lower.tail = FALSE),
                 qnorm(0.5/2,   lower.tail = FALSE),
                 qnorm(0.4/2,   lower.tail = FALSE),
                 qnorm(0.3/2,   lower.tail = FALSE),
                 qnorm(0.2/2,   lower.tail = FALSE),
                 qnorm(0.1/2,   lower.tail = FALSE),
                 qnorm(0.05/2,  lower.tail = FALSE),
                 qnorm(0.01/2,  lower.tail = FALSE),
                 qnorm(0.001/2, lower.tail = FALSE),
                 Inf)				

# p grid for plot
p.grid.plot <- qnorm(seq(1, 0, -0.005), lower.tail = FALSE)
p.grid.plot2 <-p.grid.plot[which(p.grid.plot >=0)] # this is for actual plotting of two-sided tests

# z.grid for plot
z.grid.plot  <- seq(-10.25, 10.25, 0.1025)
z.grid.plot2 <- seq(0, 10.25, 0.1025) # this is for actual plotting of absolute values


## Prepare original data for plotting and p-value table

## --------------
## micro vs macro
## --------------

p <- 1
o <- 1
zdat.micro <- list()
zdat.macro <- list()
i.micro <- NULL
i.macro <- NULL
for (i in 1:length(summary(zdat)[,1])) {
  if(zdat[[i]]$field[1] == "Macro" | zdat[[i]]$field[1] == "macro") {
    zdat.macro[[p]] <- zdat[[i]]
    i.macro[p] <- i
    p <- p+1
  }
  if(zdat[[i]]$field[1] == "micro") {
    zdat.micro[[o]] <- zdat[[i]]
    i.micro[o] <- i
    o <- o+1
  }
}

length(zdat.micro) # 131
length(zdat.macro) # 61


## ----------
## exp vs obs
## ----------

p <- 1
o <- 1
zdat.obs <- list()
zdat.exp <- list()
i.obs <- NULL
i.exp <- NULL
for (i in 1:length(summary(zdat)[,1])) {
  if(zdat[[i]]$experiment[1] == 1) { 
    zdat.exp[[p]] <- zdat[[i]]
    i.exp[p] <- i
    p <- p+1
  }
  if(zdat[[i]]$experiment[1] == 0) {
    zdat.obs[[o]] <- zdat[[i]]
    i.obs[o] <- i
    o <- o+1
  }
}

length(zdat.exp) # 30
length(zdat.obs) # 162

## ---------------------
## power: 0 vs something
## ---------------------

p <- 1
o <- 1
zdat.pow  <- list()
zdat.npow <- list()
i.pow  <- NULL
i.npow <- NULL
for (i in 1:length(summary(zdat)[,1])) {
  if(zdat[[i]]$n.waap[1] > 0) { 
    zdat.pow[[p]] <- zdat[[i]]
    i.pow[p] <- i
    p <- p+1
  }
  if(zdat[[i]]$n.waap[1] == 0) {
    zdat.npow[[o]] <- zdat[[i]]
    i.npow[o] <- i
    o <- o+1
  }
}

length(zdat.pow)  # 147
length(zdat.npow) # 45


fac <- do.call(rbind.data.frame, zdat)
facz <- abs(fac$eff / fac$se)

fac <- do.call(rbind.data.frame, zdat.micro)
facz.micro <- abs(fac$eff / fac$se)

fac <- do.call(rbind.data.frame, zdat.macro)
facz.macro <- abs(fac$eff / fac$se)

fac <- do.call(rbind.data.frame, zdat.exp)
facz.exp <- abs(fac$eff / fac$se)

fac <- do.call(rbind.data.frame, zdat.obs)
facz.obs <- abs(fac$eff / fac$se)

fac <- do.call(rbind.data.frame, zdat.pow)
facz.pow <- abs(fac$eff / fac$se)

fac <- do.call(rbind.data.frame, zdat.npow)
facz.npow <- abs(fac$eff / fac$se)

## -----------------------------------------------------------------------------

#Full
z.orig <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig[a] <- length(which(facz >= z.grid.plot2[a] & facz <= z.grid.plot2[a+1]))
}

p.orig.tab <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab[a] <- length(which(facz >= p.grid.tab2[a] & facz <= p.grid.tab2[a+1]))
}

p.orig.plot <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot[a] <- length(which(facz >= p.grid.plot2[a] & facz <= p.grid.plot2[a+1]))
}

#Micro
z.orig.micro <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.micro[a] <- length(which(facz.micro >= z.grid.plot2[a] & facz.micro <= z.grid.plot2[a+1]))
}

p.orig.tab.micro <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.micro[a] <- length(which(facz.micro >= p.grid.tab2[a] & facz.micro <= p.grid.tab2[a+1]))
}

p.orig.plot.micro <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.micro[a] <- length(which(facz.micro >= p.grid.plot2[a] & facz.micro <= p.grid.plot2[a+1]))
}

#Macro
z.orig.macro <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.macro[a] <- length(which(facz.macro >= z.grid.plot2[a] & facz.macro <= z.grid.plot2[a+1]))
}

p.orig.tab.macro <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.macro[a] <- length(which(facz.macro >= p.grid.tab2[a] & facz.macro <= p.grid.tab2[a+1]))
}

p.orig.plot.macro <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.macro[a] <- length(which(facz.macro >= p.grid.plot2[a] & facz.macro <= p.grid.plot2[a+1]))
}

#Exp
z.orig.exp <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.exp[a] <- length(which(facz.exp >= z.grid.plot2[a] & facz.exp <= z.grid.plot2[a+1]))
}

p.orig.tab.exp <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.exp[a] <- length(which(facz.exp >= p.grid.tab2[a] & facz.exp <= p.grid.tab2[a+1]))
}

p.orig.plot.exp <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.exp[a] <- length(which(facz.exp >= p.grid.plot2[a] & facz.exp <= p.grid.plot2[a+1]))
}

#Obs
z.orig.obs <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.obs[a] <- length(which(facz.obs >= z.grid.plot2[a] & facz.obs <= z.grid.plot2[a+1]))
}

p.orig.tab.obs <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.obs[a] <- length(which(facz.obs >= p.grid.tab2[a] & facz.obs <= p.grid.tab2[a+1]))
}

p.orig.plot.obs <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.obs[a] <- length(which(facz.obs >= p.grid.plot2[a] & facz.obs <= p.grid.plot2[a+1]))
}

#Pow
z.orig.pow <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.pow[a] <- length(which(facz.pow >= z.grid.plot2[a] & facz.pow <= z.grid.plot2[a+1]))
}

p.orig.tab.pow <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.pow[a] <- length(which(facz.pow >= p.grid.tab2[a] & facz.pow <= p.grid.tab2[a+1]))
}

p.orig.plot.pow <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.pow[a] <- length(which(facz.pow >= p.grid.plot2[a] & facz.pow <= p.grid.plot2[a+1]))
}

#Npow
z.orig.npow <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.npow[a] <- length(which(facz.npow >= z.grid.plot2[a] & facz.npow <= z.grid.plot2[a+1]))
}

p.orig.tab.npow <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.npow[a] <- length(which(facz.npow >= p.grid.tab2[a] & facz.npow <= p.grid.tab2[a+1]))
}

p.orig.plot.npow <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.npow[a] <- length(which(facz.npow >= p.grid.plot2[a] & facz.npow <= p.grid.plot2[a+1]))
}

## -----------------------------------------------------------------------------

## -----------
## Full sample
## -----------

cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot <- cf(dat=zdat, z.grid=z.grid.plot)
stopCluster(cl)
save(z.plot,file=paste(path_results, case, "_z_plot.RDATA",sep=""))

clu <- unique(dat.ideas$filecode) #cluster id
#itero <- 1000       # for final result
itero <- 500

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.ci <- cf.ci.cluster(dat=zdat,z.grid=z.grid.plot,iters=itero,cluster=clu)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) # about 12 hrs (for 1000 iterations)
save(z.plot.ci,file=paste(path_results,case,"_z_plot_ci.RDATA",sep=""))

## p.grid.tab (full)
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab <- cf(dat=zdat, z.grid=p.grid.tab)
stopCluster(cl)
save(p.tab, file=paste(path_results, case, "_p_tab.RDATA", sep=""))

cl <- makeCluster(7) 
registerDoParallel(cl)
s.time <- Sys.time()
p.tab.ci <- cf.ci.cluster(dat=zdat,z.grid=p.grid.tab,iters=itero,cluster=clu)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) # about 1.6 hrs.
save(p.tab.ci,file=paste(path_results, case, "_p_tab_ci.RDATA",sep=""))

## --------------
## micro vs macro
## --------------

## Micro
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.micro <- cf(dat=zdat.micro, z.grid=z.grid.plot)
stopCluster(cl)
save(z.plot.micro, file=paste(path_results,case, "_z_plot.micro.RDATA", sep=""))


s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.ci.micro <- cf.ci.cluster(dat=zdat.micro,z.grid=z.grid.plot,iters=itero,cluster=clu[i.micro])
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) # about 3.78 hrs
save(z.plot.ci.micro,file=paste(path_results,case,"_z_plot_ci.micro.RDATA", sep=""))

# p.grid.tab
#about 1.5 hours with 1000 iterations (7 cores)
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.micro <- cf(dat=zdat.micro, z.grid=p.grid.tab)
stopCluster(cl)
save(p.tab.micro, file=paste(path_results, case, "_p_tab.micro.RDATA", sep=""))

cl <- makeCluster(7) 
registerDoParallel(cl)
s.time <- Sys.time()
p.tab.ci.micro <- cf.ci.cluster(dat=zdat.micro,z.grid=p.grid.tab,iters=itero,cluster=clu[i.micro])
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) # about 35 mins.
save(p.tab.ci.micro,file=paste(path_results, case,"_p_tab_ci.micro.RDATA",sep=""))

## Macro
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.macro <- cf(dat=zdat.macro, z.grid=z.grid.plot)
stopCluster(cl)
save(z.plot.macro,file=paste(path_results,case,"_z_plot.macro.RDATA",sep=""))

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.ci.macro <- cf.ci.cluster(dat=zdat.macro, z.grid=z.grid.plot,iters=itero, cluster=clu[i.macro])
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) # about 7:00 hrs for 1000 iterations
save(z.plot.ci.macro,file=paste(path_results,case,"_z_plot_ci.macro.RDATA",sep=""))

# p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.macro <- cf(dat=zdat.macro, z.grid=p.grid.tab)
stopCluster(cl)
save(p.tab.macro,file=paste(path_results,case,"_p_tab.macro.RDATA", sep=""))

s.time <- Sys.time()
cl<-makeCluster(7) 
registerDoParallel(cl)
p.tab.ci.macro <- cf.ci.cluster(dat=zdat.macro,z.grid=p.grid.tab, iters=itero, cluster=clu[i.macro])
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) # about 1 hr for 1000 iterations
save(p.tab.ci.macro,file=paste(path_results,case,"_p_tab_ci.macro.RDATA",sep=""))


## ----------
## exp vs obs
## ----------

## Exp
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.exp <- cf(dat=zdat.exp, z.grid=z.grid.plot)
stopCluster(cl)
save(z.plot.exp,file=paste(path_results,case,"_z_plot.exp.RDATA",sep=""))

itero <- 1000
s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.ci.exp <- cf.ci.cluster(dat=zdat.exp,z.grid=z.grid.plot,iters=itero,cluster=clu[i.exp])
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) # about 20 mins
save(z.plot.ci.exp, file=paste(path_results,case, "_z_plot_ci.exp.RDATA", sep=""))


# p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.exp <- cf(dat=zdat.exp, z.grid=p.grid.tab)
stopCluster(cl)
save(p.tab.exp,file=paste(path_results,case,"_p_tab.exp.RDATA", sep=""))

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.ci.exp <- cf.ci.cluster(dat=zdat.exp,z.grid=p.grid.tab,iters=itero, cluster=clu[i.exp])
stopCluster(cl)
save(p.tab.ci.exp,file=paste(path_results,case,"_p_tab_ci.exp.RDATA",sep=""))
e.time <- Sys.time()
print(e.time - s.time) # about 4 mins

## Obs
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.obs <- cf(dat=zdat.obs, z.grid=z.grid.plot)
stopCluster(cl)
save(z.plot.obs,file=paste(path_results,case,"_z_plot.obs.RDATA",sep=""))

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.ci.obs <- cf.ci.cluster(dat=zdat.obs,z.grid=z.grid.plot,iters=itero, cluster=clu[i.obs])
stopCluster(cl)
save(z.plot.ci.obs,file=paste(path_results,case,"_z_plot_ci.obs.RDATA",sep=""))
e.time <- Sys.time()
print(e.time - s.time)

# p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.obs <- cf(dat=zdat.obs, z.grid=p.grid.tab)
stopCluster(cl)
save(p.tab.obs,file=paste(path_results,case,"_p_tab.obs.RDATA", sep=""))

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.ci.obs <- cf.ci.cluster(dat=zdat.obs,z.grid=p.grid.tab,iters=itero, cluster=clu[i.obs])
stopCluster(cl)
save(p.tab.ci.obs, file=paste(path_results, case, "_p_tab_ci.obs.RDATA", sep=""))
e.time <- Sys.time()
print(e.time - s.time) # about 1.6 hrs (7 cores)


## ---------------------
## power: 0 vs something
## ---------------------

## Pow
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.pow <- cf(dat=zdat.pow, z.grid=z.grid.plot)
stopCluster(cl)
save(z.plot.pow,file=paste(path_results,case,"_z_plot.pow.RDATA",sep=""))

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.ci.pow <- cf.ci.cluster(dat=zdat.pow,z.grid=z.grid.plot,iters=itero, cluster=clu[i.pow])
stopCluster(cl)
save(z.plot.ci.pow,file=paste(path_results,case,"_z_plot_ci.pow.RDATA",sep=""))
e.time <- Sys.time()
print(e.time - s.time) # about86 hrs (7 cores)

# p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.pow <- cf(dat=zdat.pow, z.grid=p.grid.tab)
stopCluster(cl)
save(p.tab.pow,file=paste(path_results,case,"_p_tab.pow.RDATA", sep=""))

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.ci.pow <- cf.ci.cluster(dat=zdat.pow,z.grid=p.grid.tab,iters=itero,cluster=clu[i.pow])
stopCluster(cl)
save(p.tab.ci.pow,file=paste(path_results,case,"_p_tab_ci.pow.RDATA",sep=""))
e.time <- Sys.time()
print(e.time - s.time) # about 1 hr (7 cores)


## Npow
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.npow <- cf(dat=zdat.npow, z.grid=z.grid.plot)
stopCluster(cl)
save(z.plot.npow,file=paste(path_results,case,"_z_plot.npow.RDATA",sep=""))

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.ci.npow <- cf.ci.cluster(dat=zdat.npow,z.grid=z.grid.plot,iters=itero, cluster=clu[i.npow])
stopCluster(cl)
save(z.plot.ci.npow,file=paste(path_results,case,"_z_plot_ci.npow.RDATA",sep=""))
e.time <- Sys.time()
print(e.time - s.time) # about 4 hrs (7 cores)

# p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.npow <- cf(dat=zdat.npow, z.grid=p.grid.tab)
stopCluster(cl)
save(p.tab.npow,file=paste(path_results,case,"_p_tab.npow.RDATA",sep=""))

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.ci.npow <- cf.ci.cluster(dat=zdat.npow,z.grid=p.grid.tab,iters=itero, cluster=clu[i.npow])
stopCluster(cl)
save(p.tab.ci.npow,file=paste(path_results,case,"_p_tab_ci.npow.RDATA",sep=""))
e.time <- Sys.time()
print(e.time - s.time)  # about 40 mins (7 cores)


## ----------
## Main graph
## ----------

load(paste(path_results, case, "_z_plot.RDATA", sep=""))
load(paste(path_results, case, "_z_plot_ci.RDATA", sep=""))

## Figure 1
pdf(paste(path_results,case,"_Fig1.pdf",sep=""),width=14,height=8)
# z.grid
xs <- z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2
N <- sum(p.orig.plot) #be careful z.orig goes until 10 and is smaller than the full sample size!,
# but the bootstrapping of the confidence interval is based on the entire sample
# and given in shares. So we need to divide here by the full sample to be consistent.

q025 <- apply(z.plot.ci[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(z.plot.ci[[1]], 2, quantile, probs=c(0.975))

matplot(xs, q975, t="l", lty=3, col=4, xlab="|z|-value", ylab="Frequency", xaxt="n", xlim=c(0,8))	
matplot(xs, z.plot/N, t="o", lty=1, cex=0.8, pch=20, col=4, 
        xlab="|z|-value", ylab="Frequency", xaxt="n", add=T)	
matplot(xs, q025, t="l", lty=3, col=4, add = T)	

axis(1, c(1.64, 1.96, 2.58))
axis(1, c(0, 5, 8))
#axis(1, seq(0,10,1), label = FALSE)
matplot(xs, z.orig/N , t="o", pch=20, cex=0.8, lty=2, add=T)
abline(v=qnorm(0.10/2, lower.tail = FALSE), col =2, lty=2)
abline(v=qnorm(0.05/2, lower.tail = FALSE), col =2, lty=2)
abline(v=qnorm(0.01/2, lower.tail = FALSE), col =2, lty=2)
dev.off()


## Figure 2
load(paste(path_results, case, "_z_plot.micro.RDATA", sep=""))
load(paste(path_results, case, "_z_plot_ci.micro.RDATA", sep=""))
load(paste(path_results, case, "_z_plot.macro.RDATA", sep=""))
load(paste(path_results, case, "_z_plot_ci.macro.RDATA", sep=""))
load(paste(path_results, case, "_z_plot.exp.RDATA", sep=""))
load(paste(path_results, case, "_z_plot_ci.exp.RDATA", sep=""))
load(paste(path_results, case, "_z_plot.obs.RDATA", sep=""))
load(paste(path_results, case, "_z_plot_ci.obs.RDATA", sep=""))
load(paste(path_results, case, "_z_plot.pow.RDATA", sep=""))
load(paste(path_results, case, "_z_plot_ci.pow.RDATA", sep=""))
load(paste(path_results, case, "_z_plot.npow.RDATA", sep=""))
load(paste(path_results, case, "_z_plot_ci.npow.RDATA", sep=""))

pdf(paste(path_results,case,"_Fig2.pdf",sep=""),width=14,height=12)
par(mfrow=c(3,2))
# micro
xs <- z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2
N.micro <- sum(p.orig.plot.micro) 
q025 <- apply(z.plot.ci.micro[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(z.plot.ci.micro[[1]], 2, quantile, probs=c(0.975))
matplot(xs, q975, t="l", lty=3, col=4, xlab="|z|-value", ylab="Frequency", xaxt="n", xlim=c(0,8), 
        ylim=c(0,0.08), main="Microeconomics")	
matplot(xs, z.plot.micro/N.micro, t="o", lty=1, cex=0.8, pch=20, col=4, 
        xlab="|z|-value", ylab="Frequency", xaxt="n", add=T)	
matplot(xs, q025, t="l", lty=3, col=4, add = T)	
axis(1, c(1.64, 1.96, 2.58))
axis(1, c(0, 5, 8))
matplot(xs, z.orig.micro/N.micro , t="o", pch=20, cex=0.8, lty=2, add=T)
abline(v=qnorm(0.10/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.05/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.01/2, lower.tail = FALSE), col=2, lty=2)

#macro
xs <- z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2
N.macro <- sum(p.orig.plot.macro) 
q025 <- apply(z.plot.ci.macro[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(z.plot.ci.macro[[1]], 2, quantile, probs=c(0.975))
matplot(xs, q975, t="l", lty=3, col=4, xlab="|z|-value", ylab="Frequency", xaxt="n", xlim=c(0,8), 
        ylim=c(0,0.08), main="Macroeconomics")	
matplot(xs, z.plot.macro/N.macro, t="o", lty=1, cex=0.8, pch=20, col=4, 
        xlab="|z|-value", ylab="Frequency", xaxt="n", add=T)	
matplot(xs, q025, t="l", lty=3, col=4, add = T)	
axis(1, c(1.64, 1.96, 2.58))
axis(1, c(0, 5, 8))
matplot(xs, z.orig.macro/N.macro , t="o", pch=20, cex=0.8, lty=3, add=T)
abline(v=qnorm(0.1/2, lower.tail = FALSE),  col=2, lty=2)
abline(v=qnorm(0.05/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.01/2, lower.tail = FALSE), col=2, lty=2)

# exp
xs <- z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2
N.exp <- sum(p.orig.plot.exp) 
q025 <- apply(z.plot.ci.exp[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(z.plot.ci.exp[[1]], 2, quantile, probs=c(0.975))
matplot(xs, q975, t="l", lty=3, col=4, xlab="|z|-value", ylab="Frequency", xaxt="n", xlim=c(0,8), 
        ylim=c(0,0.08), main="Experimental")	
matplot(xs, z.plot.exp/N.exp, t="o", lty=1, cex=0.8, pch=20, col=4, 
        xlab="|z|-value", ylab="Frequency", xaxt="n", add=T)	
matplot(xs, q025, t="l", lty=3, col=4, add = T)	
axis(1, c(1.64, 1.96, 2.58))
axis(1, c(0, 5, 8))
matplot(xs, z.orig.exp/N.exp , t="o", pch=20, cex=0.8, lty=2, add=T)
abline(v=qnorm(0.10/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.05/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.01/2, lower.tail = FALSE), col=2, lty=2)

# obs
xs <- z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2
N.obs <- sum(p.orig.plot.obs) 
q025 <- apply(z.plot.ci.obs[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(z.plot.ci.obs[[1]], 2, quantile, probs=c(0.975))
matplot(xs, q975, t="l", lty=3, col=4, xlab="|z|-value", ylab="Frequency", xaxt="n", xlim=c(0,8), 
        ylim=c(0,0.08), main="Observational")	
matplot(xs, z.plot.obs/N.obs, t="o", lty=1, cex=0.8, pch=20, col=4, 
        xlab="|z|-value", ylab="Frequency", xaxt="n", add=T)	
matplot(xs, q025, t="l", lty=3, col=4, add = T)	
axis(1, c(1.64, 1.96, 2.58))
axis(1, c(0, 5, 8))
matplot(xs, z.orig.obs/N.obs , t="o", pch=20, cex=0.8, lty=3, add=T)
abline(v=qnorm(0.10/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.05/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.01/2, lower.tail = FALSE), col=2, lty=2)

# power
xs <- z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2
N.pow <- sum(p.orig.plot.pow) 
q025 <- apply(z.plot.ci.pow[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(z.plot.ci.pow[[1]], 2, quantile, probs=c(0.975))
matplot(xs, q975, t="l", lty=3, col=4, xlab="|z|-value", ylab="Frequency", xaxt="n", xlim=c(0,8), 
        ylim=c(0,0.08), main="Share of APS > 0")	
matplot(xs, z.plot.pow/N.pow, t="o", lty=1, cex=0.8, pch=20, col=4, 
        xlab="|z|-value", ylab="Frequency", xaxt="n", add=T)	
matplot(xs, q025, t="l", lty=3, col=4, add = T)	
axis(1, c(1.64, 1.96, 2.58))
axis(1, c(0, 5, 8))
matplot(xs, z.orig.pow/N.pow , t="o", pch=20, cex=0.8, lty=2, add=T)
abline(v=qnorm(0.10/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.05/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.01/2, lower.tail = FALSE), col=2, lty=2)

# npow
xs <- z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2
N.npow <- sum(p.orig.plot.npow) 
q025 <- apply(z.plot.ci.npow[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(z.plot.ci.npow[[1]], 2, quantile, probs=c(0.975))
matplot(xs, q975, t="l", lty=3, col=4, xlab="|z|-value", ylab="Frequency", xaxt="n", xlim=c(0,8), 
        ylim=c(0,0.08), main="Share of APS = 0")	
matplot(xs, z.plot.npow/N.npow, t="o", lty=1, cex=0.8, pch=20, col=4, 
        xlab="|z|-value", ylab="Frequency", xaxt="n", add=T)	
matplot(xs, q025, t="l", lty=3, col=4, add = T)	
axis(1, c(1.64, 1.96, 2.58))
axis(1, c(0, 5, 8))
matplot(xs, z.orig.npow/N.npow , t="o", pch=20, cex=0.8, lty=3, add=T)
abline(v=qnorm(0.10/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.05/2, lower.tail = FALSE), col=2, lty=2)
abline(v=qnorm(0.01/2, lower.tail = FALSE), col=2, lty=2)
dev.off()

## --------------------
## Main table (Table 2) 
## --------------------

# p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab <- cf(dat=zdat, z.grid=p.grid.tab)
stopCluster(cl)
save(p.tab,file=paste(path_results,case,"_p_tab.RDATA",sep=""))

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.ci <- cf.ci.cluster(dat=zdat, z.grid=p.grid.tab,iters=itero, cluster=clu)
stopCluster(cl)
save(p.tab.ci,file=paste(path_results,case,"_p_tab_ci.RDATA",sep=""))
e.time <- Sys.time()
print(e.time - s.time) # about 1.6 hrs (7 cores)

load(paste(path_results, case, "_p_tab.RDATA", sep=""))
load(paste(path_results, case, "_p_tab_ci.RDATA", sep=""))

## p-value table (difference in frequencies between f & cf)
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "95% CI")

N <- sum(p.orig.tab)
p.table[,1] <- round((p.orig.tab - p.tab) / N,  3)
q025 <- apply(matrix(p.orig.tab/N, nrow=nrow(p.tab.ci[[1]]), ncol=ncol(p.tab.ci[[1]]), 
                     byrow=TRUE) - p.tab.ci[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab/N, nrow=nrow(p.tab.ci[[1]]), ncol=ncol(p.tab.ci[[1]]), 
                     byrow=TRUE) - p.tab.ci[[1]], 2, quantile, probs=c(0.975))
p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

# share of inflated p-values in total p-values (10% level)
infp    <- round(sum((p.orig.tab - p.tab)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab/N)[10:13]) - infp.bs, 
                            probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab/N)[10:13]) - infp.bs, 
                            probs=c(0.975)),3)
p.table <- rbind(p.table,c(infp,paste("[",infp.q025,", ",infp.q975, "]", sep="")))

# share of inflated p-values in total p-values (5% level)
infp    <- round(sum((p.orig.tab - p.tab)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab/N)[11:13]) - infp.bs, 
                            probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab/N)[11:13]) - infp.bs, 
                            probs=c(0.975)),3)
p.table <- rbind(p.table,c(infp,paste("[",infp.q025,", ",infp.q975, "]",sep="")))

# share of inflated p-values in sig. p-values (10% level)
shp.10      <- round(sum((p.orig.tab - p.tab)[10:13] / sum(p.orig.tab[10:13])), 3)
infp.bs.10  <- apply(p.tab.ci[[2]][,10:13], 1, sum)
N.sig.10    <- sum((p.orig.tab)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab/N.sig.10)[10:13]) - infp.bs.10), 
                              probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab/N.sig.10)[10:13]) - infp.bs.10), 
                              probs=c(0.975)),3)
p.table <- rbind(p.table,c(shp.10,paste("[",shp.q025.10,", ",shp.q975.10, "]",sep="")))

# share of inflated p-values in sig. p-values (5% level)
shp.5      <- round(sum((p.orig.tab - p.tab)[11:13] / sum(p.orig.tab[11:13])),  3)
infp.bs.5  <- apply(p.tab.ci[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5    <- sum((p.orig.tab)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab/N.sig.5)[11:13]) - infp.bs.5), 
                             probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab/N.sig.5)[11:13]) - infp.bs.5), 
                             probs=c(0.975)),3)
p.table <- rbind(p.table,c(shp.5,paste("[",shp.q025.5, ", ",shp.q975.5,"]", sep="")))

## sample size
p.table <- rbind(p.table, c(length(summary(zdat)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("p > 0.9", "0.9 > p > 0.8", "0.8 > p > 0.7", 
                       "0.7 > p > 0.6","0.6 > p > 0.5", 
                       "0.5 > p > 0.4","0.4 > p > 0.3", 
                       "0.3 > p > 0.2","0.2 > p > 0.1", 
                       "0.1 > p > 0.05","0.05 > p > 0.01", 
                       "0.01 > p > 0.001","0.001 > p", 
                       "IS_{0.10}","IS_{0.05}",
                       "IS_{0.10}","IS_{0.05}",
                       "No. of meta-analysis", "No. of tests")
print(p.table)
write.csv(p.table, paste(path_results, case, "_p.table.csv", sep=""))

## ----------------------
## Tables in the appendix
## ----------------------

## -----
## micro
## -----

load(paste(path_results, case, "_p_tab.micro.RDATA", sep=""))
load(paste(path_results, case, "_p_tab_ci.micro.RDATA", sep=""))

p.table.fin <- NULL
#p-value table (difference in frequencies between factual and cf)
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")
N <- sum(p.orig.tab.micro)
p.table[,1] <- round((p.orig.tab.micro - p.tab.micro) / N,  3)
q025 <- apply(matrix(p.orig.tab.micro/N, nrow=nrow(p.tab.ci.micro[[1]]), ncol=ncol(p.tab.ci.micro[[1]]), byrow=TRUE) - p.tab.ci.micro[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.micro/N, nrow=nrow(p.tab.ci.micro[[1]]), ncol=ncol(p.tab.ci.micro[[1]]), byrow=TRUE) - p.tab.ci.micro[[1]], 2, quantile, probs=c(0.975))
p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.micro - p.tab.micro)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.micro[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.micro/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.micro/N)[10:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.micro - p.tab.micro)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.micro[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.micro/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.micro/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.micro - p.tab.micro)[10:13] / sum(p.orig.tab.micro[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.micro[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.micro)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.micro/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.micro/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.micro - p.tab.micro)[11:13] / sum(p.orig.tab.micro[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.micro[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.micro)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.micro/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.micro/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(zdat.micro)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))
p.table.fin <- cbind(p.table.fin,p.table)

## -----
## macro
## -----

load(paste(path_results, case, "_p_tab.macro.RDATA", sep=""))
load(paste(path_results, case, "_p_tab_ci.macro.RDATA", sep=""))

#p-value table (difference in frequencies between factual and cf)
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")
N <- sum(p.orig.tab.macro)
p.table[,1] <- round((p.orig.tab.macro - p.tab.macro) / N,  3)
q025 <- apply(matrix(p.orig.tab.macro/N, nrow=nrow(p.tab.ci.macro[[1]]), ncol=ncol(p.tab.ci.macro[[1]]), byrow=TRUE) - p.tab.ci.macro[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.macro/N, nrow=nrow(p.tab.ci.macro[[1]]), ncol=ncol(p.tab.ci.macro[[1]]), byrow=TRUE) - p.tab.ci.macro[[1]], 2, quantile, probs=c(0.975))
p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.macro - p.tab.macro)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.macro[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.macro/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.macro/N)[10:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.macro - p.tab.macro)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.macro[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.macro/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.macro/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.macro - p.tab.macro)[10:13] / sum(p.orig.tab.macro[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.macro[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.macro)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.macro/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.macro/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.macro - p.tab.macro)[11:13] / sum(p.orig.tab.macro[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.macro[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.macro)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.macro/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.macro/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(zdat.macro)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))
p.table.fin <- cbind(p.table.fin, p.table)

rownames(p.table.fin) <- c("p > 0.9", "0.9 > p > 0.8", "0.8 > p > 0.7", 
                           "0.7 > p > 0.6","0.6 > p > 0.5", 
                           "0.5 > p > 0.4","0.4 > p > 0.3", 
                           "0.3 > p > 0.2","0.2 > p > 0.1", 
                           "0.1 > p > 0.05","0.05 > p > 0.01", 
                           "0.01 > p > 0.001","0.001 > p", 
                           "IS_{0.10}","IS_{0.05}",
                           "IS_{0.10}","IS_{0.05}",
                           "No. of meta-analysis", "No. of tests")
print(p.table.fin)
write.csv(p.table.fin, paste(path_results,case,"_p.table.micro.macro.csv", sep=""))

## ---
## exp
## ---

load(paste(path_results, case, "_p_tab.exp.RDATA", sep=""))
load(paste(path_results, case, "_p_tab_ci.exp.RDATA", sep=""))

p.table.fin <- NULL
## p-value table (difference in frequencies between factual and cf)
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")
N <- sum(p.orig.tab.exp)
p.table[,1] <- round((p.orig.tab.exp - p.tab.exp) / N,  3)
q025 <- apply(matrix(p.orig.tab.exp/N, nrow=nrow(p.tab.ci.exp[[1]]), ncol=ncol(p.tab.ci.exp[[1]]), byrow=TRUE) - p.tab.ci.exp[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.exp/N, nrow=nrow(p.tab.ci.exp[[1]]), ncol=ncol(p.tab.ci.exp[[1]]), byrow=TRUE) - p.tab.ci.exp[[1]], 2, quantile, probs=c(0.975))
p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.exp - p.tab.exp)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.exp[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.exp/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.exp/N)[10:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.exp - p.tab.exp)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.exp[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.exp/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.exp/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.exp - p.tab.exp)[10:13] / sum(p.orig.tab.exp[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.exp[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.exp)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.exp/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.exp/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.exp - p.tab.exp)[11:13] / sum(p.orig.tab.exp[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.exp[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.exp)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.exp/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.exp/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(zdat.exp)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))
p.table.fin <- cbind(p.table.fin,p.table)

## ---
## obs
## ---

load(paste(path_results, case, "_p_tab.obs.RDATA", sep=""))
load(paste(path_results, case, "_p_tab_ci.obs.RDATA", sep=""))

## p-value table (difference in frequencies between factual and cf)
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")
N <- sum(p.orig.tab.obs)
p.table[,1] <- round((p.orig.tab.obs - p.tab.obs) / N,  3)
q025 <- apply(matrix(p.orig.tab.obs/N, nrow=nrow(p.tab.ci.obs[[1]]), ncol=ncol(p.tab.ci.obs[[1]]), byrow=TRUE) - p.tab.ci.obs[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.obs/N, nrow=nrow(p.tab.ci.obs[[1]]), ncol=ncol(p.tab.ci.obs[[1]]), byrow=TRUE) - p.tab.ci.obs[[1]], 2, quantile, probs=c(0.975))
p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.obs - p.tab.obs)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.obs[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.obs/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.obs/N)[10:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.obs - p.tab.obs)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.obs[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.obs/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.obs/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.obs - p.tab.obs)[10:13] / sum(p.orig.tab.obs[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.obs[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.obs)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.obs/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.obs/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.obs - p.tab.obs)[11:13] / sum(p.orig.tab.obs[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.obs[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.obs)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.obs/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.obs/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(zdat.obs)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))
p.table.fin <- cbind(p.table.fin,p.table)

rownames(p.table.fin) <- c("p > 0.9", "0.9 > p > 0.8", "0.8 > p > 0.7", 
                           "0.7 > p > 0.6","0.6 > p > 0.5", 
                           "0.5 > p > 0.4","0.4 > p > 0.3", 
                           "0.3 > p > 0.2","0.2 > p > 0.1", 
                           "0.1 > p > 0.05","0.05 > p > 0.01", 
                           "0.01 > p > 0.001","0.001 > p", 
                           "IS_{0.10}","IS_{0.05}",
                           "IS_{0.10}","IS_{0.05}",
                           "No. of meta-analysis", "No. of tests")
print(p.table.fin)
write.csv(p.table.fin,paste(path_results,case,"_p.table.exp.obs.csv",sep=""))

## ---
## pow
## ---

load(paste(path_results, case, "_p_tab.pow.RDATA", sep=""))
load(paste(path_results, case, "_p_tab_ci.pow.RDATA", sep=""))

p.table.fin <- NULL
## p-value table (difference in frequencies between factual and cf)
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")
N <- sum(p.orig.tab.pow)
p.table[,1] <- round((p.orig.tab.pow - p.tab.pow) / N,  3)
q025 <- apply(matrix(p.orig.tab.pow/N, nrow=nrow(p.tab.ci.pow[[1]]), ncol=ncol(p.tab.ci.pow[[1]]), byrow=TRUE) - p.tab.ci.pow[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.pow/N, nrow=nrow(p.tab.ci.pow[[1]]), ncol=ncol(p.tab.ci.pow[[1]]), byrow=TRUE) - p.tab.ci.pow[[1]], 2, quantile, probs=c(0.975))
p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.pow - p.tab.pow)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.pow[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.pow/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.pow/N)[10:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.pow - p.tab.pow)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.pow[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.pow/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.pow/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.pow - p.tab.pow)[10:13] / sum(p.orig.tab.pow[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.pow[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.pow)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.pow/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.pow/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.pow - p.tab.pow)[11:13] / sum(p.orig.tab.pow[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.pow[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.pow)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.pow/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.pow/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(zdat.pow)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))
p.table.fin <- cbind(p.table.fin,p.table)

## ----
## npow
## ----

load(paste(path_results, case, "_p_tab.npow.RDATA", sep=""))
load(paste(path_results, case, "_p_tab_ci.npow.RDATA", sep=""))

## p-value table (difference in frequencies between factual and cf)
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")
N <- sum(p.orig.tab.npow)
p.table[,1] <- round((p.orig.tab.npow - p.tab.npow) / N,  3)
q025 <- apply(matrix(p.orig.tab.npow/N, nrow=nrow(p.tab.ci.npow[[1]]), ncol=ncol(p.tab.ci.npow[[1]]), byrow=TRUE) - p.tab.ci.npow[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.npow/N, nrow=nrow(p.tab.ci.npow[[1]]), ncol=ncol(p.tab.ci.npow[[1]]), byrow=TRUE) - p.tab.ci.npow[[1]], 2, quantile, probs=c(0.975))
p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.npow - p.tab.npow)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.npow[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.npow/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.npow/N)[10:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.npow - p.tab.npow)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.npow[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.npow/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.npow/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.npow - p.tab.npow)[10:13] / sum(p.orig.tab.npow[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.npow[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.npow)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.npow/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.npow/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.npow - p.tab.npow)[11:13] / sum(p.orig.tab.npow[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.npow[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.npow)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.npow/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.npow/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(zdat.npow)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))
p.table.fin <- cbind(p.table.fin,p.table)

rownames(p.table.fin) <- c("p > 0.9", "0.9 > p > 0.8", "0.8 > p > 0.7", 
                           "0.7 > p > 0.6","0.6 > p > 0.5", 
                           "0.5 > p > 0.4","0.4 > p > 0.3", 
                           "0.3 > p > 0.2","0.2 > p > 0.1", 
                           "0.1 > p > 0.05","0.05 > p > 0.01", 
                           "0.01 > p > 0.001","0.001 > p", 
                           "IS_{0.10}","IS_{0.05}",
                           "IS_{0.10}","IS_{0.05}",
                           "No. of meta-analysis", "No. of tests")
print(p.table.fin)
write.csv(p.table.fin,paste(path_results,case,"_p.table.pow.npow.csv",sep=""))



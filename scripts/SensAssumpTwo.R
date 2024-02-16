
rm(list = ls(all = TRUE))                     #clear environment

library(doParallel) ;  library(foreach)
library(metafor);  library(openxlsx)
library(readxl) 

# Paths for results and functions
path_results <- "G:/My Drive/InflatedSigInEco/results/Assumption2/"
path_fcts    <- "G:/My Drive/InflatedSigInEco/scripts/"
path_res_meta_rem <- "G:/My Drive/InflatedSigInEco/results/Assumption2/REM_Per_Meta/"

# Functions
source(paste(path_fcts,"econ_functions.r", sep=""))

## Load data
dat_final <- read_excel("G:/My Drive/InflatedSigInEco/data/data_final.xlsx")
dat_final$vi <- dat_final$se^2  #calculating variance
dim(dat_final) #74288x12


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


## -------------------------------------------------------------
## Assumption 2: Alternative meta-analytic estimators
## WLS, PET-PEESE, and random-effects model are used
## -------------------------------------------------------------

## (i) Estimating genuine effect using WLS
wlsGE <- NULL
for (i in 1:length(summary(zdat)[,1])) {
  wls <- sum(zdat[[i]]$eff / zdat[[i]]$se^2) / sum(1/zdat[[i]]$se^2)
  wlsGE <- c(wlsGE, rep(wls, length(zdat[[i]]$eff)))
}

case   <- "wlsGE"
caseGE <- wlsGE

## ------------------------------------
## Add estimated genuine effect to zdat
## ------------------------------------

dat.ideas <- do.call(rbind.data.frame, zdat)
ii <- 1
for (i in unique(dat.ideas$filecode)) {	
  zdat[[ii]]$GE <- caseGE[which(dat.ideas$filecode == i)]
  #zdat[[ii]]$n.waap <- n.waap[[ii]]
  #zdat[[ii]]$clu <- clu[[ii]]
  ii <- ii + 1
}

head(dat.ideas)
dim(dat.ideas)  # 70399 x 12


## -------------------------
## z-value and p-value grids
## -------------------------

## grid of p-values and grid of z-values
# z.grid is used for visualization and goes from 1 to 10 while
# p.grid is used for exact calculations and goes from z-values 
# of 0 to Inf (absolute values).

# p.grid for main results table
p.grid.tab <- c(-Inf,
                qnorm(0.001/2, lower.tail = TRUE),
                qnorm(0.01/2, lower.tail = TRUE),
                qnorm(0.05/2, lower.tail = TRUE),
                qnorm(0.1/2, lower.tail = TRUE),
                qnorm(0.2/2, lower.tail = TRUE),
                qnorm(0.3/2, lower.tail = TRUE),
                qnorm(0.4/2, lower.tail = TRUE),
                qnorm(0.5/2, lower.tail = TRUE),
                qnorm(0.6/2, lower.tail = TRUE),
                qnorm(0.7/2, lower.tail = TRUE),
                qnorm(0.8/2, lower.tail = TRUE),
                qnorm(0.9/2, lower.tail = TRUE),
                0,
                qnorm(0.9/2, lower.tail = FALSE),
                qnorm(0.8/2, lower.tail = FALSE),
                qnorm(0.7/2, lower.tail = FALSE),
                qnorm(0.6/2, lower.tail = FALSE),
                qnorm(0.5/2, lower.tail = FALSE),
                qnorm(0.4/2, lower.tail = FALSE),
                qnorm(0.3/2, lower.tail = FALSE),
                qnorm(0.2/2, lower.tail = FALSE),
                qnorm(0.1/2, lower.tail = FALSE),
                qnorm(0.05/2, lower.tail = FALSE),
                qnorm(0.01/2, lower.tail = FALSE),
                qnorm(0.001/2, lower.tail = FALSE),
                Inf)

# We need to calculate for the orginal z-values also the frequencies
# and we can do this directly for the two-sided test grid
p.grid.tab2 <- c(0, 
                 qnorm(0.9/2, lower.tail = FALSE),
                 qnorm(0.8/2, lower.tail = FALSE),
                 qnorm(0.7/2, lower.tail = FALSE),
                 qnorm(0.6/2, lower.tail = FALSE),
                 qnorm(0.5/2, lower.tail = FALSE),
                 qnorm(0.4/2, lower.tail = FALSE),
                 qnorm(0.3/2, lower.tail = FALSE),
                 qnorm(0.2/2, lower.tail = FALSE),
                 qnorm(0.1/2, lower.tail = FALSE),
                 qnorm(0.05/2, lower.tail = FALSE),
                 qnorm(0.01/2, lower.tail = FALSE),
                 qnorm(0.001/2, lower.tail = FALSE),
                 Inf)				

## p grid for plot
p.grid.plot <- qnorm(seq(1, 0, -0.005), lower.tail = FALSE)
p.grid.plot2 <-p.grid.plot[which(p.grid.plot >=0)] # this is for actual plotting of two-sided tests

## z.grid for plot
z.grid.plot  <- seq(-10.25, 10.25, 0.1025)
z.grid.plot2 <- seq(0, 10.25, 0.1025) # this is for actual plotting of absolute values

## ----------------------------------------------------
## Prepare original data for plotting and p-value table
## ----------------------------------------------------

fac  <- do.call(rbind.data.frame, zdat)
facz <- abs(fac$eff / fac$se)

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

## -----------------------------------------------------------------------------

## -----------------------
## Sensitivity table (WLS)
## -----------------------

## p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab <- cf(dat = zdat, z.grid = p.grid.tab)
stopCluster(cl)
save(p.tab,file=paste(path_results, case,"_p_tab.RDATA",sep=""))

clu <- unique(dat.ideas$filecode) #cluster id
itero  <- 1000       # for final result

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.ci <- cf.ci.cluster(dat = zdat, z.grid = p.grid.tab, iters = itero, cluster = clu)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  # about 8 mins
save(p.tab.ci,file=paste(path_results,case,"_p_tab_ci.RDATA",sep=""))

load(paste(path_results, case, "_p_tab.RDATA", sep=""))
load(paste(path_results, case, "_p_tab_ci.RDATA", sep=""))

#p-value table (difference in frequencies between f & cf)

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
infp.q025 <- round(quantile(sum((p.orig.tab/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab/N)[10:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

# share of inflated p-values in total p-values (5% level)
infp    <- round(sum((p.orig.tab - p.tab)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab/N)[11:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#10%
shp.10     <- round(sum((p.orig.tab - p.tab)[10:13] / sum(p.orig.tab[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab/N.sig.10)[10:13]) - infp.bs.10), probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab/N.sig.10)[10:13]) - infp.bs.10), probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#5%
shp.5 <- round(sum((p.orig.tab - p.tab)[11:13] / sum(p.orig.tab[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))


## sample size
p.table <- rbind(p.table, c(length(summary(zdat)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("$p$ $>$ 0.9", "0.9 $>$ $p$ $>$ 0.8", "0.8 $>$ $p$ $>$ 0.7", 
                       "0.7 $>$ $p$ $>$ 0.6", "0.6 $>$ $p$ $>$ 0.5", 
                       "0.5 $>$ $p$ $>$ 0.4", "0.4 $>$ $p$ $>$ 0.3", 
                       "0.3 $>$ $p$ $>$ 0.2", "0.2 $>$ $p$ $>$ 0.1", 
                       "0.1 $>$ $p$ $>$ 0.05", "0.05 $>$ $p$ $>$ 0.01", 
                       "0.01 $>$ $p$ $>$ 0.001", "0.001 $>$ $p$", 
                       "$IS_{0.10}$",
                       "$IS_{0.05}$",
                       "$IS_{0.10}$" ,
                       "$IS_{0.05}$" ,
                       "No. of meta-analysis", 
                       "No. of tests")
write.csv(p.table, paste(path_results,case,"_p.table.csv", sep=""))


## -----------------------------------------------------------------------------

## (ii) Estimating genuine effect using PET-PEESE
ppGE <- NULL
for (i in 1:length(summary(zdat)[,1])) {
  
  # PET-PEESE (PET if alpha is not significant, PEESE if it is significant.)
  if (coefficients(summary(lm(zdat[[i]]$eff ~ zdat[[i]]$se, weights=(1/zdat[[i]]$se)^2)))[1,4] > 0.1) {
    ppGE <- c(ppGE, rep(coefficients(summary(lm(zdat[[i]]$eff ~ zdat[[i]]$se, 
              weights=(1/zdat[[i]]$se)^2)))[1,1], length(zdat[[i]]$eff)))
  } else {
    se2 <- zdat[[i]]$se^2
    ppGE <- c(ppGE, rep(coefficients(summary(lm(zdat[[i]]$eff ~ se2, 
              weights=(1/zdat[[i]]$se)^2)))[1,1], length(zdat[[i]]$eff)))
  }
}

case   <- "ppGE"
caseGE <- ppGE

## ------------------------------------
## Add estimated genuine effect to zdat
## ------------------------------------

dat.ideas <- do.call(rbind.data.frame, zdat)
ii <- 1
for (i in unique(dat.ideas$filecode)) {	
  zdat[[ii]]$GE <- caseGE[which(dat.ideas$filecode == i)]
  #zdat[[ii]]$n.waap <- n.waap[[ii]]
  #zdat[[ii]]$clu <- clu[[ii]]
  ii <- ii + 1
}

head(dat.ideas)
dim(dat.ideas)  # 70399 x 12


## -------------------------
## z-value and p-value grids
## -------------------------

## grid of p-values and grid of z-values
# z.grid is used for visualization and goes from 1 to 10 while
# p.grid is used for exact calculations and goes from z-values 
# of 0 to Inf (absolute values).

# p.grid for main results table
p.grid.tab <- c(-Inf,
                qnorm(0.001/2, lower.tail = TRUE),
                qnorm(0.01/2, lower.tail = TRUE),
                qnorm(0.05/2, lower.tail = TRUE),
                qnorm(0.1/2, lower.tail = TRUE),
                qnorm(0.2/2, lower.tail = TRUE),
                qnorm(0.3/2, lower.tail = TRUE),
                qnorm(0.4/2, lower.tail = TRUE),
                qnorm(0.5/2, lower.tail = TRUE),
                qnorm(0.6/2, lower.tail = TRUE),
                qnorm(0.7/2, lower.tail = TRUE),
                qnorm(0.8/2, lower.tail = TRUE),
                qnorm(0.9/2, lower.tail = TRUE),
                0,
                qnorm(0.9/2, lower.tail = FALSE),
                qnorm(0.8/2, lower.tail = FALSE),
                qnorm(0.7/2, lower.tail = FALSE),
                qnorm(0.6/2, lower.tail = FALSE),
                qnorm(0.5/2, lower.tail = FALSE),
                qnorm(0.4/2, lower.tail = FALSE),
                qnorm(0.3/2, lower.tail = FALSE),
                qnorm(0.2/2, lower.tail = FALSE),
                qnorm(0.1/2, lower.tail = FALSE),
                qnorm(0.05/2, lower.tail = FALSE),
                qnorm(0.01/2, lower.tail = FALSE),
                qnorm(0.001/2, lower.tail = FALSE),
                Inf)

# We need to calculate for the orginal z-values also the frequencies
# and we can do this directly for the two-sided test grid
p.grid.tab2 <- c(0, 
                 qnorm(0.9/2, lower.tail = FALSE),
                 qnorm(0.8/2, lower.tail = FALSE),
                 qnorm(0.7/2, lower.tail = FALSE),
                 qnorm(0.6/2, lower.tail = FALSE),
                 qnorm(0.5/2, lower.tail = FALSE),
                 qnorm(0.4/2, lower.tail = FALSE),
                 qnorm(0.3/2, lower.tail = FALSE),
                 qnorm(0.2/2, lower.tail = FALSE),
                 qnorm(0.1/2, lower.tail = FALSE),
                 qnorm(0.05/2, lower.tail = FALSE),
                 qnorm(0.01/2, lower.tail = FALSE),
                 qnorm(0.001/2, lower.tail = FALSE),
                 Inf)				

# p grid for plot
p.grid.plot <- qnorm(seq(1, 0, -0.005), lower.tail = FALSE)
p.grid.plot2 <-p.grid.plot[which(p.grid.plot >=0)] # this is for actual plotting of two-sided tests

## z.grid for plot
z.grid.plot  <- seq(-10.25, 10.25, 0.1025)
z.grid.plot2 <- seq(0, 10.25, 0.1025) # this is for actual plotting of absolute values

## ----------------------------------------------------
## Prepare original data for plotting and p-value table
## ----------------------------------------------------

fac  <- do.call(rbind.data.frame, zdat)
facz <- abs(fac$eff / fac$se)

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

## -----------------------------------------------------------------------------

## -----------------------------
## Sensitivity table (PET-PEESE)
## -----------------------------

## p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab <- cf(dat = zdat, z.grid = p.grid.tab)
stopCluster(cl)
save(p.tab,file=paste(path_results,case,"_p_tab.RDATA",sep=""))

clu <- unique(dat.ideas$filecode) #cluster id
itero  <- 1000       # for final result
s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.ci <- cf.ci.cluster(dat = zdat, z.grid = p.grid.tab, iters = itero, cluster = clu)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  # about 2 mins
save(p.tab.ci,file=paste(path_results,case,"_p_tab_ci.RDATA",sep=""))

load(paste(path_results, case, "_p_tab.RDATA", sep=""))
load(paste(path_results, case, "_p_tab_ci.RDATA", sep=""))

#p-value table (difference in frequencies between f & cf)
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
infp.q025 <- round(quantile(sum((p.orig.tab/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab/N)[10:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

# share of inflated p-values in total p-values (5% level)
infp    <- round(sum((p.orig.tab - p.tab)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab/N)[11:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#10%
shp.10     <- round(sum((p.orig.tab - p.tab)[10:13] / sum(p.orig.tab[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab/N.sig.10)[10:13]) - infp.bs.10), probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab/N.sig.10)[10:13]) - infp.bs.10), probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#5%
shp.5 <- round(sum((p.orig.tab - p.tab)[11:13] / sum(p.orig.tab[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

## sample size
p.table <- rbind(p.table, c(length(summary(zdat)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("$p$ $>$ 0.9", "0.9 $>$ $p$ $>$ 0.8", "0.8 $>$ $p$ $>$ 0.7", 
                       "0.7 $>$ $p$ $>$ 0.6", "0.6 $>$ $p$ $>$ 0.5", 
                       "0.5 $>$ $p$ $>$ 0.4", "0.4 $>$ $p$ $>$ 0.3", 
                       "0.3 $>$ $p$ $>$ 0.2", "0.2 $>$ $p$ $>$ 0.1", 
                       "0.1 $>$ $p$ $>$ 0.05", "0.05 $>$ $p$ $>$ 0.01", 
                       "0.01 $>$ $p$ $>$ 0.001", "0.001 $>$ $p$", 
                       "$D_{0.10}$",
                       "$D_{0.05}$",
                       "$IS_{0.10}$" ,
                       "$IS_{0.05}$" ,
                       "No. of meta-analysis", 
                       "No. of tests")
write.csv(p.table, paste(path_results,case,"_p.table.csv", sep=""))

## -----------------------------------------------------------------------------

## (iii) Estimating genuine effect & between-study heterogeneity (tau2)
##       using random-effects model & add to zdat[[i]]

## Load data
dat_final <- read_excel("G:/My Drive/InflatedSigInEco/Data/data_final.xlsx")
dat_final$vi <- dat_final$se^2  #calculating variance
dim(dat_final) #74288x12


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
dat.ideas <- do.call(rbind.data.frame, zdat)

s.time <- Sys.time()
ii <- 1
mod <- list()
for (i in unique(dat.ideas$filecode)) {	
  zdat[[ii]] <- dat.ideas[dat.ideas$filecode==i, ]
  mod[[ii]] <- rma.uni(yi=zdat[[ii]]$eff, sei=zdat[[ii]]$se, method="REML")
  
  zdat[[ii]]$GE <- as.numeric(mod[[ii]]$beta) 
  zdat[[ii]]$tau2 <- mod[[ii]]$tau2
  saveRDS(zdat[[ii]], file=paste(path_res_meta_rem,"meta_", zdat[[ii]]$filecode[1],".rds",sep=""))
  ii <- ii + 1
}
e.time <- Sys.time()
print(e.time - s.time) #about 25 minutes

head(dat.ideas)
dim(dat.ideas)

## One can start from HERE!
## Reload the data after GE & tau2 are estimated
dat.ideas <- list.files("G:/My Drive/InflatedSigInEco/Results/Assumption2/REM_Per_Meta/",pattern=".rds",full.names=TRUE) %>%
  map_dfr(readRDS)
dim(dat.ideas) #72799x16
names(dat.ideas)

zdat <- list()
jj <- 1
for (i in unique(dat.ideas$filecode)){
  zdat[[jj]] <- dat.ideas[which(dat.ideas$filecode==i), ]
  jj <- jj + 1
}

## Check
mss <- NULL
for (i in 1:dim(summary(zdat))[1]) { 
  mss[i] <- length(zdat[[i]]$se)
}
length(mss) #192
sum(mss)    #70399

## -------------------------
## z-value and p-value grids
## -------------------------

## grid of p-values and grid of z-values
# z.grid is used for visualization and goes from 1 to 10 while
# p.grid is used for exact calculations and goes from z-values 
# of 0 to Inf (absolute values).

# p.grid for main results table
p.grid.tab <- c(-Inf,
                qnorm(0.001/2, lower.tail = TRUE),
                qnorm(0.01/2, lower.tail = TRUE),
                qnorm(0.05/2, lower.tail = TRUE),
                qnorm(0.1/2, lower.tail = TRUE),
                qnorm(0.2/2, lower.tail = TRUE),
                qnorm(0.3/2, lower.tail = TRUE),
                qnorm(0.4/2, lower.tail = TRUE),
                qnorm(0.5/2, lower.tail = TRUE),
                qnorm(0.6/2, lower.tail = TRUE),
                qnorm(0.7/2, lower.tail = TRUE),
                qnorm(0.8/2, lower.tail = TRUE),
                qnorm(0.9/2, lower.tail = TRUE),
                0,
                qnorm(0.9/2, lower.tail = FALSE),
                qnorm(0.8/2, lower.tail = FALSE),
                qnorm(0.7/2, lower.tail = FALSE),
                qnorm(0.6/2, lower.tail = FALSE),
                qnorm(0.5/2, lower.tail = FALSE),
                qnorm(0.4/2, lower.tail = FALSE),
                qnorm(0.3/2, lower.tail = FALSE),
                qnorm(0.2/2, lower.tail = FALSE),
                qnorm(0.1/2, lower.tail = FALSE),
                qnorm(0.05/2, lower.tail = FALSE),
                qnorm(0.01/2, lower.tail = FALSE),
                qnorm(0.001/2, lower.tail = FALSE),
                Inf)

# We need to calculate for the orginal z-values also the frequencies
# and we can do this directly for the two-sided test grid
p.grid.tab2 <- c(0, 
                 qnorm(0.9/2, lower.tail = FALSE),
                 qnorm(0.8/2, lower.tail = FALSE),
                 qnorm(0.7/2, lower.tail = FALSE),
                 qnorm(0.6/2, lower.tail = FALSE),
                 qnorm(0.5/2, lower.tail = FALSE),
                 qnorm(0.4/2, lower.tail = FALSE),
                 qnorm(0.3/2, lower.tail = FALSE),
                 qnorm(0.2/2, lower.tail = FALSE),
                 qnorm(0.1/2, lower.tail = FALSE),
                 qnorm(0.05/2, lower.tail = FALSE),
                 qnorm(0.01/2, lower.tail = FALSE),
                 qnorm(0.001/2, lower.tail = FALSE),
                 Inf)				

# p grid for plot
p.grid.plot <- qnorm(seq(1, 0, -0.005), lower.tail = FALSE)
p.grid.plot2 <-p.grid.plot[which(p.grid.plot >=0)] # this is for actual plotting of two-sided tests

## z.grid for plot
z.grid.plot  <- seq(-10.25, 10.25, 0.1025)
z.grid.plot2 <- seq(0, 10.25, 0.1025) # this is for actual plotting of absolute values

## ----------------------------------------------------
## Prepare original data for plotting and p-value table
## ----------------------------------------------------

fac  <- do.call(rbind.data.frame, zdat)
facz <- abs(fac$eff / fac$se)

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

## -----------------------------------------------------------------------------

## -----------------------------
## Sensitivity table (REM)
## -----------------------------

## p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab <- cf.rem(dat = zdat, z.grid = p.grid.tab)
stopCluster(cl)
save(p.tab,file=paste(path_results, "rem_p_tab.RDATA",sep=""))

clu <- unique(dat.ideas$filecode) #cluster id
itero  <- 1000       # for final result

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.ci <- cf.ci.cluster.rem(dat = zdat, z.grid = p.grid.tab, iters = itero, cluster = clu)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  # about 2 mins
save(p.tab.ci,file=paste(path_results,"rem_p_tab_ci.RDATA",sep=""))

load(paste(path_results, "rem_p_tab.RDATA", sep=""))
load(paste(path_results, "rem_p_tab_ci.RDATA", sep=""))

#p-value table (difference in frequencies between f & cf)
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
infp.q025 <- round(quantile(sum((p.orig.tab/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab/N)[10:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

# share of inflated p-values in total p-values (5% level)
infp    <- round(sum((p.orig.tab - p.tab)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#10%
shp.10     <- round(sum((p.orig.tab - p.tab)[10:13] / sum(p.orig.tab[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab/N.sig.10)[10:13]) - infp.bs.10), probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab/N.sig.10)[10:13]) - infp.bs.10), probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#5%
shp.5 <- round(sum((p.orig.tab - p.tab)[11:13] / sum(p.orig.tab[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

## sample size
p.table <- rbind(p.table, c(length(summary(zdat)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("$p$ $>$ 0.9", "0.9 $>$ $p$ $>$ 0.8", "0.8 $>$ $p$ $>$ 0.7", 
                       "0.7 $>$ $p$ $>$ 0.6", "0.6 $>$ $p$ $>$ 0.5", 
                       "0.5 $>$ $p$ $>$ 0.4", "0.4 $>$ $p$ $>$ 0.3", 
                       "0.3 $>$ $p$ $>$ 0.2", "0.2 $>$ $p$ $>$ 0.1", 
                       "0.1 $>$ $p$ $>$ 0.05", "0.05 $>$ $p$ $>$ 0.01", 
                       "0.01 $>$ $p$ $>$ 0.001", "0.001 $>$ $p$", 
                       "$IS_{0.10}$",
                       "$IS_{0.05}$",
                       "$IS_{0.10}$" ,
                       "$IS_{0.05}$" ,
                       "No. of meta-analysis", 
                       "No. of tests")
print(p.table)
write.csv(p.table, paste(path_results,"rem_p.table.csv", sep=""))






rm(list = ls(all = TRUE))                     #clear environment
setwd("G:/My Drive/InflatedSigInEco/R codes") #set the directory

library(doParallel)
library(foreach)
library(metafor)
library(openxlsx)
library(readxl) 

# Paths for results and functions
path_results <- "G:/My Drive/InflatedSigInEco/results/Assumption1/"
path_fcts    <- "G:/My Drive/InflatedSigInEco/R codes/"

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

## -----------------------------------------------------------
## Assumption 1: Two & three genuine effects per meta-analysis
## Each meta-analysis is splited randomly into two and three 
## groups
## -----------------------------------------------------------

## -------------------------
## z-value and p-value grids
## -------------------------

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

#### z.grid for plot
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

## ------------------------------------------
## Robustness analyses: 2 & 3 genuine effects
## ------------------------------------------

cl <- makeCluster(7) 
registerDoParallel(cl)
s.time <- Sys.time()
p.tab.2GE <- cf.mc(dat=zdat, z.grid=p.grid.tab, GE.type=1, mc.iter=500)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) # about 8 hrs for 7500 iterations (7 cores)
save(p.tab.2GE,file=paste(path_results,"wls_p_tab_2GE.RDATA",sep=""))

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.3GE <- cf.mc(dat=zdat, z.grid=p.grid.tab, GE.type=2, mc.iter=500)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) # about 8 hrs for 7500 iterations (7 cores)
save(p.tab.3GE,file=paste(path_results,"wls_p_tab_3GE.RDATA",sep=""))

## Two genuine effects 
load(paste(path_results, "wls_p_tab_2GE.RDATA", sep=""))

p.tab.2GE <- do.call(rbind.data.frame, p.tab.2GE)
p.table <- matrix(ncol=6, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Mean", "Min", "0.25", "0.50", "0.75", "Max")
N <- sum(p.orig.tab)
p.table[,1]   <- round((p.orig.tab - apply(p.tab.2GE, 2, mean)) / N,  3) 
p.table[,2:6] <- round(t(apply((t(t(p.orig.tab)) - p.tab.2GE), 2, quantile)) / N,  3)

for (pp in 1:13) {
  p.table[pp,2:6] <- round(quantile((p.orig.tab[pp] - p.tab.2GE[,pp]) / N),  3)
}

# share of inflated p-values in total p-values (10% level)
infp   <- round(mean((sum(p.orig.tab[10:13]) - apply(p.tab.2GE[,10:13], 1, sum)) / N),  3)
infp.q <- round(quantile((sum(p.orig.tab[10:13]) - apply(p.tab.2GE[,10:13], 1, sum)) / N),  3)
p.table <- rbind(p.table, c(infp, infp.q))

# share of inflated p-values in total p-values (5% level)
infp   <- round(mean((sum(p.orig.tab[11:13]) - apply(p.tab.2GE[,11:13], 1, sum)) / N),  3)
infp.q <- round(quantile((sum(p.orig.tab[11:13]) - apply(p.tab.2GE[,11:13], 1, sum)) / N),  3)
p.table <- rbind(p.table, c(infp, infp.q))

#10%
shp   <- round(mean((sum(p.orig.tab[10:13]) - apply(p.tab.2GE[,10:13], 1, sum)) / sum(p.orig.tab[10:13])),  3)
shp.q <- round(quantile((sum(p.orig.tab[10:13]) - apply(p.tab.2GE[,10:13], 1, sum)) / sum(p.orig.tab[10:13])),  3)
p.table <- rbind(p.table, c(shp, shp.q))

#5%
shp   <- round(mean((sum(p.orig.tab[11:13]) - apply(p.tab.2GE[,11:13], 1, sum)) / sum(p.orig.tab[11:13])),  3)
shp.q <- round(quantile((sum(p.orig.tab[11:13]) - apply(p.tab.2GE[,11:13], 1, sum)) / sum(p.orig.tab[11:13])),  3)
p.table <- rbind(p.table, c(shp, shp.q))


# sample size
p.table <- rbind(p.table, c(length(summary(zdat)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("$p$ $>$ 0.9", "0.9 $>$ $p$ $>$ 0.8", "0.8 $>$ $p$ $>$ 0.7", 
                       "0.7 $>$ $p$ $>$ 0.6", "0.6 $>$ $p$ $>$ 0.5", 
                       "0.5 $>$ $p$ $>$ 0.4", "0.4 $>$ $p$ $>$ 0.3", 
                       "0.3 $>$ $p$ $>$ 0.2", "0.2 $>$ $p$ $>$ 0.1", 
                       "0.1 $>$ $p$ $>$ 0.05", "0.05 $>$ $p$ $>$ 0.01", 
                       "0.01 $>$ $p$ $>$ 0.001", "0.001 $>$ $p$", 
                       "$IS_{0.10}$","$IS_{0.05}$","$IS_{0.10}$", "$IS_{0.05}$" ,
                       "No. of meta-analysis","No. of tests")
write.csv(p.table,paste(path_results,"wls_p.table_2GE.csv",sep=""))


## Three genuine effects 
load(paste(path_results, "wls_p_tab_3GE.RDATA", sep=""))

p.tab.3GE <- do.call(rbind.data.frame, p.tab.3GE)
p.table <- matrix(ncol=6, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Mean", "Min", "0.25", "0.50", "0.75", "Max")
N <- sum(p.orig.tab)
p.table[,1]   <- round((p.orig.tab - apply(p.tab.3GE, 2, mean)) / N,  3)
p.table[,2:6] <- round(t(apply((t(t(p.orig.tab)) - p.tab.3GE), 2, quantile)) / N,  3)

for (pp in 1:13) {
  p.table[pp,2:6] <- round(quantile((p.orig.tab[pp] - p.tab.3GE[,pp]) / N),  3)
}

# share of inflated p-values in total p-values (10% level)
infp   <- round(mean((sum(p.orig.tab[10:13]) - apply(p.tab.3GE[,10:13], 1, sum)) / N),  3)
infp.q <- round(quantile((sum(p.orig.tab[10:13]) - apply(p.tab.3GE[,10:13], 1, sum)) / N),  3)
p.table <- rbind(p.table, c(infp, infp.q))

# share of inflated p-values in total p-values (5% level)
infp   <- round(mean((sum(p.orig.tab[11:13]) - apply(p.tab.3GE[,11:13], 1, sum)) / N),  3)
infp.q <- round(quantile((sum(p.orig.tab[11:13]) - apply(p.tab.3GE[,11:13], 1, sum)) / N),  3)
p.table <- rbind(p.table, c(infp, infp.q))

#10%
shp   <- round(mean((sum(p.orig.tab[10:13]) - apply(p.tab.3GE[,10:13], 1, sum)) / sum(p.orig.tab[10:13])),  3)
shp.q <- round(quantile((sum(p.orig.tab[10:13]) - apply(p.tab.3GE[,10:13], 1, sum)) / sum(p.orig.tab[10:13])),  3)
p.table <- rbind(p.table, c(shp, shp.q))

#5%
shp   <- round(mean((sum(p.orig.tab[11:13]) - apply(p.tab.3GE[,11:13], 1, sum)) / sum(p.orig.tab[11:13])),  3)
shp.q <- round(quantile((sum(p.orig.tab[11:13]) - apply(p.tab.3GE[,11:13], 1, sum)) / sum(p.orig.tab[11:13])),  3)
p.table <- rbind(p.table, c(shp, shp.q))

# sample size
p.table <- rbind(p.table, c(length(summary(zdat)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("$p$ $>$ 0.9", "0.9 $>$ $p$ $>$ 0.8", "0.8 $>$ $p$ $>$ 0.7", "0.7 $>$ $p$ $>$ 0.6",
                       "0.6 $>$ $p$ $>$ 0.5", "0.5 $>$ $p$ $>$ 0.4", "0.4 $>$ $p$ $>$ 0.3", "0.3 $>$ $p$ $>$ 0.2",
                       "0.2 $>$ $p$ $>$ 0.1", "0.1 $>$ $p$ $>$ 0.05", "0.05 $>$ $p$ $>$ 0.01", "0.01 $>$ $p$ $>$ 0.001",
                       "0.001 $>$ $p$","$IS_{0.1}^{all}$", 
                       "$IS_{0.05}^{all}$","$IS_{0.1}{sig.}$","$IS_{0.05}^{sig.}$",
                       "No. of meta-analysis", "No. of tests")
write.csv(p.table,paste(path_results,"wls_p.table_3GE.csv",sep=""))





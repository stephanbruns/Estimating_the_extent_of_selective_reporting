
rm(list = ls(all = TRUE))                     #clear environment

library(doParallel)  
library(foreach)
library(metafor)
library(openxlsx)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggthemes)

## Paths for results and functions
path_results <- "G:/My Drive/InflatedSigInEco/results/AssumptionFurther2/"
path_fcts    <- "G:/My Drive/InflatedSigInEco/scripts/"

## Functions
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
dim(dat_final)  #73198x13

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

## -----------------------------------------------------------------------------
## Further assumption
## Combination of WAAP & PET-PEESE is used
## --------------------------------------------------

wlsGE  <- NULL; ppGE   <- NULL
waapGE <- NULL; n.waap <- NULL

for (i in 1:length(summary(zdat)[,1])) {	
  
  #wls
  wls   <- sum(zdat[[i]]$eff / zdat[[i]]$se^2) / sum(1/zdat[[i]]$se^2)
  wlsGE <- c(wlsGE, rep(wls, length(zdat[[i]]$eff)))
  
  #waap
  se.crit <- abs(wls) / 2.8
  pos <- which(zdat[[i]]$se <= se.crit)
  n.waap[i] <- length(pos) / length(zdat[[i]]$se) #prop. of adequately powered estimates
  
  #at least two adequately powered studies are needed
  if (length(pos) > 1) {
    waapGE <- c(waapGE, rep(sum(zdat[[i]]$eff[pos] / zdat[[i]]$se[pos]^2) / 
                      sum(1/zdat[[i]]$se[pos]^2),length(zdat[[i]]$eff)))
  } else {
    waapGE <- c(waapGE, rep(NA,length(zdat[[i]]$eff)))
  }
  
  #PET-PEESE (PET if alpha is not significant, PEESE if it is significant.)
  if (coefficients(summary(lm(zdat[[i]]$eff ~ zdat[[i]]$se, weights=(1/zdat[[i]]$se)^2)))[1,4] > 0.1) {
    ppGE <- c(ppGE, rep(coefficients(summary(lm(zdat[[i]]$eff ~ zdat[[i]]$se, 
                    weights=(1/zdat[[i]]$se)^2)))[1,1], length(zdat[[i]]$eff)))
  } else {
    se2 <- zdat[[i]]$se^2
    ppGE <-c(ppGE, rep(coefficients(summary(lm(zdat[[i]]$eff ~ se2, 
                    weights=(1/zdat[[i]]$se)^2)))[1,1], length(zdat[[i]]$eff)))
  }
}

## wappGE inserts ppGE whenever there are not at least 
## two adequately powered estimates.
wappGE <- waapGE
for (i in 1:sum(mss)) {
  if (is.na(wappGE[i])==TRUE) {
    wappGE[i] <- ppGE[i]
  }
}

## Robustness (Further assumption)
wappGE_0.1 <- 0.1*wappGE
case   <- "wappGE_0.1"
caseGE <- wappGE_0.1
zdat.orig <- zdat

wappGE_0.2 <- 0.2*wappGE
case   <- "wappGE_0.2"
caseGE <- wappGE_0.2
zdat.orig <- zdat

wappGE_0.3 <- 0.3*wappGE
case   <- "wappGE_0.3"
caseGE <- wappGE_0.3
zdat.orig <- zdat

wappGE_0.4 <- 0.4*wappGE
case   <- "wappGE_0.4"
caseGE <- wappGE_0.4
zdat.orig <- zdat

#wappGE_0.5 <- 0.5*wappGE
#case   <- "wappGE_0.5"
#caseGE <- wappGE_0.5
#zdat.orig <- zdat

wappGE_0.6 <- 0.6*wappGE
case   <- "wappGE_0.6"
caseGE <- wappGE_0.6
zdat.orig <- zdat

wappGE_0.7 <- 0.7*wappGE
case   <- "wappGE_0.7"
caseGE <- wappGE_0.7
zdat.orig <- zdat

wappGE_0.8 <- 0.8*wappGE
case   <- "wappGE_0.8"
caseGE <- wappGE_0.8
zdat.orig <- zdat

wappGE_0.9 <- 0.9*wappGE
case   <- "wappGE_0.9"
caseGE <- wappGE_0.9
zdat.orig <- zdat


## ------------------------------------
## Add estimated genuine effect to zdat
## ------------------------------------

dat.ideas <- do.call(rbind.data.frame, zdat)
ii <- 1
for (i in unique(dat.ideas$filecode)) {	
  zdat[[ii]]$GE <- caseGE[which(dat.ideas$filecode == i)]
  zdat[[ii]]$n.waap <- n.waap[[ii]]
  #zdat[[ii]]$clu <- clu[[ii]]
  ii <- ii + 1
}

head(dat.ideas)
dim(dat.ideas)  #70399 x 13

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

## z.grid for plot
z.grid.plot  <- seq(-10.25, 10.25, 0.1025)
z.grid.plot2 <- seq(0, 10.25, 0.1025) # this is for actual plotting of absolute values

## ----------------------------------------------------
## Prepare original data for plotting and p-value table
## ----------------------------------------------------

fac <- do.call(rbind.data.frame, zdat.orig)
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

## ---------------------------
## Sensitivity table (for all)
## --------------------------

## p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab <- cf(dat = zdat, z.grid = p.grid.tab)
stopCluster(cl)
save(p.tab,file=paste(path_results,case,"_p_tab.RDATA",sep=""))

#itero <- 1000    # for the final results
itero <- 500
clu <- unique(dat.ideas$filecode)

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.ci <- cf.ci.cluster(dat = zdat, z.grid = p.grid.tab,iters = itero, cluster = clu)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  # about 1.7 hrs for 1000 iterations (7 cores)
save(p.tab.ci,file=paste(path_results,case,"_p_tab_ci.RDATA",sep=""))

load(paste(path_results,case,"_p_tab.RDATA",sep=""))
load(paste(path_results,case,"_p_tab_ci.RDATA",sep=""))

# p-value table (difference in frequencies between f & cf)
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
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ", shp.q975.5, "]", sep="")))

# sample size
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
                       "$IS_{0.05}$",
                       "No. of meta-analysis", 
                       "No. of tests")
print(p.table)
write.csv(p.table, paste(path_results,case,"_p.table.csv",sep=""))


## Boxplot of IS (further assumptions)
x <- c("A","B","C","D","E","F","G","H","I","J")
y <- c(0.895,0.853,0.803,0.755,0.713,0.675,0.641,0.612,0.585,0.562)
y.rev <- rev(y)

l <- c(0.889,0.830,0.759,0.692,0.633,0.583,0.539,0.503,0.469,0.441)
l.rev <- rev(l)

u <- c(0.900,0.874,0.845,0.816,0.788,0.762,0.740,0.717,0.696,0.677)
u.rev <- rev(u)

df <- data.frame(x=x, y=y.rev, lower=l.rev, upper=u.rev)
df$x <- as.factor(df$x)
df

sens <- read_excel("G:/My Drive/InflatedSigInEco/Results/AssumptionFurther2/dataBoxPlotSens.xlsx")

pdf(paste(path_results,"SensetivityISFig.pdf",sep=""),width=4,height=2)
ggplot(sens, aes(letter, is)) +        
  geom_point(size=1.5,color="blue",shape=20,alpha=0.5) +
  geom_errorbar(aes(ymin=ll, ymax=ul),width=0.3,size=0.3,alpha=0.6) +
  ylim(0,1) +
  xlab("Genuine effects") +
  ylab("Extent of selective reporting (95% CI)") +
  theme_bw() +
  theme(axis.title=element_text(size=5),axis.text=element_text(size=4))
dev.off()

#ggsave(paste(path_results,"Boxplot_sensitivity_IS_Fig.pdf",sep=""))




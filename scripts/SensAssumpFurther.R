
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
path_results <- "G:/My Drive/InflatedSigInEco/results/AssumptionFurther/"
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
dim(dat_final)  #70399x13 #73198x13 #73738x13

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
sum(mss)    #70399 #73198 #73738


## ---------------------------------------------------
## Further assumption
## Combination of WAAP & PET-PEESE is used
## ---------------------------------------------------

n.waap <- NULL
#wlsGE_pow100 <- NULL; ppGE_pow100 <- NULL; waapGE_pow100 <- NULL
wlsGE_pow90  <- NULL; ppGE_pow90  <- NULL; waapGE_pow90  <- NULL
#wlsGE_pow80  <- NULL; ppGE_pow80  <- NULL; waapGE_pow80  <- NULL
#wlsGE_pow70  <- NULL; ppGE_pow70  <- NULL; waapGE_pow70  <- NULL
#wlsGE_pow60  <- NULL; ppGE_pow60  <- NULL; waapGE_pow60  <- NULL
#wlsGE_pow50  <- NULL; ppGE_pow50  <- NULL; waapGE_pow50  <- NULL
#wlsGE_pow40  <- NULL; ppGE_pow40  <- NULL; waapGE_pow40  <- NULL
#wlsGE_pow30  <- NULL; ppGE_pow30  <- NULL; waapGE_pow30  <- NULL
#wlsGE_pow20  <- NULL; ppGE_pow20  <- NULL; waapGE_pow20  <- NULL
#wlsGE_pow10  <- NULL; ppGE_pow10  <- NULL; waapGE_pow10 <- NULL
#wlsGE_pow0  <- NULL; ppGE_pow0  <- NULL; waapGE_pow0 <- NULL
#wlsGE_tlarger50 <- NULL; ppGE_tlarger50 <- NULL; waapGE_tlarger50 <- NULL
#wlsGE_tlarger100 <- NULL; ppGE_tlarger100 <- NULL; waapGE_tlarger100 <- NULL


for (i in 1:length(summary(zdat)[,1])) {	
  
  #wls
  wls   <- sum(zdat[[i]]$eff / zdat[[i]]$se^2) / sum(1/zdat[[i]]$se^2)
  wlsGE_pow90 <- c(wlsGE_pow90, rep(wls, length(zdat[[i]]$eff)))
  
  #waap
  se.crit <- abs(wls) / 2.8
  pos <- which(zdat[[i]]$se <= se.crit)
  n.waap[i] <- length(pos) / length(zdat[[i]]$se) #prop. of adequately powered estimates
  
  #at least two adequately powered studies are needed
  if (length(pos) > 1) {
    waapGE_pow90 <- c(waapGE_pow90, rep(sum(zdat[[i]]$eff[pos] / zdat[[i]]$se[pos]^2) / 
                      sum(1/zdat[[i]]$se[pos]^2),length(zdat[[i]]$eff)))
  } else {
    waapGE_pow90 <- c(waapGE_pow90, rep(NA,length(zdat[[i]]$eff)))
  }
  
  #PET-PEESE (PET if alpha is not significant, PEESE if it is significant.)
  if (coefficients(summary(lm(zdat[[i]]$eff ~ zdat[[i]]$se, weights=(1/zdat[[i]]$se)^2)))[1,4] > 0.1) {
    ppGE_pow90 <- c(ppGE_pow90, rep(coefficients(summary(lm(zdat[[i]]$eff ~ zdat[[i]]$se, 
                    weights=(1/zdat[[i]]$se)^2)))[1,1], length(zdat[[i]]$eff)))
  } else {
    se2 <- zdat[[i]]$se^2
    ppGE_pow90 <-c(ppGE_pow90, rep(coefficients(summary(lm(zdat[[i]]$eff ~ se2, 
                    weights=(1/zdat[[i]]$se)^2)))[1,1], length(zdat[[i]]$eff)))
  }
}


## wappGE inserts ppGE whenever there are not at least 
## two adequately powered estimates.
wappGE_pow90 <- waapGE_pow90
for (i in 1:sum(mss)) {
  if (is.na(wappGE_pow90[i])==TRUE) {
    wappGE_pow90[i] <- ppGE_pow90[i]
  }
}

## Robustness (Further assumption)
#0 Full sample
#case   <- "wappGE_pow100"
#caseGE <- wappGE_pow100
#zdat.orig <- zdat

#1 Excluding MAs with power >90% estimates
case   <- "wappGE_pow90"
caseGE <- wappGE_pow90
zdat.orig <- zdat

#2 Excluding MAs with >80% estimates
#case   <- "wappGE_pow80"
#caseGE <- wappGE_pow80
#zdat.orig <- zdat

#3 Excluding MAs with >70% estimates
#case   <- "wappGE_pow70"
#caseGE <- wappGE_pow70
#zdat.orig <- zdat

#4 Excluding MAs with >60% estimates
#case   <- "wappGE_pow60"
#caseGE <- wappGE_pow60
#zdat.orig <- zdat

#5 Excluding MAs with >50% estimates
#case   <- "wappGE_pow50"
#caseGE <- wappGE_pow50
#zdat.orig <- zdat

#6 Excluding MAs with >60% estimates
#case   <- "wappGE_pow40"
#caseGE <- wappGE_pow40
#zdat.orig <- zdat

#7 Excluding MAs with >60% estimates
#case   <- "wappGE_pow30"
#caseGE <- wappGE_pow30
#zdat.orig <- zdat

#8 Excluding MAs with >20% estimates
#case   <- "wappGE_pow20"
#caseGE <- wappGE_pow20
#zdat.orig <- zdat

#9 Excluding MAs with >10% estimates
#case   <- "wappGE_pow10"
#caseGE <- wappGE_pow10
#zdat.orig <- zdat

#10 Excluding all MAs with at least one powered estimate
#case   <- "wappGE_pow0"
#caseGE <- wappGE_pow0
#zdat.orig <- zdat

#11 tlarger50
#case   <- "wappGE_tlarger50"
#caseGE <- wappGE_tlarger50
#zdat.orig <- zdat

#12 tlarger100
#case   <- "wappGE_tlarger100"
#caseGE <- wappGE_tlarger100
#zdat.orig <- zdat

## ------------------------------------
## Add estimated genuine effect to zdat
## ------------------------------------

dat.ideas <- do.call(rbind.data.frame, zdat)
ii <- 1
for (i in unique(dat.ideas$filecode)) {	
  zdat[[ii]]$GE <- caseGE[which(dat.ideas$filecode == i)]
  zdat[[ii]]$n.waap <- n.waap[[ii]]
  zdat[[ii]]$clu <- clu[[ii]]
  ii <- ii + 1
}

View(head(dat.ideas))
dim(dat.ideas)  #70399 x 13

## -----------------------------
## Delete n.waap > 0.9, >0.8 etc
## -----------------------------

# Level of powers
xx90 <- 0.9; xx80 <- 0.8; xx70 <- 0.7
xx60 <- 0.6; xx50 <- 0.5; xx40 <- 0.4
xx30 <- 0.3; xx20 <- 0.2; xx10 <- 0.1; xx0 <- 0

zdat.temp <- list()
jjj <- 1
for (i in 1:length(summary(zdat)[,1])) {
  if(zdat[[i]]$n.waap[1] <= xx90) {
    zdat.temp[[jjj]] <- zdat[[i]]
    jjj <- jjj + 1
  }
}

clu <- NULL
zdat <- zdat.temp
clu <- dat.ideas[-which(n.waap > xx90)]
mss <- mss[-which(n.waap > xx90)]
n.waap <- n.waap[-which(n.waap > xx60)]
length(mss)  #183   #173   #166   #163   #160   #159   #149   #136   #124   #88
sum(mss)     #65796 #63572 #62501 #61697 #61500 #61428 #60538 #57265 #51050 #26622

## After the meta-analyses are dropped
case   <- "wappGE_pow60"
caseGE <- wappGE_pow60
zdat.orig <- zdat

dat.ideas <- do.call(rbind.data.frame, zdat)
ii <- 1
for (i in unique(dat.ideas$filecode)) {	
  zdat[[ii]]$GE <- caseGE[which(dat.ideas$filecode == i)]
  zdat[[ii]]$n.waap <- n.waap[[ii]]
  zdat[[ii]]$clu <- clu[[ii]]
  ii <- ii + 1
}

head(dat.ideas)
dim(dat.ideas) 


## =======================================
# To know exactly which MAs are with >=.90
# Skip this part!
dat90 <- dat.ideas %>%
           subset(n.waap >= 0.90)
head(dat90)
tail(dat90)
dim(dat90) #4613 x 13
length(unique(dat90$filecode)) #10
                               # 7  20  52  75  79  83 219 258 313 350

## -------------------------
## z-value and p-value grids
## -------------------------

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

## z.grid for plot
z.grid.plot  <- seq(-10.25, 10.25, 0.1025)
z.grid.plot2 <- seq(0, 10.25, 0.1025) # this is for actual plotting of absolute values

## ----------------------------------------------------
## Prepare original data for plotting and p-value table
## ----------------------------------------------------

  fac  <- do.call(rbind.data.frame, zdat.orig)
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

cl <- makeCluster(7)
registerDoParallel(cl)
z.plot <- cf(dat = zdat, z.grid = z.grid.plot)
stopCluster(cl)
save(z.plot,file=paste(path_results,case,"_z_plot.RDATA",sep="")) 

## About 19:00 hours with 1000 iterations (7 cores)
#itero <- 1000

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
z.plot.ci <- cf.ci.cluster(dat = zdat, z.grid = z.grid.plot, iters = itero, cluster = clu)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) # about 12 mins.

save(z.plot.ci,file=paste(path_results,case,"_z_plot_ci.RDATA",sep=""))

## --------------------------------
## Sensitivity graph (for all cases) 
## --------------------------------

load(paste(path_results,case,"_z_plot.RDATA",sep=""))
load(paste(path_results,case,"_z_plot_ci.RDATA",sep=""))

## z-plot
pdf(paste(path_results,case,"_zplot_Fig.pdf",sep=""),width=14,height=12)
#par(mfrow = c(2,1))
# z.grid
xs <- z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2
N  <- sum(p.orig.plot) #be careful z.orig goes until 10 and is smaller than the full sample size!,
# but the bootstrapping of the confidence interval is based on the entire sample
# and given in shares. So we need to divide here by the full sample to be consistent.

#means <- apply(cf.wls[[1]] / N, 1, mean)
q025 <- apply(z.plot.ci[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(z.plot.ci[[1]], 2, quantile, probs=c(0.975))

##Quantiles are too narrow using 300 bs observations
matplot(xs, q975, t="l", lty=3, col=4, xlab="z-value", 
        ylab="Frequency", main="",xaxt="n", cex.main=0.9, xlim=c(0,8),ylim=c(0,0.08))	
matplot(xs, z.plot/N, t="o", lty=1, cex=0.8, pch=20, col=4, 
        xlab="z-value", ylab="Frequency", xaxt="n", add=TRUE)	
matplot(xs, q025, t="l", lty=3, col=4, add = TRUE)	

axis(1, c(1.64, 1.96, 2.58))
#axis(1, c(0, 5, 8))
axis(1, seq(0,10,1), label = FALSE)
matplot(xs, z.orig/N,  t="o", pch=20, cex=0.8, lty=2, add=TRUE)
abline(v=qnorm(0.10/2, lower.tail = FALSE), col = 3, lty=2)
abline(v=qnorm(0.05/2, lower.tail = FALSE), col = 2, lty=2)
abline(v=qnorm(0.01/2, lower.tail = FALSE), col = 6, lty=2)
dev.off()

## ---------------------------
## Sensitivity table (for all)
## ---------------------------

## p.grid.tab
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab <- cf(dat = zdat, z.grid = p.grid.tab)
stopCluster(cl)
save(p.tab,file=paste(path_results,case,"_p_tab.RDATA",sep=""))

clu <- unique(dat.ideas$filecode)
itero  <- 500

s.time <- Sys.time()
cl <- makeCluster(7) 
registerDoParallel(cl)
p.tab.ci <- cf.ci.cluster(dat = zdat, z.grid = p.grid.tab, iters = itero, cluster = clu)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  # about 1.5 hrs (7 cores)
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

## sample size
p.table <- rbind(p.table, c(length(summary(zdat)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("$p$ $>$ 0.9", "0.9 $>$ $p$ $>$ 0.8", "0.8 $>$ $p$ $>$ 0.7", 
                       "0.7 $>$ $p$ $>$ 0.6", "0.6 $>$ $p$ $>$ 0.5", 
                       "0.5 $>$ $p$ $>$ 0.4", "0.4 $>$ $p$ $>$ 0.3", 
                       "0.3 $>$ $p$ $>$ 0.2", "0.2 $>$ $p$ $>$ 0.1", 
                       "0.1 $>$ $p$ $>$ 0.05", "0.05 $>$ $p$ $>$ 0.01", 
                       "0.01 $>$ $p$ $>$ 0.001", "0.001 $>$ $p$", 
                       "$IS_{0.10}^{all}$", "$IS_{0.05}^{all}$",
                       "$IS_{0.10}^{sig.}$","$IS_{0.05}^{sig.}$",
                       "No. of meta-analysis","No. of tests")
print(p.table)
write.csv(p.table, paste(path_results,case,"_p.table.csv",sep=""))

## Boxplot of IS (for the main report)
## Load the data
dis <- read_excel("G:/My Drive/InflatedSigInEco/Results/AssumptionFurther/data4BoxplotIS.xlsx")
head(dis)

#dis <-  data.frame(x=cat, y=is, lower=lis, upper=uis)
#dis$cat <- as.factor(dis$cat)
#dis

pdf(paste(path_results,"MainISFig.pdf", sep=""),width=4,height=2)
ggplot(dis, aes(cat, is)) +        
  geom_point(size=1.5,color="blue",shape=20,alpha=0.5) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl),width=0.4,size=0.3,alpha=0.6) +
  ylim(0,0.85) +
  xlab("Sensitivity analyses") +
  ylab("Extent of selective reporting (95% CI)") +
  theme_bw() +
  #scale_fill_few("Medium")+
  theme(axis.title=element_text(size=5),axis.text=element_text(size=4))
dev.off()



## Boxplot of IS (further assumptions)
x <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N")
#x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14) 
y <- c(0.562,0.713,0.461,0.439,0.407,0.394,0.591,0.646,0.700,0.703,0.414,0.608,0.584,0.287)
l <- c(0.517,0.633,0.295,0.415,0.386,0.263,0.537,0.563,0.652,0.634,0.255,0.553,0.534,0.026)
u <- c(0.658,0.781,0.527,0.564,0.545,0.487,0.658,0.737,0.735,0.741,0.458,0.703,0.668,0.445)

df <- data.frame(x=x, y=y, lower=l, upper=u)
df$x <- as.factor(df$x)

pdf(paste(path_results,"SensetivityISFig.pdf",sep=""),width=6,height=4)
ggplot(df, aes(x, y)) +        # ggplot2 plot with CIs
  geom_point(size=2,color="blue",shape=20) +
  geom_errorbar(aes(ymin=lower, ymax=upper),width=0.35,size=0.3,alpha=0.9) +
  ylim(0,0.8) +
  xlab("Sensitivity analyses") +
  ylab("Inflated significance") +
  theme_bw() +
  scale_fill_few("Medium")+
  theme(axis.title=element_text(size=5),axis.text=element_text(size=4))
dev.off()

#ggsave(paste(path_results,"Boxplot_sensitivity_IS_Fig.pdf",sep=""))





rm(list = ls(all = TRUE))                     #clear environment

library(doParallel)  
library(foreach)
library(metafor)  
library(openxlsx)
library(readxl) 
library(dplyr)

# Paths for results and functions
path_results <- "G:/My Drive/InflatedSigInEco/results/CaliperElliott/"
path_fcts    <- "G:/My Drive/InflatedSigInEco/scripts/"

## Functions
## Disclaimer: These functions are from Elliott et al. (2022)
source(paste(path_fcts,"elliott_functions.r",sep="")) 

## Load data
dat_final <- read_excel("G:/My Drive/InflatedSigInEco/data/data_final.xlsx")
dat_final$vi <- dat_final$se^2  #variance
dim(dat_final)

## Exclusion of large t-values
dat_final$tvs <- dat_final$eff / dat_final$se
dat_final <- dat_final %>%
  dplyr::filter(abs(tvs) <= 20)
dim(dat_final)  #70399x13

## --------------------------------------------------------------------

dat_all <- dat_final
dat_all$zval <- dat_all$eff/dat_all$se
dat_all$abs.zval <- abs(dat_all$zval)

## Find p-value for two-tailed test
dat_all$pval <- 2*pnorm(q=dat_all$abs.zval, lower.tail=FALSE)

## ------------
## Caliper test
## ------------

## caliper test function (raw data)
c.test.raw <- function(z.values, cv, calipersize) {
  
  z.values <- z.values[which(z.values < 10)]
  
  #caliper test
  overcaliper  <- length(which(z.values >  cv & z.values < cv+calipersize))
  undercaliper <- length(which(z.values <= cv & z.values > cv-calipersize))
  
  res95 <- binom.test(overcaliper, (overcaliper+undercaliper), p=0.5,
                      alternative="greater", conf.level=0.95)
  
  data.frame(p=res95$estimate, CI95.low = res95$conf.int[1], CI95.up = res95$conf.int[2], 
             pv=res95$p.value, successes=res95$statistic, trials=res95$parameter)
  
}

#critical values
cv1 <- qnorm(0.005, lower.tail = FALSE)
cv5 <- qnorm(0.025, lower.tail = FALSE)
cv10 <- qnorm(0.05, lower.tail = FALSE)

## Table S1
MAIN.results <- as.data.frame(matrix(nrow=12, ncol=7))
colnames(MAIN.results) <- c("alpha", "c", "p", "uc", "oc", "ci95.low","pv")
MAIN.results$alpha <- c(rep(0.1, 4), rep(0.05, 4), rep(0.01, 4))
MAIN.results[,"c"] <- rep(c(0.025, 0.05, 0.1, 0.2), 3)

siglevels <- c(rep(cv10,4),rep(cv5,4),rep(cv1,4))
csizes <- rep(c(0.025,0.05,0.1,0.2),3)

for (i in 1:12) {
  mc.caliper <- c.test.raw(z.values=dat_all$zval, cv=siglevels[i], calipersize=csizes[i])
  
  MAIN.results[i,"uc"] <- round(mc.caliper$trials-mc.caliper$successes, 3) 
  MAIN.results[i,"oc"] <- round(mc.caliper$successes, 3) 
  MAIN.results[i,"p"] <- round(mc.caliper$p, 3) 
  MAIN.results[i,"ci95.low"] <- round(mc.caliper$CI95.low, 3) 
  MAIN.results[i,"pv"] <- mc.caliper$pv
  
}

MAIN.results 

write.csv(MAIN.results,paste(path_results,"caliper_all_Table S1.csv",sep=""))

## ---------------------
## Elliott et al. (2022)
## ---------------------

dat_all_0.15 <- dat_all %>%  
  filter(pval <= 0.15)
P_all <- dat_all_0.15$pval
id <- 1         #no dependence
p_min <- 0
p_max <- 0.15
d_point <- 0.05 #the target cutoff for the discontinuity test
J <- 30         #use 30 bins for CS1 and CS2B tests

LCM_sup_all <- LCM(P_all, p_min, p_max)
CS_1_all <- CoxShi(P_all, id, p_min, p_max, J, 1, 0)  #Test for 1-monotonicity
CS_2B_all <- CoxShi(P_all, id, p_min, p_max, J, 2, 1) #Test for 2-monotonicity and bounds

ELLIOTT.results <- as.data.frame(matrix(nrow=1, ncol=4))
colnames(ELLIOTT.results) <- c('LCM','CS1','CS2B', 'Total obs')
ELLIOTT.results[1,1] <- round(LCM_sup_all, 3)
ELLIOTT.results[1,2] <- round(CS_1_all, 3)
ELLIOTT.results[1,3] <- round(CS_2B_all, 3)
ELLIOTT.results[1,4] <- length(P_all)

ELLIOTT.results
write.csv(ELLIOTT.results, paste(path_results,"elliott_all_Table S2.csv",sep=""))

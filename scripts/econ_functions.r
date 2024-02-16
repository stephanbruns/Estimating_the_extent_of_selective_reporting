## Function to generate counterfactual z-values
cf <- function(dat, z.grid) {

	# This is the loop over the meta-analyses indexed by i
	freqs <- foreach(i = 1:length(summary(dat)[,1])) %dopar% { 
		
		# This is the loop over the observations per meta-analyses indexed by j
		freq <- matrix(ncol = (length(z.grid)-1), nrow = length(dat[[i]]$se))
		for (j in 1:length(dat[[i]]$se)) {
			
			# This is for each distribution/observation the loop over the grid indexed by a
			for (a in 1:(length(z.grid)-1)) {
				
				# This calculates the density between the intervals of the z.grid
				# freq is a matrix with the z.grid intervals as columns and 
				# the observations/distributions per meta-analysis
				freq[j,a] <- (pnorm(z.grid[a+1], mean = dat[[i]]$GE[j] / dat[[i]]$se[j], sd = 1) - 
				 				      pnorm(z.grid[a],   mean = dat[[i]]$GE[j] / dat[[i]]$se[j], sd = 1))
			}
		}
		
		return(freq)
	}
	
	# freqs is a list with each element being the results for one meta-analysis. Each
	# of these lists contains the matrix "freq" (see above).
	freqs <- do.call(rbind.data.frame, freqs)
	freqs <- apply(freqs, 2, sum) # sums the columns for each interval.
	names(freqs) <- z.grid[-1] #some of the names change in an odd way.
	
	# This takes the sum of the first and last column and the second and the second last
	# and so on. Therefore, it is essential that the z.grid is symmetric around zero.
	abs.freqs <- NULL
	for (ii in 0:((length(freqs) / 2) -1)) {
	
		abs.freqs[(length(freqs) / 2)-ii] <- sum(freqs[1+ii], freqs[length(freqs)-ii])
	}
	names(abs.freqs) <- z.grid[which(z.grid > 0)] #some of the names change in an odd way.
	
	# abs.freqs returns a vector with the cumulated probability mass for z.grid - mirrored at zero
	# The columns names are the upper limits of the intervals given by z.grid, but only for positive
	# values as the frequencies are mirrored at zero to obtain p-values for two sided tests
	# or absolute z-values for visualization.
	return(abs.freqs)
} 

	
cf.ci.cluster <- function(dat, z.grid, iters, cluster) {
  
  # This is the loop over iterations indexed by dd
  freqs <- foreach(dd = 1:iters) %dopar% {
    
    # Within each foreach loop the random numbers are set in a reproducible way.
    set.seed(dd+12)
    
    #resampling of studies which may contain multiple meta-analyses
    bs.cluster <- sample(unique(cluster), replace=TRUE)
    
    boots <- NULL
    for (ii in 1:length(bs.cluster)) {
      boots <- c(boots,	which(cluster == bs.cluster[ii]))
    }
    
    k<-1
    # This is the loop over meta-analyses indexed by i
    freq.per.iter <- list()
    sig.1 <- NULL
    sig.5 <- NULL
    for (i in boots) {
      
      # This is the loop over the observations per meta-analyses indexed by j
      freq <- matrix(ncol = (length(z.grid)-1), nrow = length(dat[[i]]$se))
      for (j in 1:length(dat[[i]]$se)) {
        
        # This is for each distribution/observation the loop over the grid indexed by a
        for (a in 1:(length(z.grid)-1)) {
          
          # This calculates the density between the intervals of the z.grid
          # freq is a matrix with the z.grid intervals as columns and 
          # the observations/distributions per meta-analysis as rows
          freq[j,a] <- (pnorm(z.grid[a+1], mean = dat[[i]]$GE[j] / dat[[i]]$se[j], sd = 1) - 
                        pnorm(z.grid[a], mean = dat[[i]]$GE[j] / dat[[i]]$se[j], sd = 1))
        }
        
        sig.1.temp <- length(which(abs(dat[[i]]$eff / dat[[i]]$se) > abs(qnorm(0.1/2))))
        sig.5.temp <- length(which(abs(dat[[i]]$eff / dat[[i]]$se) > abs(qnorm(0.05/2))))
        
      }

      sig.1 <- c(sig.1, sig.1.temp)
      sig.5 <- c(sig.5, sig.5.temp)
      freq.per.iter[[k]] <- freq 
      k <- k+1
    }
    
    # freq.per.iter is a list with each entry being the freq matrix (see above) for one
    # meta-analysis
    freq.per.iter <- do.call(rbind.data.frame, freq.per.iter)
    freq.per.iter.all  <- apply(freq.per.iter, 2, sum) / dim(freq.per.iter)[1]# sums the columns for each interval and divides by the number of observations
    freq.per.iter.sig1 <- apply(freq.per.iter, 2, sum) / sum(sig.1)
    freq.per.iter.sig5 <- apply(freq.per.iter, 2, sum) / sum(sig.5)
    
    names(freq.per.iter.all)  <- z.grid[-1] #some of the names change in an odd way.
    names(freq.per.iter.sig1) <- z.grid[-1]
    names(freq.per.iter.sig5) <- z.grid[-1]
    
    # freq.per.iter is one row with the counterfactual frequencies per grid
    return(list(freq.per.iter.all, freq.per.iter.sig1, freq.per.iter.sig5))
  }
  
  # freqs is a list with each element being the outcome of one iteration (freq.per.iter - see above)
  freqs.all  <- NULL
  freqs.sig1 <- NULL
  freqs.sig5 <- NULL
  for (jj in 1:iters) {
    freqs.all[jj]  <- freqs[[jj]][1]
    freqs.sig1[jj] <- freqs[[jj]][2]
    freqs.sig5[jj] <- freqs[[jj]][3]
  }
  
  freqs.all  <- do.call(rbind.data.frame, freqs.all)	
  freqs.sig1 <- do.call(rbind.data.frame, freqs.sig1)
  freqs.sig5 <- do.call(rbind.data.frame, freqs.sig5)
  names(freqs.all)  <- z.grid[-1] #some of the names change in an odd way.
  names(freqs.sig1) <- z.grid[-1]
  names(freqs.sig5) <- z.grid[-1]
  
  # This takes the sum of the first and last column and the second and the second last
  # and so on. Therefore, it is essential that the z.grid is symmetric around zero.
  abs.freqs.all  <- matrix(ncol=(dim(freqs.all)[2] / 2),  nrow=iters)
  abs.freqs.sig1 <- matrix(ncol=(dim(freqs.sig1)[2] / 2), nrow=iters)
  abs.freqs.sig5 <- matrix(ncol=(dim(freqs.sig5)[2] / 2), nrow=iters)
  
  for (k in 0:((dim(freqs.all)[2] / 2) -1)) {
    abs.freqs.all[,(dim(freqs.all)[2] / 2)-k]   <- apply(matrix(c(freqs.all[,1+k],  freqs.all[,dim(freqs.all)[2]-k]),   ncol=2), 1, sum)
    abs.freqs.sig1[,(dim(freqs.sig1)[2] / 2)-k] <- apply(matrix(c(freqs.sig1[,1+k], freqs.sig1[,dim(freqs.sig1)[2]-k]), ncol=2), 1, sum)
    abs.freqs.sig5[,(dim(freqs.sig5)[2] / 2)-k] <- apply(matrix(c(freqs.sig5[,1+k], freqs.sig5[,dim(freqs.sig5)[2]-k]), ncol=2), 1, sum)
  }
  colnames(abs.freqs.all)  <- z.grid[which(z.grid > 0)] #some of the names change in an odd way.
  colnames(abs.freqs.sig1) <- z.grid[which(z.grid > 0)]
  colnames(abs.freqs.sig5) <- z.grid[which(z.grid > 0)]
  
  return(list(abs.freqs.all, abs.freqs.sig1, abs.freqs.sig5))
}	

## Functions associated to the random-effects model which takes into account the 
## heterogeneity in each meta-analysis when generating counterfactuals

cf.rem <- function(dat, z.grid) {
  
  # This is the loop over the meta-analyses indexed by i
  freqs <- foreach(i = 1:length(summary(dat)[,1])) %dopar% { 
    
    # This is the loop over the observations per meta-analyses indexed by j
    freq <- matrix(ncol = (length(z.grid)-1), nrow = length(dat[[i]]$se))
    for (j in 1:length(dat[[i]]$se)) {
      
      # This is for each distribution/observation the loop over the grid indexed by a
      for (a in 1:(length(z.grid)-1)) {
        
        # This calculates the density between the intervals of the z.grid
        # freq is a matrix with the z.grid intervals as columns and 
        # the observations/distributions per meta-analysis
        freq[j,a] <- (pnorm(z.grid[a+1], mean = dat[[i]]$GE[j] / sqrt(dat[[i]]$vi[j]+dat[[i]]$tau2[j]), sd = 1) - 
                      pnorm(z.grid[a], mean = dat[[i]]$GE[j] / sqrt(dat[[i]]$vi[j]+dat[[i]]$tau2[j]), sd = 1))
      }
    }
    
    return(freq)
  }
  
  # freqs is a list with each element being the results for one meta-analysis. Each
  # of these lists contains the matrix "freq" (see above).
  freqs <- do.call(rbind.data.frame, freqs)
  freqs <- apply(freqs, 2, sum) # sums the columns for each interval.
  names(freqs) <- z.grid[-1] #some of the names change in an odd way.
  
  # This takes the sum of the first and last column and the second and the second last
  # and so on. Therefore, it is essential that the z.grid is symmetric around zero.
  abs.freqs <- NULL
  for (ii in 0:((length(freqs) / 2) -1)) {
    
    abs.freqs[(length(freqs) / 2)-ii] <- sum(freqs[1+ii], freqs[length(freqs)-ii])
  }
  names(abs.freqs) <- z.grid[which(z.grid > 0)] #some of the names change in an odd way.
  
  # abs.freqs returns a vector with the cumulated probability mass for z.grid - mirrored at zero
  # The columns names are the upper limits of the intervals given by z.grid, but only for positive
  # values as the frequencies are mirrored at zero to obtain p-values for two sided tests
  # or absolute z-values for visualization.
  return(abs.freqs)
} 

cf.ci.cluster.rem <- function(dat, z.grid, iters, cluster) {
  
  # This is the loop over iterations indexed by dd
  freqs <- foreach(dd = 1:iters) %dopar% {
    
    # Within each foreach loop the random numbers are set in a reproducible way.
    set.seed(dd+12)
    
    #resampling of studies which may contain multiple meta-analyses
    bs.cluster <- sample(unique(cluster), replace=TRUE)
    
    boots <- NULL
    for (ii in 1:length(bs.cluster)) {
      boots <- c(boots,	which(cluster == bs.cluster[ii]))
    }
    
    k<-1
    # This is the loop over meta-analyses indexed by i
    freq.per.iter <- list()
    sig.1 <- NULL
    sig.5 <- NULL
    for (i in boots) {
      
      # This is the loop over the observations per meta-analyses indexed by j
      freq <- matrix(ncol = (length(z.grid)-1), nrow = length(dat[[i]]$se))
      for (j in 1:length(dat[[i]]$se)) {
        
        # This is for each distribution/observation the loop over the grid indexed by a
        for (a in 1:(length(z.grid)-1)) {
          
          # This calculates the density between the intervals of the z.grid
          # freq is a matrix with the z.grid intervals as columns and 
          # the observations/distributions per meta-analysis as rows
          freq[j,a] <- (pnorm(z.grid[a+1], mean = dat[[i]]$GE[j] / sqrt(dat[[i]]$vi[j]+dat[[i]]$tau2[j]), sd = 1) - 
                          pnorm(z.grid[a], mean = dat[[i]]$GE[j] / sqrt(dat[[i]]$vi[j]+dat[[i]]$tau2[j]), sd = 1))
          
        }
        
        sig.1.temp <- length(which(abs(dat[[i]]$eff / dat[[i]]$se) > abs(qnorm(0.1/2))))
        sig.5.temp <- length(which(abs(dat[[i]]$eff / dat[[i]]$se) > abs(qnorm(0.05/2))))
        
      }
      
      sig.1 <- c(sig.1, sig.1.temp)
      sig.5 <- c(sig.5, sig.5.temp)
      freq.per.iter[[k]] <- freq 
      k <- k+1
    }
    
    # freq.per.iter is a list with each entry being the freq matrix (see above) for one
    # meta-analysis
    freq.per.iter <- do.call(rbind.data.frame, freq.per.iter)
    freq.per.iter.all  <- apply(freq.per.iter, 2, sum) / dim(freq.per.iter)[1]# sums the columns for each interval and divides by the number of observations
    freq.per.iter.sig1 <- apply(freq.per.iter, 2, sum) / sum(sig.1)
    freq.per.iter.sig5 <- apply(freq.per.iter, 2, sum) / sum(sig.5)
    
    names(freq.per.iter.all)  <- z.grid[-1] #some of the names change in an odd way.
    names(freq.per.iter.sig1) <- z.grid[-1]
    names(freq.per.iter.sig5) <- z.grid[-1]
    
    # freq.per.iter is one row with the counterfactual frequencies per grid
    return(list(freq.per.iter.all, freq.per.iter.sig1, freq.per.iter.sig5))
  }
  
  # freqs is a list with each element being the outcome of one iteration (freq.per.iter - see above)
  freqs.all  <- NULL
  freqs.sig1 <- NULL
  freqs.sig5 <- NULL
  for (jj in 1:iters) {
    freqs.all[jj]  <- freqs[[jj]][1]
    freqs.sig1[jj] <- freqs[[jj]][2]
    freqs.sig5[jj] <- freqs[[jj]][3]
  }
  
  freqs.all  <- do.call(rbind.data.frame, freqs.all)	
  freqs.sig1 <- do.call(rbind.data.frame, freqs.sig1)
  freqs.sig5 <- do.call(rbind.data.frame, freqs.sig5)
  names(freqs.all)  <- z.grid[-1] #some of the names change in an odd way.
  names(freqs.sig1) <- z.grid[-1]
  names(freqs.sig5) <- z.grid[-1]
  
  # This takes the sum of the first and last column and the second and the second last
  # and so on. Therefore, it is essential that the z.grid is symmetric around zero.
  abs.freqs.all  <- matrix(ncol=(dim(freqs.all)[2] / 2),  nrow=iters)
  abs.freqs.sig1 <- matrix(ncol=(dim(freqs.sig1)[2] / 2), nrow=iters)
  abs.freqs.sig5 <- matrix(ncol=(dim(freqs.sig5)[2] / 2), nrow=iters)
  
  for (k in 0:((dim(freqs.all)[2] / 2) -1)) {
    abs.freqs.all[,(dim(freqs.all)[2] / 2)-k]   <- apply(matrix(c(freqs.all[,1+k],  freqs.all[,dim(freqs.all)[2]-k]),   ncol=2), 1, sum)
    abs.freqs.sig1[,(dim(freqs.sig1)[2] / 2)-k] <- apply(matrix(c(freqs.sig1[,1+k], freqs.sig1[,dim(freqs.sig1)[2]-k]), ncol=2), 1, sum)
    abs.freqs.sig5[,(dim(freqs.sig5)[2] / 2)-k] <- apply(matrix(c(freqs.sig5[,1+k], freqs.sig5[,dim(freqs.sig5)[2]-k]), ncol=2), 1, sum)
  }
  colnames(abs.freqs.all)  <- z.grid[which(z.grid > 0)] #some of the names change in an odd way.
  colnames(abs.freqs.sig1) <- z.grid[which(z.grid > 0)]
  colnames(abs.freqs.sig5) <- z.grid[which(z.grid > 0)]
  
  return(list(abs.freqs.all, abs.freqs.sig1, abs.freqs.sig5))
}	

## This is a separate function for the robustness check with increased standard errors. In the original function, we bootstrap shares,
## this means we divide in each bs iteration by the number of all p-values or all significant p-values in the given iteration.
## For the case with increased standard errors, it is important to divide by all significant p-values without inflating the 
## standard error. 

cf.ci.cluster.se <- function(dat, z.grid, iters, cluster, dat.orig) {
  
  # This is the loop over iterations indexed by dd
  freqs <- foreach(dd = 1:iters) %dopar% {
    
    # Within each foreach loop the random numbers are set in a reproducible way.
    set.seed(dd+12)
  
    #resampling of studies which may contain multiple meta-analyses
    bs.cluster <- sample(unique(cluster), replace=TRUE)
    
    boots <- NULL
    for (ii in 1:length(bs.cluster)) {
      boots <- c(boots,	which(cluster == bs.cluster[ii]))
    }
    
    # This is the loop over meta-analyses indexed by i
    freq.per.iter <- list()
    sig.1 <- NULL
    sig.5 <- NULL
    k <- 1
    
    for (i in boots) {
      
      freq <- matrix(ncol = (length(z.grid)-1), nrow = length(dat[[i]]$se))
      
      # This is the loop over the observations per meta-analyses indexed by j
      for (j in 1:length(dat[[i]]$se)) {
        
        # This is for each distribution/observation the loop over the grid indexed by a
        for (a in 1:(length(z.grid)-1)) {
          
          # This calculates the density between the intervals of the z.grid
          # freq is a matrix with the z.grid intervals as columns and 
          # the observations/distributions per meta-analysis as rows
          freq[j,a] <- (pnorm(z.grid[a+1], mean = dat[[i]]$GE[j] / dat[[i]]$se[j], sd = 1) - 
                        pnorm(z.grid[a],   mean = dat[[i]]$GE[j] / dat[[i]]$se[j], sd = 1))
        }
        
        sig.1.temp <- length(which(abs(dat.orig[[i]]$eff / dat.orig[[i]]$se) > abs(qnorm(0.10/2))))
        sig.5.temp <- length(which(abs(dat.orig[[i]]$eff / dat.orig[[i]]$se) > abs(qnorm(0.05/2))))
      }
      
      sig.1 <- c(sig.1, sig.1.temp)
      sig.5 <- c(sig.5, sig.5.temp)
      freq.per.iter[[k]] <- freq #for each MA
      k <- k+1
    }
    # freq.per.iter is a list with each entry being the freq matrix (see above) for one
    # meta-analysis
    freq.per.iter <- do.call(rbind.data.frame, freq.per.iter)
    # sums the columns for each interval and divides by the number of observations
    freq.per.iter.all  <- apply(freq.per.iter, 2, sum) / dim(freq.per.iter)[1]
    freq.per.iter.sig1 <- apply(freq.per.iter, 2, sum) / sum(sig.1)
    freq.per.iter.sig5 <- apply(freq.per.iter, 2, sum) / sum(sig.5)
    
    names(freq.per.iter.all)  <- z.grid[-1] #some of the names change in an odd way.
    names(freq.per.iter.sig1) <- z.grid[-1]
    names(freq.per.iter.sig5) <- z.grid[-1]
    
    # freq.per.iter is one row with the counterfactual frequencies per grid
    return(list(freq.per.iter.all, freq.per.iter.sig1, freq.per.iter.sig5))
  }
  # freqs is a list with each element being the outcome of one iteration (freq.per.iter - see above)
  freqs.all  <- NULL
  freqs.sig1 <- NULL
  freqs.sig5 <- NULL
  
  for (jj in 1:iters) {
    freqs.all[jj]  <- freqs[[jj]][1]
    freqs.sig1[jj] <- freqs[[jj]][2]
    freqs.sig5[jj] <- freqs[[jj]][3]
  }
  
  freqs.all  <- do.call(rbind.data.frame, freqs.all)	
  freqs.sig1 <- do.call(rbind.data.frame, freqs.sig1)
  freqs.sig5 <- do.call(rbind.data.frame, freqs.sig5)
  names(freqs.all)  <- z.grid[-1] #some of the names change in an odd way.
  names(freqs.sig1) <- z.grid[-1]
  names(freqs.sig5) <- z.grid[-1]
  
  # This takes the sum of the first and last column and the second and the second last
  # and so on. Therefore, it is essential that the z.grid is symmetric around zero.
  abs.freqs.all  <- matrix(ncol=(dim(freqs.all)[2] / 2),  nrow=iters)
  abs.freqs.sig1 <- matrix(ncol=(dim(freqs.sig1)[2] / 2), nrow=iters)
  abs.freqs.sig5 <- matrix(ncol=(dim(freqs.sig5)[2] / 2), nrow=iters)
  
  for (k in 0:((dim(freqs.all)[2] / 2) -1)) {
    abs.freqs.all[,(dim(freqs.all)[2] / 2)-k]   <- apply(matrix(c(freqs.all[,1+k],  freqs.all[,dim(freqs.all)[2]-k]),   ncol=2), 1, sum)
    abs.freqs.sig1[,(dim(freqs.sig1)[2] / 2)-k] <- apply(matrix(c(freqs.sig1[,1+k], freqs.sig1[,dim(freqs.sig1)[2]-k]), ncol=2), 1, sum)
    abs.freqs.sig5[,(dim(freqs.sig5)[2] / 2)-k] <- apply(matrix(c(freqs.sig5[,1+k], freqs.sig5[,dim(freqs.sig5)[2]-k]), ncol=2), 1, sum)
  }
  colnames(abs.freqs.all)  <- z.grid[which(z.grid > 0)] #some of the names change in an odd way.
  colnames(abs.freqs.sig1) <- z.grid[which(z.grid > 0)]
  colnames(abs.freqs.sig5) <- z.grid[which(z.grid > 0)]
  
  return(list(abs.freqs.all, abs.freqs.sig1, abs.freqs.sig5)) 
}	


## Functions associated to multiple genuine effects - Monte Carlo simulation

cf.mc <- function(dat, z.grid, GE.type, mc.iter) {

	# loop over MC iterations
	abs.freqs.mc <- list()
	for (eee in 1:mc.iter) {
	
		# This is the loop over the meta-analyses indexed by i
		freqs <- foreach(i = 1:length(summary(dat)[,1])) %dopar% { 
			
			set.seed(eee*i*1.1) #reproducable
			
			# random split - two genuine effects
			if (GE.type == 1) { 
		
				#split randomly and calculate genuine effect
				sub.size <- floor(length(dat[[i]]$se) / 2) 
				sub.ints <- c(seq(1, length(dat[[i]]$se), sub.size)[1:2], 
				              length(dat[[i]]$se)+1) #last +1 because sub.ids are 
				                                     #determined always -1 to the next interval border
				
				mc.draws <- sample(1:length(dat[[i]]$se))
				mc1 <- mc.draws[sub.ints[1]:(sub.ints[2]-1)]
				mc2 <- mc.draws[sub.ints[2]:(sub.ints[3]-1)]
				
				wls1 <- sum(dat[[i]]$eff[mc1] / dat[[i]]$se[mc1]^2) / sum(1/dat[[i]]$se[mc1]^2)
				wls2 <- sum(dat[[i]]$eff[mc2] / dat[[i]]$se[mc2]^2) / sum(1/dat[[i]]$se[mc2]^2)
				
				means <- c(wls1 / dat[[i]]$se[mc1], wls2 / dat[[i]]$se[mc2]) 
				
				# This is the loop over the observations per meta-analyses indexed by j
				freq <- matrix(ncol = (length(z.grid)-1), nrow = length(dat[[i]]$se))
				for (j in 1:length(dat[[i]]$se)) {
					
					# This is for each distribution/observation the loop over the grid indexed by a
					for (a in 1:(length(z.grid)-1)) {
						
						# This calculates the density between the intervals of the z.grid
						# freq is a matrix with the z.grid intervals as columns and 
						# the observations/distributions per meta-analysis
						freq[j,a] <- (pnorm(z.grid[a+1], mean = means[j], sd = 1) - 
										      pnorm(z.grid[a],   mean = means[j], sd = 1))
					}
				}
			}	
			
			# random split - three genuine effects
			if (GE.type == 2) { 
		
				#split randomly and calculate genuine effect
				sub.size <- floor(length(dat[[i]]$se) / 3) 
				sub.ints <- c(seq(1, length(dat[[i]]$se), sub.size)[1:3], length(dat[[i]]$se)+1) #last +1 because sub.ids are determined always -1 to the next interval border
				
				mc.draws <- sample(1:length(dat[[i]]$se))
				mc1 <- mc.draws[sub.ints[1]:(sub.ints[2]-1)]
				mc2 <- mc.draws[sub.ints[2]:(sub.ints[3]-1)]
				mc3 <- mc.draws[sub.ints[3]:(sub.ints[4]-1)]
				
				wls1 <- sum(dat[[i]]$eff[mc1] / dat[[i]]$se[mc1]^2) / sum(1/dat[[i]]$se[mc1]^2)
				wls2 <- sum(dat[[i]]$eff[mc2] / dat[[i]]$se[mc2]^2) / sum(1/dat[[i]]$se[mc2]^2)
				wls3 <- sum(dat[[i]]$eff[mc3] / dat[[i]]$se[mc3]^2) / sum(1/dat[[i]]$se[mc3]^2)
				
				means <- c(wls1 / dat[[i]]$se[mc1], wls2 / dat[[i]]$se[mc2],  wls3 / dat[[i]]$se[mc3]) 
				
				# This is the loop over the observations per meta-analyses indexed by j
				freq <- matrix(ncol = (length(z.grid)-1), nrow = length(dat[[i]]$se))
				for (j in 1:length(dat[[i]]$se)) {
					
					# This is for each distribution/observation the loop over the grid indexed by a
					for (a in 1:(length(z.grid)-1)) {
						
						# This calculates the density between the intervals of the z.grid
						# freq is a matrix with the z.grid intervals as columns and 
						# the observations/distributions per meta-analysis
						freq[j,a] <- (pnorm(z.grid[a+1], mean = means[j], sd = 1) - 
									      	pnorm(z.grid[a],   mean = means[j], sd = 1))
					}
				}
			}
			
			return(freq)
		}
		
		# freqs is a list with each element being the results for one meta-analysis. 
		# Each of these lists contains the matrix "freq" (see above).
		freqs <- do.call(rbind.data.frame, freqs)
		freqs <- apply(freqs, 2, sum) # sums the columns for each interval.
		names(freqs) <- z.grid[-1] #some of the names change in an odd way.
		
		# This takes the sum of the first and last column and the second and the second last
		# and so on. Therefore, it is essential that the z.grid is symmetric around zero.
		abs.freqs <- NULL
		for (ii in 0:((length(freqs) / 2) -1)) {
			abs.freqs[(length(freqs) / 2)-ii] <- sum(freqs[1+ii], freqs[length(freqs)-ii])
		}
		names(abs.freqs) <- z.grid[which(z.grid > 0)] #some of the names change in an odd way.
			
		# abs.freqs returns a vector with the cumulated probability mass for z.grid - mirrored at zero
		# The columns names are the upper limits of the intervals given by z.grid, but only for positive
		# values as the frequencies are mirrored at zero to obtain p-values for two sided tests
		# or absolute z-values for visualization.
		abs.freqs.mc[[eee]] <- abs.freqs
	}
	
	return(abs.freqs.mc)
} 










# Generate permutation distributions of t values for independent groups
## cond1 and cond2 are Nt trials x Nf time frames matrices
permtdist <- function(cond1, cond2, Nt, Nf, nboot = 2000){
  
  st2 <- Nt+1 # index to start condition 2 
  en2 <- Nt+Nt # index to end condition 2
  all <- rbind(cond1, cond2)
  
  perm.tvals <- matrix(data = 0, nrow = nboot, ncol = Nf)
  # perm.pvals <- matrix(data = 0, nrow = nboot, ncol = Nf)

  for(B in 1:nboot){
    
    permdata <- all[sample(1:en2, size = en2, replace = FALSE),]
    permcond1 <- permdata[1:Nt,]
    permcond2 <- permdata[st2:en2,]

    perm.tvals[B,] <- Rfast::ttests(permcond1, permcond2, paired=FALSE)[,1]    
    # for(F in 1:Nf){ # for each time frame:
    #   perm.tvals[B,F] <- t.test(permcond1[,F], permcond2[,F])$statistic
    # }
  }
  # list(tvals = perm.tvals, pvals = perm.pvals)
  perm.tvals
}

# Cluster correction

## Form clusters using binary vector: pvals < alpha
## clusters must be at least 2 time frames
cluster.make <- function(x){
  y <- rle(x)
  cmap <- vector(mode = "numeric", length = 0)
  nC <- length(y$values) # number of clusters
  indx <- 0 # cluster counter
  for(CL in 1:nC){
    if(y$values[CL] == 0 || y$lengths[CL] == 1){
      val <- 0
    } else {
      indx <- indx + 1
      val <- indx
    }
    cmap <- c(cmap, rep(val, y$lengths[CL]))
  }
  cmap
}

## Save sum for each cluster
# values = statistics (t values)
# cmap = output from cluster.make
cluster.sum <- function(values, cmap){
  csum <- vector(mode = "numeric", length = max(cmap))
  if(max(cmap)>0){
    for(CL in 1:max(cmap)){
      csum[CL] <- sum(values[cmap==CL])
    }
  } else {
    csum <- 0
  }
  csum
}

# Cluster test
# values = statistics (t values)
# cmap = output from cluster.make
# boot.th = bootstrap quantile of distribution of max cluster sums
cluster.test <- function(values, cmap, boot.th){
    csig <- vector(mode = "logical", length = length(cmap))
  if(max(cmap)>0){
    for(CL in 1:max(cmap)){
      csig[cmap==CL] <- sum(values[cmap==CL]) > boot.th
    }
  } else {
    csig <- FALSE
  }
  csig
}

sim.counter <- function(S, nsim, inc){
  if(S == 1){
    # print(paste(nsim,"iterations:",S))
    cat(nsim,"iterations:",S)
    beep(2)
  }
  if(S %% inc == 0){
    # print(paste("iteration",S,"/",nsim))
    cat(" /",S)
    beep(2)
  }
}
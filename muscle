# The functions related to the article
# "MUltivariate Sparse CLustering for Extremes"
# by Nicolas Meyer and Olivier Wintenberger

# rm(list=ls())
library(rlang)

vect_of_thresholds <- function(v){
  # This function gives a vector listing the ordered thresholds for which
  # a vector v has a different number of nonnull coordinates
  d <- length(v)
  u <- sort(v, TRUE)
  threshold <- rep(0,d)
  su <- cumsum(u)
  
  for(k in 1:(d-1)){
    threshold[k] <- su[k]-k*u[k+1]}
  threshold[d] <- su[d] # threshold above which the projection does not make sense
  
  return(threshold)
}

cones_function <- function(v, Gamma){
  # This function gives the subsets which contains the mass of Z
  # with Gamma by decreasing order
  d <- length(v)
  ord <- order(v)
  thres <- vect_of_thresholds(v)
  length_Gamma <- length(Gamma)
  
  J <- which(Gamma < thres[d], arr.ind = TRUE) # the coordinates of the useful Gamma
  length_J <- length(J)
  # the one that are too big provide a vector of NA
  
  cones <- matrix(NA, nrow = d, ncol = length_Gamma)
  if (length_J > 0){
    # Initialization
    j <- J[1]
    r1 <- sum(thres<Gamma[j]) + 1 # number of positive coordinates for the threshold Gamma[1]
    cones[ord[1:(d-r1)], j] <- 0 # the coordinates on which v is not projected
    cones[ord[seq2(d-r1+1, d)], j] <- 1 # the coordinates on which v is projected
    
    j <- J[2]
    r2 <- d+2 # just need to be above d+1 I think
  }
  
  if (length_J > 1){
    # the loop
    while ((r2 >= 1)&&(j<=length_Gamma)) {
      r2 <- sum(thres<Gamma[j]) + 1 # number of positive coordinates for the threshold Gamma[1]
      if (r2 <= d){
        cones[ , j] <- cones[ , j-1]
        cones[ord[seq2(d-r1+1, d-r2)], j] <- 0 # seq2 only works for increasing sequence which helps if r1 < r2
      } # we let the NA if the Gamma[1] is too huge
      j <- j + 1
      r1 <- r2
    }
  }
  return(cones)
}

occ_subsets <- function(testmtx){
  # This function takes a matrix and gives the number of occurence of each column
  nc <- ncol(testmtx)
  occ <- seq(nc)  
  for (i in seq(nc-1)) {
    dup <- colSums( abs( testmtx[,seq(i+1,nc),drop=F] - testmtx[,i] ) ) == 0
    occ[which(dup)+i] <- occ[i]
  }
  result <- rbind(as.matrix(testmtx[ ,unique(occ)]), table(occ))
  colnames(result) <- NULL
  return(result)
}

muscle_plot <- function(X, prop){
  # This function gives for each k the value of s_hat which minimizes the KL
  # and the value of this minimizer
  n <- ncol(X)
  d <- nrow(X)
  p <- length(prop) # prop should be given in an increasing order
  
  sum_norm <- apply(X, 2, sum)
  sort_sum_norm <- sort(sum_norm, decreasing = TRUE)
  thresholds <- sort_sum_norm[round(n*prop+1, 0)] # this corresponds to the list of thresholds u_n
  # we use the function 'round' since otherwise there are some approximations in n*Proportion
  X <- X[ , sum_norm > thresholds[p]] # this is done only to avoid computing too many NA's
  
  length_thresh <- length(thresholds)
  
  # Initialization: the "adjacency matrix"
  binary <- apply(X, 2, cones_function, thresholds)
  s_hat <- rep(NA, length_thresh)
  result <- matrix(NA, ncol = 4)
  colnames(result) <- c("k", "s_tilde", "hat_s", "minimizer")
  
  for (j in 1:length_thresh){
    binary_bis <- binary[( (j-1)*d + 1 ): (j*d) , !is.na(binary[(j-1)*d+1, ]) ] # we remove the columns of NA of the considered block
    k <- ncol(binary_bis) # it corresponds to k_n
    M <- occ_subsets(binary_bis)
    r <- ncol(M)
    
    # We order the matrix M in increasing order regarding the occurence of the columns
    ordered_subcones <- order(M[d+1, ], decreasing=TRUE)
    M <- as.matrix(M[ , ordered_subcones])
    T <- M[(d+1), ]
    
    # we now optimize in s_hat
    optim <- 1:(r-1) - lfactorial(k) + k*log(k) + sum(lfactorial(T)) - cumsum(T[-r]*(log(T)[-r])) - (k-cumsum(T[-r]))*log( (k-cumsum(T[-r])) / ((r-1):1) )
    s_hat <- which.min(optim)
    
    minimizer <- (optim[s_hat])/k + k/n
    result <- rbind(result, c(k, r, s_hat, minimizer))
  }

  result <- result[-1, , drop=F] # we have to remove the first row of NA's
  return(result)
}

##########################################################
muscle_clusters <- function(X, prop){
  # This function uses the previous algo, chooses the optimal k,
  # and gives the associated subsets on which the angular vector Z puts mass
  n <- ncol(X)
  d <- nrow(X)
  p <- length(prop) # prop should be given in an increasing order
  
  sum_norm <- apply(X, 2, sum)
  sort_sum_norm <- sort(sum_norm, decreasing = TRUE)
  
  thresholds <- sort_sum_norm[round(n*prop+1, 0)] # this corresponds to the list of thresholds u_n
  # we use the function 'round' since otherwise there are some approximations in n*Proportion
  X <- X[ , sum_norm > thresholds[p]] # this is done only to avoid computing too many NA's
  
  length_thresh <- length(thresholds)
  
  # Initialization: the "adjacency matrix"
  binary <- apply(X, 2, cones_function, thresholds)
  minimizer <- Inf
  s_hat <- NA
  k_hat <- NA
  
  for ( j in 1:length_thresh ){
    binary_bis <- binary[( (j-1)*d + 1 ): (j*d) , !is.na(binary[(j-1)*d+1, ]) ] # we remove the columns of NA of the considered block
    k <- ncol(binary_bis) # it corresponds to k_n
    M <- occ_subsets(binary_bis)
    r <- ncol(M)
    
    # We order the matrix M in increasing order regarding the occurence of the columns
    ordered_subcones <- order(M[d+1, ], decreasing=TRUE)
    M <- as.matrix(M[ , ordered_subcones])
    T <- M[(d+1), ]
    
    # we now optimize in s_hat
    # we now optimize in s_hat
    optim <- 1:(r-1) - lfactorial(k) + k*log(k) + sum(lfactorial(T)) - cumsum(T[-r]*(log(T)[-r])) - (k-cumsum(T[-r]))*log( (k-cumsum(T[-r])) / ((r-1):1) )
    min_loc <- optim[which.min(optim)]/k + k/n

    if (minimizer > min_loc){
      minimizer <- min_loc
      s_hat <- which.min(optim)
      k_hat <- k
      Mat <- M[, 1:s_hat]
    }
  }
  M <- as.matrix(Mat)
  u <- sort_sum_norm[round(k_hat + 1, 0)]
  weights <- M[(d+1), ]/ sum(M[(d+1), ])
  return(list(M, k_hat, u, s_hat, weights))
}


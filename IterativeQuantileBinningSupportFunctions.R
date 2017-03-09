### Helper functions for Iterative Quantile Binning Functions
# Move to this seperate file for organization that focuses on primary tasks

#--------------------------------------
## function to make a matrix for duplicating the N rows of a matrix M times each
# support funciton for IQ binning function
#!# Update later for speed 
make_stack_matrix <- function(N,M){ 
  mat <- unname(model.matrix(~as.factor(rep(1:N,each=M))-(1)))
  attributes(mat)[2:3]<-NULL
  return(mat)
} 
make_stack_matrix(3,4)

#--------------------------------------
### Helper function for checking if a vector is in a p-dimensional bin, defined by 2*p boundaries
# x = p-dimensional vector
# bin_bounds = 2*p dimensional boundary matrix (like in iq-binning definition list)
#!# need to adapt to allow bin allocations for observations outside of observed bins
bin_index_finder <- function(x, bin_bounds, nbins, strict=TRUE){ 
  p = length(x)
  b = nrow(bin_bounds)
  if(strict==TRUE) {
    xrep_mat = matrix(rep(x,b),ncol=p,byrow=TRUE)
    idx <- which(rowSums(bin_bounds[,seq(1,2*p-1,by=2)] < xrep_mat & xrep_mat <= bin_bounds[,seq(2,2*p,by=2)])==p)
    if(length(idx)==0L) idx <- NA
  }
  if(strict==FALSE) {
    #!# put process for allocating outside bins
    check_matrix <- matrix(rep(NA,b*p*2),ncol=p*2,byrow=TRUE)
    for (d in 1:p){
      blocks <- prod(nbins[1:d-1])
      blocks_n <- b/blocks
      subblocks <- prod(nbins[1:d])
      subblocks_n <- b/subblocks
      # rows with lowest bin in each strata from last dimension
      cond_lower <- rep(seq(0,b-blocks_n,by=blocks_n),each=subblocks_n) + rep(seq(1,subblocks_n,by=1),blocks)
      # rows with highest bin in each strata from last dimension
      cond_upper <- rep(seq(0,b-blocks_n,by=blocks_n),each=subblocks_n) + rep(seq(blocks_n-subblocks_n+1, blocks_n,by=1),blocks)
      
      above_lb <- bin_bounds[,d*2-1] < as.numeric(x[d]) 
      above_lb[cond_lower] <- TRUE
      check_matrix[,d*2-1] <- above_lb
      
      below_ub <- as.numeric(x[d]) <= bin_bounds[,d*2]
      below_ub[cond_upper] <- TRUE
      check_matrix[,d*2] <- below_ub
    }
    idx <- which(rowSums(check_matrix)==p*2)
    
  }
  return(idx)
} 
# myiq <- iqnn(iris, y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), nbins=c(3,5,2), jit=rep(.001,3))
# new_row <- iris[1,c("Sepal.Length","Sepal.Width","Petal.Width")]
# new_row_index <- bin_index_finder(new_row, myiq$bin_bounds, myiq$nbins)
# new_row
# myiq$bin_bounds[new_row_index,]
# myiq$bin_centers[new_row_index,]
# myiq$bin_stats[new_row_index,]
# 
# bin_index_finder(c(5,7,7), myiq$bin_bounds, myiq$nbins, strict=TRUE)
# bin_index_finder(c(5,7,7), myiq$bin_bounds, myiq$nbins, strict=FALSE)


#--------------------------------------
### function to suggesting best number of bins under contraints
# Constraint 1: nbins_1=nbins_2=...=nbins_p
# Constraint 2: each bin with <= k neighbors (conservative approx to knn)
roots_for_nbins <- function(x, p, k){
  nbin_opt <- x/k
  rep(ceiling(nbin_opt^(1/p)),p)
}
# Test with goal to mimic 10-nn with p=3 dimensions
roots_for_nbins(270,3,10)
270/prod(roots_for_nbins(270,3,10))
roots_for_nbins(122,3,10)
122/prod(roots_for_nbins(122,3,10)) # very rough approx between 81 and 270 due to x^3*10

### function to suggesting best number of bins under more appropriate contraints
# Constraint 1: nbins_1>=nbins_2>=...>=nbins_p
# Constraint 2: each bin with <= k neighbors (conservative approx to knn)
approxknn_nbins <- function(x, p, k){
  #!# work out code for repeatedly adding 1 to some nbins until tips over k per bin
}

#--------------------------------------
### CV cohort additions
# use this function to add K grouping indeces
make_cv_cohorts <- function(dat,cv_K){
  if(nrow(dat) %% cv_K == 0){ # if perfectly divisible
    cv_cohort <- sample(rep(1:cv_K, each=(nrow(dat)%/%cv_K)))
  } else { # if not perfectly divisible
    cv_cohort <- sample(c(rep(1:(nrow(dat) %% cv_K), each=(nrow(dat)%/%cv_K + 1)),
                              rep((nrow(dat) %% cv_K + 1):cv_K,each=(nrow(dat)%/%cv_K)) ) )
  }
  return(cv_cohort)
}

#--------------------------------------
### Helper function for suggesting parameters for jittering number of bins in each dimension
# based on number of ties and data resolution
#!#

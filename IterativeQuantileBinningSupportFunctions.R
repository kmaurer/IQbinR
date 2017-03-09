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

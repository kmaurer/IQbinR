### Helper functions for Iterative Quantile Binning Functions
# Move to this seperate file for organization that focuses on primary tasks

#--------------------------------------
make_bin_list <- function(bin_bounds,nbins){
  bin_dim = length(nbins)
  ### build nested list version of bin_bounds to speed up future searching for bins
  lower_level_list <- list(NULL)
  for(i in 1:nrow(bin_bounds)){
    lower_level_list[[i]] <- i
  } 
  for(d in bin_dim:1){
    # for each dimension from second lowest to highest, group up observations from lower_level_list into items in upper_level_list 
    upper_level_list <- list(NULL)
    upper_blocks <- ifelse(d==1,1,prod(nbins[1:(d-1)]))
    lower_block_size <- nbins[d]
    upper_indeces <- prod(nbins[d:length(nbins)])
    lower_indeces <- ifelse(d==bin_dim,1,prod(nbins[(d+1):bin_dim]))
    
    for(ul in 1:upper_blocks){
      # create upper level groups 
      upper_level_list[[ul]] <- list(NULL)
      for(ll in 1:lower_block_size){
        upper_level_list[[ul]][[ll]] <- lower_level_list[[(ul-1)*lower_block_size+ll]]
      }
      # upper_level_list[[ul]][[lower_block_size+1]] <- bin_bounds[(ul-1)*upper_indeces + 1:lower_block_size*lower_indeces,(d-1)*2+1:2]
      upper_level_list[[ul]][[lower_block_size+1]] <- unique(as.vector(bin_bounds[(ul-1)*upper_indeces + 1:lower_block_size*lower_indeces,(d-1)*2+1:2]))
    }
    lower_level_list <- upper_level_list
  }
  bin_list <- lower_level_list
  return(bin_list)
}
#---------------------------------------

bin_index_finder_nest <- function(x, bin_def, strict=TRUE){ 
  bin_dim = length(bin_def$nbins)
  nest_list <- bin_def$bin_list[[1]]
  x <- as.numeric(x)
    for(d in 1:bin_dim){
      nest_index <- .bincode(x[[d]], nest_list[[bin_def$nbins[d]+1]],T,T)
      if(strict == FALSE){
        if( x[[d]] < min(nest_list[[bin_def$nbins[d]+1]]) ) nest_index <- 1
        if( x[[d]] > max(nest_list[[bin_def$nbins[d]+1]]) ) nest_index <- bin_def$nbins[d]
      }
      # if(length(nest_index)==0) return(print("Observation outside of observed bins, set strict=FALSE "))
      nest_list <- nest_list[[nest_index]]
    }
    idx <- nest_list
  return(idx)
} 
# test_index <- c(1,2,51,52,101,102)
# test_data <- iris[test_index,]
# iq_def <- iterative_quant_bin(data=iris[-test_index,], bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#                               nbins=c(3,2,2), output="both")
# bin_by_iq_def(bin_def=iq_def$bin_def, new_data=test_data, output="data")
# bin_index_finder_nest(x=c(6,3,1.5),bin_def, strict=TRUE)
# bin_index_finder_nest(x=c(6,3,15),bin_def, strict=TRUE)
# bin_index_finder_nest(x=c(6,3,1.5),bin_def, strict=FALSE)
# bin_index_finder_nest(x=c(6,3,15),bin_def, strict=FALSE)

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
# bin_index_finder <- function(x, bin_bounds, nbins, strict=TRUE){ 
#   p = length(x)
#   b = nrow(bin_bounds)
#   if(strict==TRUE) {
#     xrep_mat = matrix(rep(x,b),ncol=p,byrow=TRUE)
#     idx <- which(rowSums(bin_bounds[,seq(1,2*p-1,by=2)] < xrep_mat & xrep_mat <= bin_bounds[,seq(2,2*p,by=2)])==p)
#     if(length(idx)==0L) idx <- NA
#   }else{
#     #!# put process for allocating outside bins
#     check_matrix <- matrix(rep(NA,b*p*2),ncol=p*2,byrow=TRUE)
#     for (d in 1:p){
#       blocks <- prod(nbins[1:d-1])
#       blocks_n <- b/blocks
#       subblocks <- prod(nbins[1:d])
#       subblocks_n <- b/subblocks
#       # rows with lowest bin in each strata from last dimension
#       cond_lower <- rep(seq(0,b-blocks_n,by=blocks_n),each=subblocks_n) + rep(seq(1,subblocks_n,by=1),blocks)
#       # rows with highest bin in each strata from last dimension
#       cond_upper <- rep(seq(0,b-blocks_n,by=blocks_n),each=subblocks_n) + rep(seq(blocks_n-subblocks_n+1, blocks_n,by=1),blocks)
#       
#       above_lb <- bin_bounds[,d*2-1] < as.numeric(x[d]) 
#       above_lb[cond_lower] <- TRUE
#       check_matrix[,d*2-1] <- above_lb
#       
#       below_ub <- as.numeric(x[d]) <= bin_bounds[,d*2]
#       below_ub[cond_upper] <- TRUE
#       check_matrix[,d*2] <- below_ub
#     }
#     idx <- which(rowSums(check_matrix)==p*2)
#     
#   }
#   return(idx)
# } 
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


# #--------------------------------------
# ### function to suggesting best number of bins under contraints
# # Constraint 1: nbins_1=nbins_2=...=nbins_p
# # Constraint 2: each bin with <= k neighbors (conservative approx to knn)
# roots_for_nbins <- function(x, p, k){
#   nbin_opt <- x/k
#   rep(ceiling(nbin_opt^(1/p)),p)
# }
# # Test with goal to mimic 10-nn with p=3 dimensions
# roots_for_nbins(270,3,10)
# 270/prod(roots_for_nbins(270,3,10))
# roots_for_nbins(122,3,10)
# 122/prod(roots_for_nbins(122,3,10)) # very rough approx between 81 and 270 due to x^3*10
# 
# ### function to suggesting best number of bins under more appropriate contraints
# # Constraint 1: nbins_1>=nbins_2>=...>=nbins_p
# # Constraint 2: each bin with <= k neighbors (conservative approx to knn)
# approxknn_nbins <- function(x, p, k){
#   #!# work out code for repeatedly adding 1 to some nbins until tips over k per bin
# }



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

#' Create feature pair data frame
#'
#' @description Based on https://gist.github.com/avsmith/e6f4f654451da139230b to round all numeric variables
#' @param x data frame 
#' @param digits number of digits to round
#' 
round_df <- function(x, digits=2) {
  numeric_columns <- sapply(x, class) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}


#--------------------------------------
#' Simple Majority Vote Counter
#'
#' @description Identify the maximum vote earners, then randomly pick winner if there is a tie to break
#'
#' @param votes character or factor vector
#' votes <- c("a","a","a","b","b","c")
#' majority_vote(votes)
majority_vote <- function(votes){
  top_votes <- names(which.max(table(votes))) # collect top vote earner (ties allowed)
  return(sample(top_votes,1)) # randomly select to break any ties for best
}



#--------------------------------------
#' Function to create list of nbins vectors to put into tuning iqnn 
#'
#' @description create a list of nbins vectors, use progression that increases number of bins in each dimension while always staying balanced between dimensions
#'
#' @param nbin_range positive integer vector containing lower and upper bounds on number of bins in each dimension
#' @param p number of binning dimensions
#' 
#' @return list of nbins vectors
#' @examples 
#' make_nbins_list(c(2,3),3)

make_nbins_list <- function(nbin_range, p){
  nbins_list <- list(rep(nbin_range[1],p))
  counter = 1
  for(i in 1:(nbin_range[2]-nbin_range[1])){
    for(j in 1:p){
      nbins_list[[counter+1]] <- nbins_list[[counter]]
      nbins_list[[counter+1]][j] <- nbins_list[[counter+1]][j] + 1
      counter <- counter+1
    }
  }
  return(nbins_list)
}

#--------------------------------------
### Helper function for suggesting parameters for jittering number of bins in each dimension
# based on number of ties and data resolution
#!#

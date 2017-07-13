### iterative quantile binning functions

# first source in all helper functions
# source("IterativeQuantileBinningSupportFunctions.R")

## Function List
# quant_bin_1d for quantile binning in one dimensioncv_iqnn
# iterative_quant_bin
# bin_by_IQdef
# iqnn
# predict_iqnn
# cv_pred_iqnn
# make_stack_matrix helper funtion for iterative binner to create stacking J matrices

# Note: a #!# tag will be added on items that need improved efficiency/clarity

#--------------------------------------
## Quantile 1d Binning
# used for binning the counts values by quantile
# define vector of counts and number of bin
# xs <- ggplot2::diamonds$price; nbin=4
quant_bin_1d <- function(xs, nbin, output="data",jit=0){
  if(jit > 0)  xs <- xs + runif(length(xs),-jit,jit)
  quants <- quantile(xs, seq(0, 1, by=1/(2*nbin)))
  bin_centers <- quants[seq(2,length(quants)-1, by=2)]
  bin_bounds <- quants[seq(1,length(quants)+1, by=2)]
  if(jit > 0) bin_bounds[c(1,length(bin_bounds))] <- bin_bounds[c(1,length(bin_bounds))]+c(-jit,jit)
  bin_data <- rep(bin_centers[1],length(xs))
  if(output=="definition") {
    return(list(bin_centers=bin_centers,bin_bounds=bin_bounds))
  } else{
    for (i in 2:length(bin_centers)){
      bin_data[bin_bounds[i] < xs] <- bin_centers[i]
    } 
    if(output=="data") return(bin_data)
    if(output=="both") return(list(bin_data=bin_data,bin_centers=bin_centers,bin_bounds=bin_bounds))
  }
}
# quant_bin_1d(ggplot2::diamonds$price,4,output="data")
# quant_bin_1d(ggplot2::diamonds$price,4,output="definition")
# quant_bin_1d(runif(1000,0,10),nbin=4,output="both")


#--------------------------------------
### Iterative Quantile Binning
# Input:
#   data = data frame to be binned (will coerce matrix or tibble to simple data frame)
#   bin_cols = vector of column names of variables to iteratively bin, ordered first to last
#   nbins = vector of number of bins per step of iterative binning, ordered first to last
#   jit = vector of margins for uniform jitter to each dimension to create seperability of tied obs due to finite precision
iterative_quant_bin <- function(data, bin_cols, nbins, output="data",jit = rep(0,length(bin_cols))){
  data <- as.data.frame(data)
  bin_dim <- length(bin_cols)
  bin_data <- matrix(NA,nrow=nrow(data),ncol=bin_dim, dimnames=list(row.names(data),paste(bin_cols,"binned",sep="_")))
  # Initialize with first binning step
  step_bin_info <- quant_bin_1d(data[,bin_cols[1]], nbins[1],output="both",jit[1])
  bin_bounds <- matrix(c(step_bin_info$bin_bounds[1:nbins[1]],
                         step_bin_info$bin_bounds[2:(nbins[1]+1)]),
                       nrow=nbins[1],byrow=FALSE )
  bin_centers <- matrix(step_bin_info$bin_centers, nrow=nbins[1])
  bin_data[,1] <- step_bin_info$bin_data
  # Loop over remaining variables to use quantile binning WITHIN each of previous state bins
  for(d in 2:bin_dim){
    stack_size <- nrow(bin_centers)
    stack_matrix <- make_stack_matrix(stack_size,nbins[d])
    bin_centers <- cbind(stack_matrix %*% bin_centers,matrix(rep(NA,stack_size*nbins[d]),ncol=1))
    bin_bounds <- cbind(stack_matrix %*% bin_bounds,matrix(rep(NA,2*stack_size*nbins[d]),ncol=2))
    # iterate through unique bins from prior step which are the {1,1+nbins[d],1+2*nbins[d],...} rows of the bin matrices
    for(b in seq(1,1+(stack_size-1)*nbins[d],by=nbins[d]) ){
      in_bin_b <- apply(matrix(bin_data[,1:(d-1)],ncol=(d-1)),1,identical,y=bin_centers[b,-d])
      step_bin_info <- quant_bin_1d(data[in_bin_b,bin_cols[d]], nbins[d],output="both",jit[d])
      bin_bounds[b:(b+nbins[d]-1),c(2*d-1,2*d)] <- matrix(c(step_bin_info$bin_bounds[1:nbins[d]],
                                                            step_bin_info$bin_bounds[2:(nbins[d]+1)]),
                                                          nrow=nbins[d],byrow=FALSE)
      bin_centers[b:(b+nbins[d]-1),d] <- matrix(step_bin_info$bin_centers, nrow=nbins[d])
      bin_data[in_bin_b,d] <- step_bin_info$bin_data
    }
  }
  # add bin index column to bin_data
  bin_centers_idx <- as.data.frame(bin_centers)
  names(bin_centers_idx) <- colnames(bin_data)
  bin_centers_idx$bin_index <- 1:nrow(bin_centers_idx)
  bin_data <- merge(round_df(bin_data,6),round_df(bin_centers_idx,6), all.x=TRUE)
  
  
  if(output=="data") return(bin_data)
  if(output=="definition") return(list(bin_centers=bin_centers, bin_bounds=bin_bounds,bin_cols=bin_cols, nbins=nbins, jit=jit))
  if(output=="both"){
    return(list(bin_data=bin_data, 
                bin_def=list(bin_centers=bin_centers, bin_bounds=bin_bounds, bin_cols=bin_cols, nbins=nbins, jit=jit)))
  } 
}
iterative_quant_bin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
                    nbins=c(3,5,2), output="both",jit=rep(0.001,3))

# iterative_quant_bin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#                     nbins=c(3,5,2), output="both")
# 
# iq_def <- iterative_quant_bin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#                     nbins=c(3,5,2), output="definition")



#--------------------------------------
### Adding Tolerance Buffer to outermost bins from Iterative Quantile Binning
# iq_def = iq binning definition list
# tol = vector of tolerance values to stretch each dimension for future binning
stretch_iq_bins <- function(iq_def, tolerance){
  b = nrow(iq_def$bin_def$bin_bounds)
  p = length(iq_def$bin_def$nbins)
  for (d in 1:p){
    blocks <- prod(iq_def$bin_def$nbins[1:d-1])
    blocks_n <- b/blocks
    subblocks <- prod(iq_def$bin_def$nbins[1:d])
    subblocks_n <- b/subblocks
    # strec
    iq_def$bin_def$bin_bounds[seq(1,subblocks_n,by=1),2*d-1] <- iq_def$bin_def$bin_bounds[seq(1,subblocks_n,by=1),2*d-1] - tolerance[d]
    iq_def$bin_def$bin_bounds[seq(blocks_n-subblocks_n+1, blocks_n,by=1),2*d] <- iq_def$bin_def$bin_bounds[seq(blocks_n-subblocks_n+1, blocks_n,by=1),2*d] + tolerance[d]
  }
  return(iq_def)
}
iq_def <- iterative_quant_bin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
                    nbins=c(3,5,2), output="definition")
# iq_def$bin_bounds
# 
# iq_def <- iterative_quant_bin(sim_data[-test_index,], bin_cols=c("x1","x2","x3","x4"),
#                               nbins=rep(4,4), jit=rep(0.001,4), output="both" )
# iq_stretch <- stretch_iq_bins(iq_def, tolerance=rep(5,4))
# iq_stretch$bin_def$bin_bounds


#--------------------------------------
### Iterative Quantile Binning New Data from defined bins
# iq_def= IQ bin definition list from iterative_quant_bin or iqnn
# new_data = data frame with column names matching the binned columns from bin-training data
# output matches format of iterative_quant_bin and inherets properties from iqnn if applicable
bin_by_IQdef <- function(iq_def, new_data, output="data", strict=TRUE){
  new_data <- as.data.frame(new_data)
  #!# need to introduce similar jitter to new data as in definition so "boundary" points allocated randomly
  # loop over each obs in new data, identify the bin indeces then return bin centers for associated bins
  bin_indeces <- sapply(1:nrow(new_data), function(i){
    bin_index_finder(new_data[i,iq_def$bin_cols],iq_def$bin_bounds, iq_def$nbins, strict=strict)
  })
  
  if(output=="data") return(list(data=new_data,bin_data=iq_def$bin_centers[bin_indeces,],bin_indeces=bin_indeces))
  if(output=="both"){
    return(list(bin_data=list(data=new_data,bin_data=iq_def$bin_centers[bin_indeces,],bin_indeces=bin_indeces), 
                bin_def=iq_def))
  } 
} 
# # Testing bin_by_IQdef
# test_index <- c(1,2,51,52,101,102)
# iqnn_mod <- iqnn(iris[-test_index,], y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#                  nbins=c(3,5,2), jit=rep(0.001,3))
# test_data <- iris[test_index,]
# bin_by_IQdef(iq_def=iqnn_mod, new_data=test_data, output="data")
# bin_by_IQdef(iq_def=iqnn_mod, new_data=test_data, output="data", strict=FALSE)


#--------------------------------------
### Iterative Quantile Binned Nearest Neighbors Regression
# takes in data, response column and binning parameters
iqnn <- function(data, y, bin_cols, nbins, jit = rep(0,length(bin_cols)), stretch=FALSE, tolerance = rep(0,length(bin_cols)) ){
  data <- as.data.frame(data)
  ## make bins
  iq_bin <- iterative_quant_bin(data, bin_cols, nbins, output="both",jit)
  # stretch bins if requested
  if(stretch) iq_bin <- stretch_iq_bins(iq_bin, tolerance=tolerance)
  ## For each bin, find indeces from original data where bins match, take average y value
  iq_bin$bin_def$y <- y
  total_bins = nrow(iq_bin$bin_def$bin_centers)
  iq_bin$bin_def$bin_stats <- data.frame(avg = sapply(1:total_bins, function(b) mean(data[iq_bin$bin_data$bin_index==b,y], na.rm=TRUE)),
                                         obs = sapply(1:total_bins, function(b) sum(iq_bin$bin_data$bin_index==b)) )
  ## Return bin definition with predictions added
  return(iq_bin$bin_def)
}
# iqnn(iris, y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), nbins=c(3,5,2))
myiq <- iqnn(data=iris, y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
             nbins=c(3,5,2), jit=rep(0.001,3))
myiq$bin_bounds
# 
# myiq <- iqnn(iris, y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#              nbins=c(3,5,2), jit=rep(0.001,3), stretch=TRUE, tolerance=rep(.1,3))
# myiq$bin_bounds




#--------------------------------------
### predict for new data from iqnn model
predict_iqnn <- function(iqnn_mod,test_data, type="estimate",strict=FALSE){
  test_data <- as.data.frame(test_data)
  test_bin <- bin_by_IQdef(iqnn_mod, test_data, output="data",strict=strict)
  if(type=="estimate") return(iqnn_mod$bin_stats$avg[test_bin$bin_indeces])
  if(type=="binsize") return(iqnn_mod$bin_stats$obs[test_bin$bin_indeces])
  if(type=="both") return(iqnn_mod$bin_stats[test_bin$bin_indeces,])
}
# test_index <- c(1,2,51,52,101,102)
# iqnn_mod <- iqnn(iris[-test_index,], y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#                  nbins=c(3,5,2), jit=rep(0.001,3))
# test_data <- iris[test_index,]
# predict_iqnn(iqnn_mod, test_data,strict=TRUE)
# predict_iqnn(iqnn_mod, test_data,type="binsize")
# predict_iqnn(iqnn_mod, test_data,strict=FALSE)
# predict_iqnn(iqnn_mod, test_data,type="binsize")


#--------------------------------------
### Cross Validated predictions for iqnn models
cv_pred_iqnn <- function(iqnn_mod, data, cv_method="kfold", cv_k=10, strict=FALSE){
  data <- as.data.frame(data)
  if(cv_method=="kfold") cv_cohorts <- make_cv_cohorts(data, cv_k)
  if(cv_method=="LOO") cv_cohorts <- 1:nrow(data)
  cv_preds <- rep(NA,nrow(data))
  for(fold in 1:length(unique(cv_cohorts))){
    test_index <- which(cv_cohorts==fold)
    iqnn_mod <- iqnn(data[-test_index,], y=iqnn_mod$y, bin_cols=iqnn_mod$bin_cols, 
                     nbins=iqnn_mod$nbins, jit=iqnn_mod$jit)
    cv_preds[test_index] <- predict_iqnn(iqnn_mod, data[test_index,],strict=strict)
  }
  cv_preds
}
# iqnn_mod <- iqnn(iris, y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), 
#                  nbins=c(3,5,2), jit=rep(0.001,3))
# cv_pred_iqnn(iqnn_mod,iris, cv_method="kfold", cv_k=10, strict=FALSE)
# cv_pred_iqnn(iqnn_mod,iris, cv_method="LOO", strict=FALSE)


#--------------------------------------
### Cross Validation for assessment for iqnn models
cv_iqnn <- function(iqnn_mod, data, cv_method="kfold", cv_k=10, strict=FALSE){
  data <- as.data.frame(data)
  cv_preds <- cv_pred_iqnn(iqnn_mod, data, cv_method, cv_k, strict)
  PRESS <- sum((data[,iqnn_mod$y]-cv_preds)^2)
  MSE <- PRESS/nrow(data)
  RMSE <- sqrt(MSE)
  c(PRESS=PRESS,MSE=MSE,RMSE=RMSE)
}
# iqnn_mod <- iqnn(iris, y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#                  nbins=c(3,5,2), jit=rep(0.001,3))
# cv_iqnn(iqnn_mod,iris, cv_method="kfold", cv_k=10, strict=FALSE)
# cv_iqnn(iqnn_mod,iris, cv_method="LOO", strict=FALSE)

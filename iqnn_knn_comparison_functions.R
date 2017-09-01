# Functions used in iqnn_knn_comparison.R

kdtree_nn_predict <- function(train,test,k=10){
  nearest <- nn2(data=train,query=test, k=k)
  sapply(1:nrow(test), function(x) {
    mean(train[nearest$nn.idx[x,],1])
  })
}

#--------------------------------------
### Cross Validated predictions for knn model using knn.reg from FNN package
cv_pred_knn <- function(dat, y_name, x_names, cv_method="kfold", cv_k = 10, k=5, knn_algorithm = "brute"){
  dat <- as.data.frame(dat)
  if(cv_method=="kfold") cv_cohorts <- make_cv_cohorts(dat, cv_k)
  if(cv_method=="LOO") cv_cohorts <- 1:nrow(dat)
  cv_preds <- rep(NA,nrow(dat))
  for(fold in 1:length(unique(cv_cohorts))){
    test_index <- which(cv_cohorts==fold)
    knn_mod <- knn.reg(train=dat[-test_index,x_names], test=dat[test_index,x_names], 
                       y=dat[-test_index,y_name], k = k, algorithm = knn_algorithm)
    cv_preds[test_index] <- knn_mod$pred
  }
  cv_preds
}
# cv_preds <- cv_pred_knn(dat=bb_players_st, y_name="hr", x_names=c("b2","b3","hit","ab"),
#            cv_method="kfold", cv_k = 10, k=5, knn_algorithm = "brute")
# head(cv_preds)
# cv_preds <- cv_pred_knn(dat=bb_players_st, y_name="hr", x_names=c("b2","b3","hit","ab"),
#                        cv_method="LOO", k=5, knn_algorithm = "brute")
# head(cv_preds)


#--------------------------------------
### Cross Validated predictions for knn model using knn.reg from FNN package
cv_pred_knn_class <- function(dat, y_name, x_names, cv_method="kfold", cv_k = 10, k=5, knn_algorithm = "brute"){
  dat <- as.data.frame(dat)
  if(cv_method=="kfold") cv_cohorts <- make_cv_cohorts(dat, cv_k)
  if(cv_method=="LOO") cv_cohorts <- 1:nrow(dat)
  cv_preds <- rep(NA,nrow(dat))
  for(fold in 1:length(unique(cv_cohorts))){
    test_index <- which(cv_cohorts==fold)
    knn_mod <- knn(train=dat[-test_index,x_names], test=dat[test_index,x_names], 
                   cl=dat[-test_index,y_name], k = k, algorithm = knn_algorithm)
    cv_preds[test_index] <- as.character(knn_mod)
  }
  cv_preds <- factor(cv_preds, levels=levels(dat[,y_name]))
  cv_preds
}
# cv_preds <- cv_pred_knn_class(dat=iris, y_name="Species", x_names=c("Petal.Length","Sepal.Length"),
#                        cv_method="kfold", cv_k = 10, k=5, knn_algorithm = "brute")
# cv_preds
#--------------------------------------
### Tuning function for knn regression
tune_knn_reg <- function(dat, y_name, x_names, cv_method="kfold", cv_k = 10, k_values=NULL, knn_algorithm = "brute"){
  if(!is.integer(k_values)) return(print("Please specify k_values as an integer vector of neightborhood sizes (k) to be tuned"))
  cv_results <- data.frame(k=k_values,MSE = NA)
  for(k_idx in 1:length(k_values)){
    cv_preds <- cv_pred_knn(dat, y_name, x_names, cv_method="kfold", cv_k = cv_k, k=k_values[k_idx], knn_algorithm = "brute")
    cv_results$MSE[k_idx] <- mean((dat[,y_name]-cv_preds)^2)
  }
  cv_results$RMSE <- sqrt(cv_results$MSE)
  return(cv_results)
}
# timer <- Sys.time()
# tune_knn_results_2 <- tune_knn_reg(dat=bb_players_st, y_name="hr", x_names=c("b2","b3","hit","ab"), k_values=1:50, knn_algorithm = "brute")
# Sys.time()-timer
# 
# timer <- Sys.time()
# tune_knn_results_4 <- tune_knn_reg(dat=bb_players_st, y_name="hr", x_names=c("hit","ab"), k_values=1:50, knn_algorithm = "brute")
# Sys.time()-timer
# 
# head(cv_tune_results)
# 
# ggplot() +
#   geom_point(aes(x=k,y=RMSE),data=cv_tune_results)
# k=25 --> RMSE= 24.25

#--------------------------------------
### Tuning function for knn classifier
tune_knn_class <- function(dat, y_name, x_names, cv_method="kfold", cv_k = 10, k_values=NULL, knn_algorithm = "brute"){
  if(!is.integer(k_values)) return(print("Please specify k_values as an integer vector of neightborhood sizes (k) to be tuned"))
  cv_results <- data.frame(k=k_values,error = NA)
  for(k_idx in 1:length(k_values)){
    cv_preds <- cv_pred_knn_class(dat, y_name, x_names, cv_method="kfold", cv_k = cv_k, k=k_values[k_idx], knn_algorithm = "brute")
    cv_results$error[k_idx] <- sum(cv_preds!=dat[,y_name]) / nrow(dat)
  }
  return(cv_results)
}
# timer <- Sys.time()
# tune_iris <- tune_knn_class(dat=iris, y_name="Species", x_names=c("Petal.Length","Sepal.Length"),
#                                cv_method="kfold", cv_k = nrow(iris), k_values=1:40, knn_algorithm = "brute")
# Sys.time()-timer
# tune_iris
# 
# tune_knn_class

#--------------------------------------
### knn_cv_predict with output as list containing predictions AND timing values for fitting and predicting
knn_cv_pred_timer <- function(data, y, x_names, cv_k = 10, k=5, knn_algorithm = "brute"){
  data <- as.data.frame(data)
  cv_cohorts <- make_cv_cohorts(data, cv_k)
  cv_preds <- factor(rep(NA,nrow(data)),levels=levels(data[,y]))
  time_knn <- 0
  for(fold in 1:length(unique(cv_cohorts))){
    test_index <- which(cv_cohorts==fold)
    timer <- Sys.time()
    knn_mod <- knn(train=data[-test_index,x_names], test=data[test_index,x_names], 
                   cl=data[-test_index,y], k = k, algorithm = knn_algorithm)
    time_knn <- time_knn + difftime(Sys.time(),timer,units="secs")
    cv_preds[test_index] <- knn_mod
  }
  timer <- Sys.time()
  knn_mod <- knn(train=data[-test_index,x_names], test=data[1,x_names], 
                 cl=data[-test_index,y], k = k, algorithm = knn_algorithm) 
  time_fit <- cv_k * difftime(Sys.time(),timer,units="secs")
  
  return(list(preds=cv_preds,time_fit=time_fit, pred_time=time_knn-time_fit))
}

#--------------------------------------
### iqnn_cv_predict with output as list containing predictions AND timing values for fitting and predicting
iqnn_cv_predict_timer <- function(data, y, mod_type = "reg", bin_cols, nbins, 
                                  jit = rep(0,length(bin_cols)), stretch = FALSE, 
                                  tol = rep(0, length(bin_cols)), strict = FALSE, cv_k = 10) {
  data <- as.data.frame(data)
  cv_cohorts <- make_cv_cohorts(data, cv_k)
  cv_preds <- factor(rep("NA", nrow(data)), levels(data[, y]))
  time_fit <- 0
  pred_time <- 0
  for (fold in 1:length(unique(cv_cohorts))) {
    test_index <- which(cv_cohorts == fold)
    train_data_temp <- data[-test_index, ]
    row.names(train_data_temp) <- 1:nrow(train_data_temp)
    timer <- Sys.time()
    iqnn_mod <- iqnn(train_data_temp, y = y, mod_type = mod_type, 
                     bin_cols = bin_cols, nbins = nbins, jit = jit, stretch = stretch, 
                     tol = tol)
    time_fit <- time_fit + difftime(Sys.time(),timer,units="sec")
    
    timer <- Sys.time()
    cv_preds[test_index] <- iqnn_predict(iqnn_mod, data[test_index, 
                                                        ], strict = strict, type = "estimate")
    pred_time <- pred_time + difftime(Sys.time(),timer,units="secs")
  }
  return(list(preds=cv_preds,time_fit=time_fit,pred_time=pred_time))
}

#make function for finding closest binning structure for a given data size, neighborhood size, and dimensionality
find_bin_root <- function(n,k,p){
  # nbins <- rep(1,p)
  # nbins[1] <- floor((n/k)^(1/p))
  # for(i in 2:p){
  #   nbins[i] <- floor((n/k/prod(nbins))^(1/(p-i+1)))
  # }
  # return(nbins[length(nbins):1])
  #--------------------
  # find all possible best ~semi-even nbins vectors
  nbins_list <- iqbin::make_nbins_list(c(floor((n/k)^(1/p)),ceiling((n/k)^(1/p))),p)
  # pick the one that is closest to the k size structure desired
  nbins_list[[which.min(sapply(nbins_list, function(x) abs(k-n/prod(x)) ))]]
}
# find_bin_root(6324,k=3,p=8)
# 6324/prod(find_bin_root(6324,k=3,p=8))

### Comparison of knn and iqnn for regression setting

# Load up functions and packages for iqnn and knn regression
library(devtools)
install_github(repo="kmaurer/iqbin")

library(iqbin)
help(package="iqbin")
?iqnn

library(FNN)
library(tidyverse)
library(stringr)
library(randomForest)
library(RANN)
library(mvtnorm)
library(gridExtra)
library(data.table)

setwd("~/GitHub/iqnnProject")
source("iqnn_knn_comparison_functions.R")

###-----------------------------------------------------------------------------------------------------
# Organize Classification Data sets
setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\classification")

### combine all sizes
all_sets <- c("iris","wpbc","pima","yeast","abalone","waveform","optdigits","satimage","marketing","youtube", "skin")
all_responses <- c("V5","V2","V9","V10","Sex","Class","Class","Class","Sex", "category", "V4")
all_sizes <- c(150,198,768,1484,4174,5000,5620,6435,6876, 168286,245057)

# # Clean all using helper function
# for(set in 1:length(all_sets)){
#   load(file=paste0(all_sets[set],"_raw.Rdata"))
#   # name of response variable
#   y <- all_responses[set]
#   # use helper function to standardize and drop non-numeric/constant-valued input variables
#   data <- clean_data_for_iqnn_knn(as.data.frame(data),y)
#   save(data,file=paste0(all_sets[set],"_cleaned.Rdata"))
# }

###-----------------------------------------------------------------------------------------------------
# Tune neighborhood/bin parameters for each data set using 10-fold cv then save results for running accuracy tests
setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\classification")

max_p <- 2 # max number of dimensions for inputs
cv_k <- 100 # cv folds

tuned_param_list <- list(NULL)
tuned_performance_list <- list(iqnn=list(NULL),knn=list(NULL))
for(set in 1:9){
  print(set)
  load(file=paste0(all_sets[set],"_cleaned.Rdata"))
  y <- all_responses[set]
  data[,y] <- factor(data[,y])
  
  ## Variable selection
  # Find column names in order of importance for randomForest (heuristic for doing variable selection)
  set.seed(12345)
  myforest <- randomForest(as.formula(paste0("as.factor(as.character(",y,"))~ .")) , data=sample_n(data,min(1000,nrow(data))))
  important_cols <- dimnames(importance(myforest))[[1]][order(importance(myforest),decreasing=TRUE)]
  # allow a cap to be put on number of variables considered
  p <- min(length(important_cols),max_p)
  
  ## Parameterize for binning to best match k-nn structure specified with n, k, p, and cv_k
  train_n <- floor(nrow(data)*((cv_k-1)/cv_k))
  bin_cols <- important_cols[1:p]
  
  ## Tune the iqnn shoot for no fewer than 2 per bin (otherwise problems with allocation on boundaries)
  set.seed(1234)
  tune_iqnn_out <- iqnn_tune(data=data, y=y, mod_type = "class", bin_cols=bin_cols, nbins_range=c(2,floor((train_n/2)^(1/p))),
                             jit = rep(0.0001,length(bin_cols)), stretch = FALSE,strict=FALSE,oom_search = ifelse(nrow(data)>100000, TRUE,FALSE), cv_k=cv_k)
  nbins <- tune_iqnn_out$nbins[[which.min(tune_iqnn_out$error)]]
  tuned_performance_list$iqnn[[set]] <- tune_iqnn_out
  
  ## Tune the knn over same range of neighborhood size equivalents
  set.seed(1234)
  tune_knn_out <- tune_knn_class(data=data, y_name=y, x_names=bin_cols, cv_method="kfold", cv_k = cv_k,
                                 k_values=as.integer(round(tune_iqnn_out$nn_equiv)), knn_algorithm = "brute")
  k <- tune_knn_out$k[which.min(tune_knn_out$error)]
  tuned_performance_list$knn[[set]] <- tune_knn_out
 
  tuned_param_list[[set]] <- list(bin_cols=bin_cols, nbins=nbins, k=k, n=nrow(data), cv_k=cv_k)
}
tuned_param_list
tuned_performance_list
# save(tuned_param_list,tuned_performance_list, file="tuned_param_list.Rdata")
# load(file="tuned_param_list.Rdata")

#--------------------------------------------------------------------------------------------------------
# Collecting Accuracy : use average over many 10-fold cv results


nreps <- 100 # number of times to run k-fold comparisons
# # Initialize an empty data structure to put timing/accuracy measurements into
# # Use a list with one df for each method, this is done to allow consistant storage when randomizing order of methods in each trial
results <- data.frame(data_name=rep(all_sets,each=nreps),obs = NA, nn_size = NA, cv_accuracy = NA,
                      time_fit = NA, time_pred = NA, seed = NA)
accuracy_results_list <- list(results_iqnn=results, results_knn=results,
                              results_knn_cover=results, results_knn_kd=results)

setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\classification")
big_timer <- Sys.time()

# Loop over all data sets and repetitions to record accuracies and times. 
for(set in 1:9){
  print(all_sets[set])
  
  load(file=paste0(all_sets[set],"_cleaned.Rdata"))
  y <- all_responses[set]
  data[,y] <- factor(data[,y])
  bin_cols <- tuned_param_list[[set]]$bin_cols
  nbins <- tuned_param_list[[set]]$nbins
  k <-tuned_param_list[[set]]$k
  
  ## Compare knn/iqnn method timing and accuracy with k-fold CV
  # loop over nreps for each method
  for(rep in 1:nreps){
    print(rep)
    for(method in 1:4){
      seed <- rep # set seed as rep number
      # find 10-fold CV predictions, Record time/accuracy for each
      if(method==1){
        set.seed(seed)
        pred_times <-list(preds=iqnn_cv_predict(data=data, y=y, mod_type="class", bin_cols=bin_cols,
                                                nbins=nbins, jit=rep(0.00001,length(nbins)), strict=FALSE,
                                                cv_k=cv_k),time_fit=NA,pred_time=NA)
      } else if(method==2){
        pred_times <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
                                        knn_algorithm = "brute", seed=seed)
      } else if(method==3){
        pred_times <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
                                        knn_algorithm = "cover_tree", seed=seed)
      } else {
        pred_times <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
                                        knn_algorithm = "kd_tree", seed=seed)
      }
      # store results in proper mehtod/set/rep values
      accuracy_results_list[[method]]$obs[(set-1)*nreps + rep] <- nrow(data)
      accuracy_results_list[[method]]$nn_size[(set-1)*nreps + rep] <- ifelse(method==1, nrow(data)/prod(nbins), k)
      accuracy_results_list[[method]]$cv_accuracy[(set-1)*nreps + rep] <- sum(as.character(pred_times$preds)==as.character(data[,y]))/nrow(data)
      accuracy_results_list[[method]]$time_fit[(set-1)*nreps + rep] <- pred_times$time_fit
      accuracy_results_list[[method]]$time_pred[(set-1)*nreps + rep] <- pred_times$pred_time
      accuracy_results_list[[method]]$seed[(set-1)*nreps + rep] <- seed
    }
  }
}
Sys.time() - big_timer
str(accuracy_results_list)

results_class <- data.frame(do.call("rbind", accuracy_results_list),
                          type=rep(c("iqnn","knn - brute","knn - cover tree","knn - kd tree"),each=nrow(accuracy_results_list$results_iqnn))) %>%
  group_by(data_name,obs,nn_size,type) %>%
  summarize(avg_cv_accuracy=mean(cv_accuracy,na.rm=TRUE)) %>%
  as.data.frame() %>%
  arrange(data_name, type) %>%
  group_by(data_name) %>%
  na.omit() %>%
  mutate(diff_acc = avg_cv_accuracy - avg_cv_accuracy[2]) %>% #!# only works due to knn brute being 2nd in order - danger of hard coding (don't know simple alternative)
  ungroup() %>%
  mutate(data_name = factor(data_name, levels=all_sets[order(all_sizes)]),
         diff_acc_perc = diff_acc*100) %>%
  as.data.frame()

head(results_class)
str(results_class)


# save(accuracy_results_list, results_class,tuned_param_list,tuned_performance_list, file="classification_testing.Rdata")

ggplot()+
  geom_line(aes(x=data_name,y=diff_acc,color=type,group=type),size=1,data=results_class)+
  theme_bw()


# 
# 
# 
# setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\classification")
# big_timer <- Sys.time()
# # Loop over all data sets and repetitions to record accuracies and times. 
# for(set in 1:8){
#   load(file=paste0(all_sets[set],"_cleaned.Rdata"))
#   y <- all_responses[set]
#   data[,y] <- factor(data[,y])
#   k <- ceiling(nrow(data)/81)
# 
#   ## Variable selection
#   # Find column names in order of importance for randomForest (heuristic for doing variable selection)
#   set.seed(12345)
#   myforest <- randomForest(as.formula(paste0("as.factor(as.character(",y,"))~ .")) , data=sample_n(data,min(1000,nrow(data))))
#   important_cols <- dimnames(importance(myforest))[[1]][order(importance(myforest),decreasing=TRUE)]
#   # allow a cap to be put on number of variables considered
#   p <- min(length(important_cols),max_p)
#   
#   ## Parameterize for binning to best match k-nn structure specified with n, k, p, and cv_k
#   train_n <- nrow(data)*((cv_k-1)/cv_k)
#   nbins <- find_bin_root(n=train_n,k=k,p=p)
#   bin_cols <- important_cols[1:p]
#   
#   print(paste(all_sets[set],"with",ncol(data),"columns and k=",k))
#   
#   ## Compare knn/iqnn method timing and accuracy with k-fold CV
#   # loop over nreps for each method
#   for(rep in 1:nreps){
#     # set seed for CV partitioning so that each method uses same train/test splits
#     seed <- sample(1:100000,1)
#     # pick order for methods at random
#     method_order <- sample(1:4)
#     # seed <-  12345 # fixed value to check if all reps identical predictions made **Confirmed as identical for accuracy**
#     for(method in method_order){
#       # find 10-fold CV predictions, Record time/accuracy for each
#       if(method==1){
#         # pred_times <- iqnn_cv_predict_timer(data=data, y=y, mod_type="class", bin_cols=bin_cols,
#         #                                     nbins=nbins, jit=rep(0.000001,length(nbins)), strict=FALSE, 
#         #                                     cv_k=cv_k, seed=seed)
#         set.seed(seed)
#         pred_times <-list(preds=iqnn_cv_predict(data=data, y=y, mod_type="class", bin_cols=bin_cols,
#                                             nbins=nbins, jit=rep(0.000001,length(nbins)), strict=FALSE,
#                                             cv_k=cv_k),time_fit=NA,pred_time=NA)
#       } else if(method==2){
#         pred_times <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
#                                         knn_algorithm = "brute", seed=seed)
#       } else if(method==3){
#         pred_times <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
#                                         knn_algorithm = "cover_tree", seed=seed)
#       } else {
#         # pred_times <- kdtree_nn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
#         #                                       eps=1, seed=seed)
#         pred_times <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
#                                         knn_algorithm = "kd_tree", seed=seed)
#       }
#       # store results in proper mehtod/set/rep values
#       accuracy_results_list[[method]]$obs[(set-1)*nreps + rep] <- nrow(data)
#       accuracy_results_list[[method]]$nn_size[(set-1)*nreps + rep] <- k
#       accuracy_results_list[[method]]$cv_accuracy[(set-1)*nreps + rep] <- sum(as.character(pred_times$preds)==as.character(data[,y]))/nrow(data)
#       accuracy_results_list[[method]]$time_fit[(set-1)*nreps + rep] <- pred_times$time_fit
#       accuracy_results_list[[method]]$time_pred[(set-1)*nreps + rep] <- pred_times$pred_time
#       accuracy_results_list[[method]]$seed[(set-1)*nreps + rep] <- seed
#     }
#   }
# }
# Sys.time() - big_timer
# str(accuracy_results_list)


###--------------------------------------------------------------------------------------------
# Accuracy Comparisons for Regression

setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\regression")

# all_reg_sets <- c("air_quality","casp","ccpp","laser","puma","quake","skillcraft","treasury","wankara","wpbc")
# # Clean all using helper function
# for(set in 1:length(all_reg_sets)){
#   load(file=paste0(all_reg_sets[set],"_raw.Rdata"))
#   # name of response variable
#   data$y <- scale(data$y)
#   # use helper function to standardize and drop non-numeric/constant-valued input variables
#   data <- clean_data_for_iqnn_knn(as.data.frame(data),y)
#   save(data,file=paste0(all_reg_sets[set],"_cleaned.Rdata"))
# }

# Tune neighborhood/bin parameters for each data set using 10-fold cv then save results for running accuracy tests
setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\regression")

max_p <- 2 # max number of dimensions for inputs
cv_k <- 100 # cv folds

all_reg_sets <- c("air_quality","casp","ccpp","laser","puma","quake","skillcraft","treasury","wankara","wpbc")

tuned_reg_param_list <- list(NULL)
tuned_reg_performance_list <- list(iqnn=list(NULL),knn=list(NULL))
for(set in 1:length(all_reg_sets)){
  print(set)
  load(file=paste0(all_reg_sets[set],"_cleaned.Rdata"))
  head(data)
  
  ## Variable selection
  # Find column names in order of importance for randomForest (heuristic for doing variable selection)
  set.seed(12345)
  myforest <- randomForest(y~. , data=sample_n(data,min(1000,nrow(data))))
  important_cols <- dimnames(importance(myforest))[[1]][order(importance(myforest),decreasing=TRUE)]
  # allow a cap to be put on number of variables considered
  # p <- min(length(important_cols),max_p)
  p=2
  
  ## Parameterize for binning to best match k-nn structure specified with n, k, p, and cv_k
  train_n <- floor(nrow(data)*((cv_k-1)/cv_k))
  bin_cols <- important_cols[1:p]
  
  ## Tune the iqnn shoot for no fewer than 2 per bin (otherwise problems with allocation on boundaries)
  set.seed(1234)
  tune_iqnn_out <- iqnn_tune(data=data, y="y", mod_type = "reg", bin_cols=bin_cols, nbins_range=c(2,floor((train_n/2)^(1/p))),
                             jit = rep(0.0001,length(bin_cols)), stretch = FALSE,strict=FALSE,oom_search = ifelse(nrow(data)>100000, TRUE,FALSE), cv_k=cv_k)
  nbins <- tune_iqnn_out$nbins[[which.min(tune_iqnn_out$MSE)]]
  tuned_reg_performance_list$iqnn[[set]] <- tune_iqnn_out
  
  ## Tune the knn over same range of neighborhood size equivalents
  set.seed(1234)
  tune_knn_out <- tune_knn_reg(dat=data, y_name="y", x_names=bin_cols, cv_method="kfold", cv_k = cv_k,
                                 k_values=as.integer(round(tune_iqnn_out$nn_equiv)), knn_algorithm = "brute")
  k <- tune_knn_out$k[which.min(tune_knn_out$MSE)]
  tuned_reg_performance_list$knn[[set]] <- tune_knn_out
  
  tuned_reg_param_list[[set]] <- list(bin_cols=bin_cols, nbins=nbins, k=k, n=nrow(data), cv_k=cv_k)
}
tuned_reg_param_list
tuned_reg_performance_list
# save(tuned_reg_param_list,tuned_reg_performance_list, file="tuned_reg_param_list.Rdata")
# load(file="tuned_param_list.Rdata")



#---------------------------------------------------------------------------------------------------

# write.csv(results_all,"resultsToShareKarsten.csv", row.names=FALSE)

results_all_3nn  <- results_all %>%
  select(-avg_cv_accuracy) %>%
  gather(key="metric",value="value",avg_time_fit:diff_acc) %>% 
  group_by(data_name,obs,nn_size,type,metric) %>%
  summarize(value=mean(value,na.rm=TRUE)) %>%
  as.data.frame() %>%
  mutate(data_name = factor(data_name, levels=medium_sets[order(sizes)]),
         metric_pretty = factor(metric, labels=c("Preprocess Time (sec)","Prediction Time (sec)","Test Accuracy (% diff from KNN)")),
         metric_pretty = factor(metric, labels=c("Preprocess Time (sec)","Prediction Time (sec)","Test Accuracy (% diff from KNN)"))) 

# Combine into data frame for plots
results_all <- data.frame(do.call("rbind", accuracy_results_list),
                          type=rep(c("iqnn","knn - brute","knn - cover tree","knn - kd tree"),each=nrow(accuracy_results_list$results_iqnn))) %>%
  gather(key="metric",value="value",cv_accuracy:time_pred) %>%
  group_by(data_name,obs,nn_size,type,metric) %>%
  summarize(value=mean(value,na.rm=TRUE)) %>%
  as.data.frame() %>%
  mutate(data_name = factor(data_name, levels=medium_sets[order(sizes)]),
         metric_pretty = factor(metric, labels=c("Test Accuracy Rate","Preprocess Time (sec)","Prediction Time (sec)")),
         metric_pretty = factor(metric, labels=c("Test Accuracy Rate","Preprocess Time (sec)","Prediction Time (sec)")))
head(results_all)
levels(results_all$metric_pretty)

label_data <- arrange(unique(results_all[,c("metric_pretty","data_name")]),data_name)[1:18,]
label_data$n <- NA
label_data[label_data$metric_pretty=="Prediction Time (sec)",]$n <- paste0("n=",sort(sizes))
label_data

# save(results_all_3nn,label_data,file="plot_data_3nn.Rdata")


# plot the accuracy/fit time/prediction time
ggplot()+
  geom_hline(yintercept = 0)+
  geom_line(aes(x=data_name, y=value,color=type, group=type),size=1, data=results_all_3nn)+
  facet_grid(metric_pretty ~ ., scales="free_y") +
  # geom_text(aes(x=data_name, label=n), y=.2,data=label_data)+
  theme_bw()+
  labs(title="3-NN Classifier vs IQNN Classifier (~3 per bin)",
       subtitle="Based on averages across 1000 repetitions of 10-fold CV",
       x="Data Set", y="", 
       caption="Data Source: UCI Machine Learning Data Repository")


#----------------------------------------------------------------------------------------------------
library(dplyr)
library(RWeka)
WPM("load-packages")   

## Set seed
seed <- 100

## Create a place to hold results from instance selection methods
results <- data.frame(Dataset = numeric(0), TSSMethod = numeric(0),
                      Fold = numeric(0), Size = numeric(0),
                      TrainAccuracy = numeric(0), TestAccuracy = numeric(0), 
                      ReductionTime = numeric(0), PredictionTime = numeric(0))

## Point to location of datasets
datasets <- list.files("C:\\Users\\maurerkt\\Google Drive\\AFRLSFFP\\Fall2017\\mediumDatasets\\in",
                       full.names = TRUE)

## Create filters
folder <- make_Weka_filter("weka.filters.supervised.instance.StratifiedRemoveFolds")
greedy <- make_Weka_filter("weka.filters.GreedyThreaded_SuperSpecialForKarsten")
drop3 <- make_Weka_filter("weka.filters.Drop3_Wrapper")

## Create a function to assess accuracy
RA<-function(confusion_matrix){
  row_dim<-dim(confusion_matrix)[1]
  s1<-1
  diag_sum<-0
  accuracy<-0
  while(s1<=row_dim)
  {
    s2<-1
    while(s2<=row_dim)
    {
      if(s1==s2)
      {
        diag_sum<-diag_sum+confusion_matrix[s1,s2]
      }
      s2<-s2+1
    }
    s1<-s1+1
  }
  accuracy<-diag_sum/sum(confusion_matrix)
  return(accuracy)
}

## Create a nearest neighbors classifier
knn <- RWeka::make_Weka_classifier("weka/classifiers/lazy/IBk")

## Cycle through 10-fold and fill in results
for(i in datasets){
  dat <- read.arff(i)
  name1 <- unlist(stringr::str_split(i, pattern = "/"))
  name <- name1[length(name1)]
  colnames(dat)[dim(dat)[2]] <- "Class"
  
  for(j in c(1:10)){
    train <- folder(Class ~ ., data = dat, control = Weka_control(V = TRUE, N = 10, F = j, S = seed)) 
    test <- folder(Class ~ ., data = dat, control = Weka_control(V = FALSE, N = 10, F = j, S = seed)) 
    
    ## Original (no filter)
    start <- Sys.time()
    classifier <- knn(Class ~., train, control = Weka_control(K = 3))
    timeToPrep <- Sys.time() - start

    trainPred <- predict(classifier, train[,-dim(train)[2]])
    trainAcc <- RA(table(trainPred,train$Class))
    
    start <- Sys.time()
    testPred <- predict(classifier, test[,-dim(test)[2]])
    testAcc <- RA(table(testPred,test$Class))
    timeToClassify <- Sys.time() - start
    toAdd <- data.frame(Dataset = name, TSSMethod = "None",
                        Fold = j, Size = dim(train)[1],
                        TrainAccuracy = trainAcc, TestAccuracy = testAcc, 
                        ReductionTime = timeToPrep, PredictionTime = timeToClassify)
    results <- rbind(results, toAdd)
    
    ## Greedy
    start <- Sys.time()
    selected <- greedy(Class ~ ., dat = train)
    classifier <- knn(Class ~., dat = selected, control = Weka_control(K = 3))
    timeToFilter <- Sys.time() - start
    
    trainPred <- predict(classifier, train[,-dim(train)[2]])
    trainAcc <- RA(table(trainPred,train$Class))
    
    start <- Sys.time()
    testPred <- predict(classifier, test[,-dim(test)[2]])
    testAcc <- RA(table(testPred,test$Class))
    timeToClassify <- Sys.time() - start
    toAdd <- data.frame(Dataset = name, TSSMethod = "Greedy",
                        Fold = j, Size = dim(selected)[1],
                        TrainAccuracy = trainAcc, TestAccuracy = testAcc, 
                        ReductionTime = timeToFilter, PredictionTime = timeToClassify)
    results <- rbind(results, toAdd)
    
    ## DROP3
    start <- Sys.time()
    selected <- drop3(Class ~ ., dat = train)
    classifier <- knn(Class ~., dat = selected, control = Weka_control(K = 3))
    timeToFilter <- Sys.time() - start
    
    trainPred <- predict(classifier, train[,-dim(train)[2]])
    trainAcc <- RA(table(trainPred,train$Class))
    
    start <- Sys.time()
    testPred <- predict(classifier, test[,-dim(test)[2]])
    testAcc <- RA(table(testPred,test$Class))
    timeToClassify <- Sys.time() - start
    toAdd <- data.frame(Dataset = name, TSSMethod = "DROP3",
                        Fold = j, Size = dim(selected)[1],
                        TrainAccuracy = trainAcc, TestAccuracy = testAcc, 
                        ReductionTime = timeToFilter, PredictionTime = timeToClassify)
    results <- rbind(results, toAdd)
  } 
}

# process results into similar form to those from knn/iqnn
results2 <- results %>% group_by(Dataset, TSSMethod) %>%
  summarize(Size = mean(Size),
            TrainAccuracy = mean(TrainAccuracy),
            TestAccuracy = mean(TestAccuracy),
            TimeReduce = mean(ReductionTime),
            TimePredict = mean(PredictionTime)) %>%
  ungroup() %>% as.data.frame()

results2 %>% knitr::kable()
# write.csv(results2, "resultsToShareWalter.csv", row.names=FALSE)

walter_data <- read.csv("resultsToShare.csv")
head(walter_data)
results2 <- read.csv("resultsToShareWalter.csv")
head(results2)

library(stringr)
results_instance_selection <- results2 %>%
  select(Dataset,TSSMethod,TestAccuracy,TimeReduce,TimePredict) %>%
  gather(key="metric",value="value",TestAccuracy:TimePredict) %>%
  mutate(data_name =  str_sub(Dataset,9,-6),
         metric_pretty = factor(metric, labels=c("Test Accuracy Rate","Prediction Time (sec)","Preprocess Time (sec)")),
         metric_pretty = factor(metric_pretty, levels=c("Test Accuracy Rate","Preprocess Time (sec)","Prediction Time (sec)")))
  
head(results_instance_selection)
levels(results_instance_selection$metric_pretty)

ggplot()+
  geom_hline(yintercept = 0)+
  geom_line(aes(x=data_name, y=value,color=type, group=type),size=1, data=results_all)+
  geom_line(aes(x=data_name, y=value,color=TSSMethod, group=TSSMethod),size=1, data=results_instance_selection)+
  facet_grid(metric_pretty ~ ., scales="free_y")+
  theme_bw()+
  labs(title="3-NN, Instance Selection and IQNN Classifier (~3 per bin)",
       # subtitle="Based on averages across 1000 repetitions of 10-fold CV",
       x="Data Set", y="", 
       caption="Data Source: UCI Machine Learning Data Repository")+
  geom_text(aes(x=data_name, label=n), y=.2,data=label_data)






# ------------------------------------------------------------------------------------------------------
# Baseball batting data from sean lahmann's database 
# - http://www.seanlahman.com/baseball-archive/statistics/
baseball <- read.csv("http://kmaurer.github.io/documents/SLahman_Batting2014.csv")
head(baseball)

baseball <- na.omit(baseball %>%
                      select(playerID:HR))

bb_players <- baseball %>%
  select(playerID:HR, -lgID) %>%
  mutate(hit_rate = H/G) %>%
  arrange(playerID, yearID) %>%
  group_by(playerID) %>%
  summarise(hr = sum(HR,na.rm=TRUE),
            b2 = sum(X2B,na.rm=TRUE),
            b3 = sum(X3B,na.rm=TRUE),
            hit = sum(H,na.rm=TRUE),
            ab = sum(AB,na.rm=TRUE))
bb_players <- as.data.frame(na.omit(bb_players))
head(bb_players)

# need standardized variables in knn, add to time taken for computation
bb_players_st <- bb_players %>%
  mutate(b2 = scale(b2),
         b3 = scale(b3),
         hit = scale(hit),
         ab = scale(ab))
head(bb_players_st)

## Check that we can fit models to batting career data
# iqdef <- iterative_quant_bin(dat=bb_players, bin_cols=c("b2","b3","hit","ab"),
#                     nbins=c(2,2,2,2), jit=rep(0.001,4), output="both")
# 
# iqnn_mod <- iqnn(dat=bb_players, y="hr", bin_cols=c("b2","b3","hit","ab"),
#                  nbins=c(2,2,2,2), jit=rep(0.001,4))
# cv_iqnn(iqnn_mod,bb_players, cv_method="kfold", cv_k=5, strict=FALSE)
# cv_iqnn(iqnn_mod,bb_players, cv_method="LOO", strict=FALSE)

#### knn.reg
test_index <- 1:100
knnTest <- knn.reg(train = bb_players_st[-test_index,c("b2","b3","hit","ab")],
                   test = bb_players_st[test_index,c("b2","b3","hit","ab")],
                   y = bb_players_st$hr[-test_index], k = 5, algorithm = "brute")
knnTest$pred

# from building model to predicting for new

# Comparing times with baseball data
test_index <- 1:8820
timer <- Sys.time()
knnTest <- knn.reg(train = bb_players_st[-test_index,c("b2","b3","hit","ab")],
                   test = bb_players_st[test_index,c("b2","b3","hit","ab")],
                   y = bb_players_st$hr[-test_index], k = 5, algorithm = "brute")
Sys.time() - timer

timer <- Sys.time()
iqnn_mod <- iqnn(bb_players_st[-test_index,], y="hr", bin_cols=c("b2","b3","hit","ab"),
                 nbins=c(7,7,6,6), jit=rep(0.00001,4), tol=rep(0.0001,4))
iqnn_preds <- predict_iqnn(iqnn_mod, bb_players_st[test_index,],strict=FALSE)
Sys.time() - timer


test_index <- 1:45000
timer <- Sys.time()
knnTest <- knn.reg(train = baseball[-test_index,c("X2B","H","AB")],
                   test = baseball[test_index,c("X2B","H","AB")],
                   y = baseball$HR[-test_index], k = 50, algorithm = "brute")
Sys.time() - timer

timer <- Sys.time()
iqnn_mod <- iqnn(baseball[-test_index,], y="HR", bin_cols=c("X2B","H","AB"),
                 nbins=c(9,9,9), jit=rep(0.001,3))
iqnn_preds <- predict_iqnn(iqnn_mod, baseball[test_index,],strict=FALSE)
Sys.time() - timer

#-----------------------------------------------------------------------------------------------
### Testing with cabs data
setwd("C:\\Users\\maurerkt\\Documents\\Data")
load("onePercentSample.Rdata")
library(tidyverse)
library(lubridate)

sample_size <- 200000

set.seed(12345)
taxi <- onePercentSample %>%
  select(payment_type,pickup_datetime,passenger_count,trip_distance,pickup_longitude,pickup_latitude,fare_amount,tip_amount) %>%
  na.omit()  %>%
  filter(payment_type %in% c("credit","cash")) %>%
  sample_n(sample_size) %>% 
  mutate(time = 60*60*hour(pickup_datetime) + 60*minute(pickup_datetime) + second(pickup_datetime),
         wday = wday(pickup_datetime),
         payment_type = factor(payment_type))
names(taxi)[1] <- "true_class"
head(taxi)

taxi_std <- taxi %>%
  mutate(time = scale(time),
         pickup_longitude = scale(pickup_longitude),
         pickup_latitude = scale(pickup_latitude))
head(taxi_std)

set.seed(12345)
test_index <- sample(1:sample_size,(sample_size/2))

# Compare Regression
#--------------
# unstandardized
timer <- Sys.time()
knnTest <- knn.reg(train = taxi[-test_index,c("time","pickup_longitude","pickup_latitude")],
                   test = taxi[test_index,c("time","pickup_longitude","pickup_latitude")],
                   y = taxi$fare_amount[-test_index], k = 100, algorithm = "brute")
Sys.time() - timer
sqrt(mean((taxi[test_index,"fare_amount"]-knnTest$pred)^2))

#standardized
timer <- Sys.time()
taxi_std <- taxi %>%
  mutate(time = scale(time),
         pickup_longitude = scale(pickup_longitude),
         pickup_latitude = scale(pickup_latitude))
knnTest <- knn.reg(train = taxi_std[-test_index,c("time","pickup_longitude","pickup_latitude")],
                   test = taxi_std[test_index,c("time","pickup_longitude","pickup_latitude")],
                   y = taxi_std$fare_amount[-test_index], k = 100, algorithm = "brute")
Sys.time() - timer
sqrt(mean((taxi_std[test_index,"fare_amount"]-knnTest$pred)^2))

#--------------
timer <- Sys.time()
iqbin(data=taxi[-test_index,], bin_cols=c("time","pickup_longitude","pickup_latitude"),
      nbins=c(10,10,10), jit=rep(0.001,3), output="both")
Sys.time() - timer

timer <- Sys.time()
iqnn_mod <- iqnn(taxi[-test_index,], y="fare_amount", bin_cols=c("time","pickup_longitude","pickup_latitude"),
                 nbins=c(10,10,10), jit=rep(0.001,3))
Sys.time() - timer

timer <- Sys.time()
iqnn_preds <- predict_iqnn(iqnn_mod, taxi[test_index,],strict=FALSE)
Sys.time() - timer
round(mean(iqnn_mod$bin_stats$obs)) #approx number of neightbors?
sqrt(mean((taxi[test_index,"fare_amount"]-iqnn_preds)^2))
#--------------



# Compare classification
#--------------
test_index <- sample(1:sample_size,(sample_size/2))
timer <- Sys.time()
knnTest <- knn(train = taxi[-test_index,c("time","pickup_longitude","pickup_latitude")],
                   test = taxi[test_index,c("time","pickup_longitude","pickup_latitude")],
                   cl = taxi$true_class[-test_index], k = 12, algorithm = "brute")
Sys.time() - timer
head(knnTest)
table(knnTest,taxi$true_class[test_index])
1-sum(diag(table(knnTest,taxi$true_class[test_index])))/length(test_index)
#--------------
timer <- Sys.time()
iqnn_mod <- iqnn(taxi[-test_index,], y="true_class",mod_type="class", bin_cols=c("time","pickup_longitude","pickup_latitude"),
                 nbins=c(20,20,20), jit=rep(0.001,3))
Sys.time() - timer
timer <- Sys.time()
iqnn_preds <- predict_iqnn(iqnn_mod, taxi[test_index,],strict=FALSE)
Sys.time() - timer
round(mean(iqnn_mod$bin_stats$obs)) #approx number of neightbors?
table(iqnn_preds,taxi$true_class[test_index])
1-sum(diag(table(iqnn_preds,taxi$true_class[test_index])))/length(test_index)
#--------------

# ------------------------------------------------------------------------------------------
# Cover type classification 
load(file="C:/Users/maurerkt/Documents/GitHub/BinStackedEnsemble/Data/cover_type.Rdata")

head(cover_type)

cover_type_std <- cover_type %>%
  mutate(elevation=scale(elevation),
         hori_dist_road=scale(hori_dist_road),
         hori_dist_fire=scale(hori_dist_fire))
# Compare classification


#--------------
sample_size <- nrow(cover_type)
test_index <- sample(1:sample_size,(sample_size/2))
timer <- Sys.time()
knnTest <- knn(train = cover_type_std[-test_index,c("elevation","hori_dist_road","hori_dist_fire")],
               test = cover_type_std[test_index,c("elevation","hori_dist_road","hori_dist_fire")],
               cl = cover_type_std$true_class[-test_index], k = 100, algorithm = "brute")
Sys.time() - timer
levels(knnTest)
knnTest <- factor(knnTest, levels=levels(cover_type_std$true_class))
table(knnTest,cover_type_std$true_class[test_index])
1-sum(diag(table(knnTest,cover_type_std$true_class[test_index])))/length(test_index)
# 6.9 minutes, 20.8% error rate

#--------------
timer <- Sys.time()
iqnn_mod <- iqnn(cover_type_std[-test_index,], y="true_class",mod_type="class", bin_cols=c("hori_dist_fire","elevation","hori_dist_road"),
                 nbins=c(14,14,14), jit=rep(0.001,3))
Sys.time() - timer

timer <- Sys.time()
iqnn_preds <- predict_iqnn(iqnn_mod, cover_type_std[test_index,],strict=FALSE)
Sys.time() - timer
round(mean(iqnn_mod$bin_stats$obs)) #approx number of neightbors?
iqnn_preds <- factor(iqnn_preds, levels=levels(cover_type_std$true_class))
table(iqnn_preds,cover_type_std$true_class[test_index])
1-sum(diag(table(iqnn_preds,cover_type_std$true_class[test_index])))/length(test_index)
# 2.5 minutes build, .35 minutes predict, 26.66% error rate


# ------------------------------------------------------------------------------------------
# abalone classification with cross validation
library(data.table)
abalone <- fread('https://archive.ics.uci.edu/ml/machine-learning-databases/abalone/abalone.data')
abalone<-as.data.frame(abalone)
names(abalone)[1]<-"true_class"
abalone$true_class<-as.factor(abalone$true_class)
# standardize
abalone <- data.frame(true_class=abalone$true_class,sapply(abalone[,-1], scale))
str(abalone)
timer <- Sys.time()
iqnn_preds <- cv_pred_iqnn(data=abalone, y="true_class",mod_type="class", bin_cols=c("V5","V8","V2"),
                         nbins=c(5,5,5), jit=rep(0.0000001,3),
                         strict=FALSE, cv_method="kfold", cv_k=100)
iqnn_preds <- factor(iqnn_preds, levels=levels(abalone$true_class))
table(iqnn_preds,abalone$true_class)
1-sum(diag(table(iqnn_preds,abalone$true_class)))/nrow(abalone)
Sys.time()-timer

timer <- Sys.time()
iqnn_preds <- cv_pred_knn_class(dat=abalone, y_name="true_class", x_names=c("V5","V8","V2"),
                       cv_method="kfold", cv_k=100, k=33, knn_algorithm = "brute")
iqnn_preds <- factor(iqnn_preds, levels=levels(abalone$true_class))
table(iqnn_preds,abalone$true_class)
1-sum(diag(table(iqnn_preds,abalone$true_class)))/nrow(abalone)
Sys.time()-timer



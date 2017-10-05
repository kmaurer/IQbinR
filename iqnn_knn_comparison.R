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

setwd("~/GitHub/iqnnProject")
source("iqnn_knn_comparison_functions.R")

#---------------------------------------------------------------------------

# simulate data from different numbers of dimensions, bins per dimension and neighborhood size
ps = 2:5 # number of dimensions
deltas = 2:5 # number of bins per dimension
ks = c(1,10,100,1000) # number in iq-neighborhood
ms = c(10000)
combinations <- expand.grid(ps,deltas,ks,ms)
names(combinations) <- c("p","d","k","m")
sim_times <- data.frame(combinations, n = NA,
                        knn_brute_fittime=NA,knn_brute_predtime=NA,
                        knn_cover_fittime=NA,knn_cover_predtime=NA,
                        knn_kd_fittime=NA,knn_kd_predtime=NA,
                        iqnn_fittime=NA, iqnn_predtime=NA)

# Alternative parameterization based on 2^x scaling
ns= 2^(seq(0,16,by=4))
ks= 2^(seq(0,8,by=4))
ps= 2^(1:2)
ms = c(1000)
combinations <- expand.grid(ns,ps,ks,ms)
names(combinations) <- c("n","p","k","m")
combinations <- combinations %>%
  filter(n > k) %>%
  mutate(d=(n/k)^(1/p))
head(combinations,10)
sim_times <- data.frame(combinations,
                        knn_brute_fittime=NA,knn_brute_predtime=NA,
                        knn_cover_fittime=NA,knn_cover_predtime=NA,
                        knn_kd_fittime=NA,knn_kd_predtime=NA,
                        iqnn_fittime=NA, iqnn_predtime=NA)
head(sim_times)

sim_all <- NULL
for(ntrials in 1:10){
  for(sim in 1:nrow(sim_times)){
    p = sim_times$p[sim]
    d = sim_times$d[sim]
    k = sim_times$k[sim]
    m = sim_times$m[sim]
    n = sim_times$n[sim]
    # sim data with number of observations to align number of bins with knn sizes
    # n <- k*d^p 
    sim_times$n[sim] <- n
    set.seed(1234)
    sim_data <- data.frame(sapply(1:p, function(x) runif(n+m)),
                           y=runif(n+m))
    # rebuild column names to proper dimension 
    xcols <- names(sim_data)[1:p]
    
    test_index <- 1:m
    #-------
    # time the knn predictions with brute force
    timer <- Sys.time()
    #!# need to add time taken for standardization?
    knnTest <- knn.reg(train = sim_data[-test_index,xcols],
                       test = sim_data[test_index,xcols],
                       y = sim_data$y[-test_index], k = k, algorithm = "brute")
    sim_times$knn_brute_predtime[sim] <- as.numeric(Sys.time() - timer,units="mins")
    #-------
    # time the knn predictions with cover tree
    knnTest2 <- aknn_predict(train_x = sim_data[-test_index,xcols],
                             test_x = sim_data[test_index,xcols],
                             y = sim_data[-test_index,"y"],
                             mod_type="reg", k=k,algorithm="cover_tree")
    sim_times$knn_cover_fittime[sim] <- knnTest2$fittime
    sim_times$knn_cover_predtime[sim] <- knnTest2$predtime
    #-------  
    # time the knn predictions with kd_tree
    knnTest3 <- aknn_predict(train_x = sim_data[-test_index,xcols],
                             test_x = sim_data[test_index,xcols],
                             y = sim_data[-test_index,"y"],
                             mod_type="reg", k=k,algorithm="kd_tree")
    # knnTest3 <- kdtree_nn_predict(train = sim_data[-test_index,xcols],
    #                               test = sim_data[test_index,],
    #                               mod_type="reg",
    #                               y = "y", k=k)
    sim_times$knn_kd_fittime[sim] <- knnTest3$fittime
    sim_times$knn_kd_predtime[sim] <- knnTest3$predtime
    #-------
    # time the fitting of the iq bin model
    timer <- Sys.time()
    iqnn_mod <- iqnn(sim_data[-test_index,], y="y", bin_cols=xcols,
                     nbins=rep(d,p), jit=rep(0.001,p), stretch=TRUE, tol=rep(5,p))
    sim_times$iqnn_fittime[sim] <- as.numeric(Sys.time() - timer,units="mins")
    # time the prediction using iq bin model
    timer <- Sys.time()
    iqnn_preds <- iqnn_predict(iqnn_mod, sim_data[test_index,],strict=FALSE)
    sim_times$iqnn_predtime[sim] <- as.numeric(Sys.time() - timer,units="mins")
    print(paste("completed p =",p,", d =",d,", k =",k,", n =",n,", m =",m),sep="")
  }
  sim_all <- rbind(sim_all,sim_times)
}

# write.csv(sim_times,"simulationTimesPowersOf2.csv", row.names=FALSE)

# sim_times <- read.csv("simulationTimesPowersOf2.csv")

head(sim_times)

types <-c("knn_brute","knn_cover","knn_kd","iqnn")

# pull out fittimes and move from wide to tall
sim_fit_times <- sim_times %>%
  select(-ends_with("predtime")) %>%
  gather(key="type",value="time",knn_brute_fittime:iqnn_fittime) %>%
  mutate(stage = "fitting",
         type = str_sub(type,1,-9),
         plabel = factor(paste0("p == 2^",log2(p)), levels=paste0("p == 2^",sort(unique(log2(sim_times$p))))),
         klabel = factor(paste0("k == 2^",log2(k)), levels=paste0("k == 2^",sort(unique(log2(sim_times$k))))) )
head(sim_fit_times)

ggplot()+
  geom_hline(yintercept = 0)+
  geom_line(aes(x=n, y=time, color=type),
            size=2,data=sim_fit_times) +
  facet_grid(plabel~klabel, labeller = label_parsed)+
  scale_y_continuous(trans="log10") +
  scale_x_continuous(trans="log2", breaks=2^seq(4,20,by=4), 
                     labels=paste0("2^",seq(4,20,by=4))) +
  labs(title="fitting times")


sim_pred_times <- sim_times %>%
  select(-ends_with("fittime")) %>%
  gather(key="type",value="time",knn_brute_predtime:iqnn_predtime) %>%
  mutate(stage="prediction",
         type=str_sub(type,1,-10))


sim_plot_data <- rbind(sim_fit_times,sim_pred_times)
head(sim_plot_data)
tail(sim_plot_data)



ggplot()+
  geom_hline(yintercept = 0)+
  geom_line(aes(x=k, y=time,color=type),
            size=2,data=sim_pred_times) +
  facet_grid(p~d)+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10") +
  labs(title="prediction times")


###-----------------------------------------------------------------------------------------------------
# Compare using Walter "medium data sets"

setwd("C:\\Users\\maurerkt\\Google Drive\\AFRLSFFP\\Fall2017\\mediumDatasets")

# "mediumDatasets" attributes
medium_sets <- c("abalone","banana","marketing","optdigits","satimage","waveform")
responses <- c("Sex","Class","Sex","Class","Class","Class")
sizes <- c(4174,5300,6876,5620,6435,5000)



nreps <- 1000 # number of times to run k-fold comparisons
results <- data.frame(data_name=rep(medium_sets,each=nreps),obs = NA, nn_size = NA, cv_accuracy = NA, 
                      time_fit = NA, time_pred = NA, seed = NA)
results_all_list <- list(results_iqnn=results, results_knn=results,
                         results_knn_cover=results, results_knn_kd=results)

max_p <- 2 # max number of dimensions for inputs
cv_k <- 10 # cv folds
k <- 3 # knn size




# # Loop over all data sets to clean and save to CSV for walter
# for(set in 1:6){
#   ## load and clean data in preparation for testing speed/accuracy with k-fold CV process
#   data <- RWeka::read.arff(paste0(medium_sets[set],".arff"))
#   # name of response variable
#   y <- responses[set]
#   # use helper function to standardize and drop non-numeric/constant-valued input variables
#   data <- clean_data_for_iqnn_knn(data,y)
#   
#   ## Variable selection
#   # Find column names in order of importance for randomForest (heuristic for doing variable selection)
#   myforest <- randomForest::randomForest(as.formula(paste0("as.factor(as.character(",y,"))~ .")) , data=sample_n(data,1000))
#   important_cols <- dimnames(importance(myforest))[[1]][order(importance(myforest),decreasing=TRUE)]
#   # allow a cap to be put on number of variables considered
#   p <- min(length(important_cols),max_p)
#   
#   clean_data <- data[,c(y,important_cols[1:p])]
#   # write.csv(clean_data,paste0("cleaned_",medium_sets[set],".csv"),row.names=FALSE)
# }



big_timer <- Sys.time()
# Loop over all data sets and repetitions to record accuracies and times. 
for(set in 1:6){
  ## load and clean data in preparation for testing speed/accuracy with k-fold CV process
  data <- RWeka::read.arff(paste0(medium_sets[set],".arff"))
  # name of response variable
  y <- responses[set]
  # use helper function to standardize and drop non-numeric/constant-valued input variables
  data <- clean_data_for_iqnn_knn(data,y)

  ## Variable selection
  # Find column names in order of importance for randomForest (heuristic for doing variable selection)
  myforest <- randomForest(as.formula(paste0("as.factor(as.character(",y,"))~ .")) , data=sample_n(data,1000))
  important_cols <- dimnames(importance(myforest))[[1]][order(importance(myforest),decreasing=TRUE)]
  # allow a cap to be put on number of variables considered
  p <- min(length(important_cols),max_p)
  
  ## Parameterize for binning to best match k-nn structure specified with n, k, p, and cv_k
  train_n <- nrow(data)*((cv_k-1)/cv_k)
  nbins <- find_bin_root(n=train_n,k=k,p=p)
  bin_cols <- important_cols[1:p]
  
  ## Compare knn/iqnn method timing and accuracy with k-fold CV


  # loop over nreps for each method
  for(rep in 1:nreps){
    # set seed for CV partitioning so that each method uses same train/test splits
    seed <- sample(1:100000,1)
    # pick order for methods at random
    method_order <- sample(1:4)
    # seed <-  12345 # fixed value to check if all reps identical predictions made **Confirmed as identical for accuracy**
    for(method in method_order){
      # find 10-fold CV predictions, Record time/accuracy for each
      if(method==1){
        pred_times <- iqnn_cv_predict_timer(data=data, y=y, mod_type="class", bin_cols=bin_cols,
                                            nbins=nbins, jit=rep(0.001,length(nbins)), strict=FALSE, 
                                            cv_k=cv_k, seed=seed)
      } else if(method==2){
        pred_times <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
                                        knn_algorithm = "brute", seed=seed)
      } else if(method==3){
        pred_times <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
                                        knn_algorithm = "cover_tree", seed=seed)
      } else {
        # pred_times <- kdtree_nn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
        #                                       eps=1, seed=seed)
        pred_times <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k,
                                        knn_algorithm = "kd_tree", seed=seed)
      }
      # store results in proper mehtod/set/rep values
      results_all_list[[method]]$obs[(set-1)*nreps + rep] <- nrow(data)
      results_all_list[[method]]$nn_size[(set-1)*nreps + rep] <- k
      results_all_list[[method]]$cv_accuracy[(set-1)*nreps + rep] <- sum(pred_times$preds==data[,y])/nrow(data)
      results_all_list[[method]]$time_fit[(set-1)*nreps + rep] <- pred_times$time_fit
      results_all_list[[method]]$time_pred[(set-1)*nreps + rep] <- pred_times$pred_time
      results_all_list[[method]]$seed[(set-1)*nreps + rep] <- seed
    }
  }
}
Sys.time() - big_timer
str(results_all_lsist)

# Save it
# save(results_all_list, file="iqnn_knn_comparisons.Rdata")
load(file="iqnn_knn_comparisons.Rdata")

# process for table form (average times/accuracy over each trial)
results_all <- data.frame(do.call("rbind", results_all_list),
                          type=rep(c("iqnn","knn - brute ","knn - cover tree","knn - kd tree"),each=nrow(results_all_list$results_iqnn))) %>%
  group_by(data_name,obs,nn_size,type) %>%
  summarize(avg_cv_accuracy=mean(cv_accuracy,na.rm=TRUE),
            avg_time_fit=mean(time_fit,na.rm=TRUE),
            avg_time_pred=mean(time_pred,na.rm=TRUE)) %>%
  as.data.frame()
head(results_all)
# write.csv(results_all,"resultsToShareKarsten.csv", row.names=FALSE)

# Combine into data frame for plots
results_all <- data.frame(do.call("rbind", results_all_list),
                          type=rep(c("iqnn","knn - brute ","knn - cover tree","knn - kd tree"),each=nrow(results_all_list$results_iqnn))) %>%
  gather(key="metric",value="value",cv_accuracy:time_pred) %>% 
  group_by(data_name,obs,nn_size,type,metric) %>%
  summarize(value=mean(value,na.rm=TRUE)) %>%
  as.data.frame() %>%
  mutate(data_name = factor(data_name, levels=medium_sets[order(sizes)]),
         metric_pretty = factor(metric, labels=c("Test Accuracy Rate","Preprocess Time (sec)","Prediction Time (sec)")),
         metric_pretty = factor(metric, labels=c("Test Accuracy Rate","Preprocess Time (sec)","Prediction Time (sec)"))) 
head(results_all)
levels(results_all$metric_pretty)

label_data <- arrange(unique(results_all[,c("metric_pretty","data_name")]),data_name)
label_data$n <- NA
label_data[label_data$metric_pretty=="Prediction Time (sec)",]$n <- paste0("n=",sort(sizes))
label_data

# plot the accuracy/fit time/prediction time
ggplot()+
  geom_hline(yintercept = 0)+
  geom_line(aes(x=data_name, y=value,color=type, group=type),size=1, data=results_all)+
  facet_grid(metric_pretty ~ ., scales="free_y") +
  theme_bw()+
  labs(title="3-NN Classifier vs IQNN Classifier (~3 per bin)",
       subtitle="Based on averages across 1000 repetitions of 10-fold CV",
       x="Data Set", y="", 
       caption="Data Source: UCI Machine Learning Data Repository")+
  geom_text(aes(x=data_name, label=n), y=.2,data=label_data)

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
load("~/onePercentSample.Rdata")
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



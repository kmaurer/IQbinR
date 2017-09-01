### Comparison of knn and iqnn for regression setting

# Load up functions and packages for iqnn and knn regression
library(devtools)
install_github(repo="kmaurer/iqbin")

library(iqbin)
help(package="iqbin")
?iqnn

library(FNN)
library(tidyverse)
library(randomForest)
library(RANN)
library(mvtnorm)

source("iqnn_knn_comparison_functions.R")

#---------------------------------------------------------------------------

# simulate data from different numbers of dimensions, bins per dimension and neighborhood size
ps = 2:4 # number of dimensions
bs = 11:12 # number of bins per dimension
ks = c(1,10,100) # number in iq-neighborhood
combinations <- expand.grid(ps,bs,ks)
names(combinations) <- c("p","b","k")
sim_times <- data.frame(combinations, 
                        n = NA, knntime_brute=NA,
                        knntime_cover=NA, knntime_kd=NA,iqnntime_total=NA,
                        iqfittime=NA, iqpredtime=NA)
set.seed(12345)
for(sim in 1:nrow(sim_times)){
  p = sim_times$p[sim]
  b = sim_times$b[sim]
  k = sim_times$k[sim]
  
  # sim data with number of observations to align number of bins with knn sizes
  n <- b^p*k*2 
  sim_times$n[sim] <- n
  sim_data <- data.frame(sapply(1:p, function(x) rnorm(n)),
                            y=rnorm(n,100,10))
  # rebuild column names to proper dimension 
  xcols <- names(sim_data)[1:p]
  
  test_index <- 1:n/2
  #-------
  # time the knn predictions with brute force
  timer <- Sys.time()
  #!# need to add time taken for standardization?
  knnTest <- knn.reg(train = sim_data[-test_index,xcols],
                     test = sim_data[test_index,xcols],
                     y = sim_data$y[-test_index], k = k, algorithm = "brute")
  sim_times$knntime_brute[sim] <- as.numeric(Sys.time() - timer,units="mins")
  #-------
  # time the knn predictions with cover tree
  timer <- Sys.time()
  knnTest2 <- knn.reg(train = sim_data[-test_index,xcols],
                     test = sim_data[test_index,xcols],
                     y = sim_data$y[-test_index], k = k, algorithm = "cover_tree")
  sim_times$knntime_cover[sim] <- as.numeric(Sys.time() - timer,units="mins")
  #-------  
  # time the knn predictions with kd_tree
  timer <- Sys.time()
  knnTest3 <- kdtree_nn_predict(train = sim_data[-test_index,xcols],
                                test = sim_data[test_index,xcols],
                                k=k)
  sim_times$knntime_kd[sim] <- as.numeric(Sys.time() - timer,units="mins")
  #-------
  # time the fitting of the iq bin model
  timer <- Sys.time()
  iqnn_mod <- iqnn(sim_data[-test_index,], y="y", bin_cols=xcols,
                   nbins=rep(b,p), jit=rep(0.001,p), stretch=TRUE, tol=rep(5,p))
  sim_times$iqfittime[sim] <- as.numeric(Sys.time() - timer,units="mins")
  # time the prediction using iq bin model
  timer <- Sys.time()
  iqnn_preds <- iqnn_predict(iqnn_mod, sim_data[test_index,],strict=TRUE)
  sim_times$iqpredtime[sim] <- as.numeric(Sys.time() - timer,units="mins")
  print(paste("completed p =",p,", b =",b,", k =",k),sep="")
}

sim_times$iqnntime_total <- sim_times$iqfittime + sim_times$iqpredtime

sim_times
# write.csv(sim_times,"simulationTimesBruteCoverKdIqnn.csv", row.names=FALSE)

sim_times <- read.csv("simulationTimesBruteCoverKdIqnn.csv")
library(tidyverse)
sim_plot_data <- sim_times %>%
  # select(-iqfittime,-iqpredtime) %>%
  gather(key="type",value="time",knntime_brute:iqpredtime) %>%
  filter(p==4, b>=9)
head(sim_plot_data)

ggplot()+
  geom_line(aes(x=k, y=time,color=type),
            size=2,data=sim_plot_data) +
  facet_grid(p~b)+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10")

###-----------------------------------------------------------------------------------------------------
# Compare using Walter "medium data sets"

setwd("C:\\Users\\maurerkt\\Google Drive\\AFRLSFFP\\Fall2017\\mediumDatasets")

medium_sets <- c("abalone","banana","marketing","optdigits","satimage","waveform")
responses <- c("Sex","Class","Sex","Class","Class","Class")

sizes <- c(4174,5300,6876,5620,6435,5000)
quant_cols <- c(7,2,12,64,36,40)

results <- data.frame(data_name=medium_sets,obs = NA, nn_size = NA, cv_accuracy = NA, 
                      time_fit = NA, time_pred = NA)
results_iqnn <- results
results_knn <- results
results_knn_cover <- results
results_knn_kd <- results

# P <- 3 # max number of dimensions for inputs
# B <- 10 # bins per dimension
max_p <- 2
cv_k <- 10 # cv folds
k <- 3 # knn size

# Loop over all data sets, record accuracies and times. 
for(set in 1:6){
  # load data
  data <- RWeka::read.arff(paste0(medium_sets[set],".arff"))
  # name of response variable
  y <- responses[set]
  # Fix accidental spaces before some column names
  names(data) <- stringr::str_replace_all(names(data)," ","")
  names(data) <- stringr::str_replace_all(names(data),"-","_")
  # Drop Rows with Missing Values
  data <- na.omit(data)
  # keep only numeric input columns
  keeper_cols <- sapply(data, is.numeric)
  keeper_cols[which(names(data)==y)] <- TRUE
  data <- data[,keeper_cols]
  # Find names in order of importance for randomForest (heuristic for doing variable selection)
  myforest <- randomForest(as.formula(paste0("as.factor(as.character(",y,"))~ .")) , data=sample_n(data,1000))
  important_cols <- dimnames(importance(myforest))[[1]][order(importance(myforest),decreasing=TRUE)]
  
  p <- min(length(important_cols),max_p)
  
  # nbins <- rep(B,min(ncol(data)-1,P))
  # bin_cols <- important_cols[min(ncol(data)-1,P)]
  train_n <- nrow(data)*((cv_k-1)/cv_k)
  nbins <- find_bin_root(n=train_n,k=k,p=p)
  bin_cols <- important_cols[1:p]
  
  # find 10-fold CV predictions, Record time/accuracy for each
  set.seed(12345)
  iqnnmod <- iqnn_cv_predict_timer(data=data, y=y, mod_type="class", bin_cols=bin_cols,
                                   nbins=nbins, jit=rep(0.001,length(nbins)), strict=FALSE, cv_k=cv_k)
  results_iqnn$obs[set] <- nrow(data)
  results_iqnn$nn_size[set] <- train_n/prod(nbins)
  results_iqnn$cv_accuracy[set] <- sum(iqnnmod$preds==data[,y])/nrow(data)
  results_iqnn$time_fit[set] <- iqnnmod$time_fit
  results_iqnn$time_pred[set] <- iqnnmod$pred_time
  
  # find closest equivalent number of neighbors
  # k <- round(nrow(data)*(cv_k-1)/cv_k / prod(nbins) )
  set.seed(12345)
  knnmod <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k, knn_algorithm = "brute")
  results_knn$obs[set] <- nrow(data)
  results_knn$nn_size[set] <- k
  results_knn$cv_accuracy[set] <- sum(knnmod$preds==data[,y])/nrow(data)
  results_knn$time_fit[set] <- knnmod$time_fit
  results_knn$time_pred[set] <- knnmod$pred_time
  
  set.seed(12345)
  knnmod_cover <- knn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k, knn_algorithm = "cover_tree")
  results_knn_cover$obs[set] <- nrow(data)
  results_knn_cover$nn_size[set] <- k
  results_knn_cover$cv_accuracy[set] <- sum(knnmod_cover$preds==data[,y])/nrow(data)
  results_knn_cover$time_fit[set] <- knnmod_cover$time_fit
  results_knn_cover$time_pred[set] <- knnmod_cover$pred_time
  
  set.seed(12345)
  knnmod_kd <- kdtree_nn_cv_pred_timer(data=data, y=y, x_names=bin_cols, cv_k=cv_k, k=k)
  results_knn_kd$obs[set] <- nrow(data)
  results_knn_kd$nn_size[set] <- k
  results_knn_kd$cv_accuracy[set] <- sum(knnmod_kd$preds==data[,y])/nrow(data)
  results_knn_kd$time_fit[set] <- knnmod_kd$time_fit
  results_knn_kd$time_pred[set] <- knnmod_kd$pred_time
}
results_iqnn
results_knn
results_knn_cover
results_knn_kd

results_all <- data.frame(rbind(results_iqnn, results_knn,results_knn_cover,results_knn_kd),
                          type=rep(c("iqnn","knn","knn_cover","knn_kd"),each=6)) %>%
  gather(key="metric",value="value",cv_accuracy:time_pred)
results_all

ggplot()+
  geom_hline(yintercept = 0)+
  geom_line(aes(x=data_name, y=value,color=type, group=type),size=1, data=results_all)+
  facet_grid(metric ~ ., scales="free_y") +
  theme_bw()+
  labs(title="3-NN Classifier (brute force) vs IQNN Classifier (~3 per bin)",
       subtitle="test accuracy / Preprocess Fit Time (sec) / Predict Time (sec)")



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



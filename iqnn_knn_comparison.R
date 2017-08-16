### Comparison of knn and iqnn for regression setting

# Load up functions and packages for iqnn and knn regression
source("IterativeQuantileBinning.R")
library(FNN)
library(tidyverse)



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

###-------------------------------------------------------------------------
library(mvtnorm)
help(package="mvtnorm")

P=5
B=10
ps = rep(2:P, each=(B-1))
bs = rep(2:B, (P-1)) 

k = 5

sim_times <- data.frame(knntime=NA, iqfittime=NA, iqpredtime=NA, size=NA)
for(sim in 1:length(ps)){
p= ps[sim]
b= bs[sim]
n <- b^p*k*2
sim_times[sim,"size"] <- n
# sim data to proper size
sim_data <- data.frame(x1=rnorm(n),
                       x2=rnorm(n),
                       x3=rnorm(n),
                       x4=rnorm(n),
                       x5=rnorm(n),
                       y=rnorm(n,100,10))
# rebuild column names to proper dimension 
xcols <- paste0("x",1:p)

test_index <- 1:n/2
# # time the knn predictions
# timer <- Sys.time()
#!# need to add time taken for standardization
# knnTest <- knn.reg(train = sim_data[-test_index,xcols],
#                    test = sim_data[test_index,xcols],
#                    y = sim_data$y[-test_index], k = k, algorithm = "brute")
# sim_times$knntime[sim] <- as.numeric(Sys.time() - timer,units="mins")

# time the fitting of the iq bin model 
timer <- Sys.time()
iqnn_mod <- iqnn(sim_data[-test_index,], y="y", bin_cols=xcols,
                 nbins=rep(b,p), jit=rep(0.001,p), stretch=TRUE, tolerance=rep(5,p))
sim_times$iqfittime[sim] <- as.numeric(Sys.time() - timer,units="mins")
# time the prediction using iq bin model
timer <- Sys.time()
iqnn_preds <- predict_iqnn(iqnn_mod, sim_data[test_index,],strict=TRUE)
sim_times$iqpredtime[sim] <- as.numeric(Sys.time() - timer,units="mins")
}
sim_times
#write.csv(sim_times,"simulationTimesNestedLists.csv", row.names=FALSE)

#------------------------------------------------------------------------------------------------------
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
                 nbins=c(7,7,6,6), jit=rep(0.00001,4), tolerance=rep(0.0001,4))
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



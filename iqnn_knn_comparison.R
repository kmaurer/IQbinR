### Comparison of knn and iqnn for regression setting

# Load up functions and packages for iqnn and knn regression
source("IterativeQuantileBinning.R")
library(FNN)
library(tidyverse)

## Load up data for testing 
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

# Total number of possible IQbin patterns 4 X variables: 
#   for combination of p columns there exists p! combinations of iterative binning
sum(choose(4,2:4)*factorial(2:4))

## Check that we can fit models to batting career data
# iqdef <- iterative_quant_bin(dat=bb_players, bin_cols=c("b2","b3","hit","ab"),
#                     nbins=c(2,2,2,2), jit=rep(0.001,4), output="both")
# 
# iqnn_mod <- iqnn(dat=bb_players, y="hr", bin_cols=c("b2","b3","hit","ab"),
#                  nbins=c(2,2,2,2), jit=rep(0.001,4))
# cv_iqnn(iqnn_mod,bb_players, cv_method="kfold", cv_k=5, strict=FALSE)
# cv_iqnn(iqnn_mod,bb_players, cv_method="LOO", strict=FALSE)

#### knn.reg
# need standardized variables
bb_players_st <- bb_players %>%
  mutate(b2 = scale(b2),
         b3 = scale(b3),
         hit = scale(hit),
         ab = scale(ab))
head(bb_players_st)

test_index <- 1:100
knnTest <- knn.reg(train = bb_players_st[-test_index,c("b2","b3","hit","ab")],
                   test = bb_players_st[test_index,c("b2","b3","hit","ab")],
                   y = bb_players_st$hr[-test_index], k = 5, algorithm = "brute")
knnTest$pred

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
cv_preds <- cv_pred_knn(dat=bb_players_st, y_name="hr", x_names=c("b2","b3","hit","ab"), 
           cv_method="kfold", cv_k = 10, k=5, knn_algorithm = "brute")
head(cv_preds)
# cv_preds <- cv_knn.reg(dat=bb_players_st, y_name="hr", x_names=c("b2","b3","hit","ab"), 
#                        cv_method="LOO", k=5, knn_algorithm = "brute")
# head(cv_preds)






#--------------------------------------
### Cross Validation for assessment for iqnn models
cv_knn <- function(dat, y_name, x_names, cv_method="kfold", cv_k = 10, k=5, knn_algorithm = "brute"){
  dat <- as.data.frame(dat)
  cv_preds <- cv_pred_knn(dat=dat, y_name=y_name, x_names=x_names, cv_method=cv_method, cv_k=cv_k, k=k, knn_algorithm = knn_algorithm)
  PRESS <- sum((dat[,iqnn_mod$y]-cv_preds)^2)
  MSE <- PRESS/nrow(dat)
  RMSE <- sqrt(MSE)
  c(PRESS=PRESS,MSE=MSE,RMSE=RMSE)
}
cv_knn(dat=bb_players_st, y_name="hr", x_names=c("b2","b3","hit","ab"), 
            cv_method="kfold", cv_k = 20, k=5, knn_algorithm = "brute")



#---------------------------------------------------------------------------------------------------
# Timing simulations


# from building model to predicting for new

test_index <- 1:9585
timer <- Sys.time()
knnTest <- knn.reg(train = bb_players_st[-test_index,c("b2","b3","hit","ab")],
                   test = bb_players_st[test_index,c("b2","b3","hit","ab")],
                   y = bb_players_st$hr[-test_index], k = 5, algorithm = "brute")
Sys.time() - timer

timer <- Sys.time()
iqnn_mod <- iqnn(bb_players_st[-test_index,], y="hr", bin_cols=c("b2","b3","hit","ab"),
                 nbins=c(7,7,6,6), jit=rep(0.001,4))
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



###-------------------------------------------------------------------------
library(mvtnorm)
help(package="mvtnorm")

P=4
B=3
ps = rep(2:P, each=(B-1))
bs = rep(2:B, (P-1)) 

k = 50

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
                       y=rnorm(n,100,10))
# rebuild column names to proper dimension 
xcols <- paste0("x",1:p)

test_index <- 1:n/2
# time the knn predictions
timer <- Sys.time()
knnTest <- knn.reg(train = sim_data[-test_index,xcols],
                   test = sim_data[test_index,xcols],
                   y = sim_data$y[-test_index], k = k, algorithm = "brute")
sim_times$knntime[sim] <- Sys.time() - timer
# time the fitting of the iq bin model 
iqnn_mod <- iqnn(sim_data[-test_index,], y="y", bin_cols=xcols,
                 nbins=rep(b,p), jit=rep(0.001,p), stretch=TRUE, tolerance=rep(5,p))
sim_times$iqfittime[sim] <- Sys.time() - timer
# time the prediction using iq bin model
timer <- Sys.time()
iqnn_preds <- predict_iqnn(iqnn_mod, sim_data[test_index,],strict=TRUE)
sim_times$iqpredtime[sim] <- Sys.time() - timer
}

write.csv(sim_times,"simulationTimes.csv", row.names=FALSE)

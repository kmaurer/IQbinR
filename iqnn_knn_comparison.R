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
sum(choose(4,2:4)*factortbl_dfial(2:4))

# Check that we can fit models to batting career data
iqdef <- iterative_quant_bin(dat=bb_players, bin_cols=c("b2","b3","hit","ab"),
                    nbins=c(2,2,2,2), jit=rep(0.001,4), output="both")

iqnn_mod <- iqnn(dat=bb_players, y="hr", bin_cols=c("b2","b3","hit","ab"),
                 nbins=c(2,2,2,2), jit=rep(0.001,4))
cv_iqnn(iqnn_mod,bb_players, cv_method="kfold", cv_k=5, strict=FALSE)
cv_iqnn(iqnn_mod,bb_players, cv_method="LOO", strict=FALSE)


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
cv_knn <- function(iqnn_mod, dat, cv_method="kfold", cv_k=10, strict=FALSE){
  dat <- as.data.frame(dat)
  cv_preds <- cv_knn.reg(iqnn_mod, dat, cv_method, cv_k, strict)
  PRESS <- sum((dat[,iqnn_mod$y]-cv_preds)^2)
  MSE <- PRESS/nrow(dat)
  RMSE <- sqrt(MSE)
  c(PRESS=PRESS,MSE=MSE,RMSE=RMSE)
}
# iqnn_mod <- iqnn(iris, y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#                  nbins=c(3,5,2), jit=rep(0.001,3))
# cv_iqnn(iqnn_mod,iris, cv_method="kfold", cv_k=10, strict=FALSE)
# cv_iqnn(iqnn_mod,iris, cv_method="LOO", strict=FALSE)


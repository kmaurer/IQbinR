### Simulated Data Timing Comparisons
# Load up functions and packages for iqnn and knn regression
library(devtools)
install_github(repo="kmaurer/iqbin")

library(iqbin)
help(package="iqbin")

library(FNN)
library(tidyverse)
library(stringr)
library(randomForest)
library(RANN)
library(mvtnorm)
library(gridExtra)
library(data.table)

setwd("~/GitHub/iqnnProject")
source("helper_functions_for_study.R")

#---------------------------------------------------------------------------
### Simulation-Based Timing Study

# simulate data from different numbers of dimensions, bins per dimension and neighborhood size
# parameterization based on 2^x scaling

## Set sample sizes and number of trials (dependent on sample size) and 
# generate vector of seed values for simulating data and generating random ordering of method timings

### For the 500 repetitions of small sample sizes 
ns= 2^(seq(4,14,by=2))
ntrials = 500
set.seed(12345)
trial_sim_seeds <- sample(1:100000000,1000000)

### For the 50 repetitions of big sample sizes (uncomment to run large n simulations)
# ns= 2^(seq(16,20,by=2))
# ntrials = 50
# set.seed(6789)
# trial_sim_seeds <- sample(1:100000000,1000000)

# align all neighborhood/bin parameters for comparability
ks= 2^(seq(0,14,by=2))
ps= 2^(1)
ms = c(1000)
combinations <- expand.grid(ns,ps,ks,ms)
names(combinations) <- c("n","p","k","m")
combinations <- combinations %>%
  mutate(d=(n/k)^(1/p)) %>%
  filter(n > k,
         d^p <= 2^14 )

# Initialize a list to store all complete trial timing data frames
sim_all <- list(NULL)
## Loop over all methods for each timing trial
for(trial in 1:ntrials){
  trial_timer <- Sys.time()
  # build empty data frame to store all timing measurements for new trial
  sim_times <- data.frame(combinations,
                          knn_brute_fittime=NA,knn_brute_predtime=NA,
                          knn_cover_fittime=NA,knn_cover_predtime=NA,
                          knn_kd_fittime=NA,knn_kd_predtime=NA,
                          iqnn_fittime=NA, iqnn_predtime=NA)
  for(sim in 1:nrow(sim_times)){
    # Set params for model structures
    p = sim_times$p[sim]; d = sim_times$d[sim]; k = sim_times$k[sim]; m = sim_times$m[sim]; n = sim_times$n[sim]
    # look up seed unique to sim/trial combination - trial=slow running index & sim=fast running index
    unique_seed <- trial_sim_seeds[(trial-1)*nrow(sim_times)+sim]
    # select method order unique to trial/parameterization combination (reproducable)
    set.seed(unique_seed)
    method_order <- sample(1:4)
    # set seed and sim data with number of observations to align number of bins with knn sizes (reproducable)
    set.seed(unique_seed)
    sim_data <- data.frame(sapply(1:p, function(x) runif(n+m)), y=runif(n+m))
    test_index <- 1:m
    ### Loop over 4 methods, selecting order at random to avoid favoring method in any given trial
    # Record time for fit/predict for each
    for(method in method_order){
      if(method==1){
        # time the knn - brute force
        timer <- Sys.time()
        knnTest <- knn.reg(train = sim_data[-test_index,1:p],
                           test = sim_data[test_index,1:p],
                           y = sim_data$y[-test_index], k = k, algorithm = "brute")
        sim_times$knn_brute_fittime[sim] <- 0
        sim_times$knn_brute_predtime[sim] <- as.numeric(Sys.time() - timer,units="mins")
        remove(knnTest)
      } else if(method==2){
        # time the aknn - cover tree
        knnTest2 <- aknn_predict(train_x = sim_data[-test_index,1:p],
                                 test_x = sim_data[test_index,1:p],
                                 y = sim_data[-test_index,"y"],
                                 mod_type="reg", k=k,algorithm="cover_tree")
        sim_times$knn_cover_fittime[sim] <- knnTest2$fittime
        sim_times$knn_cover_predtime[sim] <- knnTest2$predtime
        remove(knnTest2)
      } else if(method==3){
        # time the aknn - kd tree
        knnTest3 <- aknn_predict(train_x = sim_data[-test_index,1:p],
                                      test_x = sim_data[test_index,1:p],
                                      y = sim_data[-test_index,"y"],
                                      mod_type="reg", k=k)
        sim_times$knn_kd_fittime[sim] <- knnTest3$fittime
        sim_times$knn_kd_predtime[sim] <- knnTest3$predtime
        remove(knnTest3)
      } else{
        # time the fitting of the iq bin model
        timer <- Sys.time()
        iqnn_mod <- iqnn(sim_data[-test_index,], y="y", bin_cols=names(sim_data)[1:p],
                         nbins=rep(d,p), jit=rep(0.001,p), stretch=TRUE, tol=rep(5,p))
        sim_times$iqnn_fittime[sim] <- as.numeric(Sys.time() - timer,units="mins")
        # time the prediction using iq bin model
        timer <- Sys.time()
        iqnn_preds <- iqnn_predict(iqnn_mod, sim_data[test_index,],strict=FALSE)
        sim_times$iqnn_predtime[sim] <- as.numeric(Sys.time() - timer,units="mins")
        remove(iqnn_mod); remove(iqnn_preds)
      }
    }
    print(paste("completed p =",p,", d =",d,", k =",k,", n =",n,", m =",m),sep="")
  }
  remove(sim_data)
  print(paste("Trial",trial,"completed in", as.numeric(Sys.time() - trial_timer,units="mins"),"minutes"))
  sim_all[[trial]] <- sim_times
}
# save(sim_all,# file="sim_all_small_n.Rdata")  # use with small sample simulations (extra # inserted to avoid accidental overwrite)
# save(sim_all,# file="sim_all_big_n.Rdata") # use with big sample simulations (extra # inserted to avoid accidental overwrite)

#--------------------------------------------------------------------------------------------------------------------
### Timing Data cleaning 

## Load "small n" and "large n" timing lists, stacking each into a data frame with do.call/rbind
# load(file="sim_all_small_n.Rdata")
sim_times_small <- do.call("rbind", sim_all)
# load(file="sim_all_big_n.Rdata")
sim_times_big <- do.call("rbind", sim_all)
# stack both into single data frame of timing values
sim_times_all <- rbind(sim_times_big,sim_times_small)

# group by parameterizations, then average time for each combination
sim_times_agg <- sim_times_all %>% 
  group_by(n,p,k,m,d) %>%
  summarize(count= n(), 
            knn_brute_fittime = mean(knn_brute_fittime), knn_brute_predtime = mean(knn_brute_predtime),
            knn_cover_fittime = mean(knn_cover_fittime), knn_cover_predtime = mean(knn_cover_predtime),
            knn_kd_fittime = mean(knn_kd_fittime), knn_kd_predtime = mean(knn_kd_predtime),
            iqnn_fittime = mean(iqnn_fittime), iqnn_predtime = mean(iqnn_predtime)) %>%
  as.data.frame()
head(sim_times_agg)
tail(sim_times_agg)

## Save this cleaned/aggregated data set to Rdata file for use in paper (loaded in "setup" code chunk of iqnnPaper.Rnw file)
# save(sim_times_agg, file="sim_times_agg.Rdata") # save this iteration for general work

###------------------------------------------------------------------------------------------------------------------------------
# Construction of plots used in section 3.1 of manuscript
# (found in "sim_time" code chunk of iqnnPaper.Rnw file)

# pull out fittimes and move from wide to tall
sim_fit_times <- sim_times_agg %>%
  select(-ends_with("predtime")) %>%
  gather(key="type",value="time",knn_brute_fittime:iqnn_fittime) %>%
  mutate(stage = "fitting",
         type = str_sub(type,1,-9),
         plabel = factor(paste0("p == 2^",log2(p)), levels=paste0("p == 2^",sort(unique(log2(sim_times_agg$p))))),
         klabel = factor(paste0("k == 2^",log2(k)), levels=paste0("k == 2^",sort(unique(log2(sim_times_agg$k))))) )

# pull out prediction times and move from wide to tall
sim_pred_times <- sim_times_agg %>%
  select(-ends_with("fittime")) %>%
  gather(key="type",value="time",knn_brute_predtime:iqnn_predtime) %>%
  mutate(stage = "predicting",
         type = str_sub(type,1,-10),
         plabel = factor(paste0("p == 2^",log2(p)), levels=paste0("p == 2^",sort(unique(log2(sim_times_agg$p))))),
         klabel = factor(paste0("k == 2^",log2(k)), levels=paste0("k == 2^",sort(unique(log2(sim_times_agg$k))))),
         dlabel = factor(paste0("delta == 2^",log2(d)), levels=paste0("delta == 2^",sort(unique(log2(sim_times_agg$d))))) )


# define color palette
my_colors <- RColorBrewer::brewer.pal(n=4,name="Set1")[c(1,4,2,3)]

# combine and define transformation of strings to clarify labels
combined_times <- rbind(data.frame(sim_fit_times,dlabel=NA),
                        sim_pred_times)

combined_times$stage_pretty <- factor(ifelse(combined_times$stage=="fitting","Preprocessing","Prediction"),
                                      levels=c("Preprocessing","Prediction"))
combined_times$type_pretty <- factor(combined_times$type, labels=c("IQNN   ","KNN-brute   ","AKNN-cover   ", "AKNN-kd   "))

### lineplot of times with overlayed points (allows for better transition to grayscale)
ggplot()+
  geom_hline(yintercept = 0)+
  # geom_vline(xintercept = 2^14)+
  geom_line(aes(x=n, y=time, color=type_pretty, linetype=type_pretty),
            size=.6,data=combined_times) +
  geom_point(aes(x=n, y=time, color=type_pretty, shape=type_pretty),
             size=2.5,data=combined_times) +
  # geom_text(aes(x=n, y=time, label=dlabel),
  #           data=filter(sim_pred_times, type=="iqnn"), parse=TRUE) +
  facet_grid(stage_pretty~klabel, labeller = label_parsed)+
  scale_y_continuous("Time in log(mins)", trans="log", breaks=c(0.0001,.001,.01,.1,1,10,100)) +
  scale_x_continuous(trans="log2", breaks=2^seq(4,20,by=4), 
                     labels=parse(text=paste0("2^",seq(4,20,by=4)))) +
  # scale_color_brewer(palette="Set1")+
  # scale_linetype_manual("Neighborhood Algorithm:  ",values=c("solid","dotted","twodash","dotdash"))+
  scale_linetype_manual("Neighborhood Algorithm:  ",values=c(1,1,1,1))+
  scale_shape_manual("Neighborhood Algorithm:  ", values=c(2,5,0,1))+
  scale_color_manual("Neighborhood Algorithm:  ", values=my_colors)+
  labs(x="Training Size (n)")+
  theme_bw()+
  theme(legend.position = "bottom")
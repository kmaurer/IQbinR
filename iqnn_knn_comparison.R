### Comparison of knn and iqnn for regression and classification models
# Each model type uses 10 data sets from UCI and KEEL repos
# For each set, we clean and standardize to two continuous numeric features and the response value
# Then we tune each model (saving tuned parameterization)
# lastly using tuned setting we test accuracy for each.

#-------------------------------------------------------------------------------------------------------
### Preliminaries

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

# All data sets downloaded from web repos (links in Appendix B) and saved as setname_raw.Rdata

# vectors for setname, response columns, and initialize data sizes vector
all_sets <- c("iris","wpbc","pima","yeast","abalone","waveform","optdigits","satimage","marketing","seizure")
all_responses <- c("V5","V2","V9","V10","Sex","Class","Class","Class","Sex","y")
all_sizes <- rep(NA,10)

# # Clean all using helper function (commented out to avoid accidental overwrite of files used in study)
# setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\classification\\raw")
# for(set in 1:length(all_sets)){
#   load(file=paste0(all_sets[set],"_raw.Rdata"))
#   # name of response variable
#   y <- all_responses[set]
#   # use helper function to standardize and drop non-numeric/constant-valued input variables
#   data <- clean_data_for_iqnn_knn(as.data.frame(data),"y")
#   all_sizes[set] <- nrow(data)
#   save(data,file=paste0(all_sets[set],"_cleaned.Rdata"))
# }

###-----------------------------------------------------------------------------------------------------
# Tune neighborhood/bin parameters for each data set using 10-fold cv then save results for running accuracy tests
setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\classification")

max_p <- 2 # max number of dimensions for inputs
cv_k <- 10 # cv folds

tuned_param_list <- list(NULL)
tuned_performance_list <- list(iqnn=list(NULL),knn=list(NULL))
tune_reps <- 10
overall_timer <- Sys.time()
for(set in 1:10){
  print(set)
  load(file=paste0(all_sets[set],"_cleaned.Rdata"))
  y <- all_responses[set]
  data[,y] <- factor(data[,y])
  
  ## Variable selection
  # Find column names in order of importance for randomForest (approach for doing variable selection described in manuscript)
  set.seed(12345)
  myforest <- randomForest(as.formula(paste0("as.factor(as.character(",y,"))~ .")) , data=sample_n(data,min(1000,nrow(data))))
  important_cols <- dimnames(importance(myforest))[[1]][order(importance(myforest),decreasing=TRUE)]
  # allow a cap to be put on number of variables considered
  p <- min(length(important_cols),max_p)
  
  ## Parameterize for binning to best match k-nn structure specified with n, k, p, and cv_k
  train_n <- floor(nrow(data)*((cv_k-1)/cv_k))
  bin_cols <- important_cols[1:p]
  
  timer <- Sys.time()
  ## Tune the iqnn shoot for no fewer than 2 per bin (otherwise problems with allocation on boundaries)
  set.seed(1234)
  tune_iqnn_out <- list(NULL)
  for(tune_rep in 1:tune_reps){
    tune_iqnn_out[[tune_rep]] <- iqnn_tune(data=data, y=y, mod_type = "class", bin_cols=bin_cols, nbins_range=c(2,floor((train_n/2)^(1/p))),
                                           jit = rep(.0001,length(bin_cols)), stretch = FALSE,strict=FALSE,oom_search = FALSE, cv_k=cv_k)
  }
  tune_iqnn_out_all <- do.call("rbind", tune_iqnn_out) %>%
    as.data.frame() %>%
    group_by(nn_equiv, bin_dims,nbins_total) %>%
    summarize(avg_misclass = mean(error),
              se_misclass = sd(error),
              reps= n()) %>%
    as.data.frame()
  # find smallest MSE row, then figure out which model with lowest number of bins is within 1 SE[MSE] of the lowest (1 SE rule from page 214 in ISL book)
  lowest_row <- which.min(tune_iqnn_out_all$avg_misclass)
  nbins <-  as.numeric(str_split(tune_iqnn_out_all$bin_dims[lowest_row], "X")[[1]])
  # on_par_lowest <- tune_iqnn_out_all$avg_misclass <= tune_iqnn_out_all$avg_misclass[lowest_row]+tune_iqnn_out_all$se_misclass[lowest_row]
  # tuned_row <- tune_iqnn_out_all$nbins_total== min(tune_iqnn_out_all[on_par_lowest,]$nbins_total)
  # nbins <-  as.numeric(str_split(tune_iqnn_out_all$bin_dims[tuned_row], "X")[[1]])
  tuned_performance_list$iqnn[[set]] <- tune_iqnn_out_all
  print(Sys.time()-timer)
  
  ## Tune the knn over same range of neighborhood size equivalents
  timer <- Sys.time()
  set.seed(1234)
  tune_knn_out <- list(NULL)
  for(tune_rep in 1:tune_reps){
    tune_knn_out[[tune_rep]] <- tune_knn_class(data=data, y_name=y, x_names=bin_cols, cv_method="kfold", cv_k = cv_k,
                                         k_values=as.integer(round(tune_iqnn_out_all$nn_equiv)), knn_algorithm = "brute") 
  }
  tune_knn_out_all <- do.call("rbind", tune_knn_out) %>%
    as.data.frame() %>%
    group_by(k) %>%
    summarize(avg_misclass = mean(error),
              se_misclass = sd(error),
              reps= n()) %>%
    as.data.frame()
  # find smallest MSE row, then figure out which model with lowest number of bins is within 1 SE[MSE] of the lowest (1 SE rule from page 214 in ISL book)
  lowest_row_knn <- which.min(tune_knn_out_all$avg_misclass)
  on_par_lowest_knn <- tune_knn_out_all$avg_misclass <= tune_knn_out_all$avg_misclass[lowest_row_knn]+tune_knn_out_all$se_misclass[lowest_row_knn]
  k <- max(tune_knn_out_all[on_par_lowest_knn,]$k)
  tuned_performance_list$knn[[set]] <- tune_knn_out_all
  print(Sys.time()-timer)
  
  tuned_param_list[[set]] <- list(bin_cols=bin_cols, nbins=nbins, k=k, n=nrow(data), cv_k=cv_k)
}
Sys.time()-overall_timer


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
for(set in 1:10){
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

# stack all from list, group by unique simulation/neighborhood conditions, find average cross-validated accuracy, transform for relative interpretability.
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
# save(accuracy_results_list, results_class,tuned_param_list,tuned_performance_list, file="tuned_classification_testing.Rdata")

#---------------------------------------------------------------------------------------------------------------------------------------------
### Create Classification Accuracy Comparison Plots for Figure in Section 3.2.1

# load(file="tuned_classification_testing.Rdata")

#Setting order of data sets
iqnn_class_err <- filter(results_class, type=="iqnn") %>% arrange(avg_cv_accuracy)
iqnn_err_order <- as.character(iqnn_class_err$data_name)

# clean data for plots of accuracy relative to knn as baseline
#   - separate knn from others for "baseline" layers vs "points" layers in plotting
class_plot_data_knn <- results_class %>%
  filter(type == "knn - brute") %>%
  mutate(data_name = factor(data_name,levels=iqnn_err_order),
         x=as.numeric(data_name)-.5,
         xend=as.numeric(data_name)+.5)
class_plot_data_knn$type_pretty <- "KNN-brute"

shiftval=.25 # control the horizontal offset of points on plot
class_plot_data <- results_class %>%
  filter(type != "knn - brute") %>%
  mutate(data_name = factor(data_name,levels=iqnn_err_order),
         shift = (as.numeric(as.factor(as.character(type)))-2)*shiftval,
         diff_err_perc = -diff_acc*100)
class_plot_data$type_pretty <- factor(class_plot_data$type, labels=c("IQNN   ","AKNN-cover ", "AKNN-kd   ")) 

x_label_sizes <- as.numeric(sapply(levels(class_plot_data$data_name), function(x) class_plot_data$obs[class_plot_data$data_name==x][1]))
my_x_ticks <- paste0(levels(class_plot_data$data_name), "\n n=",x_label_sizes)

# plot relative to KNN
#   - use vertical lines to partition between datasets
#   - baseline KNN-brute as flat line segment
#   - consistently offset points for IWNN and AKNN on X, y=accuracy
#   - scale / theme to match previous plots
p1 <- ggplot()+
  geom_segment(aes(x=x, xend=xend,y=100-avg_cv_accuracy*100,yend=100-avg_cv_accuracy*100, linetype="KNN-brute"),
               color=RColorBrewer::brewer.pal(4,"Set1")[4], size=.8, data=class_plot_data_knn) +
  geom_vline(xintercept=seq(0.5,10.5,by=1),linetype=2,color="gray25")+
  geom_point(aes(x=as.numeric(data_name)+shift,y=100-avg_cv_accuracy*100,color=type_pretty,shape=type_pretty),
             size=3,data=class_plot_data)+
  theme_bw() +
  scale_linetype_manual("Baseline: ", values=1)+
  scale_x_continuous("", breaks=1:10,
                     labels=my_x_ticks,
                     limits=c(0.5,10.5)) +
  scale_y_continuous("Misclassification Rates\n (units of % error)",
                     breaks=seq(0,100,by=10),labels=paste0("  ",seq(0,100,by=10),"% ")) +
  scale_color_manual("Model Type: ", values=RColorBrewer::brewer.pal(4,"Set1")[c(1,2,3)])+
  scale_shape_manual("Model Type: ", values=c(2,0,1))+
  theme(panel.grid.major.x =element_blank(),
        panel.grid.minor.x =element_blank(),
        axis.ticks.x = element_blank(),
        # panel.border = element_blank(),
        legend.position = "none")


# plot relative to KNN
#   - use vertical lines to partition between datasets
#   - baseline at 0 to represent KNN-brute
#   - consistently offset points for IWNN and AKNN on X, y=accuracy
#   - scale / theme to match previous plots
p2 <- ggplot()+
  geom_segment(aes(x=.5, xend=10.5,y=0,yend=0, linetype="KNN-brute"),size=.8, color=RColorBrewer::brewer.pal(4,"Set1")[4])+
  geom_vline(xintercept=seq(0.5,10.5,by=1),linetype=2,color="gray25")+
  geom_point(aes(x=as.numeric(data_name)+shift,y=diff_err_perc,color=type_pretty,shape=type_pretty),
             size=3,data=class_plot_data)+
  theme_bw()+
  scale_linetype_manual("Baseline: ", values=1)+
  scale_x_continuous(" ", breaks=1:10,
                     labels=my_x_ticks,
                     limits=c(0.5,10.5))+
  scale_y_continuous("Difference from KNN-brute\n (units of % error)",
                     breaks=seq(-0.5,1.5,by=.5),labels=c("-0.5%","0.0%","0.5%","1.0%","1.5%"),
                     limits=c(min(-.5,min(class_plot_data$diff_err_perc)),max(1.5,max(class_plot_data$diff_err_perc)))) +
  scale_color_manual("Model Type: ", values=RColorBrewer::brewer.pal(4,"Set1")[c(1,2,3)])+
  scale_shape_manual("Model Type: ", values=c(2,0,1))+
  # annotate(geom="text",x=0,y=0,label="bold(KNN-brute)", color=RColorBrewer::brewer.pal(4,"Set1")[4], parse=T,vjust=-.2,hjust=0.4)+
  theme(panel.grid.major.x =element_blank(),
        panel.grid.minor.x =element_blank(),
        axis.ticks.x = element_blank(),
        # panel.border = element_blank(),
        legend.position = "bottom")

grid.arrange(p1,p2,nrow=2,heights=c(1,1.2))

###--------------------------------------------------------------------------------------------
# Accuracy Comparisons for Regression
###-----------------------------------------------------------------------------------------------------
# Organize Regression Data sets

# All data sets downloaded from web repos (links in Appendix B) and saved as setname_raw.Rdata

# vectors for setname, response columns, and initialize data sizes vector
setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\regression")
all_reg_sets <- c("wpbc","wankara","laser","treasury","quake","skillcraft","anacalt","puma","air_quality","ccpp")
all_reg_sizes <- c(198,321,993,1049,2178,3395,4052,8192,9471,9568)

# # Clean all using helper function
# for(set in 1:length(all_reg_sets)){
#   load(file=paste0("raw\\",all_reg_sets[set],"_raw.Rdata"))
#   # name of response variable
#   data$y <- scale(data$y)
#   # use helper function to standardize and drop non-numeric/constant-valued input variables
#   data <- clean_data_for_iqnn_knn(as.data.frame(data),y)
#   save(data,file=paste0(all_reg_sets[set],"_cleaned.Rdata"))
# }

# Tune neighborhood/bin parameters for each data set using 10-fold cv then save results for running accuracy tests
max_p <- 2 # max number of dimensions for inputs
cv_k <- 10 # cv folds
tune_reps <- 20 # number of reps for CV while tuning each parameterization

tuned_reg_param_list <- list(NULL)
tuned_reg_performance_list <- list(iqnn=list(NULL),knn=list(NULL))
for(set in 1:length(all_reg_sets)){
  print(set)
  load(file=paste0(all_reg_sets[set],"_cleaned.Rdata"))
  
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
  tune_iqnn_out <- list(NULL)
  for(tune_rep in 1:tune_reps){
  tune_iqnn_out[[tune_rep]] <- iqnn_tune(data=data, y="y", mod_type = "reg", bin_cols=bin_cols, nbins_range=c(2,floor((train_n/2)^(1/p))),
                             jit = rep(.0000000001,length(bin_cols)), stretch = FALSE,strict=FALSE,oom_search = ifelse(nrow(data)>100000, TRUE,FALSE), cv_k=cv_k)
  }
  tune_iqnn_out_all <- do.call("rbind", tune_iqnn_out) %>%
    as.data.frame() %>%
    group_by(nn_equiv, bin_dims,nbins_total) %>%
    summarize(mse = mean(MSE),
              se_mse = sd(MSE),
              reps= n()) %>%
    as.data.frame()
  # find smallest MSE row, then figure out which model with lowest number of bins is within 1 SE[MSE] of the lowest (1 SE rule from page 214 in ISL book)
  lowest_row <- which.min(tune_iqnn_out_all$mse)
  nbins <-  as.numeric(str_split(tune_iqnn_out_all$bin_dims[lowest_row], "X")[[1]])
  # on_par_lowest <- tune_iqnn_out_all$mse <= tune_iqnn_out_all$mse[lowest_row]+tune_iqnn_out_all$se_mse[lowest_row]
  # tuned_row <- tune_iqnn_out_all$nbins_total== min(tune_iqnn_out_all[on_par_lowest,]$nbins_total)
  # nbins <-  as.numeric(str_split(tune_iqnn_out_all$bin_dims[tuned_row], "X")[[1]])
  tuned_reg_performance_list$iqnn[[set]] <- tune_iqnn_out_all
  
  ## Tune the knn over same range of neighborhood size equivalents
  set.seed(1234)
  tune_knn_out <- list(NULL)
  for(tune_rep in 1:tune_reps){
    tune_knn_out[[tune_rep]] <- tune_knn_reg(data=data, y_name="y", x_names=bin_cols, cv_k = cv_k,
                                              k_values=as.integer(round(tune_iqnn_out_all$nn_equiv)), knn_algorithm = "brute")
  }
  tune_knn_out_all <- do.call("rbind", tune_knn_out) %>%
    as.data.frame() %>%
    group_by(k) %>%
    summarize(mse = mean(MSE),
              se_mse = sd(MSE),
              reps= n()) %>%
    as.data.frame()
  # find smallest MSE row, then figure out which model with lowest number of bins is within 1 SE[MSE] of the lowest (1 SE rule from page 214 in ISL book)
  lowest_row_knn <- which.min(tune_knn_out_all$mse)
  k <- tune_knn_out_all$k[lowest_row_knn]
  # on_par_lowest_knn <- tune_knn_out_all$mse <= tune_knn_out_all$mse[lowest_row_knn]+tune_knn_out_all$se_mse[lowest_row_knn]
  # k <- max(tune_knn_out_all[on_par_lowest_knn,]$k)
  
  tuned_reg_performance_list$knn[[set]] <- tune_knn_out_all
  
  tuned_reg_param_list[[set]] <- list(bin_cols=bin_cols, nbins=nbins, k=k, n=nrow(data), cv_k=cv_k)
}
tuned_reg_param_list
tuned_reg_performance_list
# save(tuned_reg_param_list,tuned_reg_performance_list, file="tuned_reg_param_list.Rdata")

#--------------------------------------------------------------------------------------------------------
# Collecting Accuracy : use average over many 10-fold cv results

# load(file="tuned_reg_param_list.Rdata")
nreps <- 100# number of times to run k-fold comparisons
# # Initialize an empty data structure to put timing/accuracy measurements into
# # Use a list with one df for each method, this is done to allow consistant storage when randomizing order of methods in each trial
results <- data.frame(data_name=rep(all_reg_sets,each=nreps),obs = NA, nn_size = NA, cv_mse = NA, seed = NA)
mse_results_list <- list(results_iqnn=results, results_knn=results,
                              results_knn_cover=results, results_knn_kd=results)

setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\regression")
big_timer <- Sys.time()

# Loop over all data sets and repetitions to record accuracies and times. 
for(set in 1:10){
  print(all_reg_sets[set])
  load(file=paste0(all_reg_sets[set],"_cleaned.Rdata"))
  bin_cols <- tuned_reg_param_list[[set]]$bin_cols
  k <-tuned_reg_param_list[[set]]$k
  n <- tuned_reg_param_list[[set]]$n
  nbins <- tuned_reg_param_list[[set]]$nbins
  # nbins <- tuned_reg_performance_list$iqnn[[set]][tuned_reg_performance_list$knn[[set]]$k ==k,"nbins"][[1]]
  
  ## Compare knn/iqnn method MSE with k-fold CV
  # loop over nreps for each method
  for(rep in 1:nreps){
    print(rep)
    for(method in 1:4){
      seed <- rep # set seed as rep number
      # find 10-fold CV predictions, Record time/accuracy for each
      if(method==1){
        set.seed(seed)
        preds <- iqnn_cv_predict(data=data, y="y", mod_type="reg", bin_cols=bin_cols,
                                                nbins=nbins, jit=rep(0.0000000001,length(nbins)), strict=FALSE,
                                                cv_k=cv_k)
      } else if(method==2){
        set.seed(seed)
        preds <- cv_pred_knn(data=data, y="y", x_names=bin_cols, cv_k=cv_k, k=k, knn_algorithm = "brute")
      } else if(method==3){
        set.seed(seed)
        preds <- cv_pred_knn(data=data, y="y", x_names=bin_cols, cv_k=cv_k, k=k,knn_algorithm = "cover_tree")
      } else {
        set.seed(seed)
        preds<- cv_pred_knn(data=data, y="y", x_names=bin_cols, cv_k=cv_k, k=k,knn_algorithm = "kd_tree")
      }
      # store results in proper mehtod/set/rep values
      mse_results_list[[method]]$obs[(set-1)*nreps + rep] <- nrow(data)
      mse_results_list[[method]]$nn_size[(set-1)*nreps + rep] <- ifelse(method==1, nrow(data)/prod(nbins), k)
      mse_results_list[[method]]$cv_mse[(set-1)*nreps + rep] <- mean((preds - data$y)^2)
      mse_results_list[[method]]$seed[(set-1)*nreps + rep] <- seed
    }
  }
}
Sys.time() - big_timer
str(mse_results_list)

results_reg <- data.frame(do.call("rbind", mse_results_list),
                            type=rep(c("iqnn","knn - brute","knn - cover tree","knn - kd tree"),each=nrow(mse_results_list$results_iqnn))) %>%
  group_by(data_name,obs,nn_size,type) %>%
  summarize(avg_cv_mse=mean(cv_mse,na.rm=TRUE)) %>%
  as.data.frame() %>%
  arrange(data_name, type) %>%
  group_by(data_name) %>%
  na.omit() %>%
  mutate(mse_ratio = avg_cv_mse/avg_cv_mse[2],
         mse_diff = avg_cv_mse[2]-avg_cv_mse,
         rmse = sqrt(avg_cv_mse),
         rmse_ratio = rmse/rmse[2],
         rmse_diff = rmse-rmse[2]) %>% #!# only works due to knn brute being 2nd in order - danger of hard coding (don't know simple alternative)
  ungroup() %>%
  mutate(data_name = factor(data_name, levels=all_reg_sets[order(all_reg_sizes)])) %>%
  as.data.frame()

head(results_reg)
tail(results_reg)
str(results_reg)
# save(mse_results_list, results_reg,tuned_reg_param_list,tuned_reg_performance_list, file="tuned_regression_testing.Rdata")

#--------------------------------------------------------------------------------------------
### Create Regression Accuracy Comparison Plots for Figure in Section 3.2.2
# load(file="tuned_regression_testing.Rdata")


#Setting order of data sets
iqnn_reg_rmse <- filter(results_reg, type=="iqnn") %>% arrange(desc(rmse))
iqnn_rmse_order <- as.character(iqnn_reg_rmse$data_name)

# clean data for plots of RMSE relative to knn as baseline
#   - separate knn from others for "baseline" layers vs "points" layers in plotting
reg_plot_data_knn <- results_reg %>%
  filter(type == "knn - brute") %>%
  mutate(data_name = factor(data_name,levels=iqnn_rmse_order),
         x=as.numeric(data_name)-.5,
         xend=as.numeric(data_name)+.5)
reg_plot_data_knn$type_pretty <- "KNN-brute"

shiftval=.25
reg_plot_data <- results_reg %>%
  filter(type != "knn - brute") %>%
  mutate(shift = (as.numeric(as.factor(as.character(type)))-2)*shiftval,
         data_name = factor(data_name,levels=iqnn_rmse_order))
reg_plot_data$type_pretty <- factor(reg_plot_data$type, labels=c("IQNN   ","AKNN-cover ", "AKNN-kd   ")) 

#!# data labels match CLASSIFIERS!
x_label_sizes_reg <- as.numeric(sapply(levels(reg_plot_data$data_name), function(x) reg_plot_data$obs[reg_plot_data$data_name==x][1]))
my_x_ticks_reg <- paste0(levels(reg_plot_data$data_name), "\n n=",x_label_sizes_reg)


# plot relative to KNN on raw RMSE scale
#   - use vertical lines to partition between datasets
#   - baseline KNN-brute as flat line segment
#   - consistently offset points for IWNN and AKNN on X, y=accuracy
#   - scale / theme to match previous plots
p3 <- ggplot()+
  geom_vline(xintercept=seq(0.5,10.5,by=1),linetype=2,color="gray25")+
  geom_segment(aes(x=x, xend=xend,y=rmse,yend=rmse, linetype="KNN-brute"),
               color=RColorBrewer::brewer.pal(4,"Set1")[4], size=.8, data=reg_plot_data_knn) +
  geom_point(aes(x=as.numeric(data_name)+shift,y=rmse,color=type_pretty,shape=type_pretty),
             size=3,data=reg_plot_data)+
  theme_bw() +
  scale_linetype_manual("Baseline: ", values=1)+
  scale_x_continuous("", breaks=1:10,
                     labels=my_x_ticks_reg,
                     limits=c(0.5,10.5)) +
  labs(y="Average RMSE\n (units of std dev)")+
  scale_color_manual("Model Type: ", values=RColorBrewer::brewer.pal(4,"Set1")[c(1,2,3)])+
  scale_shape_manual("Model Type: ", values=c(2,0,1))+
  theme(panel.grid.major.x =element_blank(),
        panel.grid.minor.x =element_blank(),
        axis.ticks.x = element_blank(),
        # panel.border = element_blank(),
        legend.position = "none")


# plot relative to KNN as difference in RMSE from KNN-brute
#   - use vertical lines to partition between datasets
#   - baseline at 0 to represent KNN-brute
#   - consistently offset points for IWNN and AKNN on X, y=accuracy
#   - scale / theme to match previous plots
p4 <- ggplot()+
  geom_segment(aes(x=.5, xend=10.5,y=0,yend=0, linetype="KNN-brute"),size=.8, color=RColorBrewer::brewer.pal(4,"Set1")[4])+
  geom_vline(xintercept=seq(0.5,10.5,by=1),linetype=2,color="gray25")+
  geom_point(aes(x=as.numeric(data_name)+shift,y=rmse_diff,color=type_pretty,shape=type_pretty),
             size=3,data=reg_plot_data)+
  theme_bw() +
  scale_x_continuous(" ", breaks=1:10,
                     labels=my_x_ticks_reg,
                     limits=c(0.5,10.5)) +
  labs(y="Difference from KNN-brute\n (units of std dev)")+
  # scale_y_continuous("Relative CV-RMSE \n (Difference from KNN-brute)", breaks=seq(-.01,.04,by=.01),limits=c(-0.01,0.04)) +
  scale_color_manual("Model Type: ", values=RColorBrewer::brewer.pal(4,"Set1")[c(1,2,3)])+
  scale_shape_manual("Model Type: ", values=c(2,0,1))+
  scale_linetype_manual("Baseline: ", values=1)+
  theme(panel.grid.major.x =element_blank(),
        panel.grid.minor.x =element_blank(),
        axis.ticks.x = element_blank(),
        # panel.border = element_blank(),
        legend.position = "bottom") +
  annotate(geom="text",x=0,y=0,label="bold(KNN-brute)", color=RColorBrewer::brewer.pal(4,"Set1")[4], parse=T,vjust=-.2,hjust=0.4)

grid.arrange(p3,p4,nrow=2,heights=c(1,1.2))


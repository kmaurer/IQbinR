### iterative quantile binning functions

## Function List
# quant_bin_1d for quantile binning in one dimension
# make_stack_matrix helper funtion for iterative binner to create stacking J matrices


# Note: a #!# tag will be added on items that need improved efficiency/clarity

#--------------------------------------
## Quantile 1d Binning
# used for binning the counts values by quantile
# define vector of counts and number of bin
# xs <- diamonds$price; nbin=4
quant_bin_1d <- function(xs, nbin, output="data",jit=0){
  if(jit > 0) xs <- xs + runif(length(xs),-jit,jit)
  quants <- quantile(xs, seq(0, 1, by=1/(2*nbin)))
  bin_centers <- quants[seq(2,length(quants)-1, by=2)]
  bin_bounds <- quants[seq(1,length(quants)+1, by=2)]
  data_bins <- rep(bin_centers[1],length(xs))
  if(output=="definition") {
    return(list(bin_centers=bin_centers,bin_bounds=bin_bounds))
  } else{
    for (i in 2:length(bin_centers)){
      data_bins[bin_bounds[i] < xs] <- bin_centers[i]
    } 
    if(output=="data") return(data_bins)
    if(output=="both") return(list(data_bins=data_bins,bin_centers=bin_centers,bin_bounds=bin_bounds))
  }
}
quant_bin_1d(ggplot2::diamonds$price,4,output="data")
quant_bin_1d(ggplot2::diamonds$price,4,output="definition")
quant_bin_1d(runif(1000,0,10),nbin=4,output="both")

#--------------------------------------
## function to make a matrix for duplicating the N rows of a matrix M times each
# support funciton for IQ binning function
#!# Update later for speed 
make_stack_matrix <- function(N,M){ 
  mat <- unname(model.matrix(~as.factor(rep(1:N,each=M))-(1)))
  attributes(mat)[2:3]<-NULL
  return(mat)
} 
make_stack_matrix(3,4)

### Iterative Quantile Binning
# Input:
#   dat = data frame to be binned (will coerce matrix or tibble to simple data frame)
#   bin_cols = vector of column names of variables to iteratively bin, ordered first to last
#   nbins = vector of number of bins per step of iterative binning, ordered first to last
#   jit = margin for uniform jitter to each dimension to create seperability of tied obs due to finite precision
iterative_quant_bin <- function(dat, bin_cols, nbins, output="data",jit=0){
  dat <- as.data.frame(dat)
  bin_dim <- length(bin_cols)
  bin_dat <- matrix(NA,nrow=nrow(dat),ncol=bin_dim, dimnames=list(row.names(dat),paste(bin_cols,"binned",sep="_")))
  # Initialize with first binning step
  step_bin_info <- quant_bin_1d(dat[,bin_cols[1]], nbins[1],output="both",jit)
  bin_bounds <- matrix(c(step_bin_info$bin_bounds[1:nbins[1]],
                         step_bin_info$bin_bounds[2:(nbins[1]+1)]),
                       nrow=nbins[1],byrow=FALSE )
  bin_centers <- matrix(step_bin_info$bin_centers, nrow=nbins[1])
  bin_dat[,1] <- step_bin_info$data_bins
  # Loop over remaining variables to use quantile binning WITHIN each of previous state bins
  for(d in 2:bin_dim){
    stack_size <- nrow(bin_centers)
    stack_matrix <- make_stack_matrix(stack_size,nbins[d])
    bin_centers <- cbind(stack_matrix %*% bin_centers,matrix(rep(NA,stack_size*nbins[d]),ncol=1))
    bin_bounds <- cbind(stack_matrix %*% bin_bounds,matrix(rep(NA,2*stack_size*nbins[d]),ncol=2))
    # iterate through unique bins from prior step which are the {1,1+nbins[d],1+2*nbins[d],...} rows of the bin matrices
    for(b in seq(1,1+(stack_size-1)*nbins[d],by=nbins[d]) ){
      in_bin_b <- apply(matrix(bin_dat[,1:(d-1)],ncol=(d-1)),1,identical,y=bin_centers[b,-d])
      step_bin_info <- quant_bin_1d(dat[in_bin_b,bin_cols[d]], nbins[d],output="both",jit)
      bin_bounds[b:(b+nbins[d]-1),c(2*d-1,2*d)] <- matrix(c(step_bin_info$bin_bounds[1:nbins[d]],
                                                            step_bin_info$bin_bounds[2:(nbins[d]+1)]),
                                                          nrow=nbins[d],byrow=FALSE)
      bin_centers[b:(b+nbins[d]-1),d] <- matrix(step_bin_info$bin_centers, nrow=nbins[d])
      bin_dat[in_bin_b,d] <- step_bin_info$data_bins
    }
  }
  if(output=="data") return(list(dat=dat,bin_dat=bin_dat))
  if(output=="definition") return(list(bin_centers=bin_centers, bin_bounds=bin_bounds,bin_cols=bin_cols, nbins=nbins))
  if(output=="both") return(list(bin_dat=list(dat=dat,bin_dat=bin_dat), bin_def=list(bin_centers=bin_centers, bin_bounds=bin_bounds,bin_cols=bin_cols, nbins=nbins)))
}
iterative_quant_bin(dat=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), nbins=c(3,5,2), output="both",jit=0.001)



### Iterative Quantile Binned Nearest Neighbors Regression
# takes in data, response column and binning parameters
iqnn <- function(dat, y, bin_cols, nbins, jit=0){
  ## make bins
  iq_bin<- iterative_quant_bin(dat, bin_cols, nbins, output="both",jit)
  ## For each bin, find indeces from original data where bins match, take average y value
  total_bins = nrow(iq_bin$bin_def$bin_centers)
  iq_bin$bin_def$bin_stats <- data.frame(avg = rep(NA,total_bins),
                                         obs = NA)
  for(b in 1:total_bins){
    match_matrix <- iq_bin$bin_dat$bin_dat == matrix(rep(iq_bin$bin_def$bin_centers[b,],nrow(dat)),ncol=3,byrow=TRUE)
    temp_data <- dat[,y][match_matrix[,length(bin_cols)]]
    iq_bin$bin_def$bin_stats[b,] <- c(mean(temp_data,na.rm=TRUE),length(temp_data))
  }
  ## Return bin definition with predictions added
  return(iq_bin$bin_def)
}
iqnn(iris, y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), nbins=c(3,5,2))
myiq <- iqnn(iris, y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), nbins=c(3,5,2), jit=.001)
myiq




### Helper function for checking if a vector is in a p-dimensional bin, defined by 2*p boundaries
# x = p-dimensional vector
# bin_bounds = 2*p dimensional boundary matrix (like in iq-binning definition list)
bin_index_finder <- function(x, bin_bounds){ 
  p = length(x)
  xrep_mat = matrix(rep(x,nrow(bin_bounds)),ncol=3,byrow=TRUE)
  which(rowSums(bin_bounds[,seq(1,2*p-1,by=2)] < xrep_mat & xrep_mat <= bin_bounds[,seq(2,2*p,by=2)])==p)
} 
myiq <- iqnn(iris, y="Petal.Length", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), nbins=c(3,5,2), jit=.001)
new_row <- iris[1,c("Sepal.Length","Sepal.Width","Petal.Width")]
new_row_index <- bin_index_finder(new_row, myiq$bin_bounds)
new_row
myiq$bin_bounds[new_row_index,]
myiq$bin_centers[new_row_index,]
myiq$bin_stats[new_row_index,]


### Iterative Quantile Binning New Data from defined bins
# iq_def= IQ bin definition list from iterative_quant_bin or iqnn
# new_data = data frame with column names matching the binned columns from bin-training data
# output matches format of iterative_quant_bin and inherets properties from iqnn if applicable
iq_def <- iterative_quant_bin(dat=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), nbins=c(3,5,2), output="definition",jit=0.001)
str(iqdef)
bin_by_IQdef <- finction(iq_def, new_data){
  total_bins = nrow(iq_def$bin_centers)
  total_cols = length(iq_def$bin_cols)
  bin_dat <- list(dat=as.data.frame(new_data), bin_dat=matrix(NA,nrow=nrow(new_data),ncol=total_cols))
  # loop over each obs in new data, identify the bin and its attributes to add to the bin_data 
  lapply(1:nrow(new_data), function(i){
    bin_index_finder(bin_dat$dat[i,iq_def$bin_cols],iq_def$bin_bounds)
  })
  #!# identify bins for observations outside of "true" ranges (problem with jitter shrinking true range of data when shift inward)
  
  for(i in 1:nrow(new_data)){
    dat_i <- bin_dat$dat[i,iq_def$bin_cols]
    
    in_bin_int(dat_i[1], )
    
    }
  }


                  
  if(!is.null(iq_def$bin_stats)){
    #!# work to add predictions goes here
  }           

  
} 



### Cross Validation function for assessing 


#-----------------------------------------------------------------------------------------
### Break from function writing to applications testing
# Test on iris and diamonds data
iterative_quant_bin(iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), nbins=c(3,5,2), output="both")
library(ggplot2)
iterative_quant_bin(ggplot2::diamonds, bin_cols=c("price","carat"), nbins=c(5,8), output="definition")

# try out plotting with them
dat=iris
set.seed(123)
dat$Sepal.Length <- dat$Sepal.Length + runif(nrow(dat), -.001,.001)
dat$Sepal.Width <- dat$Sepal.Width + runif(nrow(dat), -.001,.001)
bin_cols=c("Sepal.Length","Sepal.Width")
mybins <- iterative_quant_bin(dat, bin_cols, nbins=c(6,5), output="both")
library(tidyverse)
bin_aggs <- mybins$bin_dat %>%
  dplyr::group_by(Sepal.Length_binned,Sepal.Width_binned) %>%
  dplyr::summarize(binAvg=mean(Petal.Length),
                   binSd=sd(Petal.Length),
                   binCount=n())
bin_aggs_bounds <-cbind(as.data.frame(mybins$bin_bounds),as.data.frame(bin_aggs))

# Plot for ARA application
ggplot() +
  geom_rect(aes(xmin=V1,xmax=V2,ymin=V3,ymax=V4,fill=binAvg),color="black",size=.75,data=bin_aggs_bounds)+
  # geom_point(aes_string(x=bin_cols[1],y=bin_cols[2]),color="red",alpha=.1,data=mybins$bin_dat)+
  theme_bw() +
  xlab("Sepal Length (cm)") + ylab("Sepal Width (cm)") +
  scale_fill_gradient("Average Petal \n Length (cm)", low="#CAC2E3", high="#3A107D") +
  ggtitle("Iteratively Quantile Binned Heatmap of Average \n Petal Length within Sepal Size Regions")+
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks=seq(4.5,8, by=.5))


ggplot() +
  geom_rect(aes(xmin=V1,xmax=V2,ymin=V3,ymax=V4,fill=binSd),color="black",data=bin_aggs_bounds)+
  geom_point(aes_string(x=bin_cols[1],y=bin_cols[2]),color="red",data=mybins$bin_dat)+
  theme_bw()

p1 <- ggplot() +
  geom_rect(aes(xmin=V1,xmax=V2,ymin=V3,ymax=V4,fill=binAvg),color="black",size=.75,data=bin_aggs_bounds)+
  # geom_point(aes_string(x=bin_cols[1],y=bin_cols[2]),color="red",alpha=.1,data=mybins$bin_dat)+
  theme_bw() +
  xlab("Sepal Length (cm)") + ylab("Sepal Width (cm)") +
  scale_fill_gradient("Average Petal \n Length (cm)", low="#CAC2E3", high="#3A107D") +
  ggtitle("Iteratively Quantile Binned Heatmap of Average \n Petal Length within Sepal Size Regions")

p2 <- ggplot() +
  geom_jitter(aes_string(x=bin_cols[1],y=bin_cols[2],color="Petal.Length"),size=4,data=mybins$bin_dat)+
  theme_bw()

iqnn(dat,y="Petal.Length", bin_cols, nbins=c(6,5))

library(gridExtra)
grid.arrange(p2,p1,nrow=1)

library(devtools)
install_github("kmaurer/BinPackage")
help(package="BinPackage")
library(BinPackage)

xs <- ggplot2::diamonds$price; nbins=4; origin=min(xs); width=diff(range(xs))/nbins
rect_bin_1d(rnorm(1000,0,1),origin=-4,width=1)

rectbindata <- data.frame(x=rect_bin_1d(dat[,bin_cols[1]],4.3,.5),
           y=rect_bin_1d(dat[,bin_cols[2]],4.3,.5),
           z=dat$Petal.Length) %>%
  group_by(x,y) %>%
  dplyr::summarize(avgZ = mean(z))
p3 <- ggplot()+
  geom_tile(aes(x=x,y=y,fill=avgZ),data=rectbindata)+
  theme_bw()

library(gridExtra)
grid.arrange(p2,p1,p3,nrow=1)


dat=diamonds
bin_cols=c("carat","depth")
mybins <- iterative_quant_bin(dat, bin_cols, nbins=c(11,10), output="both")
bin_aggs <- mybins$bin_dat %>%
  group_by(carat_binned,depth_binned) %>%
  dplyr::summarize(binAvg=mean(price),
                   binCount=n())
bin_aggs_bounds <-cbind(as.data.frame(mybins$bin_bounds),bin_aggs)
head(bin_aggs_bounds)
ggplot() +
  geom_rect(aes(xmin=V1,xmax=V2,ymin=V3,ymax=V4,fill=binAvg),color="black",data=bin_aggs_bounds)+
  # geom_point(aes_string(x=bin_cols[1],y=bin_cols[2]),color="red",alpha=.05,data=mybins$bin_dat)+
  theme_bw()


### Load/Prep Data for Experiment
dataset = "phoneme"
trainSize = 5000
load(file=paste("data/",dataset,".Rdata",sep=""))

set.seed(0224)
dat1 <- sample_n(na.omit(.),size=trainSize)
dat1 <- dat1[order(as.numeric(row.names(dat1))),]
names(dat1) <- str_replace_all(names(dat1),"-","")
names(dat1) <- str_replace_all(names(dat1)," ","")
remove(.)
dat=dat1
bin_cols=c("Aa","Ao")
mybins <- iterative_quant_bin(dat, bin_cols, nbins=c(4,4), output="both")
bin_aggs <- mybins$bin_dat %>%
  group_by(Aa_binned,Ao_binned) %>%
  dplyr::summarize(binAvg=mean(Iy),
                   binCount=n())
bin_aggs_bounds <-cbind(as.data.frame(mybins$bin_bounds),bin_aggs)
head(bin_aggs_bounds)
p1 <- ggplot() +
  geom_rect(aes(xmin=V1,xmax=V2,ymin=V3,ymax=V4,fill=binCount),color="black",size=.75,data=bin_aggs_bounds)+
  # geom_point(aes_string(x=bin_cols[1],y=bin_cols[2]),color="red",alpha=.1,data=mybins$bin_dat)+
  theme_bw()

p2 <- ggplot() +
  geom_point(aes_string(x=bin_cols[1],y=bin_cols[2],color="Iy"),data=mybins$bin_dat)+
  theme_bw()

library(gridExtra)
grid.arrange(p2,p1,nrow=1)

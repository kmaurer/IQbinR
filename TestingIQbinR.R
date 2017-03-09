#Testing Iterative Quantile Binning Functions for Applications



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

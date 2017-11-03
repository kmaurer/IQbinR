
###-----------------------------------------------------------------------------------------------------
# Organize Classification Data sets
setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\classification")

## Small classifier set attributes
small_sets <- c("iris","pima","yeast")
web_link <- c("iris/iris.data","pima-indians-diabetes/pima-indians-diabetes.data","yeast/yeast.data")
small_responses <- c("V5","V9","V10")
small_sizes <- c(150,768,1484)

# # Wisconsin Prognostic Breast Cancer (wpbc) data
# data <- as.data.frame(fread("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wpbc.data"))
# data <- select(data, -V1,-V3)
# save(data, file="wpbc_raw.Rdata")

# # Loop over UCI data sets to clean and save to CSV
# for(set in 1:3){
#   data <- fread(paste0("http://archive.ics.uci.edu/ml/machine-learning-databases/",web_link[set]))
#   save(data, file=paste0(small_sets[set],"_raw.Rdata"))
# }


## "mediumDatasets" attributes
medium_sets <- c("abalone","waveform","optdigits","satimage","marketing")
medium_responses <- c("Sex","Class","Class","Class","Sex")
medium_sizes <- c(4174,5000,5620,6435,6876)

# Loop over all "medium" data sets from walter to clean and save to CSV
# for(set in 1:6){
#   setwd("C:\\Users\\maurerkt\\Google Drive\\AFRLSFFP\\Fall2017\\mediumDatasets")
#   ## load and clean data in preparation for testing speed/accuracy with k-fold CV process
#   data <- RWeka::read.arff(paste0(medium_sets[set],".arff"))
#   save(data, file=paste0(medium_sets[set],"_raw.Rdata"))
# }

### Epileptic Seizure Recognition data n=11500, y=5-group, p=10
data <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/00388/data.csv")
head(data)
save(data, file="seizure_raw.Rdata")

### Magic data n=19020, y=binary, p=10
# http://sci2s.ugr.es/keel/dataset.php?cod=102
# data <- as.data.frame(fread("magic.dat"))
# head(data)
# names(data)[which(names(data)=="V11")] <- "y"
# save(data, file="magic_raw.Rdata")


## Large classifier set attributes
large_sets <- c("youtube", "skin")
large_responses <- c("category", "V4")
large_sizes <- c(168286,245057)

# # download zipfolder containing youtube_videos.tsv from https://archive.ics.uci.edu/ml/datasets/Online+Video+Characteristics+and+Transcoding+Time+Dataset
# data <- fread("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\classification\\youtube_videos.tsv")
# save(data, file="youtube_raw.Rdata")
# # read skin segmentation data directly from web
# data <- fread("https://archive.ics.uci.edu/ml/machine-learning-databases/00229/Skin_NonSkin.txt")
# save(data, file="skin_raw.Rdata")


###-----------------------------------------------------------------------------------------------------
# Organize Regression Data sets
###-----------------------------------------------------------------------------------------------------
## UCI repo extensions (http://archive.ics.uci.edu/ml/datasets/-name here-)
# PM2.5+Data+of+Five+Chinese+Cities
# Physicochemical+Properties+of+Protein+Tertiary+Structure
# Air+Quality
# Combined+Cycle+Power+Plant
# SkillCraft1+Master+Table+Dataset
# Breast+Cancer+Wisconsin+%28Prognostic%29

setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\iqnnProject\\DataRepo\\regression")

#--------------------------
## FROM UCI Repo
#--------------------------
# # Wisconsin Prognostic Breast Cancer (wpbc) data (UCI) n=198
# data <- as.data.frame(fread("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wpbc.data"))
# data <- select(data, -V1,-V2)
# names(data)[which(names(data)=="V3")] <- "y"
# save(data, file="wpbc_raw.Rdata")
# 
# # SkillCraft1 data n=3395
# data <- as.data.frame(fread("http://archive.ics.uci.edu/ml/machine-learning-databases/00272/SkillCraft1_Dataset.csv"))
# head(data)
# data <- select(data, -GameID)
# names(data)[which(names(data)=="APM")] <- "y"
# save(data, file="skillcraft_raw.Rdata")
# 
# # Combined+Cycle+Power+Plant n=9568
# # http://archive.ics.uci.edu/ml/machine-learning-databases/00294/
# # Download CCPP.zip then save Folds5X2_pp.xlsx to csv, load from csv
# data <- as.data.frame(fread("ccpp.csv"))
# head(data)
# names(data)[which(names(data)=="PE")] <- "y"
# save(data, file="ccpp_raw.Rdata")
# 
# # Air Quality n=9471
# # http://archive.ics.uci.edu/ml/machine-learning-databases/00360/
# # Download AirQualityUCI.zip then load from AirQualityUCI.csv
# data <- as.data.frame(fread("AirQualityUCI.csv", dec=","))
# head(data)
# data <- select(data, -Time,-Date)
# names(data)[which(names(data)=="CO(GT)")] <- "y"
# save(data, file="air_quality_raw.Rdata")

#--------------------------
## FROM KEEL Repo 
# *note: all sets downloaded in zip files from http://sci2s.ugr.es/keel/category.php?cat=reg&order=ins#sub2 on 10/27/17
#--------------------------

# # Weather ankara (wankara) (KEEL) n=321
# data <- as.data.frame(fread("wankara.dat"))
# head(data)
# names(data)[which(names(data)=="V10")] <- "y"
# save(data, file="wankara_raw.Rdata")
# 
# # Laser data (KEEL) n=993
# data <- as.data.frame(fread("laser.dat"))
# head(data)
# names(data)[which(names(data)=="V5")] <- "y"
# save(data, file="laser_raw.Rdata")
#
# # Treasury (KEEL) n=1049
# data <- read.csv("treasury.csv", header=FALSE)
# head(data)
# names(data)[which(names(data)=="V16")] <- "y"
# save(data, file="treasury_raw.Rdata")
# 
# # Quake data (KEEL) n=2178
# data <- as.data.frame(fread("quake.dat"))
# head(data)
# names(data)[which(names(data)=="V4")] <- "y"
# save(data, file="quake_raw.Rdata")
# 
# # ANACALT (KEEL) n=4052
# data <- as.data.frame(fread("raw/ANACALT.dat"))
# head(data)
# names(data)[which(names(data)=="V8")] <- "y"
# save(data, file="anacalt_raw.Rdata")
# 
# # Pumadyn data n=8192
# data <- as.data.frame(read.csv("puma32h.csv", header=FALSE))
# head(data)
# names(data)[which(names(data)=="V33")] <- "y"
# save(data, file="puma_raw.Rdata")

all_reg_sets <- c("wpbc","wankara","laser","treasury","quake","skillcraft","anacalt","puma","air_quality","ccpp")
all_reg_sizes <- c(198,321,993,1049,2178,3395,4052,8192,9471,9568)

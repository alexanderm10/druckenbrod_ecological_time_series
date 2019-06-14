

setwd("/Users/tessamandra/Desktop")

#install.packages("TRADER")
#install.packages("dplR")
#install.packages("gdata")
#install.packages("R.matlab")
#install.packages("MASS")
#install.packages("robustbase")
#install.packages("minpack.lm")
#install.packages("stringr")

library("TRADER")
library("dplR")
library("gdata")
library("R.matlab")
library("MASS")
library("robustbase")
library("minpack.lm")
library("stringr")

# Load sub-functions
source("iGrowth.R")
source("bisqmean_dan.R")
source("ar_order.R")
source("outlier_clt.R")
source("backcast.R")

lct <- read.tucson("LCT.txt")

plotFirstYears(lct)
plotGrowth(lct)

growthAveragingALL(lct,
                   releases = NULL, 
                   m1 = 15, m2 = 15, 
                   buffer = 15, drawing = TRUE, 
                   criteria = 0.50, 
                   prefix = "ga", gfun = mean, 
                   length = 1, storedev = pdf)


### Find number of major releases ###

majorR <- read.csv(file = "ga_releases_Only_Major.csv", header = TRUE)
major.sum <- colSums(majorR[,!names(majorR) %in% "year"], na.rm = TRUE, dims = 1)
major.sum <- data.frame(major.sum)

### Find number of moderate releases ###

minorR <- read.csv(file = "ga_releases_Only_Moderate.csv", header = TRUE)
minor.sum <- colSums(minorR[,!names(minorR) %in% "year"], na.rm = TRUE, dims = 1)
minor.sum <- data.frame(minor.sum)

### Find dates of first and last release for each tree ###
### Make data frame with tree ID for each row and year of first and 
### last release for two of the columns. 

majorD <- data.frame(TreeID = names(majorR[,!names(majorR) %in% "year"]))

### dropping the year from majorR so that the columns and rows align
### adding the years as row names so that they can be identified in loop

majorR2 <- majorR[,!names(majorR) %in% "year"]
row.names(majorR2) <- majorR$year

for(i in 1:ncol(majorR2)){
  yr.first <- min(rownames(majorR2)[which(majorR2[,i]==1)])
  yr.last <- max(rownames(majorR2)[which(majorR2[,i]==1)])
  majorD[i, "First"] <- yr.first
  majorD[i, "Last"] <- yr.last
}

### Do the same for minorD

minorD <- data.frame(TreeID = names(minorR[,!names(minorR) %in% "year"]))

### dropping the year from minorR so that the columns and rows align
### adding the years as row names so that they can be identified in loop

minorR2 <- minorR[,!names(minorR) %in% "year"]
row.names(minorR2) <- minorR$year

for(i in 1:ncol(minorR2)){
  yr.first <- min(rownames(minorR2)[which(minorR2[,i]==1)])
  yr.last <- max(rownames(minorR2)[which(minorR2[,i]==1)])
  minorD[i, "First"] <- yr.first
  minorD[i, "Last"] <- yr.last
}

### Find first and last release date from both major and moderate releases
### Merge majorD and minorD
majminD <- merge(majorD, minorD, by = "TreeID")

### Create final dataframe for first and last release date

releaseD <- data.frame(TreeID = names(minorR[,!names(minorR) %in% "year"]))

#####################################################################
### Find min and max for first and last release
### If first or last release = infinity, then set it equal to the 
### age of the tree core (time since last disturbance)
### mhtdata is used to get the ages of the trees
#####################################################################

lctdata <- read.csv(file = "lctdata.csv", header = TRUE)

for(i in 1:nrow(majminD)){
  yr.first <- min(as.numeric(majminD[i, c("First.x", "First.y")]), na.rm = TRUE)
  yr.last <- max(as.numeric(majminD[i, c("Last.x", "Last.y")]), na.rm = TRUE)
  releaseD[i, "First"] <- yr.first
  releaseD[i, "Last"] <- yr.last
  if(releaseD[i, "First"] == "Inf" | releaseD[i, "First"] == "-Inf"){
    releaseD[i, "First"] <- 2017 - lctdata[i, "Age"]
  }
  if(releaseD[i, "Last"] == "Inf" | releaseD[i, "Last"] == "-Inf"){
    releaseD[i, "Last"] <- 2017 - lctdata[i, "Age"]
  }
}

#####################################################################
### Add column to releaseD for average radial growth rate for each tree
#####################################################################

gRate <- colMeans(lct, na.rm = TRUE, dims = 1)
gRate <- data.frame(gRate)
names(gRate) <- "avgRing"

#####################################################################
### Merge all above columns into one data frame
#####################################################################

final.data <- data.frame(TreeID = row.names(major.sum), 
                         majorRls = major.sum$major.sum, 
                         minorRls = minor.sum$minor.sum,
                         firstRelease = releaseD$First,
                         lastRelease = releaseD$Last, 
                         avgRing = gRate$avgRing)

write.csv(final.data, file = "release.csv", row.names = F)

#####################################################################
### Calculate column for years since last release
#####################################################################

for(i in 1:nrow(final.data)){
  rYears <- 2017 - final.data[i, "lastRelease"]
  final.data[i, "ySinceR"] <- rYears
}

#####################################################################
### Calculate column for total number of releases
#####################################################################

for(i in 1:nrow(final.data)){
  nRls <- final.data[i, "majorRls"] + final.data[i, "minorRls"]
  final.data[i, "allRls"] <- nRls
}

#####################################################################
### Pull in WRW spreadsheet as csv and append final.data to the end
#####################################################################

lctdata <- read.csv(file = "lctdata.csv", header = TRUE)
lctdf <- merge(lctdata, final.data)
write.csv(lctdf, file = "lctdf.csv", row.names = FALSE)

#####################################################################
### Run function to determine if each tree established in a gap or 
### no gap. 1 = gap. 2 = maybe gap. 3 = no gap. Add these numbers to 
### a new column in the csv. 
#####################################################################

direct <-setwd("/Users/tessamandra/Desktop")
fileN <- "LCT.txt"

for(i in names(lct)){
  trend <- v105pn(direct, fileN, i, fig=0, iter=8)
  ### Make new column in mhtdf, dump trend type into it for each core
  lctdf[lctdf$TreeID==i, "Gap"] <- as.numeric(str_sub(trend,-1,-1))
  ### str_sub takes off the "Trend Type: " from the output
}

write.csv(lctdf, file = "lctdf.csv", row.names = FALSE)

#####################################################################
### Add basal area increment column to df 
### First subset dataframe from df that has 2 columns: years and DBH
#####################################################################
#baisub <- df[, c("TreeID", "DBH")]
#baidf <- bai.out(rwl = wrw, diam = baisub)
#test <- wrw[order(row.names(wrw), decreasing = F),]
#ba.test <- bai.out(test, baisub)

#####################################################################
### Run PCA
#####################################################################

row.names(lctdf) <- lctdf$TreeID
test <- prcomp(lctdf[,c("Age", "DBH", "CanP", "LA", "Bam", "cIndex", "avgRing", "ySinceR", "allRls" )], scale=F)
test <- prcomp(lctdf[,c("DBH", "CanP", "cIndex", "avgRing", "ySinceR", "allRls" )], scale=F)

biplot(test, xlab="PC1 - 69%", ylab="PC2 - 27%")

summary(test)


## Calculate PCA and RDA for all the OTU tables

library(vegan)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library("ggvegan")
library(packfor)
library(geosphere)

### Load in data ###
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")

## new data "../results/Jericho_env_data_mean_imputed_transformed.csv"

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

## from 001 script 
Jericho_data <- read.csv("../results/Jericho_env_data_mean_imputed_transformed.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)


## amplicon richness calculated in earlier script.
amplicon_richness <- read.csv("../results/amplicon_richness_by_date.csv")
amplicon_richness$Date <- as.Date(amplicon_richness$Date)


### Subset data ###
## subset data that is interesting to me:
keep <- c("Date",
          "Average_viral_abundance",
          "Average_bacterial_abundance",
          "Average_chl_a",
          "Average_PO4",
          "Average_SiO2",
          "Average_NO3_NO2",
          "Temperature_YSI",
          "Salinity_ppt_YSI",
          "Dissolved_oxygen_percent",  
          "pH")
## might be useful to do some more cleaning before
## keep only certain rows. 
params_Jericho <- droplevels(Jericho_data[,keep])

params_Jericho$pH <- as.character(params_Jericho$pH)
params_Jericho$pH[13]  <- "6.5"
params_Jericho$pH <- as.numeric(params_Jericho$pH)
VC_dates <- Jericho_data$Date[match(rownames(normalized_18s_OTUs), Jericho_data$VC_number)]

## add in Daylength uisng Vancouver's latidue and the year month date
params_Jericho$day_length <- daylength(49.2827, params_Jericho$Date)
params_Jericho$month_number <- as.numeric(format(params_Jericho$Date, "%m"))
#params_Jericho_without_date <- params_Jericho[,-1]

### Do mean imputation ###
## Revisit!
summary(params_Jericho)
params_Jericho$Average_viral_abundance[is.na(params_Jericho$Average_viral_abundance)] = mean(params_Jericho$Average_viral_abundance, na.rm=TRUE)
params_Jericho$Average_bacterial_abundance[is.na(params_Jericho$Average_bacterial_abundance)] = mean(params_Jericho$Average_bacterial_abundance, na.rm=TRUE)
## median would be better here...
params_Jericho$Average_chl_a[is.na(params_Jericho$Average_chl_a)] = mean(params_Jericho$Average_chl_a, na.rm=TRUE)
params_Jericho$Average_PO4[is.na(params_Jericho$Average_PO4)] = mean(params_Jericho$Average_PO4, na.rm=TRUE)
params_Jericho$Average_SiO2[is.na(params_Jericho$Average_SiO2)] = mean(params_Jericho$Average_SiO2, na.rm=TRUE)
params_Jericho$Average_NO3_NO2[is.na(params_Jericho$Average_NO3_NO2)] = mean(params_Jericho$Average_NO3_NO2, na.rm=TRUE)
params_Jericho$Temperature_YSI[is.na(params_Jericho$Temperature_YSI)] = mean(params_Jericho$Temperature_YSI, na.rm=TRUE)
params_Jericho$Salinity_ppt_YSI[is.na(params_Jericho$Salinity_ppt_YSI)] = mean(params_Jericho$Salinity_ppt_YSI, na.rm=TRUE)
params_Jericho$Dissolved_oxygen_percent[is.na(params_Jericho$Dissolved_oxygen_percent)] = mean(params_Jericho$Dissolved_oxygen_percent, na.rm=TRUE)
params_Jericho$pH[is.na(params_Jericho$pH)] = mean(params_Jericho$pH, na.rm=TRUE)
summary(params_Jericho)


params_Jericho <- subset(params_Jericho, Date %in% VC_dates )
amplicon_richness_for_study <- subset(amplicon_richness, Date %in% VC_dates )
## temporary fix
amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0

## PCA on full dataset (correlation matrix: scale=TRUE)

## try it without po4 since it is highly correlated with no2+no3
env.pca <- rda(params_Jericho[,c(-1,-4)], scale=TRUE) # scale=TRUE calls for a stadardization of the variables
env.pca
head(summary(env.pca)) # Default scaling 2
## what is this scaling 2??
head(summary(env.pca, scaling = 1))


## Examine and plot partial results from PCA output

## Eigenvalues
ev <- env.pca$CA$eig

## Apply Kaiser-Guttman cirterio to select axes
## This computes the mean of all the eigenvalues and interprets only the axes larger than the mean
ev[ev > mean(ev)] ## give first 2 axes

## Broken stick model
## not sure what is happening here!
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n){
 bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))
}


bsm$p <- 100*bsm$p/n
bsm

## Plot eigenvalues and % of variance for each axis

barplot(ev, main="Eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red") ## average eigenvalue
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

## interpret only the values that are longer than the broken stick. 
## don't quite get why. 
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, main="%variance", col=c("bisque", 2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), pch=15, col=c("bisque", 2), bty="n")


#Plot the data:

# Two PCA biplots: scaling 1 and 2

biplot(env.pca, scaling=1, main="PCA- scaling 1")
biplot(env.pca, main="PCA- scaling 2") ## Default is scaling 2


### 5.3.3 PCA on Transformed species data. 
## PCA 

# Hellinger pre-transformation of the species data
spe.h <- decostand(normalized_18s_OTUs,"hellinger")
spe.h.pca <- rda(spe.h)
spe.h.pca
head(summary(spe.h.pca))

ev <- spe.h.pca$CA$eig

ev[ev > mean(ev)] ## give first 2 axes

## Broken stick model
## not sure what is happening here!
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n){
 bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))
}


bsm$p <- 100*bsm$p/n
bsm

## Plot eigenvalues and % of variance for each axis

barplot(ev, main="Eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red") ## average eigenvalue
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

## interpret only the values that are longer than the broken stick. 
## don't quite get why. 
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, main="%variance", col=c("bisque", 2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), pch=15, col=c("bisque", 2), bty="n")

biplot(spe.h.pca, scaling=1, main="PCA- scaling 1")
biplot(spe.h.pca, main="PCA- scaling 2") ## Default is scaling 2

autoplot(spe.h.pca)
## use fortify method to extract scores in ggplot-friendly format
scrs <- fortify(spe.h.pca, scaling = 3)
## take only site scores for this
sites <- with(scrs, scrs[Score == "sites", ])
species <- with(scrs, scrs[Score == "species", ])
## add in something to plot for the arrows

plt <- ggplot(sites, aes(x = Dim1, y = Dim2, label=Label)) + 
 geom_point() +
 coord_fixed()+
 geom_text(hjust=1.2)
plt

## add in arrows
plt <- ggplot(sites, aes(x = Dim1, y = Dim2, label=Label
                         #, colour = Moisture
)) + 
 geom_point() +
 coord_fixed()+
 geom_text(hjust=1.2)+
 geom_segment(data=species, aes(x=Dim1, y=Dim2, xend=0, yend=0), colour="pink")
plt

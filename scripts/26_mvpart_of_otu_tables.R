## try out mvpart

#install.packages("devtools")
#devtools::install_github("cran/mvpart")
library(mvpart)

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)


Jericho_data <- read.csv("../results/Jericho_env_data_mean_imputed.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)


parameters_to_exclude <- c("season",
                           "Average_NO3_NO2",
                           "Standard_error_NO3_NO2",
                           "Secchi_disk_disappears",
                           "Secchi_disk_reappears",
                           "Dissolved_oxygen_percent",
                           "VC_number",
                           "Tide_height",
                           "month_number",
                           "Raw_Bacteria_Rep_A",
                           "Raw_Bacteria_Rep_B",
                           "Raw_Viruses_rep_A",
                           "Raw_Viruses_rep_B",
                           "PO4_rep_A",
                           "PO4_rep_B",
                           "SiO2_rep_A",
                           "SiO2_rep_B",
                           "NO3_NO2_rep_A",
                           "NO3_NO2_rep_B",
                           "Average_viral_abundance",
                           "Standard_error_viral_abundance",
                           "Standard_error_bacterial_abundance",
                           "Standard_error_PO4",
                           "Standard_error_SiO2",
                           "Standard_error_chl_a")

params_Jericho <-  subset(Jericho_data, select= !(colnames(Jericho_data) %in% parameters_to_exclude))


mvpart(data.matrix(params_Jericho[,-1])~.,params_Jericho[,-1])
mvpart(data.matrix(params_Jericho[,c(-1,-5,-6)])~.,params_Jericho[,c(-1,-5,-6)], margin=0.08, cp=0, xv="pick", xval=nrow(params_Jericho), xvmult=100, which=4)

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
Jericho_MPL_env <- droplevels(subset(params_Jericho, VC_number %in% row.names(normalized_MPL_OTUs)))
rownames(Jericho_MPL_env) <- Jericho_MPL_env$Date
#normalized_MPL_OTUs$VC_number <- row.names(normalized_MPL_OTUs)

## just work with date
row.names(normalized_MPL_OTUs) <- Jericho_MPL_env$Date[match(rownames(normalized_MPL_OTUs),Jericho_MPL_env$VC_number)]
normalized_MPL_OTUs <- normalized_MPL_OTUs[sort(row.names(normalized_MPL_OTUs)),]
Jericho_MPL_env <- Jericho_MPL_env[sort(row.names(Jericho_MPL_env)),]


spe.ch.mvpart <- mvpart(data.matrix(normalized_MPL_OTUs) ~., Jericho_MPL_env[,c(-1, -5)])


spe.ch.mvpart <- mvpart(data.matrix(normalized_MPL_OTUs) ~., Jericho_MPL_env[,c(-1, -5)], size=6, xvmult=100)


spe.ch.mvpart <- mvpart(data.matrix(normalized_MPL_OTUs) ~., Jericho_MPL_env[,c(-1, -5)],margin=0.08, cp=0, xv="pick", xval=nrow(params_Jericho), xvmult=100, which=4)

spe.ch.mvpart <- mvpart(data.matrix(normalized_MPL_OTUs) ~., Jericho_MPL_env[,c(-1, -5)],margin=0.08, cp=0, xv="1se", xval=nrow(params_Jericho), xvmult=100, which=4)

summary(spe.ch.mvpart)
print(spe.ch.mvpart)
plot(spe.ch.mvpart)
text(spe.ch.mvpart)


# To obtain the path to the leaf nodes
leafnodeRows <- grepl("leaf",spe.ch.mvpart$frame$var)
nodevals <- as.numeric(rownames(spe.ch.mvpart$frame)[leafnodeRows])
rules <- path.rpart(spe.ch.mvpart,nodevals)
rulesdf <- do.call("rbind",lapply(rules,function(x)paste(x,collapse = " -AND- ")))
rulesdf <- data.frame(nodeNumber=rownames(rulesdf),rule=rulesdf[,1],stringsAsFactors=FALSE)


hist(residuals(spe.ch.mvpart), col='grey')
plot(predict(spe.ch.mvpart), residuals(spe.ch.mvpart), main="residuals vs. predicted")
abline(h=0, lty=3, col="grey")
## group composition
spe.ch.mvpart$where

## group identity
(groups.mrt <- levels(as.factor(spe.ch.mvpart$where)))

## composisition of first leaf
normalized_MPL_OTUs[which(spe.ch.mvpart$where==groups.mrt[1]),]

## environmental variables of the first leaf
Jericho_MPL_env[which(spe.ch.mvpart$where==groups.mrt[1]),]

leaft.sum <- matrix(0, length(groups.mrt), ncol(normalized_MPL_OTUs))
colnames(leaft.sum) <- colnames(normalized_MPL_OTUs)

for (i in 1:length(groups.mrt)){
 leaft.sum[i,] <- 
  apply(normalized_MPL_OTUs[which(spe.ch.mvpart$where==groups.mrt[i]),],2,sum)
}
leaft.sum

par(mfrow=c(2,2))
for(i in 1:length(groups.mrt)){
 pie(which(leaft.sum[i,]>0), radius=1, main=c("laef #", groups.mrt[i]))
}

## combine MRT with Indicator species 
library(labdsv)
## indicator species searc on MRT result
spe.ch.MRT.indval <- indval(normalized_MPL_OTUs, spe.ch.mvpart$where)
spe.ch.MRT.indval$pval


## For each significant species, find the leaf with the highest inval
spe.ch.MRT.indval$maxcls[which(spe.ch.MRT.indval$pval<=0.05)]

## Indval value in the best laef for each significant species
spe.ch.MRT.indval$indcls[which(spe.ch.MRT.indval$pval<=0.05)]


## try just with Date

spe.ch.mvpart <- mvpart(data.matrix(normalized_MPL_OTUs) ~ Temperature_YSI, Jericho_MPL_env,margin=0.08, xv="pick", xval=nrow(params_Jericho),which=4)

## well it works with VC number!!! picks the same divisions as I would by eye!!
spe.ch.seq <- mvpart(data.matrix(normalized_MPL_OTUs) ~ VC_number, Jericho_MPL_env,margin=0.08, xv="pick", xval=nrow(params_Jericho),which=4)

summary(spe.ch.seq)

## Group composition (labels of terminal nodes)

(gr <- spe.ch.seq$where)

## renumber clusters sequentially

aa <- 1
gr2 <- rep(1, length(gr))
for(i in 2:length(gr)){
 if (gr[i] !=gr[i-1]) aa <- aa+1
gr2[i] <- aa
}


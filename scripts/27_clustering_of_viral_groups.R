 

library (ade4)
library (vegan)
library(gclus)
library(cluster)
library(RColorBrewer)
library(labdsv)
library(ggplot2)
library(reshape2)

Jericho_data <- read.csv("../results/Jericho_env_data_mean_imputed.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)
parameters_to_exclude <- c("season",
                           "Average_NO3_NO2",
                           "Standard_error_NO3_NO2",
                           "Secchi_disk_disappears",
                           "Secchi_disk_reappears",
                           "Dissolved_oxygen_percent",
                           #"VC_number",
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

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
## what if I exclude the last sample?
normalized_MPL_OTUs <- subset(normalized_MPL_OTUs, !(rownames(normalized_MPL_OTUs) %in% "1244"))

Jericho_MPL_env <- droplevels(subset(params_Jericho, VC_number %in% row.names(normalized_MPL_OTUs)))
rownames(Jericho_MPL_env) <- Jericho_MPL_env$Date
#normalized_MPL_OTUs$VC_number <- row.names(normalized_MPL_OTUs)

## just work with date
row.names(normalized_MPL_OTUs) <- Jericho_MPL_env$Date[match(rownames(normalized_MPL_OTUs),Jericho_MPL_env$VC_number)]
normalized_MPL_OTUs <- normalized_MPL_OTUs[sort(row.names(normalized_MPL_OTUs)),]
Jericho_MPL_env <- Jericho_MPL_env[sort(row.names(Jericho_MPL_env)),]



spe.norm <- decostand(normalized_MPL_OTUs,"hellinger") ## normalize species abundance table. 
spe.ch <- vegdist(normalized_MPL_OTUs)

## UPGMA clustering

spe.ch.UPGMA <- hclust(spe.ch, method="average")
plot(spe.ch.UPGMA)
summary(spe.ch.UPGMA)
## see weird reversals in this. Cna be interpreted but it seems difficult. 
spe.ch.centroid<- hclust(spe.ch, method="centroid")
plot(spe.ch.centroid)

## skipping ward's

spe.ch.ward <- hclust(spe.ch, method="ward.D")
plot(spe.ch.ward)

library(ggdendro)
# basic option
ggdendrogram(spe.ch.UPGMA)
# another option
ggdendrogram(spe.ch.UPGMA, rotate = TRUE, size = 4, theme_dendro = FALSE, color = "tomato")
# Triangular lines

dendr <- dendro_data(spe.ch.UPGMA, type="triangle") 
ggplot() + 
 geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend))+ geom_text(data = label(dendr), aes(x = x, 
                                                                                                            y = y, label = label), angle = 90, lineheight = 0)


## cophenetic distance

spe.ch.UPGMA.coph <- cophenetic(spe.ch.UPGMA)
cor(spe.ch, spe.ch.UPGMA.coph) ## better correlation
spe.ch.Ward.coph <- cophenetic(spe.ch.ward)
cor(spe.ch, spe.ch.Ward.coph)

plot(spe.ch, spe.ch.UPGMA.coph, xlab="hellinger distance", ylab="Cophenetic distance", asp=1, xlim=c(0, sqrt(2)), ylim=c(0,sqrt(2)))

## can choose certain number of groups...
## what are the appropriate numbers?

k <- 5

spebc.UPGMA.g <- cutree(spe.ch.UPGMA, k)
spebc.Ward.g <- cutree(spe.ch.ward, k)

## compare
table(spebc.UPGMA.g, spebc.Ward.g)


#Optimal number of clusters according to Mantel statistic (pearson)

# Function to compute a binay distance matrix from gorups 

grpdist <- function(X)
{
 require(cluster)
 gr <- as.data.frame(as.factor(X))
 distgr <- daisy(gr, "gower")
 distgr
}

#Ran based on Ward clustering
kt <- data.frame(k=1:nrow(normalized_MPL_OTUs), r=0)

for (i in 2:(nrow(normalized_MPL_OTUs) - 1)){
 gr <- cutree(spe.ch.UPGMA, i)
 distgr <- grpdist(gr)
 mt <- cor(spe.ch, distgr, method="pearson")
 kt[i,2] <- mt
}

kt
k.best <- which.max(kt$r)
#Plot is then produced by function plot.silhouette (cluster)
plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - Ward", xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.best, paste("Optimum", k.best,spe="\n"), col="red", font=2, col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)

#Final plot

## Chooose number of clusters
k <- k.best

cutg <- cutree(spe.ch.UPGMA, k=k)
sil <- silhouette(cutg,spe.ch)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(normalized_MPL_OTUs)[attr(silo,"iOrd")]
plot(silo, main="Silhouette plot-Chord - Ward", cex.names=0.8, col=cutg+1)


spe.chwo <- reorder.hclust(spe.ch.UPGMA, spe.ch)

## plot reorder dendrograms with group labels
plot(spe.chwo, hang=-1, xlab="13 groups", sub="", ylab="Height", main="Chord-Ward(reordered)", label=cutree(spe.chwo, k=k))
rect.hclust(spe.chwo, k=k)

dend <- as.dendrogram(spe.chwo)
heatmap(as.matrix(spe.ch), Rowv=dend, symm=TRUE, margin=c(3,3))

## ordered community table
or <- vegemite(normalized_MPL_OTUs, spe.chwo, scale="log")

heatmap(t(normalized_MPL_OTUs[rev(or$species)]), Rowv=NA, Colv=dend, col=c("white", brewer.pal(5, "Greens")), scale="none", margin=c(4,4), ylab="species (weighted averages of sites)", xlab="sites")

## kmeans
spe.kmeans <- kmeans(spe.ch, centers=4, nstart=100)
table(spe.kmeans$cluster, spebc.UPGMA.g)
spe.kmeans <- kmeans(spe.ch, centers=5, nstart=100)
table(spe.kmeans$cluster, spebc.UPGMA.g)


spe.KM.cascade <- cascadeKM(spe.ch, inf.gr=2, sup.gr=10, iter=100, criterion="ssi")
plot(spe.KM.cascade, sortg=TRUE)
summary(spe.KM.cascade)
spe.KM.cascade$results
spe.KM.cascade$partition
## 8 and up separates properly...


## examine clusters by site

## Define the site clusters
## reorder site according to kmeans results
normalized_MPL_OTUs[order(spe.kmeans$cluster),]

ord.KM <- vegemite(normalized_MPL_OTUs,spe.kmeans$cluster, scale="log") ## idsplays compact ordered community tables
normalized_MPL_OTUs[ord.KM$sites, ord.KM$species]

## compare to environmental data! yes!

boxplot(Jericho_MPL_env$Temperature_YSI ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(Jericho_MPL_env$Salinity_ppt_YSI ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(log(Jericho_MPL_env$Salinity_ppt_YSI) ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(Jericho_MPL_env$Average_bacterial_abundance ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(Jericho_MPL_env$Average_chl_a ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(Jericho_MPL_env$Average_PO4 ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(Jericho_MPL_env$Average_SiO2 ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(Jericho_MPL_env$pH ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(Jericho_MPL_env$day_length ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)


## what about compared to richness in other amplicons
amplicon_richness <- read.csv("../results/amplicon_richness_by_date.csv")
amplicon_richness$Date <- as.Date(amplicon_richness$Date)
amplicon_richness_for_study <- subset(amplicon_richness, Date %in% Jericho_MPL_env$Date )
## temporary fix for missing data. Or should do mean imputation...?

#amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0
boxplot(amplicon_richness_for_study$richness.18S ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(amplicon_richness_for_study$richness.MPL ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(amplicon_richness_for_study$richness.gp23 ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)
boxplot(amplicon_richness_for_study$richness.16S ~ spe.kmeans$cluster, main="Temperature", las=1, ylab="Temp", varwidth=TRUE)

Jericho_with_cluster_from_OTUs <- cbind(Jericho_MPL_env, cluster_n=spe.kmeans$cluster, richness_18s=amplicon_richness_for_study$richness.18S)
Jericho_MPL_env_melt <- melt(Jericho_with_cluster_from_OTUs, id=c("cluster_n", "Date"))

p <- ggplot(Jericho_MPL_env_melt, aes(y=as.numeric(value), x=as.factor(cluster_n)))
p + geom_violin(adjust = .5)+geom_point()+facet_wrap(~variable, scales = "free")


p <- ggplot(subset(Jericho_MPL_env_melt, variable=="day_length"), aes(y=as.numeric(value), x=as.factor(cluster_n)))
p + geom_violin(adjust = .5)+geom_point()
p <- ggplot(subset(Jericho_MPL_env_melt, variable=="Average_chl_a"), aes(y=as.numeric(value), x=as.factor(cluster_n)))
p + geom_violin(adjust = .5)+geom_point()

p <- ggplot(subset(Jericho_MPL_env_melt, variable=="richness_18s"), aes(y=as.numeric(value), x=as.factor(cluster_n)))
p + geom_violin()+geom_point()

shapiro.test(resid(lm(Jericho_MPL_env$Temperature_YSI ~ as.factor(spe.kmeans$cluster))))
## homogeneity of variances
bartlett.test(Jericho_MPL_env$Temperature_YSI ~ as.factor(spe.kmeans$cluster))

## if looks ok I could do anova...if parametric...
## SOOOOO see if the environmental params are different between these clusters!!!
summary(aov(Jericho_MPL_env$Temperature_YSI ~ as.factor(spe.kmeans$cluster)))

## otherwise do Kruskal-wallis test for differences..
kruskal.test(Jericho_MPL_env$Temperature_YSI~ as.factor(spe.kmeans$cluster))

shapiro.test(resid(lm(log(Jericho_MPL_env$Salinity_ppt_YSI) ~ as.factor(spe.kmeans$cluster))))
## homogeneity of variances
bartlett.test(log(Jericho_MPL_env$Salinity_ppt_YSI) ~ as.factor(spe.kmeans$cluster))

## if looks ok I could do anova...if parametric...
## SOOOOO see if the environmental params are different between these clusters!!!
summary(aov(log(Jericho_MPL_env$Salinity_ppt_YSI) ~ as.factor(spe.kmeans$cluster)))

## otherwise do Kruskal-wallis test for differences..
kruskal.test(log(Jericho_MPL_env$Salinity_ppt_YSI)~ as.factor(spe.kmeans$cluster))

## try these plots with groups of 18s by family...


S18_ord_sums <- read.csv("../results/normalized_18s_summarized_by_order.csv", row.names = 1)
S18_fam_sums <- read.csv("../results/normalized_18s_summarized_by_family.csv",row.names = 1)

rownames(S18_fam_sums) <- Jericho_data$Date[match(rownames(S18_fam_sums),Jericho_data$VC_number)]

S18_fam_sum_for_Jer_Env <- droplevels(subset(S18_fam_sums, as.Date(rownames(S18_fam_sums)) %in% Jericho_MPL_env$Date))
S18_fam_sum_for_Jer_Env <- S18_fam_sum_for_Jer_Env[, colSums(S18_fam_sum_for_Jer_Env) > 100]

S18_fam_sum_for_Jer_Env$cluster_otus <- spe.kmeans$cluster[match(rownames(S18_fam_sum_for_Jer_Env ),names(spe.kmeans$cluster))]
S18_fam_sum_for_Jer_Env$Date <- row.names(S18_fam_sum_for_Jer_Env)

S18_fam_sum_for_Jer_Env_melt <- melt(S18_fam_sum_for_Jer_Env, id=c("Date", "cluster_otus"))


p <- ggplot(S18_fam_sum_for_Jer_Env_melt, aes(y=as.numeric(value), x=as.factor(cluster_otus)))
p + geom_violin(adjust = .5)+geom_point()+facet_wrap(~variable, scales = "free")


rownames(S18_ord_sums) <- Jericho_data$Date[match(rownames(S18_ord_sums),Jericho_data$VC_number)]

S18_ord_sum_for_Jer_Env <- droplevels(subset(S18_ord_sums, as.Date(rownames(S18_ord_sums)) %in% Jericho_MPL_env$Date))
S18_ord_sum_for_Jer_Env <- S18_ord_sum_for_Jer_Env[, colSums(S18_ord_sum_for_Jer_Env) > 100]

S18_ord_sum_for_Jer_Env$cluster_otus <- spe.kmeans$cluster[match(rownames(S18_ord_sum_for_Jer_Env ),names(spe.kmeans$cluster))]
S18_ord_sum_for_Jer_Env$Date <- row.names(S18_ord_sum_for_Jer_Env)

S18_ord_sum_for_Jer_Env_melt <- melt(S18_ord_sum_for_Jer_Env, id=c("Date", "cluster_otus"))


p <- ggplot(S18_ord_sum_for_Jer_Env_melt, aes(y=as.numeric(value), x=as.factor(cluster_otus)))
p + geom_violin(adjust = .5)+geom_point()+facet_wrap(~variable, scales = "free")


## See what is happening for 18s!

spe.ch <- vegdist(S18_ord_sum_for_Jer_Env[,-c(dim(S18_ord_sum_for_Jer_Env),dim(S18_ord_sum_for_Jer_Env)-1) ])

spe.kmeans <- kmeans(spe.ch, centers=5, nstart=100)
#table(spe.kmeans$cluster, spebc.UPGMA.g)


spe.KM.cascade <- cascadeKM(spe.ch, inf.gr=2, sup.gr=10, iter=100, criterion="ssi")
plot(spe.KM.cascade, sortg=TRUE)
summary(spe.KM.cascade)
spe.KM.cascade$results
spe.KM.cascade$partition
## 8 and up separates properly...


## examine clusters by site

## Define the site clusters
## reorder site according to kmeans results
S18_ord_sum_for_Jer_Env[order(spe.kmeans$cluster),]

ord.KM <- vegemite(S18_ord_sum_for_Jer_Env[,-c(dim(S18_ord_sum_for_Jer_Env),dim(S18_ord_sum_for_Jer_Env)-1) ],spe.kmeans$cluster, scale="log") ## idsplays compact ordered community tables
S18_ord_sum_for_Jer_Env[,-c(dim(S18_ord_sum_for_Jer_Env),dim(S18_ord_sum_for_Jer_Env)-1) ][ord.KM$sites, ord.KM$species]


S18_ord_sum_for_Jer_Env <- cbind(S18_ord_sum_for_Jer_Env[,-(dim(S18_ord_sum_for_Jer_Env)-1) ], cluster_n=spe.kmeans$cluster)

S18_ord_sum_for_Jer_Env_melt <- melt(S18_ord_sum_for_Jer_Env, id=c("cluster_n", "Date"))

p <- ggplot(S18_ord_sum_for_Jer_Env_melt, aes(y=as.numeric(value), x=as.factor(cluster_n)))
p + geom_violin(adjust = .5)+geom_point()+facet_wrap(~variable, scales = "free")

spe.ch.UPGMA_18s <- hclust(spe.ch, method="average")
## cophenetic between 18s and mpL
spe.ch.cop_18s <- cophenetic(spe.ch.UPGMA_18s)

normalized_MPL_OTUs_sub <- subset(normalized_MPL_OTUs, row.names(normalized_MPL_OTUs) %in% S18_ord_sum_for_Jer_Env$Date)
spe.ch_MPL <- vegdist(normalized_MPL_OTUs_sub)
spe.ch.UPGMA_MPL <- hclust(spe.ch_MPL, method="average")
## needs to have the same 
spe.ch.cop_MPL <- cophenetic(spe.ch.UPGMA_MPL)
cor(spe.ch.cop_18s, spe.ch.cop_MPL)
## very low correlation.
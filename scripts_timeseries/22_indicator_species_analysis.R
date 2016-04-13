## indicator species analysis

library(vegan)
library(labdsv)
library(indicspecies)
library (ade4)
library(gclus)
library(cluster)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dendextend)

## ok so what do I want to do? 
## Want to see if I can test the swapping of groups by something. 
## can I compare the phylogenetic groupings to the 18s groupings?


## check to see if cophenetic between 18s and MPL based on OTU clustering????

## would BioEnv be useful here? or mvpart??


#Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data <- read.csv("../results/Jericho_env_data_mean_imputed.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)
rownames(Jericho_data) <- Jericho_data$VC_number


normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

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




## remove certain columns
Jericho_subset <- Jericho_data[,c(1, 23,25,27,29,31)]
row.names(Jericho_subset) <- Jericho_subset[,1]
Jericho_subset <- Jericho_subset[,-1]
Jericho.gowerD = vegdist(Jericho_subset , "gower", na.rm=TRUE)
# Hierarchical clustering, Ward method
  
Jericho.cl = hclust(Jericho.gowerD, method="ward")
plot(Jericho.cl, hang=-1) # Dendrogram
summary(Jericho.cl)


test_out_bioenv <- function(normalized_otus, Jericho_data){
 union_vcs_A_and_B <- intersect(row.names(normalized_otus), Jericho_data$VC_number)
 
 Jericho_otus <- droplevels(subset(normalized_otus, row.names(normalized_otus) %in% union_vcs_A_and_B))
 
 Jericho_subset_env <- droplevels(subset(Jericho_data, row.names(Jericho_data) %in% union_vcs_A_and_B))
 row.names(Jericho_subset_env) <- Jericho_subset_env[,1]
 str(Jericho_subset_env)
 Jericho_subset_env <- Jericho_subset_env[,c(2,3,4,7,33,31,29,27,25,23,21)]
 
#  normalized_otus$date <- Jericho_data$Date[match(rownames(normalized_otus),Jericho_data$VC_number)]
#  rownames(normalized_otus) <- normalized_otus$date
#  normalized_otus<- normalized_otus[,-length(normalized_otus)]
 
 bioenv(normalized_otus,Jericho_subset_env)
}

test_out_bioenv(normalized_MPL_OTUs, Jericho_data)
test_out_bioenv(normalized_18s_OTUs, Jericho_data)
test_out_bioenv(normalized_16s_OTUs, Jericho_data)
test_out_bioenv(normalized_gp23_OTUs, Jericho_data)

## how many groups to cut

get_env_data_otu_table <- function (normalized_otus) {
 union_vcs_A_and_B <- intersect(row.names(normalized_otus), Jericho_data$VC_number)

 Jericho_subset_env <- droplevels(subset(Jericho_data, row.names(Jericho_data) %in% union_vcs_A_and_B))
 ## dates will be the row names
 row.names(Jericho_subset_env) <- Jericho_subset_env[,1]
 str(Jericho_subset_env)
 Jericho_subset_env <- Jericho_subset_env[,c(1,2,3,4,7,33,31,29,27,25,23,21)]
 
 return(Jericho_subset_env)
}

clustering_with_env_for_otus <- function (normalized_otus, groups_filename) {
  #normalized_otus <- normalized_18s_OTUs
 union_vcs_A_and_B <- intersect(row.names(normalized_otus), Jericho_data$VC_number)
  
 Jericho_otus <- droplevels(subset(normalized_otus, row.names(normalized_otus) %in% union_vcs_A_and_B))
  
 Jericho_subset_env <- droplevels(subset(Jericho_data, row.names(Jericho_data) %in% union_vcs_A_and_B))
 row.names(Jericho_subset_env) <- Jericho_subset_env[,1]
 str(Jericho_subset_env)
 Jericho_subset_env <- Jericho_subset_env[,c(2,3,4,7,33,31,29,27,25,23,21)]
 
 Jericho.gowerD = vegdist(Jericho_subset_env , "gower", na.rm=TRUE)
 Jericho.cl = hclust(Jericho.gowerD, method="ward")
 plot(Jericho.cl, hang=-1, main="Chord-Ward(reordered) based on Env params") # Dendrogram
 summary(Jericho.cl)
 groups = cutree(Jericho.cl , k=6)
 write.csv(groups, groups_filename)
 return(Jericho.cl)
}

clustering_otus_species_abundance <- function(normalized_otus, groups, groups_filename){
 normalized_otus$date <- Jericho_data$Date[match(rownames(normalized_otus),Jericho_data$VC_number)]
 rownames(normalized_otus) <- normalized_otus$date
 normalized_otus<- normalized_otus[,-length(normalized_otus)]
 
 spe.norm <- decostand(normalized_otus,"hellinger") ## normalize species abundance table. 
 spe.ch <- vegdist(normalized_otus)
 
 spe.ch.UPGMA <- hclust(spe.ch, method="ward")
 
 spe.chwo <- reorder.hclust(spe.ch.UPGMA, spe.ch)
 
 ## plot reorder dendrograms with group labels
 plot(spe.chwo, hang=-1, xlab=paste0(groups,"groups"), sub="", ylab="Height", main="Chord-Ward(reordered) based on OTU abun"
      #, label=cutree(spe.chwo, k=k)
 )
 str(spe.chwo)
 rect.hclust(spe.chwo, k=groups)
 groups = cutree(spe.chwo , k=6)
 write.csv(groups, groups_filename)
 dend <- as.dendrogram(spe.chwo)
 return(spe.chwo)
 
}

indicator_species_analysis <- function(normalized_otus, clustering_for_this_otu, groups){
 # Examine the dendrogram. How many groups does it contain?
 # Create a vector ‘groups’ showing the group to which each site belongs,
 # in the original order of entry of the sites.
 groups = cutree(clustering_for_this_otu, k=groups)
 
 union_vcs_A_and_B <- intersect(row.names(normalized_otus), row.names(Jericho_data))
 
 Jericho_otus <- droplevels(subset(normalized_otus, row.names(normalized_otus) %in% union_vcs_A_and_B))
 
 # Look for the indicator species
 indval <- multipatt(Jericho_otus, groups)
 summary(indval)
}

cluster_diagram_env_MPL <- clustering_with_env_for_otus(normalized_MPL_OTUs, "../results/normalized_MPL_env_clusters.csv")
cluster_diagram_otu_MPL <- clustering_otus_species_abundance(normalized_MPL_OTUs,6, "../results/normalized_MPL_otus_clusters.csv")

high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)
normalized_MPL_OTUs_without_high_res <- subset(normalized_MPL_OTUs,!(rownames(normalized_MPL_OTUs) %in% high_res_vcs))

cluster_diagram_otu_MPL_no_high_res <- clustering_otus_species_abundance(normalized_MPL_OTUs_without_high_res ,6, "../results/normalized_MPL_otus_clusters_no_high_res.csv")
cluster_diagram_env_MPL_no_high_res <- clustering_with_env_for_otus(normalized_MPL_OTUs_without_high_res, "../results/normalized_MPL_env_clusters_no_high_res.csv")

cluster_diagram_env_gp23 <- clustering_with_env_for_otus(normalized_gp23_OTUs, "../results/normalized_gp23_env_clusters.csv")
cluster_diagram_otu_gp23 <- clustering_otus_species_abundance(normalized_gp23_OTUs,6, "../results/normalized_gp23_otus_clusters.csv")

high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)
normalized_gp23_OTUs_without_high_res <- subset(normalized_gp23_OTUs,!(rownames(normalized_gp23_OTUs) %in% high_res_vcs))

cluster_diagram_otu_gp23_no_high_res <- clustering_otus_species_abundance(normalized_gp23_OTUs_without_high_res ,6, "../results/normalized_gp23_otus_clusters_no_high_res.csv")

cluster_diagram_env_gp23_no_high_res <- clustering_with_env_for_otus(normalized_gp23_OTUs_without_high_res, "../results/normalized_gp23_env_clusters_no_high_res.csv")

## how many groups to cut

#indicator_species_analysis(normalized_MPL_OTUs, cluster_diagram_env_MPL, 6) #problem because of missing high res.
indicator_species_analysis(normalized_MPL_OTUs, cluster_diagram_otu_MPL, 6)

## or kmeans



## clustering based on env data
cluster_diagram_env_18s <- clustering_with_env_for_otus(normalized_18s_OTUs, "../results/normalized_18s_env_clusters.csv")
cluster_diagram_otu_18s <- clustering_otus_species_abundance(normalized_18s_OTUs,6, "../results/normalized_18s_otus_clusters.csv")

indicator_species_analysis(normalized_18s_OTUs, cluster_diagram_env_18s, 6)



dendograms_with_same_dates <- function (normalized_otusA, normalized_otusB) {
 #normalized_otus <- normalized_18s_OTUs
 union_vcs_A_and_B <- intersect(row.names(normalized_otusA), row.names(normalized_otusB))
 
 Jericho_otusA <- droplevels(subset(normalized_otusA, row.names(normalized_otusA) %in% union_vcs_A_and_B))
 Jericho_otusB <- droplevels(subset(normalized_otusB, row.names(normalized_otusB) %in% union_vcs_A_and_B))
  Jericho.gowerA <-  vegdist(Jericho_otusA , "gower", na.rm=TRUE)
 Jericho.cl.A = hclust(Jericho.gowerA, method="ward")
 
 Jericho.gowerB <-  vegdist(Jericho_otusB , "gower", na.rm=TRUE)
 Jericho.cl.B = hclust(Jericho.gowerB, method="ward")
 
 Jericho_dendlist <- dendlist(as.dendrogram(Jericho.cl.A), as.dendrogram(Jericho.cl.B))

 return(Jericho_dendlist)
}

dendograms_env_vs_otu <- function (normalized_otusA, env_data) {
 union_vcs_A_and_B <- intersect(row.names(normalized_otusA), Jericho_data$VC_number)
 
 Jericho_otus <- droplevels(subset(normalized_otusA, row.names(normalized_otusA) %in% union_vcs_A_and_B))
 
 Jericho_subset_env <- droplevels(subset(Jericho_data, row.names(Jericho_data) %in% union_vcs_A_and_B))
 row.names(Jericho_subset_env) <- Jericho_subset_env[,1]
 Jericho_subset_env <- Jericho_subset_env[,c(2,3,4,7,33,31,29,27,25,23,21)]
 
 Jericho_otus$date <- Jericho_data$Date[match(rownames(Jericho_otus),Jericho_data$VC_number)]
 rownames(Jericho_otus) <- Jericho_otus$date
 Jericho_otus<- Jericho_otus[,-length(Jericho_otus)]
 
 Jericho.gower_env = vegdist(Jericho_subset_env , "gower", na.rm=TRUE)
 Jericho.cl.env = hclust(Jericho.gower_env, method="ward")
str(Jericho.cl.env)
 Jericho.gowerA <-  vegdist(Jericho_otus , "gower", na.rm=TRUE)
 Jericho.cl.A = hclust(Jericho.gowerA, method="ward")

 
 Jericho_dendlist <- dendlist(as.dendrogram(Jericho.cl.env), as.dendrogram(Jericho.cl.A))
 
 return(Jericho_dendlist)
}

test_dend <- dendograms_with_same_dates(normalized_MPL_OTUs, normalized_18s_OTUs)

test_dend <- dendograms_env_vs_otu(normalized_MPL_OTUs, Jericho_data)
dend_diff(test_dend )
tanglegram(test_dend , sort=TRUE)
tanglegram(test_dend , common_subtrees_color_branches = TRUE, sort=TRUE)
entanglement(test_dend)

test_dend %>% untangle(method = "step1side") %>% 
 tanglegram(common_subtrees_color_branches = TRUE)

x <- test_dend
x %>% plot(main = paste("entanglement =", round(entanglement(x), 2)))


test_dend <- dendograms_with_same_dates(normalized_16s_OTUs, normalized_18s_OTUs)

dend_diff(test_dend )
tanglegram(test_dend , sort=TRUE)
tanglegram(test_dend , common_subtrees_color_branches = TRUE, sort=TRUE)
entanglement(test_dend)

test_dend %>% untangle(method = "step1side") %>% 
 tanglegram(common_subtrees_color_branches = TRUE)

x <- test_dend
x %>% plot(main = paste("entanglement =", round(entanglement(x), 2)))

test_dend <- dendograms_with_same_dates(normalized_16s_OTUs, normalized_gp23_OTUs)

dend_diff(test_dend )
tanglegram(test_dend , sort=TRUE)
tanglegram(test_dend , common_subtrees_color_branches = TRUE, sort=TRUE)
entanglement(test_dend)

test_dend %>% untangle(method = "step1side") %>% 
 tanglegram(common_subtrees_color_branches = TRUE)

x <- test_dend
x %>% plot(main = paste("entanglement =", round(entanglement(x), 2)))

MPL_env_data <- get_env_data_otu_table(normalized_MPL_OTUs)
### look at groups split up by env data

groups = cutree(cluster_diagram_otu_MPL, k=6)
Jericho_with_cluster_from_OTUs <- cbind(MPL_env_data, cluster_n=groups )
Jericho_MPL_env_melt <- melt(Jericho_with_cluster_from_OTUs, id=c("cluster_n", "Date"))

p <- ggplot(Jericho_MPL_env_melt, aes(y=as.numeric(value), x=as.factor(cluster_n)))
p + geom_violin(adjust = .5)+geom_point()+facet_wrap(~variable, scales = "free")


gp23_env_data <- get_env_data_otu_table(normalized_gp23_OTUs)
### look at groups split up by env data

groups = cutree(cluster_diagram_otu_gp23, k=6)
Jericho_with_cluster_from_OTUs <- cbind(gp23_env_data, cluster_n=groups )
Jericho_gp23_env_melt <- melt(Jericho_with_cluster_from_OTUs, id=c("cluster_n", "Date"))

p <- ggplot(Jericho_gp23_env_melt, aes(y=as.numeric(value), x=as.factor(cluster_n)))
p + geom_violin(adjust = .5)+geom_point()+facet_wrap(~variable, scales = "free")



S18_env_data <- get_env_data_otu_table(normalized_18s_OTUs)
### look at groups split up by env data

groups = cutree(cluster_diagram_env_18s, k=6)
Jericho_with_cluster_from_OTUs <- cbind(S18_env_data, cluster_n=groups )
Jericho_18s_env_melt <- melt(Jericho_with_cluster_from_OTUs, id=c("cluster_n", "Date"))

p <- ggplot(Jericho_18s_env_melt, aes(y=as.numeric(value), x=as.factor(cluster_n)))
p + geom_violin(adjust = .5)+geom_point()+facet_wrap(~variable, scales = "free")



shapiro.test(resid(lm(S18_env_data$Temperature_YSI ~ as.factor(groups))))
## homogeneity of variances
bartlett.test(S18_env_data$Temperature_YSI ~ as.factor(groups))

## if looks ok I could do anova...if parametric...
## SOOOOO see if the environmental params are different between these clusters!!!
summary(aov(S18_env_data$Temperature_YSI ~ as.factor(groups)))

## otherwise do Kruskal-wallis test for differences..
kruskal.test(S18_env_data$Temperature_YSI~ as.factor(groups))

shapiro.test(resid(lm(log(S18_env_data$Temperature_YSI) ~ as.factor(groups))))
## homogeneity of variances
bartlett.test(log(S18_env_data$Temperature_YSI) ~ as.factor(groups))


### what if I used the clustering in MPL to look at indicator species in 18s

## need same dates...

indicator_species_using_clustering_of_another <- function (normalized_otusA, normalized_otusB) {
 #normalized_otus <- normalized_18s_OTUs
 union_vcs_A_and_B <- intersect(row.names(normalized_otusA), row.names(normalized_otusB))
 
 Jericho_otusA <- droplevels(subset(normalized_otusA, row.names(normalized_otusA) %in% union_vcs_A_and_B))
 Jericho_otusB <- droplevels(subset(normalized_otusB, row.names(normalized_otusB) %in% union_vcs_A_and_B))
 Jericho.gowerA <-  vegdist(Jericho_otusA , "gower", na.rm=TRUE)
 Jericho.cl.A = hclust(Jericho.gowerA, method="ward")
 
 groups = cutree(Jericho.cl.A , k=6)
 spe.ch.indval  <- multipatt(Jericho_otusB, groups
                             #,control = how(nperm=999)
                             )
 summary(spe.ch.indval)

 return(spe.ch.indval)
}


test_ind <- indicator_species_using_clustering_of_another(normalized_MPL_OTUs, normalized_18s_OTUs)

test_ind <- indicator_species_using_clustering_of_another(normalized_16s_OTUs, normalized_18s_OTUs)


## would it make sense to plot cluster membership underneath the bars in script with phylogenetic groups?? 

#then look if there are any indicator species pulled from 18s and plot their relative abundances??? 


















#### MVPart regression based tree

## check out just the env data
mvpart(data.matrix(params_Jericho[,-1])~.,params_Jericho[,-1])
mvpart(data.matrix(params_Jericho[,c(-1,-5,-6)])~.,params_Jericho[,c(-1,-5,-6)], margin=0.08, cp=0, xv="pick", xval=nrow(params_Jericho), xvmult=100, which=4)

Jericho_MPL_env <- droplevels(subset(params_Jericho, VC_number %in% row.names(normalized_MPL_OTUs)))
rownames(Jericho_MPL_env) <- Jericho_MPL_env$Date
#normalized_MPL_OTUs$VC_number <- row.names(normalized_MPL_OTUs)

## just work with date
row.names(normalized_MPL_OTUs) <- Jericho_MPL_env$Date[match(rownames(normalized_MPL_OTUs),Jericho_MPL_env$VC_number)]
normalized_MPL_OTUs <- normalized_MPL_OTUs[sort(row.names(normalized_MPL_OTUs)),]
Jericho_MPL_env <- Jericho_MPL_env[sort(row.names(Jericho_MPL_env)),]


spe.ch.mvpart <- mvpart(data.matrix(normalized_MPL_OTUs) ~., Jericho_MPL_env[,c(-1, -5)])


spe.ch.mvpart <- mvpart(data.matrix(normalized_MPL_OTUs) ~., Jericho_MPL_env[,c(-1, -5)], size=6, xvmult=100)

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

## compare this to the clustering results??

## can i put this on the clustering tree???


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

## combine MRT with Indicator species 
library(labdsv)
## indicator species searc on MRT result
spe.ch.MRT.indval <- indval(normalized_MPL_OTUs, spe.ch.mvpart$where)
spe.ch.MRT.indval$pval


spe.ch.MRT.indval  <- multipatt(normalized_MPL_OTUs, spe.ch.mvpart$where)
summary(spe.ch.MRT.indval)


## ok the results are really different between the clusters generated by MRT and by cluster analysis. 



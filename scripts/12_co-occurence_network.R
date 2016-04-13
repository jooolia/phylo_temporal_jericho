
library(vegan)
library(phyloseq)
library(igraph)
library(Hmisc)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(corrgram)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

#### want site as rows.

### still working on this.!!!
# data(sipoo)
# ## Matrix temperature
# out <- nestedtemp(sipoo)
# out
# plot(out)
# plot(out, kind="incid")
# ## Use oecosimu to assess the nonrandomness of checker board units
# nestedchecker(sipoo)
# oecosimu(sipoo, nestedchecker, "quasiswap")
# ## Another Null model and standardized checkerboard score
# oecosimu(sipoo, nestedchecker, "r00", statistic = "C.score")
# 
 out <- nestedtemp(normalized_16s_OTUs)
 out
 plot(out)
# plot(out, kind="incid")
 nestedchecker(normalized_16s_OTUs)
 oecosimu(normalized_16s_OTUs, nestedchecker, "quasiswap")
## would it make sense to only take those with 5 reads or more each?

generate_network_via_spearman_from_otu_table <- function (normalized_OTU_table, name) {
  normalized_OTUs_more_than_10 <- normalized_OTU_table[, colSums(normalized_OTU_table) > 10]  
  
  normalized_OTUs_more_than_100 <- normalized_OTU_table[, colSums(normalized_OTU_table) > 100]  
  
  test <- rcorr(as.matrix(normalized_OTUs_more_than_100),type="spearman") 
  correlations <- test$r
  p_values <- test$P
  

  pdf(paste("../figures/corrgram", name,".pdf", sep=""),width = 11, height = 11)
   cor_colour <- colorRampPalette(brewer.pal(5, "Purples"))
  corrgram(correlations, order=TRUE, lower.panel=panel.shade,
           upper.panel=NULL, text.panel=panel.txt,
           main="sorted", 
           label.pos=c(0.5,0.5),
           cex.label=0.25,
           col.regions=cor_colour)
  dev.off()
  ## get rid of the upper section by turning to NA. 
  correlations[upper.tri(correlations )] <- NA
  p_values[upper.tri(p_values)] <- NA
  
   
  melted_cor <- melt(correlations)
  melted_p <- melt(p_values)
  
  melted_together <- cbind(melted_p$value, melted_cor)
  melted_together <- na.omit(melted_together) ## gets rid of the leftover diagonals
  names(melted_together)
  names(melted_together)[1] <- "p_value"
  #filter for correlation about 0.6 and p value less than 0.01
  filtered_data <- filter(melted_together, p_value <= 0.01  & value > 0.6)
  
  names(filtered_data)
  names(filtered_data)[4] <- "weight"
  
  filtered_data <- subset(filtered_data, select=c(Var1, Var2, weight))
  
  
  graph <- graph.data.frame(filtered_data,
                            directed=FALSE)
  write.graph(graph,file=paste("../figures/network_cooccurence_spearmand", name,".graphml", sep=""),format="graphml")
  return(graph)
}

## Generate network plot ####
MPL_network <- generate_network_via_spearman_from_otu_table(normalized_MPL_OTUs, "MPL")

plot(MPL_network)
gp23_network <- generate_network_via_spearman_from_otu_table(normalized_gp23_OTUs, "gp23")

plot(gp23_network)

#AVS_network <- generate_network_via_spearman_from_otu_table(normalized_AVS_OTUs, "AVS")

#plot(AVS_network)

S16_network <- generate_network_via_spearman_from_otu_table(normalized_16s_OTUs, "16S")

plot(S16_network)

S18_network <- generate_network_via_spearman_from_otu_table(normalized_18s_OTUs, "18s")

plot(S18_network)

### get taxonomic info to add...####
taxnomic_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)


### what if I use only those that are more than 10% of populations ####

## first change all the OTUs to something with descriptive name
ided_normalized_MPL_OTUs <- normalized_MPL_OTUs
names(ided_normalized_MPL_OTUs) <- gsub("OTU", "MPL_", names(ided_normalized_MPL_OTUs))

proportional_MPL <- prop.table(as.matrix(ided_normalized_MPL_OTUs),  margin=1)
ided_normalized_MPL_OTUs <- as.matrix(ided_normalized_MPL_OTUs)
## get rid of those below 10% of relative abundance
ided_normalized_MPL_OTUs[which(proportional_MPL<0.01)] <- 0


ided_normalized_gp23_OTUs <- normalized_gp23_OTUs
names(ided_normalized_gp23_OTUs) <- gsub("OTU", "gp23_", names(ided_normalized_gp23_OTUs))

proportional_gp23 <- prop.table(as.matrix(ided_normalized_gp23_OTUs),  margin=1)
ided_normalized_gp23_OTUs <- as.matrix(ided_normalized_gp23_OTUs)
## get rid of those below 10% of relative abundance
ided_normalized_gp23_OTUs[which(proportional_gp23 <0.01)] <- 0


ided_normalized_AVS_OTUs <- normalized_AVS_OTUs
names(ided_normalized_AVS_OTUs) <- gsub("OTU", "AVS_", names(ided_normalized_AVS_OTUs))

proportional_AVS <- prop.table(as.matrix(ided_normalized_AVS_OTUs),  margin=1)
ided_normalized_AVS_OTUs <- as.matrix(ided_normalized_AVS_OTUs)
## get rid of those below 10% of relative abundance
ided_normalized_AVS_OTUs[which(proportional_AVS <0.01)] <- 0


ided_normalized_16S_OTUs <- normalized_16s_OTUs
names(ided_normalized_16S_OTUs) <- gsub("OTU", "S16_", names(ided_normalized_16S_OTUs))

proportional_16S <- prop.table(as.matrix(ided_normalized_16S_OTUs),  margin=1)
ided_normalized_16S_OTUs <- as.matrix(ided_normalized_16S_OTUs)
## get rid of those below 10% of relative abundance
ided_normalized_16S_OTUs[which(proportional_16S <0.01)] <- 0


ided_normalized_18S_OTUs <- normalized_18s_OTUs
names(ided_normalized_18S_OTUs) <- gsub("OTU", "S18_", names(ided_normalized_18S_OTUs))

proportional_18S <- prop.table(as.matrix(ided_normalized_18S_OTUs),  margin=1)
ided_normalized_18S_OTUs <- as.matrix(ided_normalized_18S_OTUs)
## get rid of those below 10% of relative abundance
ided_normalized_18S_OTUs[which(proportional_18S <0.01)] <- 0




merged_gp23_and_16s <- merge(ided_normalized_gp23_OTUs, ided_normalized_16S_OTUs,by=0,all=TRUE)
row.names(merged_gp23_and_16s) <- merged_gp23_and_16s$Row.names
merged_gp23_and_16s  <- subset(merged_gp23_and_16s, select=-c(Row.names))
merged_gp23_and_16s[is.na(merged_gp23_and_16s)] <- 0


merged_18s_and_AVS <- merge(ided_normalized_18S_OTUs, ided_normalized_AVS_OTUs,by=0,all=TRUE,  incomparables = 0)
row.names(merged_18s_and_AVS) <- merged_18s_and_AVS$Row.names
merged_18s_and_AVS  <- subset(merged_18s_and_AVS, select=-c(Row.names))
merged_18s_and_AVS[is.na(merged_18s_and_AVS)] <- 0


merged_18s_and_AVS_MPL <- merge(merged_18s_and_AVS, ided_normalized_MPL_OTUs,by="row.names",all=TRUE)
row.names(merged_18s_and_AVS_MPL) <- merged_18s_and_AVS_MPL$Row.names
merged_18s_and_AVS_MPL  <- subset(merged_18s_and_AVS_MPL, select=-c(Row.names))
merged_18s_and_AVS_MPL[is.na(merged_18s_and_AVS_MPL)] <- 0


merged_16s_and_18s <- merge(ided_normalized_16S_OTUs, ided_normalized_18S_OTUs,by="row.names",all=TRUE)
row.names(merged_16s_and_18s) <- merged_16s_and_18s$Row.names
merged_16s_and_18s  <- subset(merged_16s_and_18s, select=-c(Row.names))
merged_16s_and_18s[is.na(merged_16s_and_18s)] <- 0


merged_all_together <- merge(merged_18s_and_AVS_MPL, merged_gp23_and_16s,by="row.names",all=TRUE)
row.names(merged_all_together) <- merged_all_together$Row.names
merged_all_together  <- subset(merged_all_together, select=-c(Row.names))
merged_all_together[is.na(merged_all_together)] <- 0
##get rid of garbage factors
str(merged_all_together)
 merged_all_together_numeric <- apply(as.matrix(merged_all_together), 2,as.numeric)
row.names(merged_all_together_numeric) <- row.names(merged_all_together) 
str(merged_all_together_numeric)

## get rid of columns only 0
merged_for_correlation <- merged_all_together_numeric[,colSums(merged_all_together_numeric)>0]
merged_18s_and_AVS_MPL <- merged_18s_and_AVS_MPL[,colSums(merged_18s_and_AVS_MPL)>0]
merged_gp23_and_16s <- merged_gp23_and_16s[,colSums(merged_gp23_and_16s)>0]

### Make the networks:

rm(ided_normalized_16S_OTUs,ided_normalized_18S_OTUs,ided_normalized_AVS_OTUs,ided_normalized_gp23_OTUs,ided_normalized_MPL_OTUs)

spearman_correlation_and_generate_network <- function (merged_identified_OTU_table, names) {
  spearman <- rcorr(as.matrix(merged_identified_OTU_table),type="spearman") 
  correlations <- spearman$r
  p_values <- spearman$P
  
  
  ## get rid of the upper section by turning to NA. 
  correlations[upper.tri(correlations )] <- NA
  p_values[upper.tri(p_values)] <- NA
  
  
  melted_cor <- melt(correlations)
  melted_p <- melt(p_values)
  
  melted_together <- cbind(melted_p$value, melted_cor)
  melted_together <- na.omit(melted_together) ## gets rid of the leftover diagonals
  names(melted_together)
  names(melted_together)[1] <- "p_value"
  
  #filter for correlation about 0.6 and p value less than 0.01
  filtered_data <- filter(melted_together, p_value <= 0.01  & value > 0.6)
  
  names(filtered_data)
  names(filtered_data)[4] <- "weight"
  
  filtered_data <- subset(filtered_data, select=c(Var1, Var2, weight))
  graph <- graph.data.frame(filtered_data,
                            directed=FALSE)
  write.graph(graph,file=paste("../figures/network_cooccurence_spearman_", names,"_all_types_proportion_above0_01.graphml",sep=""),format="graphml")
  return(graph)
}

# Network_all <- spearman_correlation_and_generate_network(merged_for_correlation, "merged_all")
# 
# Network_18s_AVS_MPL <- spearman_correlation_and_generate_network(merged_18s_and_AVS_MPL, "merged_18s_and_AVS_MPL")
# Network_16s_gp23 <- spearman_correlation_and_generate_network(merged_gp23_and_16s, "merged_16s_and_gp23")
# Network_16s_and_18s <- spearman_correlation_and_generate_network(merged_16s_and_18s, "merged_16s_and_18s")

## open in cytoscape. Involves manual curation of the style, the placement algorithm and the taxonomic data. Not sure how to get it to be better and how to do less clicking. 
##4 combine the 16s and 18s and export and then import into cytoscape


## for not I am getting them so that I can save it as a network and just reediting it. Not ideal. The OTU names for the 16s and 18s need to be hand-edited and it is potentially error prone. 

## want to look over time!

library(vegan)
library(phyloseq)
library(igraph)
library(Hmisc)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(reshape2)
library(corrgram)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)
taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
## Reformat date
# %d is day as a number, %b is abreviated month in words, %y is 2digit year
Jericho_data$Date <- as.Date(Jericho_data$Date)


### what if I use only those that are more than 10% of populations ####
## or try with only those that occur in more than 3 sites.
min_number_sites <- 10

## first change all the OTUs to something with descriptive name
ided_normalized_MPL_OTUs <- normalized_MPL_OTUs
names(ided_normalized_MPL_OTUs) <- gsub("OTU", "MPL_", names(ided_normalized_MPL_OTUs))

## get columns where OTUs are in more than 3 sites
nonzero <- function(x) sum(x != 0)
count_sites_MPL <- numcolwise(nonzero)(ided_normalized_MPL_OTUs)
## so want OTUs that occur in 3 or more samples 
count_sites_MPL <- count_sites_MPL[,count_sites_MPL> min_number_sites]
## subset the table by the otus that are there 3 or more times. 
ided_normalized_MPL_OTUs <- subset(ided_normalized_MPL_OTUs,select = colnames(count_sites_MPL))

proportional_MPL <- prop.table(as.matrix(ided_normalized_MPL_OTUs),  margin=1)
ided_normalized_MPL_OTUs <- as.matrix(ided_normalized_MPL_OTUs)
## get rid of those below 10% of relative abundance
## test it with them still there. 
ided_normalized_MPL_OTUs[which(proportional_MPL<0.005)] <- 0

ided_normalized_gp23_OTUs <- normalized_gp23_OTUs
names(ided_normalized_gp23_OTUs) <- gsub("OTU", "gp23_", names(ided_normalized_gp23_OTUs))

count_sites_gp23 <- numcolwise(nonzero)(ided_normalized_gp23_OTUs)
## so want OTUs that occur in 3 or more samples 
count_sites_gp23 <- count_sites_gp23[,count_sites_gp23> min_number_sites]
## subset the table by the otus that are there 3 or more times. 
ided_normalized_gp23_OTUs <- subset(ided_normalized_gp23_OTUs,select = colnames(count_sites_gp23))

proportional_gp23 <- prop.table(as.matrix(ided_normalized_gp23_OTUs),  margin=1)
ided_normalized_gp23_OTUs <- as.matrix(ided_normalized_gp23_OTUs)
## get rid of those below 10% of relative abundance
ided_normalized_gp23_OTUs[which(proportional_gp23 <0.005)] <- 0

ided_normalized_AVS_OTUs <- normalized_AVS_OTUs
names(ided_normalized_AVS_OTUs) <- gsub("OTU", "AVS_", names(ided_normalized_AVS_OTUs))

count_sites_AVS <- numcolwise(nonzero)(ided_normalized_AVS_OTUs)
## so want OTUs that occur in 3 or more samples 
count_sites_AVS <- count_sites_AVS[,count_sites_AVS> min_number_sites]
## subset the table by the otus that are there 3 or more times. 
ided_normalized_AVS_OTUs <- subset(ided_normalized_AVS_OTUs,select = colnames(count_sites_AVS))

proportional_AVS <- prop.table(as.matrix(ided_normalized_AVS_OTUs),  margin=1)
ided_normalized_AVS_OTUs <- as.matrix(ided_normalized_AVS_OTUs)
## get rid of those below 10% of relative abundance
ided_normalized_AVS_OTUs[which(proportional_AVS <0.005)] <- 0


ided_normalized_16S_OTUs <- normalized_16s_OTUs
names(ided_normalized_16S_OTUs) <- gsub("OTU", "S16", names(ided_normalized_16S_OTUs))
names(ided_normalized_16S_OTUs) <- gsub(".size.*.", "", names(ided_normalized_16S_OTUs))

count_sites_16s <- numcolwise(nonzero)(ided_normalized_16S_OTUs)
## so want OTUs that occur in 3 or more samples 
count_sites_16s <- count_sites_16s[,count_sites_16s> min_number_sites]
## subset the table by the otus that are there 3 or more times. 
ided_normalized_16S_OTUs <- subset(ided_normalized_16S_OTUs,select = colnames(count_sites_16s))

proportional_16S <- prop.table(as.matrix(ided_normalized_16S_OTUs),  margin=1)
ided_normalized_16S_OTUs <- as.matrix(ided_normalized_16S_OTUs)
## get rid of those below 10% of relative abundance
ided_normalized_16S_OTUs[which(proportional_16S <0.005)] <- 0

ided_normalized_18S_OTUs <- normalized_18s_OTUs
names(ided_normalized_18S_OTUs) <- gsub(".size.*.", "", names(ided_normalized_18S_OTUs))

## want to remove the opistokonts
# create new row called Phylum
ided_normalized_18S_OTUs["Phylum",] <- taxonomy_18s$Phylum[match(colnames(ided_normalized_18S_OTUs),taxonomy_18s$otu_number)]
id_normalized_18s_extra_row_no_opist <- ided_normalized_18S_OTUs[, grep("Opisthokonta", ided_normalized_18S_OTUs["Phylum",], invert=TRUE)]

new_normalized_18S_OTUs <- (subset(id_normalized_18s_extra_row_no_opist, rownames(id_normalized_18s_extra_row_no_opist) != "Phylum"))
names(new_normalized_18S_OTUs) <- gsub("OTU", "S18", names(new_normalized_18S_OTUs))
keep_rownames <- row.names(new_normalized_18S_OTUs)

new_normalized_18S_OTUs <- as.data.frame(sapply(new_normalized_18S_OTUs, as.numeric))
#long_18s_with_taxonomy <- droplevels(subset(long_18s_with_taxonomy, Domain == "Eukaryota"))
row.names(new_normalized_18S_OTUs) <- keep_rownames
ided_normalized_18S_OTUs <- new_normalized_18S_OTUs

count_sites_18s <- numcolwise(nonzero)(ided_normalized_18S_OTUs)
## so want OTUs that occur in 3 or more samples 
count_sites_18s <- count_sites_18s[,count_sites_18s> min_number_sites]
## subset the table by the otus that are there 3 or more times. 
ided_normalized_18S_OTUs <- subset(ided_normalized_18S_OTUs,select = colnames(count_sites_18s))

proportional_18S <- prop.table(as.matrix(ided_normalized_18S_OTUs),  margin=1)
ided_normalized_18S_OTUs <- as.matrix(ided_normalized_18S_OTUs)
## get rid of those below 10% of relative abundance
ided_normalized_18S_OTUs[which(proportional_18S <0.005)] <- 0

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
str(merged_16s_and_18s)

merged_all_together <- merge(merged_18s_and_AVS_MPL, merged_gp23_and_16s,by="row.names",all=TRUE)
row.names(merged_all_together) <- merged_all_together$Row.names
merged_all_together  <- subset(merged_all_together, select=-c(Row.names))
merged_all_together[is.na(merged_all_together)] <- 0
##get rid of garbage factors
str(merged_all_together)
#merged_all_together_numeric <- apply(as.matrix(merged_all_together), 2,as.numeric)
#row.names(merged_all_together_numeric) <- row.names(merged_all_together) 
#str(merged_all_together_numeric)


## get rid of columns only 0
#merged_for_correlation <- merged_all_together_numeric[,colSums(merged_all_together_numeric)>0]
merged_for_correlation <- merged_all_together[,colSums(merged_all_together)>0]
merged_16s_and_18s <- merged_16s_and_18s[,colSums(merged_16s_and_18s) > 0]
merged_18s_and_AVS_MPL <- merged_18s_and_AVS_MPL[,colSums(merged_18s_and_AVS_MPL)>0]
merged_gp23_and_16s <- merged_gp23_and_16s[,colSums(merged_gp23_and_16s)>0]

### Make the networks:

#rm(ided_normalized_16S_OTUs,ided_normalized_18S_OTUs,ided_normalized_AVS_OTUs,ided_normalized_gp23_OTUs,ided_normalized_MPL_OTUs)


taxonomic_16_and_18s <- read.csv(  "../results/cleaned_up_taxonomies_with_renamed_OTUs.csv")

spearman_correlation_and_generate_network_with_temporal_presence <- function (merged_identified_OTU_table, 
                                                                              names, 
                                                                              first_amplicon,
                                                                              second_amplicon, 
                                                                              third_amplicon, 
                                                                              fourth_amplicon,
                                                                              fifth_amplicon) {
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
 
 #filter for correlation about 0.6 (use absolute to make it positive or negative) and p value less than 0.01
 filtered_data <- filter(melted_together, p_value <= 0.01  & abs(value) > 0.6)
 
 names(filtered_data)
 names(filtered_data)[4] <- "weight"
 
 filtered_data <- subset(filtered_data, select=c(Var1, Var2, weight))
 
 filtered_data$cor <- filtered_data$weight
 filtered_data$weight <- abs(filtered_data$weight)
write.csv(filtered_data, paste("../results/significant_strong_correlation",names, ".csv",sep="" ))
 ## so from here I would make the time series networks to see how the connections would change..
 
 for (row in seq(nrow(merged_identified_OTU_table))){
  is_there <- c()
  merged_original_otus_1_row <- droplevels(merged_identified_OTU_table[row,])
  melted_original_otus_1_row <- melt(merged_original_otus_1_row)
  filtered_original_otus_1_row <- filter(melted_original_otus_1_row , value > 0)
  
  for (row in seq(nrow(filtered_data))){
   row_of_interest <- filtered_data[row,]
  # print(row_of_interest)
   ## if these significantly correlated otus are found at one time point give add the date found to a data frame
   if (row_of_interest$Var2  %in% filtered_original_otus_1_row$variable & row_of_interest$Var1 %in% filtered_original_otus_1_row$variable){
   # print("yes")
    is_there <- c(is_there,"1")
   }
   else{
    is_there <- c(is_there,"0")
   }
  }
  filtered_data <- cbind(filtered_data, is_there)
  names(filtered_data)[ncol(filtered_data)] <- row.names(merged_original_otus_1_row)
 }
 
 graph <- graph.data.frame( filtered_data,directed=FALSE)
 
 for (row in seq(length(V(graph)$name))){
  if (grepl(first_amplicon,  V(graph)$name[row])){
   V(graph)$type[row] <- first_amplicon
  }
  else if (grepl(second_amplicon,  V(graph)$name[row])){
   V(graph)$type[row] <- second_amplicon
  }
  else if (grepl(third_amplicon,  V(graph)$name[row])){
   V(graph)$type[row] <- third_amplicon
  }
  else if (grepl(fourth_amplicon,  V(graph)$name[row])){
   V(graph)$type[row] <- fourth_amplicon
  }
  else {
   V(graph)$type[row] <- fifth_amplicon 
  }
 }

 ## Add in the taxonomic information for each node if it is available. 
 ## Could also add in lowest classification....
 V( graph)$Domain <- as.character(taxonomic_16_and_18s$Domain[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 V( graph)$Phylum <- as.character(taxonomic_16_and_18s$Phylum[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 V( graph)$Class <- as.character(taxonomic_16_and_18s$Class[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 V( graph)$Order <- as.character(taxonomic_16_and_18s$Order[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 V( graph)$Family <- as.character(taxonomic_16_and_18s$Family[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 V( graph)$Genus <- as.character(taxonomic_16_and_18s$Genus[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 
 ## would like to check on some of the weird weights...
 summary(E( graph)$weight)
 ## could I add degree to the node attributes too so that it is more automated in cytoscape
 V(graph)$degree <- degree(graph)
 ## Many networks consist of modules which are densely connected themselves but sparsely connected to other modules.
 wc <- cluster_edge_betweenness(graph)
 V(graph)$membership <- wc$membership
 ## was giving an error in the 16s and 18s sections
# V(graph)$betweenness <- betweenness(graph)
 
 ## like this when I import into Cytoscape the edge attributes are included in the plot
 write.graph(graph,file=paste("../figures/network_cooccurence_spearman_", names,"_proportion_above0_01_temporal_data.graphml",sep=""),format="graphml")
 return(graph)
}


# Network_16s_gp23 <- spearman_correlation_and_generate_network_with_temporal_presence(merged_gp23_and_16s, "merged_16s_and_gp23", "gp23", "S16", "", "","")
 
Network_18s_AVS_MPL <- spearman_correlation_and_generate_network_with_temporal_presence(merged_18s_and_AVS_MPL, "merged_18s_and_AVS_MPL","AVS", "MPL", "S18", "","" )
# 
# Network_16s_and_18s <- spearman_correlation_and_generate_network_with_temporal_presence(merged_16s_and_18s, "merged_16s_and_18s","S16", "S18", "" ,"","")
# 
 Network_all <- spearman_correlation_and_generate_network_with_temporal_presence(merged_for_correlation, "merged_all","AVS", "MPL", "S18", "gp23", "S16" )


### make just single ones too. ided_normalized_18S_OTUs
## have to turn from a matrix into a data frame
## temporarily commented out
# ided_normalized_18S_OTUs_df <- as.data.frame(ided_normalized_18S_OTUs)
# Network_18S <- spearman_correlation_and_generate_network_with_temporal_presence(ided_normalized_18S_OTUs_df, "18s", "S18", "", "", "","")
# 
# ided_normalized_16S_OTUs_df <- as.data.frame(ided_normalized_16S_OTUs)
# Network_16S <- spearman_correlation_and_generate_network_with_temporal_presence(ided_normalized_16S_OTUs_df, "16s", "S16", "", "", "","")
# 
# ided_normalized_AVS_OTUs_df <- as.data.frame(ided_normalized_AVS_OTUs)
# Network_AVS <- spearman_correlation_and_generate_network_with_temporal_presence(ided_normalized_AVS_OTUs_df, "AVS", "AVS", "", "", "","")
# 
# ided_normalized_gp23_OTUs_df <- as.data.frame(ided_normalized_gp23_OTUs)
# Network_gp23 <- spearman_correlation_and_generate_network_with_temporal_presence(ided_normalized_gp23_OTUs_df, "gp23", "gp23", "", "", "","")
# 
# ided_normalized_MPL_OTUs_df <- as.data.frame(ided_normalized_MPL_OTUs)
# Network_MPL <- spearman_correlation_and_generate_network_with_temporal_presence(ided_normalized_MPL_OTUs_df, "MPL", "MPL", "", "", "","")


spearman_correlation_and_network_with_temporal_only_edges_btw_unlike <- function (merged_identified_OTU_table,
                                                                                  names, 
                                                                                  first_amplicon,
                                                                                  second_amplicon,
                                                                                  third_amplicon,
                                                                                  fourth_amplicon,
                                                                                  fifth_amplicon) {
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
## filter out those that have the same type
library(stringr)
melted_together <- filter(melted_together,str_split_fixed(melted_together$Var1, "_", n=2)[,1]!= str_split_fixed(melted_together$Var2, "_", n=2)[,1])

#filter for correlation about 0.6 (use absolute to make it positive or negative) and p value less than 0.01
 filtered_data <- filter(melted_together, p_value <= 0.01  & abs(value) > 0.6)
 
 names(filtered_data)
 names(filtered_data)[4] <- "weight"
 
 filtered_data <- subset(filtered_data, select=c(Var1, Var2, weight))
 ## so from here I would make the time series networks to see how the connections would change..
 ## doing this so that the 
 filtered_data$cor <- filtered_data$weight
 filtered_data$weight <- abs(filtered_data$weight)
 
 
 for (row in seq(nrow(merged_identified_OTU_table))){
  is_there <- c()
  merged_original_otus_1_row <- droplevels(merged_identified_OTU_table[row,])
  melted_original_otus_1_row <- melt(merged_original_otus_1_row)
  filtered_original_otus_1_row <- filter(melted_original_otus_1_row , value > 0)
  
  for (row in seq(nrow(filtered_data))){
   row_of_interest <- filtered_data[row,]
   # print(row_of_interest)
   ## if these significantly correlated otus are found at one time point give add the date found to a data frame
   if (row_of_interest$Var2  %in% filtered_original_otus_1_row$variable & row_of_interest$Var1 %in% filtered_original_otus_1_row$variable){
    # print("yes")
    is_there <- c(is_there,"1")
   }
   else{
    is_there <- c(is_there,"0")
   }
  }
  filtered_data <- cbind(filtered_data, is_there)
  names(filtered_data)[ncol(filtered_data)] <- row.names(merged_original_otus_1_row)
 }
 
 graph <- graph.data.frame( filtered_data,directed=FALSE)
 
 for (row in seq(length(V(graph)$name))){
  if (grepl(first_amplicon,  V(graph)$name[row])){
   V(graph)$type[row] <- first_amplicon
  }
  else if (grepl(second_amplicon,  V(graph)$name[row])){
   V(graph)$type[row] <- second_amplicon
  }
  else if (grepl(third_amplicon,  V(graph)$name[row])){
   V(graph)$type[row] <- third_amplicon
  }
  else if (grepl(fourth_amplicon,  V(graph)$name[row])){
   V(graph)$type[row] <- fourth_amplicon
  }
  else {
   V(graph)$type[row] <- fifth_amplicon 
  }
 }
 
 ## Add in the taxonomic information for each node if it is available. 
 ## Could also add in lowest classification....
 V( graph)$Domain <- as.character(taxonomic_16_and_18s$Domain[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 V( graph)$Phylum <- as.character(taxonomic_16_and_18s$Phylum[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 V( graph)$Class <- as.character(taxonomic_16_and_18s$Class[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 V( graph)$Order <- as.character(taxonomic_16_and_18s$Order[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 V( graph)$Family <- as.character(taxonomic_16_and_18s$Family[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 V( graph)$Genus <- as.character(taxonomic_16_and_18s$Genus[match(V( graph)$name,taxonomic_16_and_18s$otu_number)])
 wc <- cluster_edge_betweenness(graph)
 V(graph)$membership <- wc$membership
 ## would like to check on some of the weird weights...
 summary(E( graph)$weight)
 ## could I add degree to the node attributes too so that it is more automated in cytoscape
 V(graph)$degree <- degree(graph)
 
 wc <- cluster_edge_betweenness(graph)
 V(graph)$membership <- wc$membership
 
 ## was giving an error in the 16s and 18s sections
 # V(graph)$betweenness <- betweenness(graph)
 
 ## like this when I import into Cytoscape the edge attributes are included in the plot
 write.graph(graph,file=paste("../figures/network_cooccurence_spearman_", names,"__temporal_data_only_edges_btwn_unlike.graphml",sep=""),format="graphml")
 return(graph)
}

Network_16s_gp23_no_self_links <- spearman_correlation_and_network_with_temporal_only_edges_btw_unlike(merged_gp23_and_16s, "merged_16s_and_gp23", "gp23", "S16", "", "","")

Network_18s_AVS_MPL_no_self_links <- spearman_correlation_and_network_with_temporal_only_edges_btw_unlike(merged_18s_and_AVS_MPL, "merged_18s_and_AVS_MPL","AVS", "MPL", "S18", "","" )

#Network_16s_and_18s_no_self_links <- spearman_correlation_and_network_with_temporal_only_edges_btw_unlike(merged_16s_and_18s, "merged_16s_and_18s","S16", "S18", "" ,"","")

Network_all_no_self_links <- spearman_correlation_and_network_with_temporal_only_edges_btw_unlike(merged_for_correlation, "merged_all","AVS", "MPL", "S18", "gp23", "S16" )

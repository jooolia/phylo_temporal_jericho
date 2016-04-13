## Print out taxonomy from the files read in. 
## Use this so I don't have to keep processing the taxonomy info
library(reshape2)
normalized_18s_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_18S.tsv",
                                  row.names="VC_number")

taxonomy_18s <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ngnon_chimeras_ref97.00.nr_v119.wang.taxonomy",header= FALSE)

taxonomy_18s$V1 <- gsub(".size.*.", "", taxonomy_18s$V1)
names(taxonomy_18s) <- c("otu_number", "silva_taxonomy")

Jericho_OTUs <- unique(names(normalized_18s_OTUs))
Jericho_OTUs <- gsub(".size.*.", "", Jericho_OTUs)

taxonomy_only_Jericho <- droplevels(subset(taxonomy_18s, otu_number %in% Jericho_OTUs, select=c(otu_number, silva_taxonomy)))

## so just need OTU number and then split up...

split_up_taxonomy <- colsplit(taxonomy_only_Jericho$silva_taxonomy, ";", c("Domain", "Phylum", "Class", "Order", "Family","Genus"))
taxonomy_only_Jericho <- cbind(taxonomy_only_Jericho,split_up_taxonomy)
# add in the phylum, class, order and family data

taxonomy_only_Jericho$Domain <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Domain)
taxonomy_only_Jericho$Phylum <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Phylum)
taxonomy_only_Jericho$Class <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Class)
taxonomy_only_Jericho$Order <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Order)
taxonomy_only_Jericho$Family <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Family)
taxonomy_only_Jericho$Genus <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Genus)
taxonomy_only_Jericho$Genus <- gsub(";", "", taxonomy_only_Jericho$Genus)

clean_18s_taxonomy <- taxonomy_only_Jericho[,c("otu_number","Domain", "Phylum", "Class", "Order", "Family","Genus" )]

write.csv(clean_18s_taxonomy, "../results/cleaned_up_18s_taxonomy_Jericho.csv")

test <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)

## want unique classifications
classifications <- test[,2:7]
unique_18s_classifications <- unique(classifications)
write.csv(unique_18s_classifications, "../results/unique_18s_taxonomy_Jericho.csv")

## do it for the whole data set ####

normalized_18s_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_18S.tsv",
                                  row.names="VC_number")

Jericho_OTUs <- unique(names(normalized_18s_OTUs))
Jericho_OTUs <- gsub(".size.*.", "", Jericho_OTUs)

taxonomy_only_Jericho <- droplevels(subset(taxonomy_18s, otu_number %in% Jericho_OTUs, select=c(otu_number, silva_taxonomy)))

## so just need OTU number and then split up...

split_up_taxonomy <- colsplit(taxonomy_only_Jericho$silva_taxonomy, ";", c("Domain", "Phylum", "Class", "Order", "Family","Genus"))
taxonomy_only_Jericho <- cbind(taxonomy_only_Jericho,split_up_taxonomy)
# add in the phylum, class, order and family data

taxonomy_only_Jericho$Domain <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Domain)
taxonomy_only_Jericho$Phylum <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Phylum)
taxonomy_only_Jericho$Class <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Class)
taxonomy_only_Jericho$Order <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Order)
taxonomy_only_Jericho$Family <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Family)
taxonomy_only_Jericho$Genus <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Genus)
taxonomy_only_Jericho$Genus <- gsub(";", "", taxonomy_only_Jericho$Genus)

clean_18s_taxonomy <- taxonomy_only_Jericho[,c("otu_number","Domain", "Phylum", "Class", "Order", "Family","Genus" )]

write.csv(clean_18s_taxonomy, "../results/cleaned_up_18s_taxonomy_all_time_series_Jericho.csv")

## want unique classifications
classifications <- test[,2:7]
unique_18s_classifications <- unique(classifications)
write.csv(unique_18s_classifications, "../results/unique_18s_taxonomy_all_time_series_Jericho.csv")


## Do the 16s taxonomny now =====

taxonomy_16s <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_forwardgood.good.pcr.ngnon_chimeras_ref97.00.nr_v119.wang.taxonomy",header= FALSE)

normalized_16s_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_16S_R1.tsv", row.names="VC_number") 


taxonomy_16s$V1 <- gsub(".size.*.", "", taxonomy_16s$V1)
names(taxonomy_16s) <- c("otu_number", "silva_taxonomy")

Jericho_OTUs <- unique(names(normalized_16s_OTUs))
Jericho_OTUs <- gsub(".size.*.", "", Jericho_OTUs) # because "size" attached to OTUs

taxonomy_only_Jericho <- droplevels(subset(taxonomy_16s, otu_number %in% Jericho_OTUs, select=c(otu_number, silva_taxonomy)))

## so just need OTU number and then split up...

split_up_taxonomy <- colsplit(taxonomy_only_Jericho$silva_taxonomy, ";", c("Domain", "Phylum", "Class", "Order", "Family","Genus"))
taxonomy_only_Jericho <- cbind(taxonomy_only_Jericho,split_up_taxonomy)
# add in the phylum, class, order and family data

taxonomy_only_Jericho$Domain <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Domain)
taxonomy_only_Jericho$Phylum <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Phylum)
taxonomy_only_Jericho$Class <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Class)
taxonomy_only_Jericho$Order <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Order)
taxonomy_only_Jericho$Family <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Family)
taxonomy_only_Jericho$Genus <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Genus)
taxonomy_only_Jericho$Genus <- gsub(";", "", taxonomy_only_Jericho$Genus)

clean_16s_taxonomy <- taxonomy_only_Jericho[,c("otu_number","Domain", "Phylum", "Class", "Order", "Family","Genus" )]

write.csv(clean_16s_taxonomy, "../results/cleaned_up_16s_taxonomy_Jericho.csv")
test <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)

## want unique classifications
classifications <- test[,2:7]
unique_16s_classifications <- unique(classifications)
write.csv(unique_16s_classifications, "../results/unique_16s_taxonomy_Jericho.csv")

#### Do it for the whole time series


normalized_16s_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

Jericho_OTUs <- unique(names(normalized_16s_OTUs))
Jericho_OTUs <- gsub(".size.*.", "", Jericho_OTUs) # because "size" attached to OTUs

taxonomy_only_Jericho <- droplevels(subset(taxonomy_16s, otu_number %in% Jericho_OTUs, select=c(otu_number, silva_taxonomy)))

## so just need OTU number and then split up...

split_up_taxonomy <- colsplit(taxonomy_only_Jericho$silva_taxonomy, ";", c("Domain", "Phylum", "Class", "Order", "Family","Genus"))
taxonomy_only_Jericho <- cbind(taxonomy_only_Jericho,split_up_taxonomy)
# add in the phylum, class, order and family data

taxonomy_only_Jericho$Domain <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Domain)
taxonomy_only_Jericho$Phylum <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Phylum)
taxonomy_only_Jericho$Class <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Class)
taxonomy_only_Jericho$Order <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Order)
taxonomy_only_Jericho$Family <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Family)
taxonomy_only_Jericho$Genus <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Genus)
taxonomy_only_Jericho$Genus <- gsub(";", "", taxonomy_only_Jericho$Genus)

clean_16s_taxonomy <- taxonomy_only_Jericho[,c("otu_number","Domain", "Phylum", "Class", "Order", "Family","Genus" )]

write.csv(clean_16s_taxonomy, "../results/cleaned_up_16s_taxonomy_all_time_series_Jericho.csv")
test <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)

## want unique classifications
classifications <- test[,2:7]
unique_16s_classifications <- unique(classifications)
write.csv(unique_16s_classifications, "../results/unique_16s_taxonomy_Jericho.csv")



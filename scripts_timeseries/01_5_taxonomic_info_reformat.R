## Print out taxonomy from the files read in. 
## Pull in the taxonomy for 16s and 18s specific for this project
## Author: Julia Gustavsen

library(reshape2)
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv",
                                  row.names="VC_number")

taxonomy_18s <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ngnon_chimeras_ref97.00.nr_v119.wang.pick.taxonomy",header= FALSE)
taxonomy_18s$V1 <- gsub(".size.*.", "", taxonomy_18s$V1)
names(taxonomy_18s) <- c("otu_number", "silva_taxonomy")

Jericho_OTUs <- unique(names(normalized_18s_OTUs))
Jericho_OTUs <- gsub(".size.*.", "", Jericho_OTUs)

taxonomy_only_Jericho <- droplevels(subset(taxonomy_18s, otu_number %in% Jericho_OTUs, select=c(otu_number, silva_taxonomy)))

split_up_taxonomy <- colsplit(taxonomy_only_Jericho$silva_taxonomy, ";", c("Domain", "Phylum", "Class", "Order", "Family","Genus"))
taxonomy_only_Jericho <- cbind(taxonomy_only_Jericho,split_up_taxonomy)
# add in the phylum, class, order and family data

## This regular expression gets rid of the confidence levels in the taxonomy file. 
## Could use the confidence to filter the classifications. 
taxonomy_only_Jericho$Domain <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Domain)
taxonomy_only_Jericho$Phylum <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Phylum)
taxonomy_only_Jericho$Class <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Class)
taxonomy_only_Jericho$Order <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Order)
taxonomy_only_Jericho$Family <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Family)
taxonomy_only_Jericho$Genus <- sub("\\(.{1,4}\\)", "", taxonomy_only_Jericho$Genus)
taxonomy_only_Jericho$Genus <- gsub(";", "", taxonomy_only_Jericho$Genus)

clean_18s_taxonomy <- taxonomy_only_Jericho[,c("otu_number","Domain", "Phylum", "Class", "Order", "Family","Genus" )]

write.csv(clean_18s_taxonomy, "../results/cleaned_up_18s_taxonomy_Jericho.csv")

reread_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)

## want unique classifications
classifications <- reread_18s[,2:7]
unique_18s_classifications <- unique(classifications)
write.csv(unique_18s_classifications, "../results/unique_18s_taxonomy_Jericho.csv")

## Do the 16s taxonomny now =====

taxonomy_16s <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_forwardgood.good.pcr.ngnon_chimeras_ref97.00.nr_v119.wang.pick.taxonomy",header= FALSE)

taxonomy_16s$V1 <- gsub(".size.*.", "", taxonomy_16s$V1)
names(taxonomy_16s) <- c("otu_number", "silva_taxonomy")

normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

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

## for use with networks.
# taxonomic_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)
# taxonomic_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)

taxonomic_16s_reformat_OTU_names <- clean_16s_taxonomy
taxonomic_16s_reformat_OTU_names$otu_number <- gsub("OTU", "S16",taxonomic_16s_reformat_OTU_names$otu_number)

taxonomic_18s_reformat_OTU_names <- clean_18s_taxonomy
taxonomic_18s_reformat_OTU_names$otu_number <- gsub("OTU", "S18",taxonomic_18s_reformat_OTU_names$otu_number)

all_otus_taxonomic_16s_and_18s <- rbind(taxonomic_16s_reformat_OTU_names, taxonomic_18s_reformat_OTU_names)

write.csv(all_otus_taxonomic_16s_and_18s, "../results/cleaned_up_taxonomies_with_renamed_OTUs.csv", row.names=FALSE)




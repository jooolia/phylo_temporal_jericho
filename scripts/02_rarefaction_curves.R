## Rarefaction

library(vegan)
library(reshape2)
library(ggplot2)
library(broom)
source("rarefaction_function.R")


# args <- commandArgs(TRUE)
# inputFile <- args[1]
# 
# ## test to see if input file is given, so I can decide whether to use this argument or the orginal one. 
# if (!file_test("-f", inputFile)) {
#  print("input theme not defined, using orginal one for manuscript.")
#  source("../../JAG_manuscript_figure.R")
# } else {
#  print("Cool you passed a nice theme file to this script")
#  source(inputFile)
# }


original_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_AVS1.tsv", row.names="VC_number")
original_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_MPL.tsv", row.names="VC_number")
original_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_gp23.tsv", row.names="VC_number")
original_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_18S.tsv", row.names="VC_number")
original_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_16S_R1.tsv", row.names="VC_number") 


#normalized_AVS_R1_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R2.tsv", row.names="VC_number")


## Collector's curves
# plot(specaccum(original_16s_OTUs), xlab = "# of samples", ylab = "# of species")
# plot(specaccum(original_18s_OTUs), xlab = "# of samples", ylab = "# of species")
# plot(specaccum(original_gp23_OTUs), xlab = "# of samples", ylab = "# of species")
# plot(specaccum(original_AVS_OTUs), xlab = "# of samples", ylab = "# of species")
# plot(specaccum(original_MPL_OTUs), xlab = "# of samples", ylab = "# of species")


#Community data should be as "species" for columns and sites as rows. 

#Rarefy the dataset

rarefaction_curves_by_amplicon <- function (original_OTU_table, title) {
 #original_OTU_table <- t(original_OTU_table)
 original_OTU_table[is.na(original_OTU_table)] <- 0
 
 #  OTU_table.rare <- rarefaction(original_OTU_table, col=F)
 #   
 #   tidied_OTU_rarefied_richness <- tidy(OTU_table.rare)
 #   melt_rarefied <- melt(tidied_OTU_rarefied_richness, id="subsample")
 #   
 #   #get rid of Standard error for now
 #   new_melt <- droplevels(subset(melt_rarefied, grepl("richness", variable)))
 #   
 #   p <- ggplot(new_melt, aes(x=subsample, y=value,group=variable)) +geom_line(aes(colour=variable))+ggtitle(title)
 pdf(paste0("../figures/Rarefaction_curves_",title,"%03d.pdf"), width = 11, height = 15, onefile=FALSE)
 plot(specaccum(original_OTU_table), xlab = "# of samples", ylab = "# of species")
 #Adding in the sum of all the OTUs
 Overall_AllOTUs <- colSums(original_OTU_table)
 AllOTUs_transposed <- rbind(original_OTU_table, Overall_AllOTUs)
 rarecurve(AllOTUs_transposed, sample=1)
 # print(p)
 dev.off()
}

## Rarefaction curves for individual amplicons. 
AVS_rarefaction_curves <- rarefaction_curves_by_amplicon(original_AVS_OTUs, "AVS1_original")

MPL_rarefaction_curves <- rarefaction_curves_by_amplicon(original_MPL_OTUs, "MPL_original")

gp23_rarefaction_curves <- rarefaction_curves_by_amplicon(original_gp23_OTUs, "gp23_original")

S18_rarefaction_curves <- rarefaction_curves_by_amplicon(original_18s_OTUs, "18S_original")

S16_rarefaction_curves <- rarefaction_curves_by_amplicon(original_16s_OTUs, "16S_R1_original")

## Print to file 
 pdf("../figures/Rarefaction_curves_original_data%03d.pdf", width = 11, height = 15, onefile=FALSE)
AVS_rarefaction_curves
MPL_rarefaction_curves
gp23_rarefaction_curves
S18_rarefaction_curves
S16_rarefaction_curves
dev.off()
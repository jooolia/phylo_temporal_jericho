## make table of all of the normalized reads
## useful for getting an idea about the normalized reads 

library(reshape2)
library(plyr)
library(dplyr)
library(gtools)
library(xtable)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)
## Summarise 

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv",
                                  row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv", 
                                   row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", 
                                  row.names="VC_number")
#normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
normalized_AVS_R1_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R2.tsv", row.names="VC_number")
#normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S.tsv", row.names="VC_number") 
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

count_16s_otus <- adply(t(normalized_16s_OTUs), 2,function(x)sum(x>0))
count_18s_otus <-adply(t(normalized_18s_OTUs), 2, function(x)sum(x>0))
count_gp23_otus <-adply(t(normalized_gp23_OTUs), 2, function(x)sum(x>0))
count_MPL_otus <-adply(t(normalized_MPL_OTUs), 2, function(x)sum(x>0))
count_AVS1_otus <-adply(t(normalized_AVS_R1_OTUs), 2, function(x)sum(x>0))
#count_AVS2_otus <-adply(t(normalized_AVS_R2_OTUs), 2, function(x)sum(x>0))
#count_AVS_concat_otus <-adply(t(normalized_AVS_OTUs), 2, function(x)sum(x>0))


all_data <- Reduce(function(x, y) merge(x, y, all=TRUE, by="X1"), list(count_16s_otus, 
                                                                       count_18s_otus, 
                                                                       count_gp23_otus,
                                                                       count_MPL_otus,
                                                                       count_AVS1_otus
                                                                       #,
                                                                       #count_AVS2_otus,
                                                                       #count_AVS_concat_otus
))
colnames(all_data) <- c("VC","OTUs_16S", "OTUS_18S", "OTUS_gp23","OTUs_MPL", "OTUs_AVS1")
# "OTUs_AVS2", "OTUs_AVS_concat"
all_data <- all_data[mixedorder(all_data$VC),]
## get rid of ones with no associated VC. 
## maybe should revisit since the MPL actually had quite a few. 
all_data <- subset(all_data, VC!="")
## add in date

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

all_data$Date <- Jericho_data$Date[match(all_data$VC, Jericho_data$VC_number)]

## xtable doesn't play with dates well
all_data$Date <- as.character(all_data$Date)
tab_all_otus <- xtable(all_data)

print(tab_all_otus, type="html", file="../results/Jericho_table_otu_counts_by_VC.html",include.rownames=FALSE)

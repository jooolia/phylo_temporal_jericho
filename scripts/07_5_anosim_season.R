
## NMDS plots of the communities of the OTUs. 

library(ggplot2)
library(gridExtra)
library(vegan)
library(reshape2)
library(phyloseq)
library(dplyr)
library(ape)
library(cowplot)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv",                                  row.names=1)


### get relevant VC numbers
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)


MPL.dist <- vegdist(normalized_MPL_OTUs)
MPL.season <- Jericho_data$season[match(as.numeric(rownames(normalized_MPL_OTUs)),
                                        Jericho_data$VC_number)]
MPL.ano <- anosim(MPL.dist, MPL.season)
summary(MPL.ano)
plot(MPL.ano)
MPL.ado <- adonis(normalized_MPL_OTUs ~ MPL.season,  permutations=990)


gp23.dist <- vegdist(normalized_gp23_OTUs)
gp23.season <- Jericho_data$season[match(as.numeric(rownames(normalized_gp23_OTUs)),
                                        Jericho_data$VC_number)]
gp23.ano <- anosim(gp23.dist, gp23.season)
summary(gp23.ano)
#plot(gp23.ano)
gp23.ado <- adonis(normalized_gp23_OTUs ~ gp23.season,  permutations=999)


s16.dist <- vegdist(normalized_16s_OTUs)
s16.season <- Jericho_data$season[match(as.numeric(rownames(normalized_16s_OTUs)),
                                         Jericho_data$VC_number)]
s16.ano <- anosim(s16.dist, s16.season)
summary(s16.ano)
s16.ado <- adonis(normalized_16s_OTUs ~ s16.season,  permutations=999)


s18.dist <- vegdist(normalized_18s_OTUs)
s18.season <- Jericho_data$season[match(as.numeric(rownames(normalized_18s_OTUs)),
                                        Jericho_data$VC_number)]
s18.ano <- anosim(s18.dist, s18.season)
summary(s18.ano)
#plot(s18.ano)
s18.ado <- adonis(normalized_18s_OTUs ~ s18.season,  permutations=999)

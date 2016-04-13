library(lattice)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(MASS)
library(GGally)



normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv", row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S.tsv", row.names="VC_number") 

# source("../../JAG_black_presentation.R")
source("../../JAG_manuscript_figure.R")

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

## don't need the standard error:
enviro_keep <- c("VC_number",
                 "Date",
                 "season",
                 "pH",
                 "Temperature_YSI",
                 "Salinity_ppt_YSI",
                 "Dissolved_oxygen_percent",
                 "Average_viral_abundance",
                 "Average_bacterial_abundance",
                 "Average_chl_a",
                 "Average_PO4",
                 "Average_SiO2",
                 "Average_NO3_NO2"
) 

Jericho_data <- Jericho_data[ , names(Jericho_data) %in% enviro_keep]

## add in the seq facility by amplicon



do_parallel_coordinates_on_otu_table <- function (normalized_OTUs,Jericho_data) {
  ## melt as matrix to capture the row names. 
  melted_OTUs <- melt(as.matrix(normalized_OTUs))
  names(melted_OTUs) <- c("VC_number", "OTU_number", "Reads")
  
  melted_OTUs_with_environmental <- merge(melted_OTUs, Jericho_data, by="VC_number", all=FALSE,)
  melted_OTUs_with_environmental <- filter(melted_OTUs_with_environmental, Reads > 1)
  
 print( ggparcoord(melted_OTUs_with_environmental,  columns = c(2,5,7, 10:15)))
  parcoord(na.omit(melted_OTUs_with_environmental[,c(1,5,7, 10:15)]), col=rainbow(length(melted_OTUs_with_environmental[,1])), var.label=TRUE)
  
  ## isn't really giving me what I want. 
  
  # what if I averaged by OTU number
  summarized_melted_OTUs_with_environmental <- ddply(melted_OTUs_with_environmental,
                                                    .(OTU_number),
                                                    summarize, 
                                                    mean_reads = mean(Reads),
                                                    mean_temp=mean(Temperature_YSI, na.rm=TRUE),
                                                    mean_sal=mean(Salinity_ppt_YSI, na.rm=TRUE), 
                                                    mean_DO = mean(Dissolved_oxygen_percent, na.rm=TRUE),
                                                    mean_pH=mean(pH, na.rm=TRUE),
                                                    mean_va = mean(Average_viral_abundance, na.rm=TRUE),
                                                    mean_ba=mean(Average_bacterial_abundance, na.rm=TRUE),
                                                    mean_chla = mean(Average_chl_a, na.rm=TRUE),
                                                    mean_PO4 = mean(Average_PO4, na.rm=TRUE),
                                                    mean_SIO2 = mean(Average_SiO2, na.rm=TRUE),
                                                    mean_NO3_NO2 = mean(Average_NO3_NO2, na.rm=TRUE)
  )
  
  print(ggparcoord(summarized_melted_OTUs_with_environmental,  columns = c(5,3,7,8,4,10,12,11,9), groupColumn=1, scale="center"))
  parcoord(na.omit(summarized_melted_OTUs_with_environmental[, c(5,3,7,8,4,10,12,11,9)]), col=rainbow(length(summarized_melted_OTUs_with_environmental[,1])), var.label=TRUE)
}

parallel_AVS <- do_parallel_coordinates_on_otu_table(normalized_AVS_OTUs, Jericho_data)
parallel_16s <- do_parallel_coordinates_on_otu_table(normalized_16s_OTUs, Jericho_data)
parallel_18s <- do_parallel_coordinates_on_otu_table(normalized_18s_OTUs, Jericho_data)
parallel_gp23 <- do_parallel_coordinates_on_otu_table(normalized_gp23_OTUs, Jericho_data)
parallel_MPL <- do_parallel_coordinates_on_otu_table(normalized_MPL_OTUs, Jericho_data)

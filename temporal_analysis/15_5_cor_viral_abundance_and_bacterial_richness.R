### check correlations between viral abundance and bacterial richness

library(Hmisc)
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)


richness <- read.csv("../results/amplicon_richness_by_date.csv")
subset_richness <- richness[!(is.na(richness$richness.16S)),]

Dates_in_common <- intersect(as.Date(subset_richness$Date), as.Date(Jericho_data$Date))

subset_Jericho <- subset(Jericho_data, Date %in% Dates_in_common)
subset_richness_16s <- subset(subset_richness, as.Date(Date) %in% Dates_in_common)


rcorr(subset_Jericho$Average_viral_abundance, subset_richness_16s$richness.16S)
## weak negative correlation, but not significant. 

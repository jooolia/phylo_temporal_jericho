## map data in the library to the VC number to get the bloom samples...


library(vegan)
library(plyr)
library(dplyr)


###### Bloom sample analysis
#source("../../Miseq_Initial_Run_Processing/scripts/functions.R")
# Read in environmental data ---------------------------------------
## subset so that it is just bloom data
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)

Jericho_data$Date <- as.Date(Jericho_data$Date)

## Use library sequencing info to identify which pool corresponds to which VC
Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=83)


S18_data <- c()
S18_data <- as.character(Jericho_data$VC_number)

S18_data <- c(S18_data, "1255_large")

S18_data <- as.data.frame(S18_data)
names(S18_data)[1] <- "VC_number"

S18_data$Pool_number <-Library_metadata$PCR_pool_number[match(S18_data$VC_number,
                                              as.character(Library_metadata$X18s_VC_number))]
S18_data$Pool_number <- as.numeric(as.character(S18_data$Pool_number))

write.csv(S18_data, "../results/S18_VC_number_with_pool.csv")


S16_data <- c()
S16_data <- as.character(Jericho_data$VC_number)

S16_data <- c(S16_data, "1255_large")

S16_data <- as.data.frame(S16_data)
names(S16_data)[1] <- "VC_number"

S16_data$Pool_number <-Library_metadata$PCR_pool_number[match(S16_data$VC_number,
                                                              as.character(Library_metadata$X16s_VC_number))]

write.csv(S16_data, "../results/S16_VC_number_with_pool.csv")



MPL_data <- c()
MPL_data <- as.character(Jericho_data$VC_number)


MPL_data <- as.data.frame(MPL_data)
names(MPL_data)[1] <- "VC_number"

MPL_data$Pool_number <-Library_metadata$PCR_pool_number[match(MPL_data$VC_number,
                                                              as.character(Library_metadata$MPL_VC_number))]
write.csv(MPL_data, "../results/MPL_VC_number_with_pool.csv")


gp23_data <- c()
gp23_data <- as.character(Jericho_data$VC_number)

gp23_data <- c(gp23_data, "1256_ctab")

gp23_data <- as.data.frame(gp23_data)
names(gp23_data)[1] <- "VC_number"

gp23_data$Pool_number <-Library_metadata$PCR_pool_number[match(gp23_data$VC_number,
                                                              as.character(Library_metadata$gp23_VC_number))]

write.csv(gp23_data, "../results/gp23_VC_number_with_pool.csv")


## map data in the library to the VC number to get the bloom samples...


### Do it with all of the data  ####

## subset so that it is all the data. 
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis_all_time_series.csv", row.names=1)

Jericho_data$Date <- as.Date(Jericho_data$Date)

## Use library sequencing info to identify which pool corresponds to which VC
Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=83)


S18_data <- c()
S18_data <- as.character(Jericho_data$VC_number)

S18_data <- c(S18_data, "1255_large","1256_ctab","1223_redo","1225_redo")

S18_data <- as.data.frame(S18_data)
names(S18_data)[1] <- "VC_number"

S18_data$Pool_number <-Library_metadata$PCR_pool_number[match(S18_data$VC_number,
                                                              as.character(Library_metadata$X18s_VC_number))]
S18_data$Pool_number <- as.numeric(as.character(S18_data$Pool_number))

write.csv(S18_data, "../results/S18_VC_number_with_pool_with_all_times.csv")


S16_data <- c()
S16_data <- as.character(Jericho_data$VC_number)

S16_data <- c(S16_data,"1255_large","1256_ctab","1223_redo","1225_redo")

S16_data <- as.data.frame(S16_data)
names(S16_data)[1] <- "VC_number"

S16_data$Pool_number <-Library_metadata$PCR_pool_number[match(S16_data$VC_number,
                                                              as.character(Library_metadata$X16s_VC_number))]

write.csv(S16_data, "../results/S16_VC_number_with_pool_with_all_times.csv")



MPL_data <- c()
MPL_data <- as.character(Jericho_data$VC_number)
MPL_data <-c(MPL_data, "1256_ctab", "1255_large")


MPL_data <- as.data.frame(MPL_data)
names(MPL_data)[1] <- "VC_number"

MPL_data$Pool_number <-Library_metadata$PCR_pool_number[match(MPL_data$VC_number,
                                                              as.character(Library_metadata$MPL_VC_number))]
write.csv(MPL_data, "../results/MPL_VC_number_with_pool_with_all_times.csv")


gp23_data <- c()
gp23_data <- as.character(Jericho_data$VC_number)

gp23_data <- c(gp23_data, "1256_ctab", "1255_large")

gp23_data <- as.data.frame(gp23_data)
names(gp23_data)[1] <- "VC_number"

gp23_data$Pool_number <-Library_metadata$PCR_pool_number[match(gp23_data$VC_number,
                                                               as.character(Library_metadata$gp23_VC_number))]

write.csv(gp23_data, "../results/gp23_VC_number_with_pool_with_all_times.csv")



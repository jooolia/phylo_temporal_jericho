library(vegan)
library(plyr)
library(dplyr)


###### Bloom sample analysis but with all of the data.
# Read in environmental data ---------------------------------------
## subset so that it is just bloom data
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis_all_time_series.csv", row.names=1)

Jericho_data$Date <- as.Date(Jericho_data$Date)

## Use library sequencing info to identify which pool corresponds to which VC
Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=83)

## Process these data and then print out an OTU table from the project and a norma

normalize_and_get_OTUs_occuring_more_than_1x_from_OTU_table <- function (OTU_table_from_project){
 #Community data should be as "species" for columns and sites as rows. 
 OTU_table_without_VC_column <- OTU_table_from_project[,-which(names(OTU_table_from_project)=="VC_number")]
 #Calculate the site with the least amount of reads
 #could do this more automatically by summing counts of each OTU.
 min_reads <- min(apply(OTU_table_without_VC_column,1,sum))
 min_reads
 # Normalize data to the library with the lowest number of reads ---------
 AllOTUs_rarefied_iteration <- lapply(as.list(1:3), 
                                      function(x) rrarefy( OTU_table_without_VC_column,
                                                           sample=min_reads)) 
 
 # takes the median value of all of the iterations. 
 AllOTUs_rarefied_by_lowest_with_zeros <- aaply(laply(AllOTUs_rarefied_iteration,
                                                      as.matrix),
                                                c(2,3),
                                                median) 
 OTU_in_time_series_more_than_1_reads <- AllOTUs_rarefied_by_lowest_with_zeros[, colSums(AllOTUs_rarefied_by_lowest_with_zeros) > 1]  
 
 #if I transpose then I can add the date and nutrients to all of them. 
 OTUs_by_site <- data.frame(OTU_in_time_series_more_than_1_reads)
 OTUs_by_site$VC_number <- OTU_table_from_project$VC_number
 
 return(OTUs_by_site)
}

get_project_samples <- function (OTU_table,
                                 Library_metadata,
                                 Library_VC,
                                 Project_data){
#  OTU_table <- OTU_table_gp23
#  Library_VC <- Library_metadata$gp23_VC_number
#  Project_data <- gp23_data
 ## Get OTUs just belonging to the specific project
 ## Will normalize afterwards
 colnames(OTU_table) <- gsub("pool", 
                             "", 
                             colnames(OTU_table))  
 ## transpose data frame so that rows are sites and OTUs are columns
 transposed_table <- as.data.frame(t(OTU_table)) 
 transposed_table$VC_number <-Library_VC[match(rownames(transposed_table),
                                               as.character(Library_metadata$PCR_pool_number))]
 ##Data frame row names are silently dropped. To preserve, convert to an explicit variable.
 transposed_table_project_only <- filter(transposed_table,
                                         transposed_table$VC_number %in% as.factor(Project_data$VC_number)) 
 ## so since there are no row names keep VC number
 
 ## Get rid of OTU columns that no longer have OTUs in them. 
 transposed_table_removed_OTUs_no_longer_present <- transposed_table_project_only[,
                                                                                  colSums(transposed_table_project_only[,
                                                                                                                        -which(names(transposed_table_project_only)=="VC_number")]) > 1]

 ## Since I have added in some redos this checks to see if there are duplicate
 ## VCs and makes them all different. 
 transposed_table_removed_OTUs_no_longer_present$VC_number <-  make.unique(as.character(transposed_table_removed_OTUs_no_longer_present$VC_number))
 
 return(transposed_table_removed_OTUs_no_longer_present)
}


normalize_otu_table_by_project <-  function (amplicon_time_series,
                                             #  OTU_table,
                                             #                                              Library_metadata,
                                             #                                              Library_VC,
                                             #                                              Project_data, 
                                             amplicon){
 ## Function to pull out the data that is specific to this project 
 ## from the overall OTU table.  It is matched by VC number and 
 ## then pulled into a new table and printed to a file.
 ## The OTU table is also normalized and then printed to a file
 ## so that I don't have to keep re-running these long processes.
 ## first just writes the raw data from the project
 
 #amplicon_time_series <- get_project_samples(OTU_table,
 #                                           Library_metadata,
 #                                           Library_VC,
 #                                           Project_data)
 
 write.table(amplicon_time_series,
             paste("../data/OTU_table_bloom_with_all_time_series_",
                   amplicon, 
                   ".tsv",
                   sep=""),
             sep="\t",
             row.names=FALSE)
 # pulls out the normalized amplicon data. 
 OTU_table_time_series_normalized_amplicon <- normalize_and_get_OTUs_occuring_more_than_1x_from_OTU_table(amplicon_time_series)
 
 write.table(OTU_table_time_series_normalized_amplicon, 
             paste("../data/OTU_table_bloom_with_all_time_series_normalized_", 
                   amplicon, 
                   ".tsv", 
                   sep=""), 
             sep="\t",
             row.names=FALSE )
 
 return(OTU_table_time_series_normalized_amplicon)
}


#### gp23 ####

OTU_table_gp23 <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_gp23_pear.assembled_annotated_quality_length_filt_otus_95.00global_OTU_table.tsv", row.names = 1, header=TRUE)


no_libs_below <- 200
filtered_OTU_table_gp23 <- OTU_table_gp23[,colSums(OTU_table_gp23) > no_libs_below]
# ## might want to think about getting rid of VC1150 since it didn't really work and 1203. Go through and look at how these all could be filtered. If they are not removed they could bias the normalization

## adding in some special VCs
gp23_data <- c()
gp23_data <- as.character(Jericho_data$VC_number)

gp23_data <- c(gp23_data, "1256_ctab", "1255_large")
gp23_data <- as.data.frame(gp23_data)
names(gp23_data)[1] <- "VC_number"

get_project_data_gp23 <- get_project_samples(filtered_OTU_table_gp23,
                                             Library_metadata,
                                             Library_metadata$gp23_VC_number,
                                             gp23_data)

get_project_data_gp23$VC_number
dim(get_project_data_gp23)[2]
rowSums(get_project_data_gp23[,-dim(get_project_data_gp23)[2]])


get_project_data_gp23$VC_number

### deal with duplicated VCs
## so see that I need to decide between 1168 and 1168.1, 1256_ctab 
## so see that I need to decide between 1168 and 1168.1, 1256_ctab 
rowSums(get_project_data_gp23[get_project_data_gp23$VC_number=="1168",-dim(get_project_data_gp23)[2]])
min(get_project_data_gp23[get_project_data_gp23$VC_number=="1168",-dim(get_project_data_gp23)[2]])
max(get_project_data_gp23[get_project_data_gp23$VC_number=="1168",-dim(get_project_data_gp23)[2]])

rowSums(get_project_data_gp23[get_project_data_gp23$VC_number=="1168.1",-dim(get_project_data_gp23)[2]])
min(get_project_data_gp23[get_project_data_gp23$VC_number=="1168.1",-dim(get_project_data_gp23)[2]])
max(get_project_data_gp23[get_project_data_gp23$VC_number=="1168.1",-dim(get_project_data_gp23)[2]])

exclude <- "1168.1"
get_project_data_gp23 <- droplevels(subset(get_project_data_gp23, !(VC_number %in% exclude)))
get_project_data_gp23$VC_number

rowSums(get_project_data_gp23[get_project_data_gp23$VC_number=="1256_ctab",-dim(get_project_data_gp23)[2]])
min(get_project_data_gp23[get_project_data_gp23$VC_number=="1256_ctab",-dim(get_project_data_gp23)[2]])
max(get_project_data_gp23[get_project_data_gp23$VC_number=="1256_ctab",-dim(get_project_data_gp23)[2]])

rowSums(get_project_data_gp23[get_project_data_gp23$VC_number=="1256",-dim(get_project_data_gp23)[2]])
min(get_project_data_gp23[get_project_data_gp23$VC_number=="1256",-dim(get_project_data_gp23)[2]])
max(get_project_data_gp23[get_project_data_gp23$VC_number=="1256",-dim(get_project_data_gp23)[2]])

exclude <- "1256"
get_project_data_gp23 <- droplevels(subset(get_project_data_gp23, !(VC_number %in% exclude)))
get_project_data_gp23$VC_number


get_project_data_gp23$VC_number <- sub("1256_ctab", "1256", get_project_data_gp23$VC_number)


read_and_normalize_gp23 <- normalize_otu_table_by_project(get_project_data_gp23,
                                                          "gp23")



#### MPL ####

OTU_table_MPL <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_MPL_pear.assembled_lines_annotated_quality_length_filt_otus_95.00global_OTU_table.tsv", row.names = 1, header=TRUE)

no_libs_below <- 1000
filtered_OTU_table_MPL <- OTU_table_MPL[,colSums(OTU_table_MPL) > no_libs_below]
# 
get_project_data_MPL <- get_project_samples(filtered_OTU_table_MPL,
                                            Library_metadata, 
                                            Library_metadata$MPL_VC_number, 
                                            Jericho_data)


MPL_data <- c()
MPL_data <- as.character(Jericho_data$VC_number)

MPL_data <- c(MPL_data, "1256_ctab", "1255_large")
MPL_data <- as.data.frame(MPL_data)
names(MPL_data)[1] <- "VC_number"


get_project_data_MPL <- get_project_samples(filtered_OTU_table_MPL,
                                            Library_metadata, 
                                            Library_metadata$MPL_VC_number, 
                                            MPL_data)



rowSums(get_project_data_MPL[,-dim(get_project_data_MPL)[2]])
head(colSums(get_project_data_MPL[,-dim(get_project_data_MPL)[2]]))


get_project_data_MPL$VC_number
## nothing to do. 

rowSums(get_project_data_MPL[get_project_data_MPL$VC_number==1193,-dim(get_project_data_MPL)[2]])
min(get_project_data_MPL[get_project_data_MPL$VC_number==1193,-dim(get_project_data_MPL)[2]])
max(get_project_data_MPL[get_project_data_MPL$VC_number==1193,-dim(get_project_data_MPL)[2]])

rowSums(get_project_data_MPL[get_project_data_MPL$VC_number==1193.1,-dim(get_project_data_MPL)[2]])
min(get_project_data_MPL[get_project_data_MPL$VC_number==1193.1,-dim(get_project_data_MPL)[2]])
max(get_project_data_MPL[get_project_data_MPL$VC_number==1193.1,-dim(get_project_data_MPL)[2]])

exclude <- "1193.1"
get_project_data_MPL <- droplevels(subset(get_project_data_MPL, !(VC_number %in% exclude)))
get_project_data_MPL$VC_number

#test_normalize <- normalize_and_get_OTUs_occuring_more_than_1x_from_OTU_table(get_project_data_MPL)
# ## need to check with data sheets, but figuring out which VCs should be excluded. 
#MPL_not_part <- c("1252", "1259","1203", "1234", "1115", "1185", "1186", "1223")
# 
# get_project_data_MPL <- get_project_samples(OTU_table_MPL, Library_metadata, Library_metadata$MPL_VC_number, Jericho_data[!(Jericho_data$VC_number %in% MPL_not_part),])
# rowSums(get_project_data_MPL[,-dim(get_project_data_MPL)[2]])
# get_project_data_MPL$VC_number

read_and_normalize_MPL <- normalize_otu_table_by_project(get_project_data_MPL, 
                                                         "MPL")



### AVS ###

# OTU_table_AVS <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_AVS_cleaned_up_concat_clustal_aligned_edited_otus_90.00global_OTU_table.tsv", row.names = 1, header=TRUE)
# 
# no_libs_below <- 20
# filtered_OTU_table_AVS <- OTU_table_AVS[,colSums(OTU_table_AVS) > no_libs_below]
# 
# AVS_data <- c()
# AVS_data <- as.character(Jericho_data$VC_number)
# 
# #AVS_data <- append(AVS_data, "1256_ctab")
# AVS_data <- as.data.frame(AVS_data)
# names(AVS_data)[1] <- "VC_number"
# 
# 
# get_project_data_AVS <- get_project_samples(filtered_OTU_table_AVS,
#                                             Library_metadata,
#                                             Library_metadata$AVS_VC_number,
#                                             AVS_data)
# get_project_data_AVS$VC_number
# 
# frame_cols <- dim(get_project_data_AVS)[2]
# rowSums(get_project_data_AVS[,-frame_cols ])

# get_project_data_AVS1$VC_number
# ## need to check with data sheets, but figuring out which VCs should be excluded. 
# AVS_not_part <- c( "1159","1211", "1259","1235", "1176","1244","1252", "1115","1101", "1223","1193","1098" )
# 
# read_and_normalize_AVS <- normalize_otu_table_by_project(get_project_data_AVS,
#                                                          "AVS_concat")
# 


#  OTU_table_AVS2 <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_AVS2_cleaned_up_otus_90.00global_OTU_table.tsv", row.names = 1, header=TRUE)
# # 
#  no_libs_below <- 200
#  filtered_OTU_table_AVS2 <- OTU_table_AVS2[,colSums(OTU_table_AVS2) > no_libs_below]
# # 
# # 
#  colSums(filtered_OTU_table_AVS2 )
#   get_project_data_AVS2 <- get_project_samples(filtered_OTU_table_AVS2, Library_metadata, Library_metadata$AVS_VC_number, Jericho_data)
# 
#  read_and_normalize_AVS2 <- normalize_otu_table_by_project(filtered_OTU_table_AVS2, 
#                                                            Library_metadata,
#                                                            Library_metadata$AVS_VC_number, 
#                                                            Jericho_data,
#                                                            "AVS2")




OTU_table_AVS1 <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_AVS1_cleaned_up_clustal_aligned_edited_otus_90.00global_OTU_table.tsv", row.names = 1, header=TRUE)

no_libs_below <- 50
filtered_OTU_table_AVS1 <- OTU_table_AVS1[,colSums(OTU_table_AVS1) > no_libs_below]
colSums(filtered_OTU_table_AVS1)

AVS1_data <- c()
AVS1_data <- as.character(Jericho_data$VC_number)

AVS1_data <- append(AVS1_data, "1256_ctab")
AVS1_data <- as.data.frame(AVS1_data)
names(AVS1_data)[1] <- "VC_number"

get_project_data_AVS1 <- get_project_samples(filtered_OTU_table_AVS1,
                                             Library_metadata,
                                             Library_metadata$AVS_VC_number,
                                             AVS1_data)

rowSums(get_project_data_AVS1[,-dim(get_project_data_AVS1)[2]])
get_project_data_AVS1$VC_number





#AVS_not_part <- c( "1157", "1235","1266")
#get_project_data_AVS1$VC_number <- sub("1256_ctab", "1256", get_project_data_AVS1$VC_number)
#get_project_data_AVS1$VC_number

read_and_normalize_AVS1 <- normalize_otu_table_by_project(get_project_data_AVS1,
                                                          "AVS1")

#### 18S ####
# 

OTU_table_18s <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ng_otus_97.00global_OTU_table.tsv", row.names = 1, header=TRUE)
colSums(OTU_table_18s)

no_libs_below <- 2000
filtered_OTU_table_18s <- OTU_table_18s[,colSums(OTU_table_18s) > no_libs_below]

S18_data <- c()
S18_data <- as.character(Jericho_data$VC_number)

S18_data <- c(S18_data, "1256_ctab", "1255_large","1223_redo","1225_redo")

S18_data <- as.data.frame(S18_data)
names(S18_data)[1] <- "VC_number"

get_project_data_18s <- get_project_samples(filtered_OTU_table_18s,
                                            Library_metadata,
                                            Library_metadata$X18s_VC_number,
                                            S18_data)

frame_cols <- dim(get_project_data_18s)[2]
rowSums(get_project_data_18s[,-frame_cols ])
get_project_data_18s$VC_number

### Dealing with the duplicated VCs #######################

get_project_data_18s$VC_number
## so see that I need to decide between 1244 and 1244.1, 1225_redo and 1225 and change 1255_large to 1255

rowSums(get_project_data_18s[get_project_data_18s$VC_number=="1244",-dim(get_project_data_18s)[2]])
min(get_project_data_18s[get_project_data_18s$VC_number=="1244",-dim(get_project_data_18s)[2]])
max(get_project_data_18s[get_project_data_18s$VC_number=="1244",-dim(get_project_data_18s)[2]])

rowSums(get_project_data_18s[get_project_data_18s$VC_number=="1244.1",-dim(get_project_data_18s)[2]])
min(get_project_data_18s[get_project_data_18s$VC_number=="1244.1",-dim(get_project_data_18s)[2]])
max(get_project_data_18s[get_project_data_18s$VC_number=="1244.1",-dim(get_project_data_18s)[2]])

## they look very similar so let's go with 1244 and delete 1244.1
exclude <- "1244"
get_project_data_18s <- subset(get_project_data_18s, !(VC_number %in% exclude))
get_project_data_18s$VC_number

## rename 1244.1 to 1244
#levels(get_project_data_18s$VC_number) <- sub("1244.1", "1244", levels(get_project_data_18s$VC_number))
get_project_data_18s$VC_number <- sub("1244.1", "1244", get_project_data_18s$VC_number)
get_project_data_18s$VC_number

rowSums(get_project_data_18s[get_project_data_18s$VC_number=="1225",-dim(get_project_data_18s)[2]])
min(get_project_data_18s[get_project_data_18s$VC_number=="1225",-dim(get_project_data_18s)[2]])
max(get_project_data_18s[get_project_data_18s$VC_number=="1225",-dim(get_project_data_18s)[2]])

rowSums(get_project_data_18s[get_project_data_18s$VC_number=="1225_redo",-dim(get_project_data_18s)[2]])
min(get_project_data_18s[get_project_data_18s$VC_number=="1225_redo",-dim(get_project_data_18s)[2]])
max(get_project_data_18s[get_project_data_18s$VC_number=="1225_redo",-dim(get_project_data_18s)[2]])

## also look similar so I will stick with 1225.
exclude <- "1225"
get_project_data_18s <- droplevels(subset(get_project_data_18s, !(VC_number %in% exclude)))
get_project_data_18s$VC_number
## rename 1225_redo
#levels(get_project_data_18s$VC_number) <- sub("1225_redo", "1225", levels(get_project_data_18s$VC_number))
get_project_data_18s$VC_number <- sub("1225_redo", "1225", get_project_data_18s$VC_number)
get_project_data_18s$VC_number


## rename 1255
get_project_data_18s$VC_number <- sub("1255_large", "1255", get_project_data_18s$VC_number)
get_project_data_18s$VC_number

read_and_normalize_18s <- normalize_otu_table_by_project(get_project_data_18s,
                                                         "18S")



#### 16S ####

# OTU_table_16s <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_forwardgood.good.pcr.ng_otus_97.00global_OTU_table.tsv", row.names = 1, header=TRUE)
# 
# 
# colSums(OTU_table_16s)
# 
# no_libs_below <- 2000
# filtered_OTU_table_16s <- OTU_table_16s[,colSums(OTU_table_16s) > no_libs_below]
# 
# S16_data <- c()
# S16_data <- as.character(Jericho_data$VC_number)
# 
# S16_data <- c(S16_data,"1255_large","1256_ctab" )
# S16_data <- as.data.frame(S16_data)
# names(S16_data)[1] <- "VC_number"
# 
# get_project_data_16s <- get_project_samples(filtered_OTU_table_16s,
#                                             Library_metadata,
#                                             Library_metadata$X16s_VC_number,
#                                             S16_data)
# 
# frame_cols <- dim(get_project_data_16s)[2]
# rowSums(get_project_data_16s[,-frame_cols ])
# get_project_data_16s$VC_number
# 
# read_and_normalize_16s <- normalize_otu_table_by_project(get_project_data_16s,
#                                                          "16S_concat")

### 16s R1

OTU_table_16s_R1 <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_forwardgood.good.pcr.ng_otus_97.00global_OTU_table.tsv", row.names = 1, header=TRUE)

colSums(OTU_table_16s_R1)

no_libs_below <- 2000
filtered_OTU_table_16s_R1 <- OTU_table_16s_R1[,colSums(OTU_table_16s_R1) > no_libs_below]

S16_data <- c()
S16_data <- as.character(Jericho_data$VC_number)
S16_data <- c(S16_data,"1255_large","1256_ctab","1223_redo","1225_redo" )

S16_data <- as.data.frame(S16_data)
names(S16_data)[1] <- "VC_number"


get_project_data_16s_R1 <- get_project_samples(filtered_OTU_table_16s_R1,
                                               Library_metadata,
                                               Library_metadata$X16s_VC_number,
                                               S16_data)

frame_cols <- dim(get_project_data_16s_R1)[2]
rowSums(get_project_data_16s_R1[,-frame_cols ])



get_project_data_16s_R1$VC_number

## so see that I need to decide between 1244 and 1244.1, 1225_redo and 1225 and change 1255_large to 1255
rowSums(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1211",-dim(get_project_data_16s_R1)[2]])
min(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1211",-dim(get_project_data_16s_R1)[2]])
max(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1211",-dim(get_project_data_16s_R1)[2]])

rowSums(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1211.1",-dim(get_project_data_16s_R1)[2]])
min(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1211.1",-dim(get_project_data_16s_R1)[2]])
max(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1211.1",-dim(get_project_data_16s_R1)[2]])

## they look very similar so let's go with 1244 and delete 1244.1
exclude <- "1211.1"
get_project_data_16s_R1 <- subset(get_project_data_16s_R1, !(VC_number %in% exclude))
get_project_data_16s_R1$VC_number

rowSums(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1255",-dim(get_project_data_16s_R1)[2]])
min(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1255",-dim(get_project_data_16s_R1)[2]])
max(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1255",-dim(get_project_data_16s_R1)[2]])

rowSums(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1255_large",-dim(get_project_data_16s_R1)[2]])
min(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1255_large",-dim(get_project_data_16s_R1)[2]])
max(get_project_data_16s_R1[get_project_data_16s_R1$VC_number=="1255_large",-dim(get_project_data_16s_R1)[2]])

## 
exclude <- "1255_large"
get_project_data_16s_R1 <- subset(get_project_data_16s_R1, !(VC_number %in% exclude))
get_project_data_16s_R1$VC_number

read_and_normalize_16s_R1 <- normalize_otu_table_by_project(get_project_data_16s_R1,
                                                            "16S_R1")

# 
# OTU_table_16s_R2 <- read.delim("../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_R2_annotated_trimmed_paired.good.good.ng_otus_97.00global_OTU_table.tsv", row.names = 1, header=TRUE)
# 
# colSums(OTU_table_16s_R2)
# 
# no_libs_below <- 2000
# filtered_OTU_table_16s_R2 <- OTU_table_16s_R2[,colSums(OTU_table_16s_R2) > no_libs_below]
# 
# S16_data <- c()
# S16_data <- as.character(Jericho_data$VC_number)
# 
# non_numeric_VCs <- c("1223_redo","1225_redo","1255_large")
# S16_data <- append(S16_data,non_numeric_VCs )
# S16_data <- as.data.frame(S16_data)
# names(S16_data)[1] <- "VC_number"
# 
# get_project_data_16s_R2 <- get_project_samples(filtered_OTU_table_16s_R2,
#                                                Library_metadata,
#                                                Library_metadata$X16s_VC_number,
#                                                Jericho_data)
# 
# frame_cols <- dim(get_project_data_16s_R2)[2]
# rowSums(get_project_data_16s_R2[,-frame_cols ])
# 
# 
# read_and_normalize_16s_R2 <- normalize_otu_table_by_project(filtered_OTU_table_16s_R2,
#                                                             Library_metadata,
#                                                             Library_metadata$X16s_VC_number,
#                                                             S16_data,
#                                                             "16S_R2")








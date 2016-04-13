## get top 100 and maybe even top 500 otus
library(reshape2)
library(plyr)
## while reprocessing it is missing the VC number
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv",                                  row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv", 
                                   row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", 
                                  row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
#normalized_AVS_R1_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R2.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

find_top100_OTUs <- function (normalized_otus) {
 long_otus <- melt(normalized_otus)
 summarised_by_OTU <- ddply(long_otus,~variable,summarise,relative_abundance=sum(value))
 OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
 vector_top_100_OTUs <- droplevels(head(OTU_sorted_by_rel_abun, n=100))
 top_100_OTUs <- droplevels(subset(normalized_otus, select=colnames(normalized_otus) %in% vector_top_100_OTUs$variable))
 return(top_100_OTUs)
}


find_top10_OTUs <- function (normalized_otus) {
 long_otus <- melt(normalized_otus)
 summarised_by_OTU <- ddply(long_otus,~variable,summarise,relative_abundance=sum(value))
 OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
 vector_top_10_OTUs <- droplevels(head(OTU_sorted_by_rel_abun, n=10))
 top_10_OTUs <- droplevels(subset(normalized_otus, select=colnames(normalized_otus) %in% vector_top_10_OTUs$variable))
 return(top_10_OTUs)
} 


MPL_top_100 <- find_top100_OTUs(normalized_MPL_OTUs)
write.table(MPL_top_100,"../data/OTU_table_Jericho_time_series_MPL_normalized_top_100.tsv",
            sep="\t",
            col.names=NA)

AVS_top_100 <- find_top100_OTUs(normalized_AVS_OTUs)
write.table(AVS_top_100,"../data/OTU_table_Jericho_time_series_AVS1_normalized_top_100.tsv",
            sep="\t",
            col.names=NA)
gp23_top_100 <- find_top100_OTUs(normalized_gp23_OTUs)
write.table(gp23_top_100,"../data/OTU_table_Jericho_time_series_gp23_normalized_top_100.tsv",
            sep="\t",
            col.names=NA)
S18_top_100 <- find_top100_OTUs(normalized_18s_OTUs)
write.table(S18_top_100,"../data/OTU_table_Jericho_time_series_18s_normalized_top_100.tsv",
            sep="\t",
            col.names=NA)
S16_top_100 <- find_top100_OTUs(normalized_16s_OTUs)
write.table(S16_top_100,"../data/OTU_table_Jericho_time_series_16s_R1_normalized_top_100.tsv",
            sep="\t",
            col.names=NA)


## top 10


MPL_top_10 <- find_top10_OTUs(normalized_MPL_OTUs)
write.table(MPL_top_10,"../data/OTU_table_Jericho_time_series_MPL_normalized_top_10.tsv",
            sep="\t",
            col.names=NA)

AVS_top_10 <- find_top10_OTUs(normalized_AVS_OTUs)
write.table(AVS_top_10,"../data/OTU_table_Jericho_time_series_AVS1_normalized_top_10.tsv",
            sep="\t",
            col.names=NA)
gp23_top_10 <- find_top10_OTUs(normalized_gp23_OTUs)
write.table(gp23_top_10,"../data/OTU_table_Jericho_time_series_gp23_normalized_top_10.tsv",
            sep="\t",
            col.names=NA)
S18_top_10 <- find_top10_OTUs(normalized_18s_OTUs)
write.table(S18_top_10,"../data/OTU_table_Jericho_time_series_18s_normalized_top_10.tsv",
            sep="\t",
            col.names=NA)
S16_top_10 <- find_top10_OTUs(normalized_16s_OTUs)
write.table(S16_top_10,"../data/OTU_table_Jericho_time_series_16s_R1_normalized_top_10.tsv",
            sep="\t",
            col.names=NA)


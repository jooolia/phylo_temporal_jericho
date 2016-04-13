## Betadiversity ana;lyses

library(vegan)
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv", row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
#normalized_AVS_R1_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R2.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S.tsv", row.names="VC_number") 


## Betadiversity

## ok maybe that is it right now. 
## 
shared_16s <- betadiver(normalized_16s_OTUs, triangular = FALSE)
#this is what I wanted....how many are shared between pools!
#so now I need to sum and then average over time....need some sort aaply here....eep.
shared_16s$a #these are the shared OTUs between each pool
str(shared_16s)
#so this is overall. Would I want to do this in linear time? yes
mean(shared_16s$a)
min(shared_16s$a)
max(shared_16s$a)

shared_18s <- betadiver(normalized_18s_OTUs, triangular = FALSE)

shared_18s$a #these are the shared OTUs between each pool
str(shared_18s)
mean(shared_18s$a)
min(shared_18s$a)
max(shared_18s$a)

shared_AVS_R1 <- betadiver(normalized_AVS_R1_OTUs, triangular = FALSE)
shared_AVS_R1$a #these are the shared OTUs between each pool
str(shared_AVS_R1)
mean(shared_AVS_R1$a)
min(shared_AVS_R1$a)
max(shared_AVS_R1$a)


shared_MPL <- betadiver(normalized_MPL_OTUs, triangular = FALSE)

shared_MPL$a #these are the shared OTUs between each pool
str(shared_MPL$a)
#so this is overall. Would I want to do this in linear time? yes
mean(shared_MPL$a)
min(shared_MPL$a)
max(shared_MPL$a)


shared_gp23 <- betadiver(normalized_gp23_OTUs, triangular = FALSE)
shared_gp23$a #these are the shared OTUs between each pool
str(shared_gp23$a)
#so this is overall. Would I want to do this in linear time? yes
mean(shared_gp23$a)
min(shared_gp23$a)
max(shared_gp23$a)
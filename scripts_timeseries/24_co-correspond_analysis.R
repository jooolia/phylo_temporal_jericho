## perform co-correspondence analysis

library(cocorresp)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

## symmetric CoCA
data(beetles)
## log transform the bettle data
beetles <- log(beetles + 1)
data(plants)
## fit the model
bp.sym <- coca(beetles ~ ., data = plants, method = "symmetric")
bp.sym
summary(bp.sym)
plot(bp.sym)

## try with 18s and MPL
union_vcs_MPL_and_18s <- intersect(row.names(normalized_MPL_OTUs), row.names(normalized_18s_OTUs))

Jericho_MPL <- droplevels(subset(normalized_MPL_OTUs, row.names(normalized_MPL_OTUs) %in% union_vcs_MPL_and_18s))

Jericho_18s <- droplevels(subset(normalized_18s_OTUs, row.names(normalized_18s_OTUs) %in% union_vcs_MPL_and_18s))

MPL_18s.sym <- coca(Jericho_MPL ~ ., data = Jericho_18s, method = "symmetric")
summary(MPL_18s.sym )
plot(MPL_18s.sym )

## try with AVS and 18s
union_vcs_AVS_and_18s <- intersect(row.names(normalized_AVS_OTUs), row.names(normalized_18s_OTUs))

Jericho_AVS <- droplevels(subset(normalized_AVS_OTUs, row.names(normalized_AVS_OTUs) %in% union_vcs_AVS_and_18s))
Jericho_18s <- droplevels(subset(normalized_18s_OTUs, row.names(normalized_18s_OTUs) %in% union_vcs_AVS_and_18s))

AVS_18s.sym <- coca(Jericho_AVS ~ ., data = Jericho_18s, method = "symmetric")
head(summary(AVS_18s.sym ))
plot(AVS_18s.sym )



## try with gp23 and 16s
union_vcs_gp23_and_16s <- intersect(row.names(normalized_gp23_OTUs), row.names(normalized_16s_OTUs))

Jericho_gp23 <- droplevels(subset(normalized_gp23_OTUs, row.names(normalized_gp23_OTUs) %in% union_vcs_gp23_and_16s))

Jericho_16s <- droplevels(subset(normalized_16s_OTUs, row.names(normalized_16s_OTUs) %in% union_vcs_gp23_and_16s))

gp23_16s.sym <- coca(Jericho_gp23 ~ ., data = Jericho_16s, method = "symmetric")
summary(gp23_16s.sym )
plot(gp23_16s.sym )


union_vcs_gp23_and_18s <- intersect(row.names(normalized_gp23_OTUs), row.names(normalized_18s_OTUs))

Jericho_gp23 <- droplevels(subset(normalized_gp23_OTUs, row.names(normalized_gp23_OTUs) %in% union_vcs_gp23_and_18s))

Jericho_18s <- droplevels(subset(normalized_18s_OTUs, row.names(normalized_18s_OTUs) %in% union_vcs_gp23_and_18s))


gp23_18s.sym <- coca(Jericho_gp23 ~ ., data = Jericho_18s, method = "symmetric")
summary(gp23_18s.sym )
plot(gp23_18s.sym )

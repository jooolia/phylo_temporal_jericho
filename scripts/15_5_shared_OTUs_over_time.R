## shared OTUs over time

library(ggplot2)
library(vegan)
library(scales)
library(plyr)
library(reshape2)
library(dplyr)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)
## Summarise 

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv",row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",  row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
normalized_AVS_R1_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS2.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
## Reformat date
# %d is day as a number, %b is abreviated month in words, %y is 2digit year
Jericho_data$Date <- as.Date(Jericho_data$Date)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)


args <- commandArgs(TRUE)
inputFile <- args[1]

## test to see if input file is given, so I can decide whether to use this argument or the orginal one. 
if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
 path_colour <- "black"
 line_colour <- "black"
 point_colour <- "black"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
  point_colour <- "white"
 }
}

# adding in the seasons and the spring bloom for the ggplots
season_line <- geom_vline(xintercept = as.numeric(c(as.Date("2010-03-22"), 
                                                    as.Date("2010-06-22"), 
                                                    as.Date("2010-09-22"),
                                                    as.Date("2010-12-22"),
                                                    as.Date("2011-03-22"),
                                                    as.Date("2011-06-22"))),
                          colour="grey",
                          size=1.5)
spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=1)
date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")))

count_16s_otus <- adply(t(normalized_16s_OTUs), 2,function(x)sum(x>0))
count_18s_otus <-adply(t(normalized_18s_OTUs), 2, function(x)sum(x>0))
count_gp23_otus <-adply(t(normalized_gp23_OTUs), 2, function(x)sum(x>0))
count_MPL_otus <-adply(t(normalized_MPL_OTUs), 2, function(x)sum(x>0))
count_AVS1_otus <-adply(t(normalized_AVS_R1_OTUs), 2, function(x)sum(x>0))
count_AVS2_otus <-adply(t(normalized_AVS_R2_OTUs), 2, function(x)sum(x>0))
count_AVS_concat_otus <-adply(t(normalized_AVS_OTUs), 2, function(x)sum(x>0))


shared_16s <- betadiver(normalized_16s_OTUs, triangular = FALSE)
#this is what I wanted....how many are shared between pools!
#so now I need to sum and then average over time....need some sort aaply here....eep.
#shared_16s$a #these are the shared OTUs between each pool
shared_16s_melt <- melt(as.matrix(shared_16s$a), value.name="shared_16s")
each_time_16s_melt <-  melt(as.matrix(shared_16s$a + shared_16s$b +shared_16s$c), value.name="total_each_time_16s")

shared_18s <- betadiver(normalized_18s_OTUs, triangular = FALSE)

shared_18s_melt <- melt(as.matrix(shared_18s$a), value.name="shared_18s")
each_time_18s_melt <-  melt(as.matrix(shared_18s$a + shared_18s$b +shared_18s$c), value.name="total_each_time_18s")

# shared_AVS_R1 <- betadiver(normalized_AVS_R1_OTUs, triangular = FALSE)
#shared_AVS_melt <- melt(as.matrix(shared_AVS$a), value.name="shared_AVS")
#each_time_AVS_melt <-  melt(as.matrix(shared_AVS$a + shared_AVs$b +shared_AVS$c), value.name="total_each_time_AVS")

shared_MPL <- betadiver(normalized_MPL_OTUs, triangular = FALSE)
#shared_MPL$a #these are the shared OTUs between each pool
shared_MPL_melt <- melt(as.matrix(shared_MPL$a), value.name="shared_MPL")
each_time_MPL_melt <-  melt(as.matrix(shared_MPL$a + shared_MPL$b +shared_MPL$c), value.name="total_each_time_MPL")

shared_gp23 <- betadiver(normalized_gp23_OTUs, triangular = FALSE)
#shared_gp23$a #these are the shared OTUs between each pool
shared_gp23_melt <- melt(as.matrix(shared_gp23$a), value.name="shared_gp23")
each_time_gp23_melt <-  melt(as.matrix(shared_gp23$a + shared_gp23$b +shared_gp23$c), value.name="total_each_time_gp23")

shared_otus <- Reduce(function(x, y) merge(x, y, 
                                           by=c("Var1", "Var2"),
                                           all.x=TRUE), list(shared_16s_melt,
                                                             each_time_16s_melt,
                                                             shared_18s_melt,
                                                             each_time_18s_melt,
                                                             shared_gp23_melt,
                                                             each_time_gp23_melt,
                                                             shared_MPL_melt,
                                                             each_time_MPL_melt
                                                             #,
                                                             #count_AVS1_otus,
                                                             #count_AVS2_otus,
                                                             #count_AVS_concat_otus
                                           ))
## want to remove the rows where it is the same comparions
shared_otus <- filter(shared_otus,(Var1!=Var2))


shared_otus$Date1 <- Jericho_data$Date[match(shared_otus$Var1, Jericho_data$VC)]
shared_otus$Date2 <- Jericho_data$Date[match(shared_otus$Var2, Jericho_data$VC)]


community_over_time <- data.frame(Date1=as.Date(character()),
                                  Date2=as.Date(character()),
                                  shared_16s=numeric(),
                                  total_each_time_16s=numeric(),
                                  shared_18s=numeric(),
                                  total_each_time_18s=numeric(),
                                  shared_gp23=numeric(),
                                  total_each_time_gp23=numeric(),
                                  shared_MPL=numeric(),
                                  total_each_time_MPL=numeric()
                                  # "OTUs_AVSshared_between_dates"
)

## or just a loop through all the dates using the indexing to get at them
for (sample_date in as.character(unique(sort(shared_otus$Date1)))){
 print(sample_date)
 date_time_a <- as.Date(sample_date)
 index_time_a <- which(as.character(unique(sort(shared_otus$Date1)))==sample_date)
 
 ## use the next row in the sorted matrix!
 index_time_b <- index_time_a + 1
 date_time_b <- as.Date(as.character(unique(sort(shared_otus$Date1)))[index_time_b])
 ## at time 2 how does it compare to the last time
 filter(shared_otus, Date1==sample_date & Date2==date_time_b)
 otus_shared_between_dates <- filter(shared_otus, Date1==sample_date & Date2==date_time_b)
 
 to_add <- otus_shared_between_dates[,c("Date1",
                                        "Date2",
                                        "shared_16s",
                                        "total_each_time_16s",
                                        "shared_18s",
                                        "total_each_time_18s",
                                        "shared_gp23",
                                        "total_each_time_gp23",
                                        "shared_MPL",
                                        "total_each_time_MPL")]
 community_over_time <- rbind(community_over_time,to_add  )
}


plot(community_over_time$Date1, community_over_time$shared_18s)
plot(community_over_time$Date1, community_over_time$shared_16s)

melted_community_over_time <- melt(community_over_time, id.vars=c("Date1", "Date2"))
ggplot(melted_community_over_time,
       aes(x=Date1, y=value))+
 geom_line(data=subset(melted_community_over_time, grepl("shared", variable)),aes(group=variable,
                                                                                  colour=variable),
           size=2)+
 geom_line(data=subset(melted_community_over_time, grepl("total", variable)),aes(group=variable,
                                                                                 colour=variable),
           size=2)

ggplot(melted_community_over_time,
       aes(x=Date1, y=value))+
 geom_point()+
 facet_wrap(~variable)



ggplot(subset(melted_community_over_time, grepl("18s", variable)),
       aes(x=Date1, y=value))+
 geom_line(aes(group=variable,
               colour=variable),
           size=2)

ggplot(subset(melted_community_over_time, grepl("16s", variable)),
       aes(x=Date1, y=value))+
 geom_line(aes(group=variable,
               colour=variable),
           size=2)

ggplot(subset(melted_community_over_time, grepl("MPL", variable)),
       aes(x=Date1, y=value))+
 geom_line(aes(group=variable,
               colour=variable),
           size=2) 
ggplot(subset(melted_community_over_time, grepl("gp23", variable)),
       aes(x=Date1, y=value))+
 geom_line(aes(group=variable,
               colour=variable),
           size=2) 



community_percent_over_time <- transmute(community_over_time, 
                                         Date1,
                                         Date2,
                                         percent_16s_shared= round(shared_16s/total_each_time_16s, digits=2) *100,
                                         percent_18s_shared= round(shared_18s/total_each_time_18s, digits=2) *100,
                                         percent_gp23_shared= round(shared_gp23/total_each_time_gp23, digits=2) *100,
                                         #percent_AVS_shared= round(shared_AVS/total_each_time_AVS,
                                         percent_MPL_shared= round(shared_MPL/total_each_time_MPL, digits=2) *100)

#kable(community_percent_over_time)
#pandoc.table(community_percent_over_time, style="multiline", justify="left")

plot(community_percent_over_time$Date1,community_percent_over_time$percent_16s_shared)
melted_community_percent_over_time <- melt(community_percent_over_time, id.vars=c("Date1", "Date2"))

pdf("../figures/shared_community_percent_over_time.pdf")
ggplot(melted_community_percent_over_time,
       aes(x=Date1, y=value))+
 season_line +
 spring_bloom_line+
 date_scaling+
 theme_JAG_presentation()+
 geom_line(aes(group=variable,
               colour=variable),
           size=2)
dev.off()


## WOrk with time series OTUs

## To do: add in Chao1 estimator?

## @knitr read_in_otu_tables_with_plots
## @knitr gp_23_all
library(ggplot2)
library(plyr)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(scales)
library(vegan)
library(phyloseq)
library(cowplot)

args <- commandArgs(TRUE)
inputFile <- args[1]

## test to see if input file is given, so I can decide whether to use this argument or the orginal one
if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
 path_colour <- "black"
 line_colour <- "black"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
 }
}

season_line <- geom_vline(xintercept = as.numeric(c(as.Date("2010-03-22"), 
                                                    as.Date("2010-06-22"), 
                                                    as.Date("2010-09-22"),
                                                    as.Date("2010-12-22"),
                                                    as.Date("2011-03-22"),
                                                    as.Date("2011-06-22"))),
                          colour="grey90",
                          size=1.5)
spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=1)
HAKA_bloom <- annotate("rect",xmin = as.Date("2011-06-20"),
                       xmax = as.Date("2011-07-06"), 
                       ymin = -Inf,
                       ymax = Inf,
                       alpha = 0.4,
                       fill = "burlywood3",
                       colour=NA)

date_scaling <-   scale_x_date(limits = c(as.Date("2011-06-20"), as.Date("2011-07-06")))

## Need to edit and make sure the font, type are consistent
## To-do:

## read in normalized OTU tables

## while reprocessing it is missing the VC number
normalized_18s_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_18S.tsv",
                                  row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_gp23.tsv", 
                                   row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_MPL.tsv", 
                                  row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
#normalized_AVS_R1_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_AVS_R1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_AVS_R2.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_16S_R1.tsv", row.names="VC_number") 



### Get relevant VC numbers
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)

## first rename the column names in the OTU table to something better
add_date_to_OTU_table_rtrn_long <- function (otu_table, 
                                             Jericho_data) {
 otu_table$VC_number = rownames(otu_table)
 long_otus <- melt(otu_table, 
                   id="VC_number", 
                   variable.name="OTUid")
 ## Add in dates:
 long_otus$Date <- Jericho_data$Date[match(long_otus$VC_number, 
                                           Jericho_data$VC_number)]
 return(long_otus)
}

## Read in OTU table and get long format with Date
long_gp23_otus <- add_date_to_OTU_table_rtrn_long(normalized_gp23_OTUs,
                                                  Jericho_data)

long_18s_otus <- add_date_to_OTU_table_rtrn_long(normalized_18s_OTUs,
                                                 Jericho_data)

long_AVS_otus <- add_date_to_OTU_table_rtrn_long(normalized_AVS_OTUs,
                                                 Jericho_data)

long_MPL_otus <- add_date_to_OTU_table_rtrn_long(normalized_MPL_OTUs,
                                                 Jericho_data)

long_16s_otus <- add_date_to_OTU_table_rtrn_long(normalized_16s_OTUs,
                                                 Jericho_data)


## how many OTUs total?
Number_of_OTUs_gp23 <- n_distinct(long_gp23_otus$OTUid)
Number_of_OTUs_18s <- n_distinct(long_18s_otus$OTUid)
Number_of_OTUs_MPL <- n_distinct(long_MPL_otus$OTUid)
Number_of_OTUs_AVS <- n_distinct(long_AVS_otus$OTUid)
Number_of_OTUs_16s <- n_distinct(long_16s_otus$OTUid)

## pull out the counts of the number of OTUs-richness over time

richness_of_otus_over_time <- function (long_otus) {
 summarise_richness_otus <- long_otus %>%
  filter(value > 0) %>%
  ddply(~Date,summarise,richness=length(OTUid))
}


evenness_of_otus_over_time <- function (long_otus) {
 #long_otus <- long_gp23_otus
 summarise_richness_otus <- long_otus %>%
  filter(value > 0) %>%
  ddply(~Date,summarise,richness=length(OTUid))
 
 OTUs_over_time <- dcast(long_otus, Date~OTUid)
 rownames(OTUs_over_time) <- OTUs_over_time$Date
 OTUs_over_time <- OTUs_over_time[-1] # remove date after made rownames
 
 H <- diversity(OTUs_over_time)
 ## Species richness (S) and Pielou's evenness (J):
 J <- H/log(summarise_richness_otus$richness)
 J <- as.data.frame(J)
 names(J)
 J$Date <- as.Date(row.names(J))
 return(J)
}


par(mar=c(0,0,0,0))
plot_richness_otus_over_time <- function (richness_over_time_data, total_count_of_OTUs, primer, clustering_method) {
 
 sampling_points <- richness_over_time_data$Date
 p <- ggplot(richness_over_time_data, aes(x=Date, y=richness))
 p + 
  HAKA_bloom +
  season_line +
  spring_bloom_line+
  geom_line(colour=line_colour )+
  date_scaling +
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16))+
  ggtitle(primer)
}


plot_evenness_otus_over_time <- function (evenness_over_time_data, total_count_of_OTUs, primer, clustering_method) {
 
 sampling_points <- evenness_over_time_data$Date
 p <- ggplot(evenness_over_time_data, aes(x=Date, y=J))
 p + 
  HAKA_bloom +
  season_line +
  spring_bloom_line+
  geom_line(colour=line_colour )+
  date_scaling +
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  ylim(c(-0.02, 1.0))+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16))+
  ggtitle(primer)
}



plot_gp23_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_gp23_otus),
                                                             Number_of_OTUs_gp23,
                                                             "gp23",
                                                             "usearch")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)

plot_gp23_richness_over_time


plot_18S_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_18s_otus),
                                                            n_distinct(long_18s_otus$OTUid),
                                                            "18S",
                                                            "usearch")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_18S_richness_over_time

plot_16S_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_16s_otus),
                                                            n_distinct(long_16s_otus$OTUid),
                                                            "16S",
                                                            "usearch")
plot_16S_richness_over_time

plot_AVS_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_AVS_otus),
                                                            n_distinct(long_AVS_R1_otus$OTUid),
                                                            "AVS R1",
                                                            "usearch")+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_AVS_richness_over_time

plot_MPL_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_MPL_otus),
                                                            n_distinct(long_MPL_otus$OTUid),
                                                            "MPL",
                                                            "usearch")+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_MPL_richness_over_time


pdf("../figures/time_series_richness_all_amplicons_line_graphs%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(plot_MPL_richness_over_time,
             plot_AVS_richness_over_time,
             plot_gp23_richness_over_time,
             plot_18S_richness_over_time, 
             plot_16S_richness_over_time,
             ncol=1)
dev.off()

## Test out density plotting ####

## pick some colours to use:
library(RColorBrewer)
colours_for_plots <- brewer.pal(n = 5, name = 'Dark2')

plot_richness_density_over_time <- function (richness_over_time_data, total_count_of_OTUs, primer, clustering_method, fillcolour) {
 sampling_points <- richness_over_time_data$Date
 p <- ggplot(richness_over_time_data, aes(x=Date, y=richness))
 p +
  HAKA_bloom +
  season_line +
  spring_bloom_line+
  geom_density(stat="identity", colour=line_colour, fill=fillcolour)+
  date_scaling +
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16))+
  ggtitle(primer)
}


plot_evenness_density_over_time <- function (evenness_over_time_data, total_count_of_OTUs, primer, clustering_method, fillcolour) {
 
 sampling_points <- evenness_over_time_data$Date
 p <- ggplot(evenness_over_time_data, aes(x=Date, y=J))
 p + 
  HAKA_bloom +
  season_line +
  spring_bloom_line+
  geom_density(stat="identity", colour=line_colour, fill=fillcolour)+
  date_scaling +
  ylim(c(-0.02, 1.0))+
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16))+
  ggtitle(primer)
}


plot_MPL_richness_over_time <- plot_richness_density_over_time(richness_of_otus_over_time(long_MPL_otus),
                                                               n_distinct(long_MPL_otus$OTUid),
                                                               "MPL",
                                                               "usearch",
                                                               colours_for_plots[1])+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_MPL_richness_over_time


plot_gp23_richness_over_time <- plot_richness_density_over_time(richness_of_otus_over_time(long_gp23_otus),
                                                                Number_of_OTUs_gp23,
                                                                "gp23",
                                                                "usearch",
                                                                colours_for_plots[2])+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)

plot_gp23_richness_over_time


plot_18S_richness_over_time <- plot_richness_density_over_time(richness_of_otus_over_time(long_18s_otus),
                                                               n_distinct(long_18s_otus$OTUid),
                                                               "18S",
                                                               "usearch",
                                                               colours_for_plots[3])+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_18S_richness_over_time

plot_16S_richness_over_time <- plot_richness_density_over_time(richness_of_otus_over_time(long_16s_otus),
                                                               n_distinct(long_16s_otus$OTUid),
                                                               "16S",
                                                               "usearch",
                                                               colours_for_plots[4])
plot_16S_richness_over_time

plot_AVS_richness_over_time <- plot_richness_density_over_time(richness_of_otus_over_time(long_AVS_otus),
                                                               n_distinct(long_AVS_R1_otus$OTUid),
                                                               "AVS R1",
                                                               "usearch",
                                                               colours_for_plots[5])+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_AVS_richness_over_time


pdf("../figures/time_series_richness_all_amplicons_density%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(
 #plot_MPL_richness_over_time,
 #           plot_AVS_richness_over_time,
 plot_gp23_richness_over_time,
 plot_18S_richness_over_time, 
 plot_16S_richness_over_time,
 ncol=1)
dev.off()



plot_MPL_evenness_over_time <- plot_evenness_density_over_time(evenness_of_otus_over_time(long_MPL_otus),
                                                               n_distinct(long_MPL_otus$OTUid),
                                                               "MPL",
                                                               "usearch",
                                                               colours_for_plots[1])+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_MPL_evenness_over_time


plot_gp23_evenness_over_time <- plot_evenness_density_over_time(evenness_of_otus_over_time(long_gp23_otus),
                                                                Number_of_OTUs_gp23,
                                                                "gp23",
                                                                "usearch",
                                                                colours_for_plots[2])+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)

plot_gp23_evenness_over_time


plot_18S_evenness_over_time <- plot_evenness_density_over_time(evenness_of_otus_over_time(long_18s_otus),
                                                               n_distinct(long_18s_otus$OTUid),
                                                               "18S",
                                                               "usearch",
                                                               colours_for_plots[3])+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_18S_evenness_over_time

plot_16S_evenness_over_time <- plot_evenness_density_over_time(evenness_of_otus_over_time(long_16s_otus),
                                                               n_distinct(long_16s_otus$OTUid),
                                                               "16S",
                                                               "usearch",
                                                               colours_for_plots[4])
plot_16S_evenness_over_time

plot_AVS_evenness_over_time <- plot_evenness_density_over_time(evenness_of_otus_over_time(long_AVS_otus),
                                                               n_distinct(long_AVS_R1_otus$OTUid),
                                                               "AVS R1",
                                                               "usearch",
                                                               colours_for_plots[5])+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_AVS_evenness_over_time


pdf("../figures/time_series_evenness_all_amplicons_density%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(
 #plot_MPL_evenness_over_time,
 #           plot_AVS_evenness_over_time,
 plot_gp23_evenness_over_time,
 plot_18S_evenness_over_time, 
 plot_16S_evenness_over_time,
 ncol=1)
dev.off()



add_date_to_OTU_table_rtrn_long <- function (otu_table, Jericho_and_SOG_data, Jericho_data) {
 long_otus <- melt(as.matrix(otu_table), id=row.names(otu_table), variable.name="OTUid")
 names(long_otus)[1] <- "VC_number"
 names(long_otus)[2] <- "OTUid"
 ## Add in dates:
 long_otus$Date <- Jericho_data$Date[match(long_otus$VC_number, Jericho_data$VC_number)]
 return(long_otus)
}



## testing out stuff from Quince et al from Stamps course ####

row.names(Jericho_data) <- Jericho_data$VC_number
str(Jericho_data)
Jericho_data <- subset(Jericho_data, select=-VC_number)
str(Jericho_data)

##test with MPL
N <- rowSums(normalized_MPL_OTUs)
S <- specnumber(normalized_MPL_OTUs)
chla <- subset(Jericho_data, select=Average_chl_a)
chla <- subset(chla, row.names(chla) %in% row.names(normalized_MPL_OTUs))

DateJ <- subset(Jericho_data, select=Date)
DateJ <- subset(DateJ, row.names(DateJ) %in% row.names(normalized_MPL_OTUs))

Salinity <- subset(Jericho_data, select=Salinity_ppt_YSI)
Salinity <- subset(Salinity, row.names(Salinity) %in% row.names(normalized_MPL_OTUs))
plot(chla$Average_chl_a, S)
#cor.test(chla$Average_chl_a, S)
## not significant

#S.lm <- lm(S ~ chla$Average_chl_a + Salinity$Salinity_ppt_YSI + DateJ$Date, #data=Jericho_data)
#summary(S.lm)
# library(ape)
# library(mgcv)
# chla_ave <- chla$Average_chl_a
# test.gam <- gam(S~s(chla_ave))
# summary(test.gam)
# plot(test.gam)


## print out the richness estimates

AVS_richness <- richness_of_otus_over_time(long_AVS_otus)
MPL_richness <- richness_of_otus_over_time(long_MPL_otus)
S16_richness <- richness_of_otus_over_time(long_16s_otus)
S18_richness <- richness_of_otus_over_time(long_18s_otus)
gp23_richness <- richness_of_otus_over_time(long_gp23_otus)

all_richness <- merge(AVS_richness, MPL_richness, by="Date",all=TRUE, suffixes = c(".AVS",".MPL"))
all_richness <- merge(all_richness, S16_richness, by="Date", all=TRUE, suffixes = c("",".16S"))
all_richness <- merge(all_richness, S18_richness, by="Date", all=TRUE, suffixes = c("",".18S"))
all_richness <- merge(all_richness, gp23_richness, by="Date", all=TRUE, suffixes = c("",".gp23"))
names(all_richness)[4] <- "richness.16S"
write.csv(all_richness, "../results/amplicon_richness_by_date.csv", row.names=FALSE)

################################################
#### Do with the whole data set ######
##############################################

## read in normalized OTU tables

## while reprocessing it is missing the VC number
normalized_18s_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_18S.tsv",
                                  row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_gp23.tsv", 
                                   row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_MPL.tsv", 
                                  row.names="VC_number")
#normalized_AVS_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
#normalized_AVS_R1_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_AVS_R1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_AVS_R2.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")))
# HAKA_bloom <- geom_rect(aes(xmin = as.Date("2011-06-20"),
#                             xmax = as.Date("2011-07-06"), 
#                             ymin = -Inf,
#                             ymax = Inf),
#                         alpha = 0.2,
#                         fill = "burlywood3",
#                         colour=NA) 

 

### Get relevant VC numbers
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis_all_time_series.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)

## first rename the column names in the OTU table to something better
add_date_to_OTU_table_rtrn_long <- function (otu_table, 
                                             Jericho_data) {
 otu_table$VC_number = rownames(otu_table)
 long_otus <- melt(otu_table, 
                   id="VC_number", 
                   variable.name="OTUid")
 ## Add in dates:
 long_otus$Date <- Jericho_data$Date[match(long_otus$VC_number, 
                                           Jericho_data$VC_number)]
 return(long_otus)
}

## Read in OTU table and get long format with Date
long_gp23_otus <- add_date_to_OTU_table_rtrn_long(normalized_gp23_OTUs,
                                                  Jericho_data)

long_18s_otus <- add_date_to_OTU_table_rtrn_long(normalized_18s_OTUs,
                                                 Jericho_data)

# long_AVS_otus <- add_date_to_OTU_table_rtrn_long(normalized_AVS_OTUs,
#                                                  Jericho_data)

long_MPL_otus <- add_date_to_OTU_table_rtrn_long(normalized_MPL_OTUs,
                                                 Jericho_data)

long_16s_otus <- add_date_to_OTU_table_rtrn_long(normalized_16s_OTUs,
                                                 Jericho_data)


## how many OTUs total?
Number_of_OTUs_gp23 <- n_distinct(long_gp23_otus$OTUid)
Number_of_OTUs_18s <- n_distinct(long_18s_otus$OTUid)
Number_of_OTUs_MPL <- n_distinct(long_MPL_otus$OTUid)
# Number_of_OTUs_AVS <- n_distinct(long_AVS_otus$OTUid)
Number_of_OTUs_16s <- n_distinct(long_16s_otus$OTUid)

## pull out the counts of the number of OTUs-richness over time

richness_of_otus_over_time <- function (long_otus) {
 summarise_richness_otus <- long_otus %>%
  filter(value > 0) %>%
  ddply(~Date,summarise,richness=length(OTUid))
}


par(mar=c(0,0,0,0))
plot_richness_otus_over_time <- function (richness_over_time_data, total_count_of_OTUs, primer, clustering_method) {
 sampling_points <- richness_over_time_data$Date
 p <- ggplot(richness_over_time_data, aes(x=Date, y=richness))
 p + 
  season_line +
  HAKA_bloom +
  spring_bloom_line+
  geom_line(colour=line_colour )+
  date_scaling +
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16))+
  ggtitle(primer)
}

plot_gp23_richness_over_time_all_times <- plot_richness_otus_over_time(richness_of_otus_over_time(long_gp23_otus),
                                                                       Number_of_OTUs_gp23,
                                                                       "gp23",
                                                                       "usearch")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)



plot_18S_richness_over_time_all_times <- plot_richness_otus_over_time(richness_of_otus_over_time(long_18s_otus),
                                                                      n_distinct(long_18s_otus$OTUid),
                                                                      "18S",
                                                                      "usearch")



plot_16S_richness_over_time_all_times <- plot_richness_otus_over_time(richness_of_otus_over_time(long_16s_otus),
                                                                      n_distinct(long_16s_otus$OTUid),
                                                                      "16S",
                                                                      "usearch")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

# plot_AVS_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_AVS_otus),
#                                                             n_distinct(long_AVS_R1_otus$OTUid),
#                                                             "AVS R1",
#                                                             "usearch")+
#  theme(axis.ticks.x = element_blank(),
#        axis.text.x = element_blank(),
#        axis.title.x=element_blank())+
#  xlab(NULL)

# plot_AVS_richness_over_time

plot_MPL_richness_over_time_all_times <- plot_richness_otus_over_time(richness_of_otus_over_time(long_MPL_otus),
                                                                      n_distinct(long_MPL_otus$OTUid),
                                                                      "MPL",
                                                                      "usearch")+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_MPL_richness_over_time


pdf("../figures/time_series_with_all_time_richness_all_amplicons_line_graphs%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(plot_MPL_richness_over_time_all_times,
             # plot_AVS_richness_over_time,
             plot_gp23_richness_over_time_all_times,
             plot_18S_richness_over_time_all_times, 
             plot_16S_richness_over_time_all_times,
             ncol=1)
dev.off()

## Test out density plotting ####

## pick some colours to use:
library(RColorBrewer)
colours_for_plots <- brewer.pal(n = 5, name = 'Dark2')

plot_richness_density_over_time <- function (richness_over_time_data, total_count_of_OTUs, primer, clustering_method, fillcolour) {
 sampling_points <- richness_over_time_data$Date
 p <- ggplot(richness_over_time_data, aes(x=Date, y=richness))
 p +
HAKA_bloom +
    season_line +
  spring_bloom_line+
  geom_density(stat="identity", colour=line_colour, fill=fillcolour)+
  date_scaling +
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16))
 #+
 # ggtitle(primer)
}

plot_MPL_richness_over_time_all_times <- plot_richness_density_over_time(richness_of_otus_over_time(long_MPL_otus),
                                                                         n_distinct(long_MPL_otus$OTUid),
                                                                         "MPL",
                                                                         "usearch",
                                                                         colours_for_plots[1])+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_MPL_richness_over_time_all_times


plot_gp23_richness_over_time_all_times <- plot_richness_density_over_time(richness_of_otus_over_time(long_gp23_otus),
                                                                          Number_of_OTUs_gp23,
                                                                          "gp23",
                                                                          "usearch",
                                                                          colours_for_plots[2])+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)

plot_gp23_richness_over_time_all_times


plot_18S_richness_over_time_all_times <- plot_richness_density_over_time(richness_of_otus_over_time(long_18s_otus),
                                                                         n_distinct(long_18s_otus$OTUid),
                                                                         "18S",
                                                                         "usearch",
                                                                         colours_for_plots[3])

plot_18S_richness_over_time_all_times

plot_16S_richness_over_time_all_times <- plot_richness_density_over_time(richness_of_otus_over_time(long_16s_otus),
                                                                         n_distinct(long_16s_otus$OTUid),
                                                                         "16S",
                                                                         "usearch",
                                                                         colours_for_plots[4])+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_16S_richness_over_time_all_times

# plot_AVS_richness_over_time <- plot_richness_density_over_time(richness_of_otus_over_time(long_AVS_otus),
#                                                                n_distinct(long_AVS_R1_otus$OTUid),
#                                                                "AVS R1",
#                                                                "usearch",
#                                                                colours_for_plots[5])+
#  theme(axis.ticks.x = element_blank(),
#        axis.text.x = element_blank(),
#        axis.title.x=element_blank())+
#  xlab(NULL)
# 
# plot_AVS_richness_over_time


pdf("../figures/time_series_with_all_time_richness_all_amplicons_density%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(
 #plot_MPL_richness_over_time,
 #           plot_AVS_richness_over_time,
 plot_gp23_richness_over_time_all_times,
 plot_18S_richness_over_time_all_times, 
 plot_16S_richness_over_time_all_times,
 ncol=1)
dev.off()


pdf("../figures/time_series_both_time_and_bloom_richness_all_amplicons_density.pdf", width = 20, height = 15, onefile = FALSE)
plot_grid(
 #plot_MPL_richness_over_time,
 #           plot_AVS_richness_over_time,
 plot_gp23_richness_over_time_all_times+
  scale_x_date(limits = c(as.Date("2011-06-20"),
                          as.Date("2011-07-06"))),
 plot_gp23_richness_over_time_all_times+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()
  ),
 plot_16S_richness_over_time_all_times+
  scale_x_date(limits = c(as.Date("2011-06-20"),
                          as.Date("2011-07-06"))),
 plot_16S_richness_over_time_all_times+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()
  ),
 plot_18S_richness_over_time_all_times+
  scale_x_date(limits = c(as.Date("2011-06-20"),
                          as.Date("2011-07-06"))),
 plot_18S_richness_over_time_all_times+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()
  ), 

 ncol=2, labels=c("A", "B", "C", "D", "E", "F"), label_size = 24)
dev.off()


plot_evenness_density_over_time <- function (evenness_over_time_data, total_count_of_OTUs, primer, clustering_method, fillcolour) {
 sampling_points <- evenness_over_time_data$Date
 p <- ggplot(evenness_over_time_data, aes(x=Date, y=J))
 p +
  HAKA_bloom +
   season_line +
  spring_bloom_line+
  geom_density(stat="identity", colour=line_colour, fill=fillcolour)+
  date_scaling +
  ylim(c(-0.02, 1.0))+
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16))
 #+
 # ggtitle(primer)
}

plot_MPL_evenness_over_time_all_times <- plot_evenness_density_over_time(evenness_of_otus_over_time(long_MPL_otus),
                                                                         n_distinct(long_MPL_otus$OTUid),
                                                                         "MPL",
                                                                         "usearch",
                                                                         colours_for_plots[1])+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_MPL_evenness_over_time_all_times


plot_gp23_evenness_over_time_all_times <- plot_evenness_density_over_time(evenness_of_otus_over_time(long_gp23_otus),
                                                                          Number_of_OTUs_gp23,
                                                                          "gp23",
                                                                          "usearch",
                                                                          colours_for_plots[2])+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)

plot_gp23_evenness_over_time_all_times


plot_18S_evenness_over_time_all_times <- plot_evenness_density_over_time(evenness_of_otus_over_time(long_18s_otus),
                                                                         n_distinct(long_18s_otus$OTUid),
                                                                         "18S",
                                                                         "usearch",
                                                                         colours_for_plots[3])

plot_18S_evenness_over_time_all_times

plot_16S_evenness_over_time_all_times <- plot_evenness_density_over_time(evenness_of_otus_over_time(long_16s_otus),
                                                                         n_distinct(long_16s_otus$OTUid),
                                                                         "16S",
                                                                         "usearch",
                                                                         colours_for_plots[4])+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)
plot_16S_evenness_over_time_all_times

# plot_AVS_evenness_over_time <- plot_evenness_density_over_time(evenness_of_otus_over_time(long_AVS_otus),
#                                                                n_distinct(long_AVS_R1_otus$OTUid),
#                                                                "AVS R1",
#                                                                "usearch",
#                                                                colours_for_plots[5])+
#  theme(axis.ticks.x = element_blank(),
#        axis.text.x = element_blank(),
#        axis.title.x=element_blank())+
#  xlab(NULL)
# 
# plot_AVS_evenness_over_time


pdf("../figures/time_series_with_all_time_evenness_all_amplicons_density%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(
 #plot_MPL_evenness_over_time,
 #           plot_AVS_evenness_over_time,
 plot_gp23_evenness_over_time_all_times,
 plot_18S_evenness_over_time_all_times, 
 plot_16S_evenness_over_time_all_times,
 ncol=1)
dev.off()


pdf("../figures/time_series_both_time_and_bloom_evenness_all_amplicons_density.pdf", width = 20, height = 15, onefile = FALSE)
plot_grid(
 #plot_MPL_evenness_over_time,
 #           plot_AVS_evenness_over_time,
 plot_gp23_evenness_over_time_all_times+
  scale_x_date(limits = c(as.Date("2011-06-20"),
                          as.Date("2011-07-06"))),
 plot_gp23_evenness_over_time_all_times+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()
  ),
 plot_16S_evenness_over_time_all_times+
  scale_x_date(limits = c(as.Date("2011-06-20"),
                          as.Date("2011-07-06"))),
 plot_16S_evenness_over_time_all_times+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()
  ),
 plot_18S_evenness_over_time_all_times+
  scale_x_date(limits = c(as.Date("2011-06-20"),
                          as.Date("2011-07-06"))),
 plot_18S_evenness_over_time_all_times+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()
  ), 
 ncol=2,
 labels=c("A", "B", "C", "D", "E", "F"), label_size = 24)
dev.off()





add_date_to_OTU_table_rtrn_long <- function (otu_table, Jericho_and_SOG_data, Jericho_data) {
 long_otus <- melt(as.matrix(otu_table), id=row.names(otu_table), variable.name="OTUid")
 names(long_otus)[1] <- "VC_number"
 names(long_otus)[2] <- "OTUid"
 ## Add in dates:
 long_otus$Date <- Jericho_data$Date[match(long_otus$VC_number, Jericho_data$VC_number)]
 return(long_otus)
}

## testing out stuff from Quince et al from Stamps course ####

row.names(Jericho_data) <- Jericho_data$VC_number
str(Jericho_data)
Jericho_data <- subset(Jericho_data, select=-VC_number)
str(Jericho_data)

##test with MPL
N <- rowSums(normalized_MPL_OTUs)
S <- specnumber(normalized_MPL_OTUs)
chla <- subset(Jericho_data, select=Average_chl_a)
chla <- subset(chla, row.names(chla) %in% row.names(normalized_MPL_OTUs))

DateJ <- subset(Jericho_data, select=Date)
DateJ <- subset(DateJ, row.names(DateJ) %in% row.names(normalized_MPL_OTUs))

Salinity <- subset(Jericho_data, select=Salinity_ppt_YSI)
Salinity <- subset(Salinity, row.names(Salinity) %in% row.names(normalized_MPL_OTUs))
plot(chla$Average_chl_a, S)
#cor.test(chla$Average_chl_a, S)
## not significant

#S.lm <- lm(S ~ chla$Average_chl_a + Salinity$Salinity_ppt_YSI + DateJ$Date, #data=Jericho_data)
#summary(S.lm)
# library(ape)
# library(mgcv)
# chla_ave <- chla$Average_chl_a
# test.gam <- gam(S~s(chla_ave))
# summary(test.gam)
# plot(test.gam)


## print out the richness estimates

#AVS_richness <- richness_of_otus_over_time(long_AVS_otus)
MPL_richness <- richness_of_otus_over_time(long_MPL_otus)
S16_richness <- richness_of_otus_over_time(long_16s_otus)
S18_richness <- richness_of_otus_over_time(long_18s_otus)
gp23_richness <- richness_of_otus_over_time(long_gp23_otus)

all_richness <- merge(MPL_richness, S16_richness, by="Date", all=TRUE, suffixes = c(".MPL",".16S"))
all_richness <- merge(all_richness, S18_richness, by="Date", all=TRUE, suffixes = c("",".18S"))
all_richness <- merge(all_richness, gp23_richness, by="Date", all=TRUE, suffixes = c("",".gp23"))
names(all_richness)[4] <- "richness.16S"
write.csv(all_richness, "../results/amplicon_richness_with_all_time_by_date.csv", row.names=FALSE)




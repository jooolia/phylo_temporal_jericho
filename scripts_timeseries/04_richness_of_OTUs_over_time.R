## Examine richness of OTUs over time
## Also examine the correlations of the richness over time with the environmental parameters. 

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
library(psych)


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



## Need to edit and make sure the font, type are consistent
## To-do:

## read in normalized OTU tables

## while reprocessing it is missing the VC number
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv",                                  row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv", 
                                   row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", 
                                  row.names="VC_number")
#normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
normalized_AVS_R1_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R2.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 


normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv",                                  row.names=1)




args <- commandArgs(TRUE)
inputFile <- args[1]

if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
 path_colour <- "black"
 line_colour <- "black"
 point_colour <- "black"
 figures_dir <- "../figures/"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
  point_colour <- "white"
  figures_dir <- "../figures_pres/"
 }
}

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

long_18s_otus_phytos <- add_date_to_OTU_table_rtrn_long(normalized_18s_OTUs_phytos,
                                                 Jericho_data)
long_18s_otus_hetero <- add_date_to_OTU_table_rtrn_long(normalized_18s_OTUs_hetero,
                                                 Jericho_data)

long_AVS_otus <- add_date_to_OTU_table_rtrn_long(normalized_AVS_R1_OTUs,
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
  ddply(~Date,
        summarise,
        richness=length(OTUid))
}

par(mar=c(0,0,0,0))
plot_richness_otus_over_time <- function (richness_over_time_data,
                                          total_count_of_OTUs,
                                          primer,
                                          clustering_method,
                                          title) {
 p <- ggplot(richness_over_time_data, aes(x=Date, y=richness))+
  theme_JAG_presentation()+
 season_line +
 spring_bloom_line+
  geom_line(colour=line_colour)+
  geom_point(colour=point_colour)+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16),
        panel.grid.major = element_blank())+
  ggtitle(NULL)+
  date_scaling+
  ggtitle(title)
}

plot_gp23_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_gp23_otus),
                                                             Number_of_OTUs_gp23,
                                                             "gp23",
                                                             "usearch", 
                                                             "gp23")

plot_gp23_richness_over_time


plot_18S_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_18s_otus),
                                                             n_distinct(long_18s_otus$OTUid),
                                                             "18S",
                                                             "usearch", "18S")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_18S_richness_over_time

#phyto
plot_18S_phytos_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_18s_otus_phytos),
                                                            n_distinct(long_18s_otus$OTUid),
                                                            "18S",
                                                            "usearch", "18S phyto")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_18S_phytos_richness_over_time

#het

plot_18S_hetero_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_18s_otus_hetero),
                                                            n_distinct(long_18s_otus$OTUid),
                                                            "18S",
                                                            "usearch", "18S hetero")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_18S_hetero_richness_over_time


plot_16S_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_16s_otus),
                                                            n_distinct(long_16s_otus$OTUid),
                                                            "16S",
                                                            "usearch","16S")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)

plot_16S_richness_over_time

plot_AVS_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_AVS_otus),
                                                            n_distinct(long_AVS_otus$OTUid),
                                                            "AVS R1",
                                                            "usearch", "AVS_R1")+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_AVS_richness_over_time

plot_MPL_richness_over_time <- plot_richness_otus_over_time(richness_of_otus_over_time(long_MPL_otus),
                                                               n_distinct(long_MPL_otus$OTUid),
                                                               "MPL",
                                                               "usearch","MPL")+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

plot_MPL_richness_over_time

pdf(paste0(figures_dir,"time_series_richness_all_amplicons_line_graphs%03d.pdf"), width = 20, height = 15, onefile = FALSE)
grid.arrange(plot_MPL_richness_over_time,
             plot_AVS_richness_over_time,
             plot_gp23_richness_over_time,
             plot_18S_richness_over_time, 
             plot_16S_richness_over_time,
             ncol=1)
dev.off()


pdf(paste0(figures_dir,"time_series_richness_all_amplicons_NO_AVS_line_graphs%03d.pdf"), width = 20, height = 15, onefile = FALSE)
grid.arrange(plot_MPL_richness_over_time,
             plot_18S_richness_over_time, 
             #plot_AVS_richness_over_time,
             plot_16S_richness_over_time,
             plot_gp23_richness_over_time,
             ncol=1)
dev.off()

pdf(paste0(figures_dir,"time_series_richness_all_amplicons_NO_AVS_with_phytos_and_het_line_graphs%03d.pdf"), width = 20, height = 15, onefile = FALSE)
grid.arrange(plot_MPL_richness_over_time,
             plot_18S_richness_over_time, 
             plot_18S_phytos_richness_over_time, 
             plot_18S_hetero_richness_over_time, 
             #plot_AVS_richness_over_time,
             plot_16S_richness_over_time,
             plot_gp23_richness_over_time,
             ncol=1)
dev.off()


### look at evenness

evenness_of_otus_over_time <- function (normalized_otus) {
 #normalized_otus <- normalized_gp23_OTUs
 H <- diversity(normalized_otus)
 ## Species richness (S) and Pielou's evenness (J):
 S <- specnumber(normalized_otus) 
 J <- as.data.frame(H/log(S))
 
 J$Date <- Jericho_data$Date[match(rownames(normalized_otus), 
                                           Jericho_data$VC_number)]
 names(J)[1] <- "evenness" 
 return(J)
}


plot_evenness_otus_over_time <- function (evenness_over_time_data,
                                          primer,
                                          title) {
 p <- ggplot(evenness_over_time_data, aes(x=Date, y=evenness))+
  theme_JAG_presentation()+
  season_line +
  spring_bloom_line+
  geom_line(colour=line_colour)+
  geom_point(colour=point_colour)+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16),
        panel.grid.major = element_blank())+
  ggtitle(NULL)+
  date_scaling+
  ggtitle(title)
 # print(p)
}

plot_gp23_evenness_over_time <- plot_evenness_otus_over_time(evenness_of_otus_over_time(normalized_gp23_OTUs),
                                                             "gp23",
                                                             "gp23")

plot_gp23_evenness_over_time


plot_18S_evenness_over_time <- plot_evenness_otus_over_time(evenness_of_otus_over_time(normalized_18s_OTUs),
                                                            "18S","18S")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())

plot_18S_evenness_over_time

#phyto
plot_18S_phytos_evenness_over_time <- plot_evenness_otus_over_time(evenness_of_otus_over_time(normalized_18s_OTUs_phytos),                                                                   "18S", "18S phyto")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())

plot_18S_phytos_evenness_over_time

#het

plot_18S_hetero_evenness_over_time <- plot_evenness_otus_over_time(evenness_of_otus_over_time(normalized_18s_OTUs_hetero),
                                                                   "18S","18S hetero")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(),
       axis.title.x=element_blank())

plot_18S_hetero_evenness_over_time


plot_16S_evenness_over_time <- plot_evenness_otus_over_time(evenness_of_otus_over_time(normalized_16s_OTUs),
                                                            "16S","16S")+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())
plot_16S_evenness_over_time

plot_AVS_evenness_over_time <- plot_evenness_otus_over_time(evenness_of_otus_over_time(normalized_AVS_R1_OTUs),
                                                            "AVS R1", "AVS_R1")+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())

plot_AVS_evenness_over_time

plot_MPL_evenness_over_time <- plot_evenness_otus_over_time(evenness_of_otus_over_time(normalized_MPL_OTUs),
                                                            "MPL","MPL")+ 
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())

plot_MPL_evenness_over_time

pdf(paste0(figures_dir,"time_series_evenness_all_amplicons_line_graphs%03d.pdf"), width = 20, height = 15, onefile = FALSE)
grid.arrange(plot_MPL_evenness_over_time,
             #plot_AVS_evenness_over_time,
             plot_gp23_evenness_over_time,
             plot_18S_evenness_over_time, 
             plot_16S_evenness_over_time,
             ncol=1)
dev.off()


pdf(paste0(figures_dir,"time_series_evenness_all_amplicons_NO_AVS_line_graphs%03d.pdf"), width = 20, height = 15, onefile = FALSE)
grid.arrange(plot_MPL_evenness_over_time,
             plot_18S_evenness_over_time, 
             #plot_AVS_evenness_over_time,
             plot_16S_evenness_over_time,
             plot_gp23_evenness_over_time,
             ncol=1)
dev.off()

pdf(paste0(figures_dir,"time_series_evenness_all_amplicons_NO_AVS_with_phytos_and_het_line_graphs%03d.pdf"), width = 20, height = 15, onefile = FALSE)
grid.arrange(plot_MPL_evenness_over_time,
             plot_18S_evenness_over_time, 
             plot_18S_phytos_evenness_over_time, 
             plot_18S_hetero_evenness_over_time, 
             #plot_AVS_evenness_over_time,
             plot_16S_evenness_over_time,
             plot_gp23_evenness_over_time,
             ncol=1)
dev.off()



## Test out density plotting ####

## pick some colours to use:
library(RColorBrewer)
colours_for_plots <- brewer.pal(n = 5, name = 'Dark2')

plot_richness_density_over_time <- function (richness_over_time_data, 
                                             total_count_of_OTUs, 
                                             primer, 
                                             clustering_method, 
                                             fillcolour) {
 p <- ggplot(richness_over_time_data, aes(x=Date, y=richness))+
  season_line +
  spring_bloom_line
 
 p2 <- p + geom_density(stat="identity", colour="black", fill=fillcolour)+
  geom_point()+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16),
        panel.grid.major = element_blank())+
  ggtitle(NULL)+ 
  date_scaling
 return(p2)
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


pdf(paste0(figures_dir,"time_series_richness_all_amplicons_density%03d.pdf"), width = 20, height = 15, onefile = FALSE)
grid.arrange(plot_MPL_richness_over_time,
             plot_AVS_richness_over_time,
             plot_gp23_richness_over_time,
             plot_18S_richness_over_time, 
             plot_16S_richness_over_time,
             ncol=1)
dev.off()

add_date_to_OTU_table_rtrn_long <- function (otu_table, 
                                             Jericho_and_SOG_data,
                                             Jericho_data) {
 long_otus <- melt(as.matrix(otu_table), 
                   id=row.names(otu_table),
                   variable.name="OTUid")
 names(long_otus)[1] <- "VC_number"
 names(long_otus)[2] <- "OTUid"
 ## Add in dates:
 long_otus$Date <- Jericho_data$Date[match(long_otus$VC_number, Jericho_data$VC_number)]
 return(long_otus)
}

add_environmental_data_to_original_data <- function (normalized_OTU_data, Jericho_data ){
 ################ Add  enviro data ----------------- 
 # Add the mean of the environmental data
 normalized_OTU_data$PO4 <- Jericho_data$Average_PO4[match(normalized_OTU_data$VC_number, Jericho_data$VC_number)]
 
 normalized_OTU_data$SiO2 <- Jericho_data$Average_SiO2[match(normalized_OTU_data$VC_number, Jericho_data$VC_number)]
 
 normalized_OTU_data$NO3_NO2 <- Jericho_data$Average_NO3_NO2[match(normalized_OTU_data$VC_number, Jericho_data$VC_number)]
 
 normalized_OTU_data$Chl_a <- Jericho_data$Average_chl_a[match(normalized_OTU_data$VC_number, Jericho_data$VC_number)]
 
 normalized_OTU_data$Dissolved_oxygen_percent <- Jericho_data$Dissolved_oxygen_percent[match(normalized_OTU_data$VC_number, Jericho_data$VC_number)]
 
 normalized_OTU_data$Salinity_ppt_YSI <- Jericho_data$Salinity_ppt_YSI[match(normalized_OTU_data$VC_number, Jericho_data$VC_number)]
 
 normalized_OTU_data$Temperature_YSI <- Jericho_data$Temperature_YSI[match(normalized_OTU_data$VC_number, Jericho_data$VC_number)]
 normalized_OTU_data$season <- Jericho_data$season[match(normalized_OTU_data$VC_number, Jericho_data$VC_number)]
 
 return(normalized_OTU_data)  
}


richness_by_environ_parem_print <- function (normalized_otu_table, 
                                             Library_metadata,
                                             Jericho_data,
                                             name) {
 OTU_long <- add_date_to_OTU_table_rtrn_long(normalized_otu_table,
                                             Library_metadata,
                                             Jericho_data)

 OTU_with_environ <- add_environmental_data_to_original_data(OTU_long, Jericho_data)
 OTU_with_environ <- subset(OTU_with_environ, value > 0)
## could summarise by date, number of OTUs found...
OTU_with_environ_count <- OTU_with_environ %>%
 group_by(Date, PO4, SiO2, NO3_NO2, Chl_a, Dissolved_oxygen_percent, Salinity_ppt_YSI, Temperature_YSI, season) %>%
 select(Date, PO4, SiO2, NO3_NO2, Chl_a, Dissolved_oxygen_percent, Salinity_ppt_YSI, Temperature_YSI, season, value) %>%
 summarize(count = n()) 
 
 OTU_temp <- ggplot( OTU_with_environ_count, aes(x=Temperature_YSI, y=count)) +geom_point(position="jitter")
 OTU_season <- ggplot( OTU_with_environ_count, aes(x=season,y=count)) +geom_point(position="jitter")
 OTU_Chla <- ggplot( OTU_with_environ_count, aes(x=Chl_a,y=count)) +geom_point(position="jitter")
 OTU_PO4 <- ggplot( OTU_with_environ_count, aes(x=PO4,y=count)) +geom_point(position="jitter")
 OTU_Sil <- ggplot( OTU_with_environ_count, aes(x=SiO2,y=count)) +geom_point(position="jitter")
 OTU_NO <- ggplot( OTU_with_environ_count, aes(x=NO3_NO2,y=count)) +geom_point(position="jitter")
 OTU_OX <- ggplot( OTU_with_environ_count, aes(x=Dissolved_oxygen_percent,y=count)) +geom_point(position="jitter")
 OTU_sal <- ggplot( OTU_with_environ_count, aes(x=Salinity_ppt_YSI,y=count)) +geom_point(position="jitter")
 pdf( paste(figures_dir,name, "_richness_histograms%03d.pdf", sep=""),width = 11, height = 17,onefile = FALSE)
 
 sidebysideplot <- grid.arrange(OTU_Chla,OTU_Sil,OTU_PO4,OTU_NO,OTU_sal,OTU_temp,OTU_season,OTU_OX,ncol=2)
 sidebysideplot
 dev.off()
}

richness_by_environ_parem_print(normalized_16s_OTUs, Library_metadata, Jericho_data, "16S")
richness_by_environ_parem_print(normalized_18s_OTUs, Library_metadata, Jericho_data, "18S")
richness_by_environ_parem_print(normalized_18s_OTUs_phytos, Library_metadata, Jericho_data, "18S")
richness_by_environ_parem_print(normalized_18s_OTUs_hetero, Library_metadata, Jericho_data, "18S")
richness_by_environ_parem_print(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23")
richness_by_environ_parem_print(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL")
richness_by_environ_parem_print(normalized_AVS_R1_OTUs, Library_metadata, Jericho_data, "AVS_R1")

## testing out stuff from Quince et al from Stamps course ####
## actually this stuff could be worth a second look. How to interpret the GAM....look over those notes again?
row.names(Jericho_data) <- Jericho_data$VC_number
str(Jericho_data)
Jericho_data <- subset(Jericho_data, select=-VC_number)
str(Jericho_data)

##test with MPL
# N <- rowSums(normalized_MPL_OTUs)
# S <- specnumber(normalized_MPL_OTUs)
# chla <- subset(Jericho_data, select=Average_chl_a)
# chla <- subset(chla, row.names(chla) %in% row.names(normalized_MPL_OTUs))
# 
# DateJ <- subset(Jericho_data, select=Date)
# DateJ <- subset(DateJ, row.names(DateJ) %in% row.names(normalized_MPL_OTUs))
# 
# Salinity <- subset(Jericho_data, select=Salinity_ppt_YSI)
# Salinity <- subset(Salinity, row.names(Salinity) %in% row.names(normalized_MPL_OTUs))
# plot(chla$Average_chl_a, S)
# cor.test(chla$Average_chl_a, S)
# ## not significant
# 
# S.lm <- lm(S ~ chla$Average_chl_a + Salinity$Salinity_ppt_YSI + DateJ$Date, data=Jericho_data)
# summary(S.lm)
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

parameters_to_exclude <- c("season",
                           # "Average_NO3_NO2",
                           "Standard_error_NO3_NO2",
                           "Secchi_disk_disappears",
                           "Secchi_disk_reappears",
                           "VC_number",
                           "Raw_Bacteria_Rep_A",
                           "Raw_Bacteria_Rep_B",
                           "Raw_Viruses_rep_A",
                           "Raw_Viruses_rep_B",
                           "PO4_rep_A",
                           "PO4_rep_B",
                           "SiO2_rep_A",
                           "SiO2_rep_B",
                           "NO3_NO2_rep_A",
                           "NO3_NO2_rep_B",
                           "Standard_error_viral_abundance",
                           "Standard_error_bacterial_abundance",
                           "Standard_error_PO4",
                           "Standard_error_SiO2",
                           "Standard_error_chl_a")
params_Jericho <-  subset(Jericho_data, select= !(colnames(Jericho_data) %in% parameters_to_exclude))


### do correlations of the richness across the different environmental parameterss
# amplicon_richness_env <- corr.test(all_richness[,-1], params_Jericho[,-1], method="spearman",use = "pairwise")
# 
# correlations <- corr.test(all_richness[,-1], params_Jericho[,-1], method="spearman",use = "pairwise")$r
# p <- corr.test(all_richness[,-1], params_Jericho[,-1], method="spearman",use = "pairwise")$p
# 
# melted_cor <- melt(correlations)
# melted_p <- melt(p)
# mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
# melted_mystars <- melt(mystars)
# 
# melted_together <- cbind( melted_cor,p_value=melted_p$value, stars=melted_mystars$value)
# melted_together <- na.omit(melted_together) ## gets rid of the leftover diagonals
# names(melted_together)
# 
# melted_together$cor_with_stars <- paste(round(melted_together$value, digits=3), " ", melted_together$stars)
# 
# significant_cor <- dcast(melted_together, Var1 ~Var2, value.var="cor_with_stars")
# write.csv(significant_cor, file="../results/cor_amplicon_richness_with_env.csv")

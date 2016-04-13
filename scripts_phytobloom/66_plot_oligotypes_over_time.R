## graphing oligotypes

## need to select appropriate oligotypes..
## so first check appropriate ones and then make plots



## would really be in ../results/
library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

# percent_oligo_MPL_5 <- read.delim("../results/oligotypingMPL_OTU_5like_reads/MATRIX-PERCENT.txt")
# percent_oligo_MPL_3 <- read.delim("../results/oligotypingMPL_OTU_3like_reads/MATRIX-PERCENT.txt")
# percent_oligo_MPL_6 <- read.delim("../results/oligotypingMPL_OTU_6like_reads/MATRIX-PERCENT.txt")
 percent_oligo_MPL_1 <- read.delim("../results/oligotypingMPL_OTU_1like_reads/MATRIX-PERCENT.txt")
# percent_oligo_MPL_2 <- read.delim("../results/oligotypingMPL_OTU_2like_reads/MATRIX-PERCENT.txt")
# percent_oligo_MPL_4 <- read.delim("../results/oligotypingMPL_OTU_4like_reads/MATRIX-PERCENT.txt")
# percent_oligo_MPL_7 <- read.delim("../results/oligotypingMPL_OTU_7like_reads/MATRIX-PERCENT.txt")
# percent_oligo_MPL_9 <- read.delim("../results/oligotypingMPL_OTU_9like_reads/MATRIX-PERCENT.txt")
# percent_oligo_MPL_11 <- read.delim("../results/oligotypingMPL_OTU_11like_reads/MATRIX-PERCENT.txt")
# percent_oligo_MPL_332 <- read.delim("../results/oligotypingMPL_OTU_332like_reads/MATRIX-PERCENT.txt")

### similar to 1
percent_oligo_MPL_283 <- read.delim("../results/oligotypingMPL_OTU_283like_reads/MATRIX-PERCENT.txt")
#percent_oligo_MPL_544 <- read.delim("../results/oligotypingMPL_OTU_544like_reads/MATRIX-PERCENT.txt")
percent_oligo_MPL_1588 <- read.delim("../results/oligotypingMPL_OTU_1588like_reads/MATRIX-PERCENT.txt")





### Get relevant VC numbers
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv")


args <- commandArgs(TRUE)
inputFile <- args[1]
## test to see if input file is given, so I can decide whether to use this argument or the orginal one. 

if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
 path_colour <- "black"
 line_colour <- "black"
 point_colour <- "black"
 heatmap_fill <- "white"
 low_heatmap <- "yellow"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
  point_colour <- "white"
  heatmap_fill <- "black"
  low_heatmap <- "grey"
 }
}

# adding in the seasons and the spring bloom for the ggplots
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
date_scaling <-      scale_x_date(limits = c(as.Date("2011-06-20"), as.Date("2011-07-06")))



## get pool and add date

## first rename the column names in the OTU table to something better
add_date_to_oligotype_table <- function (oligotypes_percent,
                                         Library_VC) {
  oligotypes_percent$samples <- gsub("Lib-pool", 
                              "", 
                              oligotypes_percent$samples) 

  oligotypes_percent$VC_number <-Library_VC[match(oligotypes_percent$samples,
                                                as.character(Library_metadata$PCR_pool_number))]
  oligotypes_percent$VC_number <- gsub("_large", "", oligotypes_percent$VC_number) 
  oligotypes_percent$Date <- Jericho_data$Date[match(oligotypes_percent$VC_number, 
                                   Jericho_data$VC_number)]
  print(oligotypes_percent)
 long_otus <- melt(oligotypes_percent, 
                   id=c("VC_number", "samples", "Date"))
 return(long_otus)
}

plot_oligotypes_over_time_percent <- function (percent_oligo_MPL,
                                               Library_VC,
                                               name) {
  MPL_oligo <- add_date_to_oligotype_table(percent_oligo_MPL,
                                           Library_VC)
  plot_over_time <- ggplot(MPL_oligo, aes(x=Date, y=value, group=variable))+
   geom_bar(aes(fill=variable), stat="identity")+
   season_line+
   spring_bloom_line+
   date_scaling+
   theme_JAG_presentation()+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.title=element_blank())
  
  print(plot_over_time)
  pdf(paste0("../figures/", name, ".pdf"), width = 15, height = 11)
  print(plot_over_time)
  dev.off()
  return(plot_over_time)
}


plot_oligotypes_over_time_percent_lines <- function (percent_oligo_MPL,
                                               Library_VC,
                                               name) {

 MPL_oligo <- add_date_to_oligotype_table(percent_oligo_MPL,
                                          Library_VC)
 number_oligotypes <- length(unique(MPL_oligo$variable))
 sampling_points <- MPL_oligo$Date
 print(number_oligotypes)
 plot_over_time <- ggplot(MPL_oligo, aes(x=Date, y=value, group=variable))+
  season_line+
  spring_bloom_line+
    geom_line(aes(colour=variable), size=3)+
  date_scaling+
  scale_colour_discrete(name="Oligotype",
                      labels=c(seq_along(1:number_oligotypes)))+
  theme_JAG_presentation()+
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
 print(plot_over_time)
 pdf(paste0("../figures/", name, ".pdf"), width = 15, height = 11)
 print(plot_over_time)
 dev.off()
 return(plot_over_time)
}

plot_oligotypes_over_time_percent_lines_log_y <- function (percent_oligo_MPL,
                                                     Library_VC,
                                                     name) {
 MPL_oligo <- add_date_to_oligotype_table(percent_oligo_MPL,
                                          Library_VC)
 plot_over_time <- ggplot(MPL_oligo, aes(x=Date, y=value, group=variable))+ scale_y_log10()+
  geom_line(aes(colour=variable), size=3)+
  season_line+
  spring_bloom_line+
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank())
 print(plot_over_time)
 pdf(paste0("../figures/", name, ".pdf"), width = 15, height = 11)
 print(plot_over_time)
 dev.off()
 return(plot_over_time)
}


# plot_oligotypes_over_time_percent(percent_oligo_MPL_2,
#                                   Library_metadata$MPL_VC_number,
#                                   "MPL_OTU_2_oligo_percent_over_time") 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_3,
#                                   Library_metadata$MPL_VC_number,
#                                   "MPL_OTU_3_oligo_percent_over_time") 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_4,
#                                   Library_metadata$MPL_VC_number,
#                                   "MPL_OTU_4_oligo_percent_over_time") 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_5,
#                                   Library_metadata$MPL_VC_number,
#                                   "MPL_OTU_5_oligo_percent_over_time") 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_6,
#                                   Library_metadata$MPL_VC_number,
#                                   "MPL_OTU_6_oligo_percent_over_time") 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_7,
#                                   Library_metadata$MPL_VC_number,
#                                   "MPL_OTU_7_oligo_percent_over_time") 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_9,
#                                   Library_metadata$MPL_VC_number,
#                                   "MPL_OTU_9_oligo_percent_over_time") 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_11,
#                                   Library_metadata$MPL_VC_number,
#                                   "MPL_OTU_11_oligo_percent_over_time") 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_332,
#                                   Library_metadata$MPL_VC_number,
#                                   "MPL_OTU_332_oligo_percent_over_time") 


## OTU 1 and friends

plot_oligotypes_over_time_percent(percent_oligo_MPL_1,
                                  Library_metadata$MPL_VC_number,
                                  "MPL_OTU_1_oligo_percent_over_time") 
plot_oligotypes_over_time_percent(percent_oligo_MPL_283,
                                  Library_metadata$MPL_VC_number,
                                  "MPL_OTU_283_oligo_percent_over_time") 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_544,
#                                   Library_metadata$MPL_VC_number,
#                                   "MPL_OTU_544_oligo_percent_over_time") 
plot_oligotypes_over_time_percent(percent_oligo_MPL_1588,
                                  Library_metadata$MPL_VC_number,
                                  "MPL_OTU_1588_oligo_percent_over_time") 



### pulled out group A from the nucleotide stuff:

# percent_oligo_MPL_groupA <- read.delim("../results/oligotyping_MPL_group_A_from_nuc_like_reads/MATRIX-PERCENT.txt")
# 
# count_oligo_MPL_groupA <- read.delim("../results/oligotyping_MPL_group_A_from_nuc_like_reads/MATRIX-COUNT.txt")
# 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_groupA,
#                                   Library_metadata$MPL_VC_number,
#                                   "percent_oligo_MPL_groupA") 
# 
# plot_oligotypes_over_time_percent_lines(percent_oligo_MPL_groupA,
#                                   Library_metadata$MPL_VC_number,
#                                   "percent_oligo_MPL_groupA_lines") 
# 
# plot_oligotypes_over_time_percent_lines(count_oligo_MPL_groupA,
#                                         Library_metadata$MPL_VC_number,
#                                         "count_oligo_MPL_groupA_lines")  
# plot_oligotypes_over_time_percent_lines_log_y(count_oligo_MPL_groupA,
#                                         Library_metadata$MPL_VC_number,
#                                         "count_oligo_MPL_groupA_lines_log_y")

#### raphidophyte
percent_oligo_raphidopyte <- read.delim("../results/oligotyping_raphidophyte_OTUs/MATRIX-PERCENT.txt")

count_oligo_raphidopyte <- read.delim("../results/oligotyping_raphidophyte_OTUs/MATRIX-COUNT.txt")

plot_oligotypes_over_time_percent(percent_oligo_raphidopyte ,
                                  Library_metadata$X18s_VC_number,
                                  "percent_oligo_raphidophyte") 

plot_oligotypes_over_time_percent_lines(percent_oligo_raphidopyte ,
                                        Library_metadata$X18s_VC_number,
                                        "percent_oligo_raphidophyte_lines") 

count_raph_bloom <- plot_oligotypes_over_time_percent_lines(count_oligo_raphidopyte ,
                                        Library_metadata$X18s_VC_number,
                                        "count_oligo_raphidophyte") 

## wanted to remove the jne ones so that I could see more of the other variation. 
plot_oligotypes_over_time_percent_lines_log_y(count_oligo_raphidopyte ,
                                        Library_metadata$X18s_VC_number,
                                        "count_oligo_raphidophyte_log_y") 


## just haka

percent_oligo_haka <- read.delim("../results/haka_only_bloom/MATRIX-PERCENT.txt")

count_oligo_haka <- read.delim("../results/haka_only_bloom/MATRIX-COUNT.txt")

plot_oligotypes_over_time_percent(percent_oligo_haka ,
                                  Library_metadata$X18s_VC_number,
                                  "percent_oligo_haka") 

plot_oligotypes_over_time_percent_lines(percent_oligo_haka  ,
                                        Library_metadata$X18s_VC_number,
                                        "percent_oligo_haka_lines") 

plot_oligotypes_over_time_percent_lines(count_oligo_haka  ,
                                        Library_metadata$X18s_VC_number,
                                        "count_oligo_haka") 

## wanted to remove the jne ones so that I could see more of the other variation. 
plot_oligotypes_over_time_percent_lines_log_y(count_oligo_haka  ,
                                              Library_metadata$X18s_VC_number,
                                              "count_oligo_haka_log_y") 


#### cyano
percent_oligo_cyano <- read.delim("../results/oligotyping_cyano_otus/MATRIX-PERCENT.txt")

count_oligo_cyano <- read.delim("../results/oligotyping_cyano_otus/MATRIX-COUNT.txt")

plot_oligotypes_over_time_percent(percent_oligo_cyano,
                                  Library_metadata$X16s_VC_number,
                                  "percent_oligo_cyano") 

plot_oligotypes_over_time_percent_lines(percent_oligo_cyano ,
                                        Library_metadata$X16s_VC_number,
                                        "percent_oligo_cyano_lines") 

plot_oligotypes_over_time_percent_lines(count_oligo_cyano ,
                                        Library_metadata$X16s_VC_number,
                                        "count_oligo_cyano") 

plot_oligotypes_over_time_percent_lines_log_y(count_oligo_cyano ,
                                              Library_metadata$X16s_VC_number,
                                              "count_oligo_cyano_log_y") 

### Dino over all time series 

percent_oligo_dino <- read.delim("../results/dino_only_bloom/MATRIX-PERCENT.txt")

count_oligo_dino <- read.delim("../results/dino_only_bloom/MATRIX-COUNT.txt")

plot_oligotypes_over_time_percent(percent_oligo_dino ,
                                  Library_metadata$X18s_VC_number,
                                  "percent_oligo_dino") 

plot_oligotypes_over_time_percent_lines(percent_oligo_dino ,
                                        Library_metadata$X18s_VC_number,
                                        "percent_oligo_dino_lines") 

count_dino1_bloom <- plot_oligotypes_over_time_percent_lines(count_oligo_dino ,
                                        Library_metadata$X18s_VC_number,
                                        "count_oligo_dino") 

## wanted to remove the jne ones so that I could see more of the other variation. 
plot_oligotypes_over_time_percent_lines_log_y(count_oligo_dino ,
                                              Library_metadata$X18s_VC_number,
                                              "count_oligo_dino_log_y") 

## Dino 2 

percent_oligo_dino2 <- read.delim("../results/dino2_only_bloom/MATRIX-PERCENT.txt")

count_oligo_dino2 <- read.delim("../results/dino2_only_bloom/MATRIX-COUNT.txt")

plot_oligotypes_over_time_percent(percent_oligo_dino2 ,
                                  Library_metadata$X18s_VC_number,
                                  "percent_oligo_dino2") 

plot_oligotypes_over_time_percent_lines(percent_oligo_dino2 ,
                                        Library_metadata$X18s_VC_number,
                                        "percent_oligo_dino2_lines") 

count_dino2_bloom <- plot_oligotypes_over_time_percent_lines(count_oligo_dino2 ,
                                        Library_metadata$X18s_VC_number,
                                        "count_oligo_dino2") 

## wanted to remove the jne ones so that I could see more of the other variation. 
plot_oligotypes_over_time_percent_lines_log_y(count_oligo_dino2 ,
                                              Library_metadata$X18s_VC_number,
                                              "count_oligo_dino2_log_y") 


######################################################
##### Over whole time series



### Get relevant VC numbers
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis_all_time_series.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv")

date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")))
# HAKA_bloom <- geom_rect(aes(xmin = as.Date("2011-06-20"),
#                             xmax = as.Date("2011-07-06"), 
#                             ymin = -Inf,
#                             ymax = Inf),
#                         alpha = 0.8,
#                         fill = "burlywood3",
#                         colour=NA) 

HAKA_bloom <- annotate("rect",xmin = as.Date("2011-06-20"),
                            xmax = as.Date("2011-07-06"), 
                            ymin = -Inf,
                            ymax = Inf,
                        alpha = 0.4,
                        fill = "burlywood3",
                        colour=NA) 


## first rename the column names in the OTU table to something better
add_date_to_oligotype_table <- function (oligotypes_percent,
                                         Library_VC) {
 oligotypes_percent$samples <- gsub("Lib-pool", 
                                    "", 
                                    oligotypes_percent$samples) 
 
 oligotypes_percent$VC_number <-Library_VC[match(oligotypes_percent$samples,
                                                 as.character(Library_metadata$PCR_pool_number))]
 oligotypes_percent$VC_number <- gsub("_large", "", oligotypes_percent$VC_number) 
 
 oligotypes_percent$VC_number <- gsub("_redo", "", oligotypes_percent$VC_number) 
 
 oligotypes_percent$Date <- Jericho_data$Date[match(oligotypes_percent$VC_number, 
                                                    Jericho_data$VC_number)]
 print(oligotypes_percent)
 long_otus <- melt(oligotypes_percent, 
                   id=c("VC_number", "samples", "Date"))
 return(long_otus)
}

plot_oligotypes_over_time_percent <- function (percent_oligo_MPL,
                                               Library_VC,
                                               name) {
 MPL_oligo <- add_date_to_oligotype_table(percent_oligo_MPL,
                                          Library_VC)
 plot_over_time <- ggplot(MPL_oligo, aes(x=Date, y=value, group=variable))+
  season_line+
  HAKA_bloom+
  spring_bloom_line+
  geom_bar(aes(fill=variable), stat="identity")+
  ylab("Reads")+
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank())
 
 print(plot_over_time)
 pdf(paste0("../figures/", name, "all_times.pdf"), width = 15, height = 11)
 print(plot_over_time)
 dev.off()
 return(plot_over_time)
}

plot_oligotypes_over_time_percent_lines <- function (percent_oligo_MPL,
                                                     Library_VC,
                                                     name) {
 MPL_oligo <- add_date_to_oligotype_table(percent_oligo_MPL,
                                          Library_VC)
 number_oligotypes <- length(unique(MPL_oligo$variable))
 sampling_points <- MPL_oligo$Date
 plot_over_time <- ggplot(MPL_oligo, aes(x=Date, y=value, group=variable))+
  season_line+
  HAKA_bloom+
  spring_bloom_line+
  geom_line(aes(colour=variable), size=2)+
  ylab("Reads")+
  scale_colour_discrete(name="Oligotype",
                        labels=c(LETTERS[1:number_oligotypes]))+
  date_scaling+
  theme_JAG_presentation()+
  annotate("point", 
           x = as.Date(sampling_points),
           y = -30,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank())
 print(plot_over_time)
 pdf(paste0("../figures/", name, "all_times.pdf"), width = 15, height = 11)
 print(plot_over_time)
 dev.off()
 return(plot_over_time)
}

plot_oligotypes_over_time_percent_area <- function (percent_oligo_MPL,
                                                     Library_VC,
                                                     name) {
 MPL_oligo <- add_date_to_oligotype_table(percent_oligo_MPL,
                                          Library_VC)
 plot_over_time <- ggplot(MPL_oligo, aes(x=Date, y=value, group=variable))+
  season_line+
  HAKA_bloom+
  spring_bloom_line+
  geom_area(aes(fill=variable))+
  ylab("Reads")+
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank())
 print(plot_over_time)
 pdf(paste0("../figures/", name, "all_times.pdf"), width = 15, height = 11)
 print(plot_over_time)
 dev.off()
 return(plot_over_time)
}

plot_oligotypes_over_time_percent_lines_log_y <- function (percent_oligo_MPL,
                                                           Library_VC,
                                                           name) {
 MPL_oligo <- add_date_to_oligotype_table(percent_oligo_MPL,
                                          Library_VC)
 plot_over_time <- ggplot(MPL_oligo, aes(x=Date, y=value, group=variable))+
  scale_y_log10()+
  season_line+
  HAKA_bloom+
  spring_bloom_line+
  geom_line(aes(colour=variable), size=2)+
  ylab("Reads")+
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank())
 print(plot_over_time)
 pdf(paste0("../figures/", name, "_all_times.pdf"), width = 15, height = 11)
 print(plot_over_time)
 dev.off()
 return(plot_over_time)
}


percent_oligo_raphidopyte <- read.delim("../results/oligotyping_raphidophyte_all_times_OTUs/MATRIX-PERCENT.txt")

count_oligo_raphidopyte <- read.delim("../results/oligotyping_raphidophyte_all_times_OTUs/MATRIX-COUNT.txt")


plot_oligotypes_over_time_percent(percent_oligo_raphidopyte ,
                                  Library_metadata$X18s_VC_number,
                                  "percent_oligo_raphidophyte") 

plot_oligotypes_over_time_percent_lines(percent_oligo_raphidopyte ,
                                        Library_metadata$X18s_VC_number,
                                        "percent_oligo_raphidophyte_lines") 

count_raph_all <- plot_oligotypes_over_time_percent_lines(count_oligo_raphidopyte ,
                                        Library_metadata$X18s_VC_number,
                                        "count_oligo_raphidophyte")

plot_oligotypes_over_time_percent_lines_log_y(count_oligo_raphidopyte ,
                                              Library_metadata$X18s_VC_number,
                                              "count_oligo_raphidophyte_log_y") 


#### cyano
# percent_oligo_cyano <- read.delim("../results/oligotyping_cyano_otus_all_times/MATRIX-PERCENT.txt")
# 
# count_oligo_cyano <- read.delim("../results/oligotyping_cyano_otus_all_times/MATRIX-COUNT.txt")
# 
# plot_oligotypes_over_time_percent(percent_oligo_cyano,
#                                   Library_metadata$X16s_VC_number,
#                                   "percent_oligo_cyano") 
# 
# plot_oligotypes_over_time_percent_lines(percent_oligo_cyano ,
#                                         Library_metadata$X16s_VC_number,
#                                         "percent_oligo_cyano_lines") 
# 
# plot_oligotypes_over_time_percent_lines(count_oligo_cyano ,
#                                         Library_metadata$X16s_VC_number,
#                                         "count_oligo_cyano") 
# 
# plot_oligotypes_over_time_percent_lines_log_y(count_oligo_cyano ,
#                                               Library_metadata$X16s_VC_number,
#                                               "count_oligo_cyano_log_y") 

### Dino over all time series 

percent_oligo_dino <- read.delim("../results/dino_all_times/MATRIX-PERCENT.txt")

count_oligo_dino <- read.delim("../results/dino_all_times/MATRIX-COUNT.txt")

plot_oligotypes_over_time_percent(percent_oligo_dino ,
                                  Library_metadata$X18s_VC_number,
                                  "percent_oligo_dino") 

plot_oligotypes_over_time_percent_lines(percent_oligo_dino ,
                                        Library_metadata$X18s_VC_number,
                                        "percent_oligo_dino_lines") 

count_dino1_all <- plot_oligotypes_over_time_percent_lines(count_oligo_dino ,
                                        Library_metadata$X18s_VC_number,
                                        "count_oligo_dino") 

## wanted to remove the jne ones so that I could see more of the other variation. 
plot_oligotypes_over_time_percent_lines_log_y(count_oligo_dino ,
                                              Library_metadata$X18s_VC_number,
                                              "count_oligo_dino_log_y") 


## Dino 2 

percent_oligo_dino2 <- read.delim("../results/dino2_all_times/MATRIX-PERCENT.txt")

count_oligo_dino2 <- read.delim("../results/dino2_all_times/MATRIX-COUNT.txt")

plot_oligotypes_over_time_percent(percent_oligo_dino2 ,
                                  Library_metadata$X18s_VC_number,
                                  "percent_oligo_dino2") 

plot_oligotypes_over_time_percent_lines(percent_oligo_dino2 ,
                                        Library_metadata$X18s_VC_number,
                                        "percent_oligo_dino2_lines") 

count_dino2_all <- plot_oligotypes_over_time_percent_lines(count_oligo_dino2 ,
                                        Library_metadata$X18s_VC_number,
                                        "count_oligo_dino2") 

## wanted to remove the jne ones so that I could see more of the other variation. 
plot_oligotypes_over_time_percent_lines_log_y(count_oligo_dino2 ,
                                              Library_metadata$X18s_VC_number,
                                              "count_oligo_dino2_log_y") 


##### MPL OTU 3
# 
percent_oligo_MPL_OTU3 <- read.delim("../results/oligotyping_MPL_OTU_3_nuc_all_times/MATRIX-PERCENT.txt")

count_oligo_MPL_OTU3  <- read.delim("../results/oligotyping_MPL_OTU_3_nuc_all_times/MATRIX-COUNT.txt")


plot_oligotypes_over_time_percent(percent_oligo_MPL_OTU3 ,
                                  Library_metadata$MPL_VC_number,
                                  "percent_oligo_MPL3") 

plot_oligotypes_over_time_percent_lines(percent_oligo_MPL_OTU3 ,
                                        Library_metadata$MPL_VC_number,
                                        "percent_oligo_MPL3_lines") 

all_time_MPL_3 <- plot_oligotypes_over_time_percent_lines(count_oligo_MPL_OTU3 ,
                                        Library_metadata$MPL_VC_number,
                                        "count_oligo_MPL3") 

plot_oligotypes_over_time_percent_area(count_oligo_MPL_OTU3,
                                       Library_metadata$MPL_VC_number,
                                       "count_oligo_MPL3_area")

# ## wanted to remove the jne ones so that I could see more of the other variation. 
all_time_MPL_3_log_y <- plot_oligotypes_over_time_percent_lines_log_y(count_oligo_MPL_OTU3 ,
                                              Library_metadata$MPL_VC_number,
                                              "count_oligo_MPL3_log_y") 

# 
# ##### MPL OTU 4
# 
# percent_oligo_MPL_OTU4 <- read.delim("../results/oligotyping_MPL_OTU_4_nuc_all_times/MATRIX-PERCENT.txt")
# 
# count_oligo_MPL_OTU4  <- read.delim("../results/oligotyping_MPL_OTU_4_nuc_all_times/MATRIX-COUNT.txt")
# 
# 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_OTU4 ,
#                                   Library_metadata$MPL_VC_number,
#                                   "percent_oligo_MPL4") 
# 
# plot_oligotypes_over_time_percent_lines(percent_oligo_MPL_OTU4 ,
#                                         Library_metadata$MPL_VC_number,
#                                         "percent_oligo_MPL4_lines") 
# 
# plot_oligotypes_over_time_percent_lines(count_oligo_MPL_OTU4 ,
#                                         Library_metadata$MPL_VC_number,
#                                         "count_oligo_MPL4") 
# 
# plot_oligotypes_over_time_percent_area(count_oligo_MPL_OTU4,
#                                        Library_metadata$MPL_VC_number,
#                                        "count_oligo_MPL4_area")
# 
# ## wanted to remove the jne ones so that I could see more of the other variation. 
# plot_oligotypes_over_time_percent_lines_log_y(count_oligo_MPL_OTU4 ,
#                                               Library_metadata$MPL_VC_number,
#                                               "count_oligo_MPL4_log_y") 
# 
# ##### MPL OTU 5
# 
# percent_oligo_MPL_OTU5 <- read.delim("../results/oligotyping_MPL_OTU_5_nuc_all_times/MATRIX-PERCENT.txt")
# 
# count_oligo_MPL_OTU5  <- read.delim("../results/oligotyping_MPL_OTU_5_nuc_all_times/MATRIX-COUNT.txt")
# 
# 
# plot_oligotypes_over_time_percent(percent_oligo_MPL_OTU5 ,
#                                   Library_metadata$MPL_VC_number,
#                                   "percent_oligo_MPL5") 
# 
# plot_oligotypes_over_time_percent_lines(percent_oligo_MPL_OTU5 ,
#                                         Library_metadata$MPL_VC_number,
#                                         "percent_oligo_MPL5_lines") 
# 
# plot_oligotypes_over_time_percent_lines(count_oligo_MPL_OTU5 ,
#                                         Library_metadata$MPL_VC_number,
#                                         "count_oligo_MPL5") 
# 
# plot_oligotypes_over_time_percent_area(count_oligo_MPL_OTU5,
#                                        Library_metadata$MPL_VC_number,
#                                        "count_oligo_MPL5_area")
# 
# ## wanted to remove the jne ones so that I could see more of the other variation. 
# plot_oligotypes_over_time_percent_lines_log_y(count_oligo_MPL_OTU5 ,
#                                               Library_metadata$MPL_VC_number,
#                                               "count_oligo_MPL5_log_y") 
# 
# 
##### MPL OTU 288

percent_oligo_MPL_OTU288 <- read.delim("../results/oligotyping_MPL_OTU_288_nuc_all_times/MATRIX-PERCENT.txt")

count_oligo_MPL_OTU288  <- read.delim("../results/oligotyping_MPL_OTU_288_nuc_all_times/MATRIX-COUNT.txt")


plot_oligotypes_over_time_percent(percent_oligo_MPL_OTU288 ,
                                  Library_metadata$MPL_VC_number,
                                  "percent_oligo_MPL288") 

plot_oligotypes_over_time_percent_lines(percent_oligo_MPL_OTU288 ,
                                        Library_metadata$MPL_VC_number,
                                        "percent_oligo_MPL288_lines") 

all_time_MPL_288 <- plot_oligotypes_over_time_percent_lines(count_oligo_MPL_OTU288 ,
                                        Library_metadata$MPL_VC_number,
                                        "count_oligo_MPL288") 

plot_oligotypes_over_time_percent_area(count_oligo_MPL_OTU288,
                                       Library_metadata$MPL_VC_number,
                                       "count_oligo_MPL288_area")

## wanted to remove the jne ones so that I could see more of the other variation. 
all_time_MPL_288_log_y <- plot_oligotypes_over_time_percent_lines_log_y(count_oligo_MPL_OTU288 ,
                                              Library_metadata$MPL_VC_number,
                                              "count_oligo_MPL288_log_y") 



# pdf("../figures/oligotypes_bloom_and_overall_euks%01d.pdf", width=20, height = 20, onefile=FALSE)
# plot_grid(count_raph_all+xlab(NULL)+
#            theme(axis.ticks.x = element_blank(),
#                  axis.text.x = element_blank(),
#                  axis.title.x=element_blank()),
# count_raph_bloom+xlab(NULL)+
#  theme(axis.ticks.x = element_blank(),
#        axis.text.x = element_blank(),
#        axis.title.x=element_blank()),
# 
# count_dino1_all+xlab(NULL)+
#  theme(axis.ticks.x = element_blank(),
#        axis.text.x = element_blank(),
#        axis.title.x=element_blank()),
# count_dino1_bloom+xlab(NULL)+
#  theme(axis.ticks.x = element_blank(),
#        axis.text.x = element_blank(),
#        axis.title.x=element_blank()),
# count_dino2_all,
# count_dino2_bloom,
# ncol=2, 
# labels = c("A", "B", "C", "D", "E", "F"))
# dev.off()


pdf("../figures/oligotypes_bloom_and_overall_euks%01d.pdf", width=20, height = 20, onefile=FALSE)
plot_grid(
          count_raph_all+
           xlab(NULL)+
           scale_x_date(breaks = date_breaks("2 days"),
                        date_labels = "%b %d"
           )+
           coord_cartesian(xlim = c(as.Date("2011-06-20"),
                                    as.Date("2011-07-06")))+
           theme(legend.position="none")+
           theme(axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank()),
          count_raph_all+
           xlab(NULL)+
           theme(axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank())+
           theme(legend.position="left"),
          count_dino1_all+
           xlab(NULL)+
           scale_x_date(breaks = date_breaks("2 days"),
                        date_labels = "%b %d"
           )+
           coord_cartesian(xlim = c(as.Date("2011-06-20"),
                                    as.Date("2011-07-06")))+
           theme(legend.position="none")+
           theme(axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank()),
          count_dino1_all+
           xlab(NULL)+
           theme(axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank())+
           theme(legend.position="left"),
          count_dino2_all+
           theme(legend.position="none")+
           scale_x_date(breaks = date_breaks("2 days"),
                        date_labels = "%b %d"
           )+
           coord_cartesian(xlim = c(as.Date("2011-06-20"),
                                    as.Date("2011-07-06"))),
          count_dino2_all+
           theme(legend.position="left"),
          ncol=2,
          align="hv",
          labels = c("A", "B", "C", "D", "E", "F"),
          label_size = 24)
dev.off()



pdf("../figures/oligotypes_bloom_and_overall_MPL%01d.pdf", width=20, height = 13, onefile=FALSE)
plot_grid(all_time_MPL_288+
           scale_x_date(breaks = date_breaks("2 days"),
                        date_labels = "%b %d"
                        )+
           coord_cartesian(xlim = c(as.Date("2011-06-20"),
                                   as.Date("2011-07-06")))+
           theme(legend.position="none"),        
          
          
          all_time_MPL_288+
           theme(legend.position="left"),
          all_time_MPL_3+
           scale_x_date(breaks = date_breaks("2 days"),
                        date_labels = "%b %d"
           )+
           coord_cartesian(xlim = c(as.Date("2011-06-20"),
                                    as.Date("2011-07-06")))+
           theme(legend.position="none"),       
          all_time_MPL_3+
           theme(legend.position="left"),
          ncol=2,
          align="hv",
          rel_widths=c(1,1.5),
          labels = c("A", "B", "C", "D"),
          label_size = 24)
dev.off()


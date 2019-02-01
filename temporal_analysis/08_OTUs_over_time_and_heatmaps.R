## Make plots of individual OTUs and heatmaps of all otus over time. 
## author: Julia Gustavsen

library(ggplot2)
library(plyr)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(scales)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv",
                                  row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",
                                   row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv",
                                  row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv",
                                  row.names="VC_number") 

args <- commandArgs(TRUE)
inputFile <- args[1]

## test to see if input file is given, so I can decide whether to use this argument or the orginal one. 
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

### Get relevant VC numbers ====
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
## Reformat date
## %d is day as a number,
## %b is abreviated month in words,
## %y is 2digit year
Jericho_data$Date <- as.Date(Jericho_data$Date)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv",
                             nrows=61)

high_res_vcs <- c(1198,
                  1199,
                  1200,
                  1201,
                  1202)

## want Jericho_dates that are not in MPL table and high res
Jericho_no_high_res_vcs <- Jericho_data$VC_number[!(Jericho_data$VC_number %in% high_res_vcs)]

## want to annotate missing data days. 
## want to annotate those days using star
missing_MPL_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_MPL_OTUs))]
missing_MPL_dates <- Jericho_data$Date[match(missing_MPL_samples,
                                             Jericho_data$VC_number)]

missing_gp23_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_gp23_OTUs))]
missing_gp23_dates <- Jericho_data$Date[match(missing_gp23_samples,
                                              Jericho_data$VC_number)]

missing_18s_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_18s_OTUs))]
missing_18s_dates <- Jericho_data$Date[match(missing_18s_samples,
                                             Jericho_data$VC_number)]

missing_16s_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_16s_OTUs))]
missing_16s_dates <- Jericho_data$Date[match(missing_16s_samples,
                                             Jericho_data$VC_number)]



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
long_MPL_otus <- add_date_to_OTU_table_rtrn_long(normalized_MPL_OTUs,
                                                 Jericho_data)
long_16s_otus <- add_date_to_OTU_table_rtrn_long(normalized_16s_OTUs,
                                                 Jericho_data)

## how many OTUs?
Number_of_OTUs_gp23 <- n_distinct(long_gp23_otus$OTUid)
Number_of_OTUs_18s <- n_distinct(long_18s_otus$OTUid)
Number_of_OTUs_MPL <- n_distinct(long_MPL_otus$OTUid)
Number_of_OTUs_16s <- n_distinct(long_16s_otus$OTUid)

## pull out the counts of the number of OTUs-richness over time

richness_of_otus_over_time <- function (long_otus) {
  summarise_richness_otus <- long_otus %>%
    filter(value > 0) %>%
    ddply(~Date,
          summarise,
          richness = length(OTUid))
}

find_top20_OTUs <- function (long_otus) {
  summarised_by_OTU <- ddply(long_otus,
                             ~OTUid,
                             summarise,
                             relative_abundance = sum(value))
  OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
  top_20 <- droplevels(head(OTU_sorted_by_rel_abun,
                            n=20))
  return(top_20$OTUid)
} 

barplot_from_long_otus <- function(long_OTUs,
                                   title){
  number_of_otus <- length(unique(long_OTUs$OTUid))
  bar_colour <- colorRampPalette(brewer.pal(n =12,
                                            name = "Set3"))(number_of_otus)
  long_OTUs$OTUid <- gsub(".size.*.",
                          "",
                          long_OTUs$OTUid)
  OTU_barplot <- ggplot(long_OTUs,
                        aes(x=Date,
                            y=value,
                            group=OTUid))+ 
    season_line +
  #  spring_bloom_line+
    geom_bar(aes(fill=OTUid),
             stat="identity", 
             position="stack",
             colour="black")+
    scale_fill_manual(values=bar_colour)+
    theme_JAG_presentation()+
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_blank())+
    date_scaling 
  # +
  #  ggtitle(title)
  return(OTU_barplot)
}


barplot_of_top_20 <- function(otu_table,
                              Jericho_and_SOG_data,
                              Jericho_data,
                              title){
  proportional_otus <- as.data.frame(prop.table(as.matrix(otu_table),
                                                margin=1))
  long_otus <- add_date_to_OTU_table_rtrn_long(proportional_otus,
                                               Jericho_data)
  long_otus$OTUid <- gsub(".size.*.",
                                        "",
                                        long_otus$OTUid)
  vector_top_20_OTUs <- find_top20_OTUs(long_otus)
  top_20_long_OTUs <- droplevels(subset(long_otus,
                                        OTUid %in% vector_top_20_OTUs))
  heatmap_top_20_OTUs <- barplot_from_long_otus(top_20_long_OTUs,
                                                title)
}

MPL_barplot <- barplot_of_top_20(normalized_MPL_OTUs,
                                 Library_metadata,
                                 Jericho_data,
                                 "MPL")+
  annotate("text",
           x = missing_MPL_dates,
           y = 0.00,
           label="x",
           colour = "grey",
           size=14)
gp23_barplot <- barplot_of_top_20(normalized_gp23_OTUs,
                                  Library_metadata,
                                  Jericho_data,
                                  "gp23")+
  annotate("text",
           x = missing_gp23_dates,
           y = 0.00,
           label="X",
           size=7,
           colour = "grey")

S16_barplot <- barplot_of_top_20(normalized_16s_OTUs,
                                 Library_metadata,
                                 Jericho_data,
                                 "16S")+
  annotate("text",
           x = missing_16s_dates,
           y = 0.00,
           label="X",
           size=7,
           colour = "grey")


S18_barplot <- barplot_of_top_20(normalized_18s_OTUs,
                                 Library_metadata,
                                 Jericho_data,
                                 "18S")+
  annotate("text",
           x = missing_18s_dates,
           y = 0.00,
           label="X",
           size=7,
           colour = "grey")

pdf(paste0(figures_dir,
           "time_series_top_20_OTUs_barplots_all_amplicons%03d.pdf"),
    width = 20,
    height = 15,
    onefile = FALSE)
gp23_barplot
MPL_barplot
S16_barplot
S18_barplot
dev.off()

## no high res
high_res_vcs <- c(1198,
                  1199,
                  1200,
                  1201, 
                  1202)
MPL_OTUs_no_high_res <- subset(normalized_MPL_OTUs,
                               !(rownames(normalized_MPL_OTUs) %in% high_res_vcs))
gp23_OTUs_no_high_res <- subset(normalized_gp23_OTUs,
                                !(rownames(normalized_gp23_OTUs) %in% high_res_vcs))
S18_OTUs_no_high_res <- subset(normalized_18s_OTUs,
                               !(rownames(normalized_18s_OTUs) %in% high_res_vcs))
S16_OTUs_no_high_res <- subset(normalized_16s_OTUs,
                               !(rownames(normalized_16s_OTUs) %in% high_res_vcs))


MPL_barplot <- barplot_of_top_20(MPL_OTUs_no_high_res,
                                 Library_metadata,
                                 Jericho_data,
                                 "MPL")+
  annotate("text",
           x = missing_MPL_dates,
           y = 0.00,
           label="X",
           size=7,
           colour = "grey")

gp23_barplot <- barplot_of_top_20(gp23_OTUs_no_high_res,
                                  Library_metadata,
                                  Jericho_data,
                                  "gp23")+
  annotate("text",
           x = missing_gp23_dates,
           y = 0.00,
           label="X",
           size=7,
           colour = "grey")

S16_barplot <- barplot_of_top_20(S16_OTUs_no_high_res,
                                 Library_metadata,
                                 Jericho_data,
                                 "16S")+
  annotate("text",
           x = missing_16s_dates,
           y = 0.00,
           label="X",
           size=7,
           colour = "grey")

S18_barplot <- barplot_of_top_20(S18_OTUs_no_high_res,
                                 Library_metadata,
                                 Jericho_data,
                                 "18S")+
  annotate("text",
           x = missing_18s_dates,
           y = 0.00,
           label="X",
           size=7,
           colour = "grey")
pdf(paste0(figures_dir,
           "time_series_top_20_OTUs_barplots_all_amplicons_no_high_res%03d.pdf"),
    width = 20,
    height = 15,
    onefile = FALSE)
gp23_barplot
MPL_barplot
S16_barplot
S18_barplot
dev.off()

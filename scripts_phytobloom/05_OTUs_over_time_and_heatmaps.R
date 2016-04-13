## WOrk with time series OTUs
## Make plots of individual OTUs and heatmaps of all otus over time. 

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
library(cowplot)

## Need to edit and make sure the font, type are consistent
## To-do:
## make it easier to change what the 18s and 16s are facetted by. 


normalized_18s_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_gp23.tsv", row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_MPL.tsv", row.names="VC_number")
#normalized_AVS_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
# normalized_AVS_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_AVS1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_AVS_R2.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

colnames(normalized_16s_OTUs) <- gsub(".size.*.", "", colnames(normalized_16s_OTUs))

colnames(normalized_18s_OTUs) <- gsub(".size.*.", "", colnames(normalized_18s_OTUs))

## not sure how I generated these!!
#normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_bloom_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
#normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_bloom_time_series_18s_normalized_Heterotrophs.tsv",
#                                         row.names=1)

total_number_reads <- sum(normalized_18s_OTUs, normalized_MPL_OTUs, 
                          #normalized_AVS_OTUs, 
                          normalized_gp23_OTUs, normalized_16s_OTUs)

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

load(file="../../phyla_palette_euks.txt")
load(file="../../order_palette_euks.txt")
load(file="../../family_palette_euks.txt")
load(file="../../phyla_palette_bac.txt")
load(file="../../order_palette_bac.txt")
load(file="../../family_palette_bac.txt")

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
date_scaling <-   scale_x_date(limits = c(as.Date("2011-06-20"), as.Date("2011-07-06")))


### Get relevant VC numbers ====
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
## Reformat date
# %d is day as a number, %b is abreviated month in words, %y is 2digit year
Jericho_data$Date <- as.Date(Jericho_data$Date)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)

## first rename the column names in the OTU table to something better
add_date_to_OTU_table_rtrn_long <- function (otu_table, 
                                             Jericho_data) {
 otu_table$VC_number = rownames(otu_table)
 long_otus <- melt(otu_table, id="VC_number", variable.name="OTUid")
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
#long_18s_otus_phytos <- add_date_to_OTU_table_rtrn_long(normalized_18s_OTUs_phytos,
#                                                        Jericho_data)
#long_18s_otus_hetero <- add_date_to_OTU_table_rtrn_long(normalized_18s_OTUs_hetero,
#                                                        Jericho_data)
#long_AVS_otus <- add_date_to_OTU_table_rtrn_long(normalized_AVS_OTUs,
 #                                                Jericho_data)
long_MPL_otus <- add_date_to_OTU_table_rtrn_long(normalized_MPL_OTUs,
                                                 Jericho_data)
long_16s_otus <- add_date_to_OTU_table_rtrn_long(normalized_16s_OTUs,
                                                 Jericho_data)


## Plot single OTUs over time
plot_individual_otus_over_time <- function (long_otus, title) {
 colourCount  <- length(unique(long_otus$OTUid))
 getPalette  <-  colorRampPalette(brewer.pal(11, "Spectral"))
 p <- ggplot(droplevels(long_otus), aes(x=Date, y=value, group=OTUid))+
  season_line +
  spring_bloom_line+
  geom_area(aes(fill=OTUid), position="stack", colour="black")+
  date_scaling +
  scale_fill_manual(values = getPalette(colourCount))+
  theme_JAG_presentation()+
  ggtitle(title)
}

## how many OTUs?
Number_of_OTUs_gp23 <- n_distinct(long_gp23_otus$OTUid)
Number_of_OTUs_18s <- n_distinct(long_18s_otus$OTUid)
Number_of_OTUs_MPL <- n_distinct(long_MPL_otus$OTUid)
#Number_of_OTUs_AVS <- n_distinct(long_AVS_otus$OTUid)
Number_of_OTUs_16s <- n_distinct(long_16s_otus$OTUid)

## pull out the counts of the number of OTUs-richness over time

richness_of_otus_over_time <- function (long_otus) {
 summarise_richness_otus <- long_otus %>%
  filter(value > 0) %>%
  ddply(~Date,summarise,richness=length(OTUid))
}

find_top10_OTUs <- function (long_otus) {
 summarised_by_OTU <- ddply(long_otus,~OTUid,summarise,relative_abundance=sum(value))
 OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
 top_10 <- droplevels(head(OTU_sorted_by_rel_abun, n=10))
 return(top_10$OTUid)
} 

find_top20_OTUs <- function (long_otus) {
 summarised_by_OTU <- ddply(long_otus,~OTUid,summarise,relative_abundance=sum(value))
 OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
 top_20 <- droplevels(head(OTU_sorted_by_rel_abun, n=20))
 return(top_20$OTUid)
}

heatmap_from_long_otus <- function(long_OTUs, high_colour, title){
 OTU_heatmap <- ggplot(long_OTUs, aes(x=Date,y=OTUid)) + 
  season_line +
  spring_bloom_line+
  geom_tile(aes(fill = value),
            colour = "grey10") +
  scale_fill_gradient(low = line_colour,
                      high = high_colour,
                      na.value = line_colour,
                      trans = "log"
  )+ 
  theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  date_scaling +
  ggtitle(title)
 
 return(OTU_heatmap)
}


function_to_do_all_prop <- function(otu_table, Jericho_and_SOG_data,
                                    Jericho_data,
                                    primer, 
                                    clustering_method, title){
 proportional_otus <- as.data.frame(prop.table(as.matrix(otu_table),  margin=1))
 
 long_otus <- add_date_to_OTU_table_rtrn_long(proportional_otus,
                                              Jericho_data) 
 total_count_of_OTUs <- n_distinct(long_otus$OTUid)
 plot_over_time <- plot_individual_otus_over_time(long_otus, paste(title, total_count_of_OTUs, sep="")) 
 pdf(paste("../figures/", primer, "_individual_OTUs_time_series_proportional%03d.pdf", sep=""), width = 15, height = 11, onefile = FALSE)
 print(plot_over_time)
 dev.off()
}

gp23_all <- function_to_do_all_prop(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23", "Usearch", "gp23 ")

MPL_all <- function_to_do_all_prop(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL", "Usearch", "MPL ")

# AVS_all <- function_to_do_all_prop(normalized_AVS_OTUs, Library_metadata, Jericho_data, "AVS", "Usearch", "AVS ")

S16_all <- function_to_do_all_prop(normalized_16s_OTUs, Library_metadata, Jericho_data, "16S", "Usearch", "16S ")

S18s_all <- function_to_do_all_prop(normalized_18s_OTUs, Library_metadata, Jericho_data, "18S", "Usearch", "18S ")

#S18s_phytos <- function_to_do_all_prop(normalized_18s_OTUs_phytos, Library_metadata, Jericho_data, "18S_phytos", "Usearch", "18S ")
#S18s_hetero <- function_to_do_all_prop(normalized_18s_OTUs_hetero, Library_metadata, Jericho_data, "18S", "Usearch", "18S_hetero ")

### print all heatmaps together
heatmaps_of_top_10 <- function(otu_table,Jericho_and_SOG_data, Jericho_data, title){
 long_otus <- add_date_to_OTU_table_rtrn_long(otu_table,
                                              Jericho_data)
 vector_top_10_OTUs <- find_top10_OTUs(long_otus)
 top_10_long_OTUs <- droplevels(subset(long_otus, OTUid %in% vector_top_10_OTUs))
 heatmap_top_10_OTUs <- heatmap_from_long_otus(top_10_long_OTUs, "chartreuse", title)
}

MPL_heatmap <- heatmaps_of_top_10(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL")

# AVS_heatmap <- heatmaps_of_top_10(normalized_AVS_OTUs, Library_metadata, Jericho_data, "AVS")

gp23_heatmap <- heatmaps_of_top_10(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23")

S16_heatmap <- heatmaps_of_top_10(normalized_16s_OTUs, Library_metadata,Jericho_data, "16S")

S18_heatmap <- heatmaps_of_top_10(normalized_18s_OTUs, Library_metadata,Jericho_data, "18S")

pdf("../figures/time_series_top_10_OTUs_heatmaps_all_amplicons%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(
 #AVS_heatmap, 
 gp23_heatmap,MPL_heatmap,S16_heatmap,S18_heatmap,ncol=1)
dev.off()

barplot_from_long_otus <- function(long_OTUs, title){
 number_of_otus <- length(unique(long_OTUs$OTUid))
 bar_colour <- colorRampPalette(brewer.pal(n =12 , name = "Set3"))(number_of_otus)
 OTU_barplot <- ggplot(long_OTUs, aes(x=Date,y=value, group=OTUid))+ 
  season_line +
  spring_bloom_line+
  geom_bar(aes(fill=OTUid),
           stat="identity", 
           position="stack",
           colour="black")+
  scale_fill_manual(values=bar_colour)+
  theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  date_scaling +
  ggtitle(title)
 return(OTU_barplot)
}


barplot_of_top_10 <- function(otu_table,Jericho_and_SOG_data, Jericho_data, title){
 proportional_otus <- as.data.frame(prop.table(as.matrix(otu_table),  margin=1))
 long_otus <- add_date_to_OTU_table_rtrn_long(proportional_otus,
                                              Jericho_data) 
 vector_top_10_OTUs <- find_top20_OTUs(long_otus)
 top_10_long_OTUs <- droplevels(subset(long_otus, OTUid %in% vector_top_10_OTUs))
 heatmap_top_10_OTUs <- barplot_from_long_otus(top_10_long_OTUs, title)
}

MPL_barplot <- barplot_of_top_10(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL")
MPL_barplot

# AVS_barplot <- barplot_of_top_10(normalized_AVS_OTUs, Library_metadata, Jericho_data, "AVS")

gp23_barplot <- barplot_of_top_10(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23")

S16_barplot <- barplot_of_top_10(normalized_16s_OTUs, Library_metadata,Jericho_data, "16S")

S18_barplot <- barplot_of_top_10(normalized_18s_OTUs, Library_metadata,Jericho_data, "18S")

#S18_barplot_phyto <- barplot_of_top_10(normalized_18s_OTUs_phytos, Library_metadata,Jericho_data, "18S_phytos")

#S18_barplot_hetero <- barplot_of_top_10(normalized_18s_OTUs_hetero, Library_metadata,Jericho_data, "18S_hetero")

pdf("../figures/time_series_top_10_OTUs_barplots_all_amplicons%03d.pdf", width = 20, height = 15, onefile = FALSE)
#AVS_barplot
gp23_barplot
MPL_barplot
S16_barplot
S18_barplot
#S18_barplot_phyto
#S18_barplot_hetero
dev.off()

areaplot_from_long_otus <- function(long_OTUs, title){
 number_of_otus <- length(unique(long_OTUs$OTUid))
 bar_colour <- colorRampPalette(brewer.pal(n =12 , name = "Set3"))(number_of_otus)
 sampling_points <- unique(long_OTUs$Date)
 OTU_areaplot <- ggplot(long_OTUs, aes(x=Date,y=value, group=OTUid)) + 
  season_line +
  spring_bloom_line+
  geom_area(aes(fill=OTUid), position="stack", colour="black")+
  scale_fill_manual(values=bar_colour)+
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           solid=TRUE,
           size = 5)+
  theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  date_scaling +
  ggtitle(title)
 
 return(OTU_areaplot)
}

areaplot_of_top_10 <- function(otu_table,Jericho_and_SOG_data, Jericho_data, title){
 proportional_otus <- as.data.frame(prop.table(as.matrix(otu_table),  margin=1))
 long_otus <- add_date_to_OTU_table_rtrn_long(proportional_otus,
                                              Jericho_data) 
 vector_top_10_OTUs <- find_top20_OTUs(long_otus)
 top_10_long_OTUs <- droplevels(subset(long_otus, OTUid %in% vector_top_10_OTUs))
 heatmap_top_10_OTUs <- areaplot_from_long_otus(top_10_long_OTUs, title)
}

MPL_areaplot <- areaplot_of_top_10(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL")

# AVS_areaplot <- areaplot_of_top_10(normalized_AVS_OTUs, Library_metadata, Jericho_data, "AVS")

gp23_areaplot <- areaplot_of_top_10(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23")

S16_areaplot <- areaplot_of_top_10(normalized_16s_OTUs, Library_metadata,Jericho_data, "16S")

S18_areaplot <- areaplot_of_top_10(normalized_18s_OTUs, Library_metadata,Jericho_data, "18S")

pdf("../figures/time_series_top_10_OTUs_areaplots_all_amplicons%03d.pdf", width = 20, height = 15, onefile = FALSE)
#AVS_areaplot
gp23_areaplot
MPL_areaplot
S16_areaplot
S18_areaplot
dev.off()





### look at the top 20 OTus but do it individually and then line up

areaplots_of_each_of_top_10 <- function(otu_table,Jericho_and_SOG_data, Jericho_data, title){
 
 areaplot_from_long_single_otus <- function(long_OTUs, title, colour_area){
  number_of_otus <- length(unique(long_OTUs$OTUid))
  sampling_points <- unique(long_OTUs$Date)
  OTU_areaplot <- ggplot(long_OTUs, aes(x=Date,y=value, group=OTUid)) + 
   season_line +
   spring_bloom_line+
   geom_area(fill=colour_area, position="stack", colour="black")+
   annotate("point", 
            x = as.Date(sampling_points),
            y = -0.02,
            colour = "grey",
            shape= 17,
            solid=TRUE,
            size = 5)+
  # scale_fill_manual(values=colour_area)+
   theme_JAG_presentation()+
   theme(panel.grid.major.y = element_blank(),
         panel.grid.major.x = element_blank())+
   date_scaling +
   ggtitle(title)
  
  return(OTU_areaplot)
 }
 
 proportional_otus <- as.data.frame(prop.table(as.matrix(otu_table),  margin=1))
 long_otus <- add_date_to_OTU_table_rtrn_long(proportional_otus,
                                              Jericho_data) 
 vector_top_10_OTUs <- find_top10_OTUs(long_otus)
 top_10_long_OTUs <- droplevels(subset(long_otus, OTUid %in% vector_top_10_OTUs))
 bar_colour <- colorRampPalette(brewer.pal(n =12 , name = "Set3"))(length(vector_top_10_OTUs))
 list_of_plots <- list()
 ## try with top 10
 #for (otu in seq(1:10)){
  for (otu in seq(length(vector_top_10_OTUs))){
  print(vector_top_10_OTUs[otu])
  current_otu_name <- vector_top_10_OTUs[otu]
  current_otu_long_OTUs <- droplevels(subset(long_otus, OTUid %in% current_otu_name))
  print(bar_colour[otu])
  heatmap_one_of_top_10_OTUs <- areaplot_from_long_single_otus(current_otu_long_OTUs, current_otu_name, bar_colour[otu])
  heatmap_one_of_top_10_OTUs
  list_of_plots <- c(list_of_plots, list(heatmap_one_of_top_10_OTUs))
  
 }
return(list_of_plots)
}

list_area_MPL <- areaplots_of_each_of_top_10(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL")
list_area_gp23 <- areaplots_of_each_of_top_10(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23")
list_area_16s <- areaplots_of_each_of_top_10(normalized_16s_OTUs, Library_metadata, Jericho_data, "16s")
list_area_18s <- areaplots_of_each_of_top_10(normalized_18s_OTUs, Library_metadata, Jericho_data, "18s")

pdf("../figures/time_series_top_10_each_OTUs_areaplots_all_amplicons%03d.pdf", width = 20, height = 20, onefile = FALSE)
#AVS_areaplot
plot_grid(plotlist = list_area_gp23, ncol=2)
plot_grid(plotlist = list_area_MPL, ncol=2)
plot_grid(plotlist = list_area_16s, ncol=2)
plot_grid(plotlist = list_area_18s, ncol=2)
dev.off()



###################### 
### For everyting ####
######################

normalized_18s_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_gp23.tsv", row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_MPL.tsv", row.names="VC_number")
#normalized_AVS_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
# normalized_AVS_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_AVS1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_AVS_R2.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

colnames(normalized_18s_OTUs) <- gsub(".size.*.", "", colnames(normalized_18s_OTUs))

colnames(normalized_16s_OTUs) <- gsub(".size.*.", "", colnames(normalized_16s_OTUs))

## not sure how I generated these!!
#normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_bloom_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
#normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_bloom_time_series_18s_normalized_Heterotrophs.tsv",
#                                         row.names=1)

total_number_reads <- sum(normalized_18s_OTUs, normalized_MPL_OTUs, 
                          #normalized_AVS_OTUs, 
                          normalized_gp23_OTUs, normalized_16s_OTUs)

args <- commandArgs(TRUE)
inputFile <- args[1]


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

HAKA_bloom <- annotate("rect",xmin = as.Date("2011-06-20"),
                       xmax = as.Date("2011-07-06"), 
                       ymin = -Inf,
                       ymax = Inf,
                       alpha = 0.4,
                       fill = "burlywood3",
                       colour=NA) 

### Get relevant VC numbers ====
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis_all_time_series.csv", row.names=1)
## Reformat date
# %d is day as a number, %b is abreviated month in words, %y is 2digit year
Jericho_data$Date <- as.Date(Jericho_data$Date)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)

## first rename the column names in the OTU table to something better
add_date_to_OTU_table_rtrn_long <- function (otu_table, 
                                             Jericho_data) {
 otu_table$VC_number = rownames(otu_table)
 long_otus <- melt(otu_table, id="VC_number", variable.name="OTUid")
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
#long_18s_otus_phytos <- add_date_to_OTU_table_rtrn_long(normalized_18s_OTUs_phytos,
#                                                        Jericho_data)
#long_18s_otus_hetero <- add_date_to_OTU_table_rtrn_long(normalized_18s_OTUs_hetero,
#                                                        Jericho_data)
#long_AVS_otus <- add_date_to_OTU_table_rtrn_long(normalized_AVS_OTUs,
#                                                Jericho_data)
long_MPL_otus <- add_date_to_OTU_table_rtrn_long(normalized_MPL_OTUs,
                                                 Jericho_data)
long_16s_otus <- add_date_to_OTU_table_rtrn_long(normalized_16s_OTUs,
                                                 Jericho_data)


## Plot single OTUs over time
plot_individual_otus_over_time <- function (long_otus, title) {
 colourCount  <- length(unique(long_otus$OTUid))
 getPalette  <-  colorRampPalette(brewer.pal(11, "Spectral"))
 sampling_points <- unique(long_otus$Date)
 p <- ggplot(droplevels(long_otus), aes(x=Date, y=value, group=OTUid))+
  season_line +
  HAKA_bloom +
  spring_bloom_line+
  geom_area(aes(fill=OTUid), position="stack", colour="black")+
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           solid=TRUE,
           size = 5)+
  date_scaling +
  scale_fill_manual(values = getPalette(colourCount))+
  theme_JAG_presentation()+
  ggtitle(title)
}

## how many OTUs?
Number_of_OTUs_gp23 <- n_distinct(long_gp23_otus$OTUid)
Number_of_OTUs_18s <- n_distinct(long_18s_otus$OTUid)
Number_of_OTUs_MPL <- n_distinct(long_MPL_otus$OTUid)
#Number_of_OTUs_AVS <- n_distinct(long_AVS_otus$OTUid)
Number_of_OTUs_16s <- n_distinct(long_16s_otus$OTUid)

## pull out the counts of the number of OTUs-richness over time

richness_of_otus_over_time <- function (long_otus) {
 summarise_richness_otus <- long_otus %>%
  filter(value > 0) %>%
  ddply(~Date,summarise,richness=length(OTUid))
}

find_top10_OTUs <- function (long_otus) {
 summarised_by_OTU <- ddply(long_otus,~OTUid,summarise,relative_abundance=sum(value))
 OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
 top_10 <- droplevels(head(OTU_sorted_by_rel_abun, n=10))
 return(top_10$OTUid)
} 

find_top20_OTUs <- function (long_otus) {
 summarised_by_OTU <- ddply(long_otus,~OTUid,summarise,relative_abundance=sum(value))
 OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
 top_20 <- droplevels(head(OTU_sorted_by_rel_abun, n=20))
 return(top_20$OTUid)
}

heatmap_from_long_otus <- function(long_OTUs, high_colour, title){
 OTU_heatmap <- ggplot(long_OTUs, aes(x=Date,y=OTUid))+ 
  geom_tile(aes(fill = value),
            colour = "grey10") +
  scale_fill_gradient(low = line_colour,
                      high = high_colour,
                      na.value = line_colour,
                      trans = "log"
  )+ 
  theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  #date_scaling +
  ggtitle(title)
 
 return(OTU_heatmap)
}


function_to_do_all_prop <- function(otu_table, Jericho_and_SOG_data,
                                    Jericho_data,
                                    primer, 
                                    clustering_method, title){
 proportional_otus <- as.data.frame(prop.table(as.matrix(otu_table),  margin=1))
 
 long_otus <- add_date_to_OTU_table_rtrn_long(proportional_otus,
                                              Jericho_data) 
 total_count_of_OTUs <- n_distinct(long_otus$OTUid)
 plot_over_time <- plot_individual_otus_over_time(long_otus, paste(title, total_count_of_OTUs, sep="")) 
 pdf(paste("../figures/", primer, "_individual_OTUs_with_all_time_series_proportional%03d.pdf", sep=""), width = 15, height = 11, onefile = FALSE)
 print(plot_over_time)
 dev.off()
}

gp23_all <- function_to_do_all_prop(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23", "Usearch", "gp23 ")

MPL_all <- function_to_do_all_prop(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL", "Usearch", "MPL ")

# AVS_all <- function_to_do_all_prop(normalized_AVS_OTUs, Library_metadata, Jericho_data, "AVS", "Usearch", "AVS ")

S16_all <- function_to_do_all_prop(normalized_16s_OTUs, Library_metadata, Jericho_data, "16S", "Usearch", "16S ")

S18s_all <- function_to_do_all_prop(normalized_18s_OTUs, Library_metadata, Jericho_data, "18S", "Usearch", "18S ")

#S18s_phytos <- function_to_do_all_prop(normalized_18s_OTUs_phytos, Library_metadata, Jericho_data, "18S_phytos", "Usearch", "18S ")
#S18s_hetero <- function_to_do_all_prop(normalized_18s_OTUs_hetero, Library_metadata, Jericho_data, "18S", "Usearch", "18S_hetero ")

### print all heatmaps together
heatmaps_of_top_10 <- function(otu_table,Jericho_and_SOG_data, Jericho_data, title){
 long_otus <- add_date_to_OTU_table_rtrn_long(otu_table,
                                              Jericho_data)
 vector_top_10_OTUs <- find_top10_OTUs(long_otus)
 top_10_long_OTUs <- droplevels(subset(long_otus, OTUid %in% vector_top_10_OTUs))
 heatmap_top_10_OTUs <- heatmap_from_long_otus(top_10_long_OTUs, "chartreuse", title)
}

# 
# MPL_heatmap <- heatmaps_of_top_10(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL")
# 
# # AVS_heatmap <- heatmaps_of_top_10(normalized_AVS_OTUs, Library_metadata, Jericho_data, "AVS")
# 
# gp23_heatmap <- heatmaps_of_top_10(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23")
# 
# S16_heatmap <- heatmaps_of_top_10(normalized_16s_OTUs, Library_metadata,Jericho_data, "16S")
# 
# S18_heatmap <- heatmaps_of_top_10(normalized_18s_OTUs, Library_metadata,Jericho_data, "18S")
# 
# pdf("../figures/time_series_with_all_top_10_OTUs_heatmaps_all_amplicons%03d.pdf", width = 20, height = 15, onefile = FALSE)
# grid.arrange(
#  #AVS_heatmap, 
#  gp23_heatmap,MPL_heatmap,S16_heatmap,S18_heatmap,ncol=1)
# dev.off()

barplot_from_long_otus <- function(long_OTUs, title){
 number_of_otus <- length(unique(long_OTUs$OTUid))
 bar_colour <- colorRampPalette(brewer.pal(n =12 , name = "Set3"))(number_of_otus)
 OTU_barplot <- ggplot(long_OTUs, aes(x=Date,y=value, group=OTUid))+ 
  season_line +
  HAKA_bloom +
  spring_bloom_line+
  geom_bar(aes(fill=OTUid),
           stat="identity", 
           position="stack",
           colour="black")+
  scale_fill_manual(values=bar_colour)+
  theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  date_scaling +
  ggtitle(title)
 return(OTU_barplot)
}


barplot_of_top_10 <- function(otu_table,Jericho_and_SOG_data, Jericho_data, title){
 proportional_otus <- as.data.frame(prop.table(as.matrix(otu_table),  margin=1))
 long_otus <- add_date_to_OTU_table_rtrn_long(proportional_otus,
                                              Jericho_data) 
 vector_top_10_OTUs <- find_top20_OTUs(long_otus)
 top_10_long_OTUs <- droplevels(subset(long_otus, OTUid %in% vector_top_10_OTUs))
 heatmap_top_10_OTUs <- barplot_from_long_otus(top_10_long_OTUs, title)
}

MPL_barplot <- barplot_of_top_10(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL")
MPL_barplot

# AVS_barplot <- barplot_of_top_10(normalized_AVS_OTUs, Library_metadata, Jericho_data, "AVS")

gp23_barplot <- barplot_of_top_10(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23")

S16_barplot <- barplot_of_top_10(normalized_16s_OTUs, Library_metadata,Jericho_data, "16S")

S18_barplot <- barplot_of_top_10(normalized_18s_OTUs, Library_metadata,Jericho_data, "18S")

#S18_barplot_phyto <- barplot_of_top_10(normalized_18s_OTUs_phytos, Library_metadata,Jericho_data, "18S_phytos")

#S18_barplot_hetero <- barplot_of_top_10(normalized_18s_OTUs_hetero, Library_metadata,Jericho_data, "18S_hetero")

pdf("../figures/time_series_with_all_top_10_OTUs_barplots_all_amplicons%03d.pdf", width = 20, height = 15, onefile = FALSE)
#AVS_barplot
gp23_barplot
MPL_barplot
S16_barplot
S18_barplot
#S18_barplot_phyto
#S18_barplot_hetero
dev.off()

areaplot_from_long_otus <- function(long_OTUs, title){
 number_of_otus <- length(unique(long_OTUs$OTUid))
 bar_colour <- colorRampPalette(brewer.pal(n =12 , name = "Set3"))(number_of_otus)
 sampling_points <- unique(long_OTUs$Date)
 OTU_areaplot <- ggplot(long_OTUs, aes(x=Date,y=value, group=OTUid)) + 
  season_line +
  HAKA_bloom +
  spring_bloom_line+
  geom_area(aes(fill=OTUid), position="stack", colour="black")+
  scale_fill_manual(values=bar_colour)+
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           solid=TRUE,
           size = 5)+
  theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  date_scaling +
  # ggtitle(title)+
  ylab("Relative abundance")
 
 return(OTU_areaplot)
}

areaplot_of_top_10 <- function(otu_table,Jericho_and_SOG_data, Jericho_data, title){
 proportional_otus <- as.data.frame(prop.table(as.matrix(otu_table),  margin=1))
 long_otus <- add_date_to_OTU_table_rtrn_long(proportional_otus,
                                              Jericho_data) 
 vector_top_10_OTUs <- find_top20_OTUs(long_otus)
 top_10_long_OTUs <- droplevels(subset(long_otus, OTUid %in% vector_top_10_OTUs))
 heatmap_top_10_OTUs <- areaplot_from_long_otus(top_10_long_OTUs, title)
}

MPL_areaplot_all <- areaplot_of_top_10(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL")

# AVS_areaplot <- areaplot_of_top_10(normalized_AVS_OTUs, Library_metadata, Jericho_data, "AVS")

gp23_areaplot_all <- areaplot_of_top_10(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23")

S16_areaplot_all <- areaplot_of_top_10(normalized_16s_OTUs, Library_metadata,Jericho_data, "16S")

S18_areaplot_all <- areaplot_of_top_10(normalized_18s_OTUs, Library_metadata,Jericho_data, "18S")

pdf("../figures/time_series_with_all_top_10_OTUs_areaplots_all_amplicons%03d.pdf", width = 20, height = 15, onefile = FALSE)
#AVS_areaplot
gp23_areaplot_all
MPL_areaplot_all
S16_areaplot_all
S18_areaplot_all
dev.off()


pdf("../figures/time_series_with_all_top_10_OTUs_areaplots_all_amplicons_together.pdf", width = 20, height = 18, onefile = FALSE)
plot_grid(gp23_areaplot_all+
           theme(legend.position="none",
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank())+
           scale_x_date(limits = c(as.Date("2011-06-20"),
                                   as.Date("2011-07-06")))+
           xlab(NULL),
          gp23_areaplot_all+ 
           theme(legend.position="none",
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank()
                 )+
           xlab(NULL),
          S16_areaplot_all+
           theme(legend.position="none",
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank())+
           scale_x_date(limits = c(as.Date("2011-06-20"),
                                   as.Date("2011-07-06")))+
           xlab(NULL),
          
          S16_areaplot_all+
           theme(legend.position="none",
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank()
                 )+
           xlab(NULL),
          S18_areaplot_all+
           theme(legend.position="none",
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank())+
           scale_x_date(limits = c(as.Date("2011-06-20"),
                                   as.Date("2011-07-06")))+
           xlab(NULL),         
          S18_areaplot_all+
           theme(legend.position="none",
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank())+
           xlab(NULL),
          MPL_areaplot_all+
           scale_x_date(limits = c(as.Date("2011-06-20"),
                                   as.Date("2011-07-06")))+
           theme(legend.position="none"), 
          MPL_areaplot_all+
           theme(legend.position="none",
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank()), 
          ncol=2,
          align="hv", labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          label_size = 24)
dev.off()


### can I try just looking at March-August

gp23_areaplot+scale_x_date(limits = c(as.Date("2011-03-20"), as.Date("2011-08-01")))
S18_areaplot+scale_x_date(limits = c(as.Date("2011-03-20"), as.Date("2011-08-01")))
S16_areaplot+scale_x_date(limits = c(as.Date("2011-03-20"), as.Date("2011-08-01")))
MPL_areaplot+scale_x_date(limits = c(as.Date("2011-03-20"), as.Date("2011-08-01")))

areaplots_of_each_of_top_10 <- function(otu_table,Jericho_and_SOG_data, Jericho_data, title){
 
 areaplot_from_long_single_otus <- function(long_OTUs, title, colour_area){
  number_of_otus <- length(unique(long_OTUs$OTUid))
  sampling_points <- unique(long_otus$Date)
  OTU_areaplot <- ggplot(long_OTUs, aes(x=Date,y=value, group=OTUid)) + 
   season_line +
   spring_bloom_line+
   HAKA_bloom+
   geom_area(fill=colour_area, position="stack", colour="black")+
   # scale_fill_manual(values=colour_area)+
   annotate("point", 
            x = as.Date(sampling_points),
            y = -0.02,
            colour = "grey",
            shape= 17,
            solid=TRUE,
            size = 5)+
   theme_JAG_presentation()+
   theme(panel.grid.major.y = element_blank(),
         panel.grid.major.x = element_blank())+
   date_scaling +
   ggtitle(title)
  
  return(OTU_areaplot)
 }
 
 proportional_otus <- as.data.frame(prop.table(as.matrix(otu_table),  margin=1))
 long_otus <- add_date_to_OTU_table_rtrn_long(proportional_otus,
                                              Jericho_data) 
 vector_top_10_OTUs <- find_top10_OTUs(long_otus)
 top_10_long_OTUs <- droplevels(subset(long_otus, OTUid %in% vector_top_10_OTUs))
 bar_colour <- colorRampPalette(brewer.pal(n =12 , name = "Set3"))(length(vector_top_10_OTUs))
 list_of_plots <- list()
 ## try with top 10
 #for (otu in seq(1:10)){
 for (otu in seq(length(vector_top_10_OTUs))){
  print(vector_top_10_OTUs[otu])
  current_otu_name <- vector_top_10_OTUs[otu]
  current_otu_long_OTUs <- droplevels(subset(long_otus, OTUid %in% current_otu_name))
  print(bar_colour[otu])
  heatmap_one_of_top_10_OTUs <- areaplot_from_long_single_otus(current_otu_long_OTUs, current_otu_name, bar_colour[otu])
  heatmap_one_of_top_10_OTUs
  list_of_plots <- c(list_of_plots, list(heatmap_one_of_top_10_OTUs))
  
 }
 return(list_of_plots)
}

list_area_MPL <- areaplots_of_each_of_top_10(normalized_MPL_OTUs, Library_metadata, Jericho_data, "MPL")
list_area_gp23 <- areaplots_of_each_of_top_10(normalized_gp23_OTUs, Library_metadata, Jericho_data, "gp23")
list_area_16s <- areaplots_of_each_of_top_10(normalized_16s_OTUs, Library_metadata, Jericho_data, "16s")
list_area_18s <- areaplots_of_each_of_top_10(normalized_18s_OTUs, Library_metadata, Jericho_data, "18s")

pdf("../figures/time_series_top_10_each_OTUs_areaplots_all_amplicons_all_dates%03d.pdf", width = 20, height = 20, onefile = FALSE)
#AVS_areaplot
plot_grid(plotlist = list_area_gp23, ncol=2)
plot_grid(plotlist = list_area_MPL, ncol=2)
plot_grid(plotlist = list_area_16s, ncol=2)
plot_grid(plotlist = list_area_18s, ncol=2)
dev.off()





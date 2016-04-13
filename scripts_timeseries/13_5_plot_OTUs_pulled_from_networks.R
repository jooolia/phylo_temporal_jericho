
## plot certain OTUs pulled from networkss



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


## Need to edit and make sure the font, type are consistent
## To-do:
## make it easier to change what the 18s and 16s are facetted by. 

## read in normalized OTU tables

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv", row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
#normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R2.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

total_number_reads <- sum(normalized_18s_OTUs, normalized_MPL_OTUs, normalized_AVS_OTUs, normalized_gp23_OTUs, normalized_16s_OTUs)

args <- commandArgs(TRUE)
inputFile <- args[1]

## test to see if input file is given, so I can decide whether to use this argument or the orginal one. 
if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
}


## a palette I made from iwanthue
## for 18s there are 5 unique phyla and unknown-grey light =====
phyla_palette_euks <- c("#585671",
                        "#E99449",
                        "#4A9F89",
                        "#023D8D",
                        "#D16477",
                        "#0AC757",
                        "#6CC5CE",
                        "#DE2BA2",
                        "gray23")

## total phyla for both bac and euks18
## get the hex colour for grey 23 for use in other programs:
col2rgb("gray23")
rgb(59, 59, 59, maxColorValue=255)

names(phyla_palette_euks) <- c("Amoebozoa",
                               "Archaeplastida",
                               "Crenarchaeota", 
                               "Cryptophyceae",
                               "Excavata",
                               "Haptophyta",
                               "Opisthokonta",
                               "SAR",
                               "unclassified")


## for 16s there are 10 unique and 1  unique phyla -grey-dark gray66 =====

phyla_palette_bac <- c("#193D32",
                       "#7488C0",
                       "#D34641",
                       "#589857",
                       "#939533",
                       "#8B2761",
                       "#D6491F",
                       "#C57FB8",
                       "#C59B33",
                       "#864B30", 
                       "gray66")

names(phyla_palette_bac) <- c("Actinobacteria",
                              "Bacteroidetes",
                              "Chloroflexi",
                              "Cyanobacteria",
                              "Deferribacteres",
                              "Firmicutes",
                              "Fusobacteria",
                              "Planctomycetes",
                              "Proteobacteria",
                              "Tenericutes",
                              "unclassified")



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
 long_otus$Date <- Jericho_data$Date[match(long_otus$VC_number, Jericho_data$VC_number)]
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

## Plot single OTUs over time
plot_individual_otus_over_time <- function (long_otus) {
 p <- ggplot(long_otus, aes(x=Date, y=value, group=OTUid))+ 
  geom_line(aes(color=OTUid))+
  scale_x_date(breaks = date_breaks("month"),
               labels = date_format("%b"),
               limits = c(as.Date("2010-06-01"),
                          as.Date("2011-07-30")))+
  theme_JAG_presentation()
}

## Plot single OTUs over time _new version
plot_individual_otus_over_time <- function (long_otus, title) {
 #cols <- colorRampPalette(brewer.pal(8,"Dark2"))(nrow(long_otus))
 colourCount  <- length(unique(long_otus$OTUid))
 getPalette  <-  colorRampPalette(brewer.pal(11, "Spectral"))
 p <- ggplot(droplevels(long_otus), aes(x=Date, y=value, group=OTUid))+ 
  geom_area(aes(fill=OTUid))+
  scale_x_date(breaks = date_breaks("month"),
               labels = date_format("%b"),
               limits = c(as.Date("2010-06-01"),
                          as.Date("2011-07-30")))+
  scale_fill_manual(values = getPalette(colourCount))+
  theme_JAG_presentation()+
  ggtitle(title)
}

## want to filter OTU by ones that would be interesting over time based on network analyses

## read in table for analysis
nodes_of_interest <- read.csv("../results/cytoscape_network_group_of_gp23_and_16s_nodes.csv")
## interested in the name parameter here

names_of_interest <- nodes_of_interest$name


long_gp23_otus$OTUid <- gsub("OTU_", "gp23__", long_gp23_otus$OTUid)
long_16s_otus$OTUid <- gsub("OTU_", "S16__", long_16s_otus$OTUid)

## subset by those found in network plot

gp23_of_interest <- subset(long_gp23_otus, OTUid %in% names_of_interest)
S16_of_interest <- subset(long_16s_otus, OTUid %in% names_of_interest)

#ind_gp23_otus_over_time <- plot_individual_otus_over_time(long_gp23_otus)
#ind_gp23_otus_over_time

gp23_and_16s_of_interest <- rbind(gp23_of_interest, S16_of_interest)

test_network <- plot_individual_otus_over_time(gp23_of_interest)
test_network
test_network <- plot_individual_otus_over_time(S16_of_interest)
test_network

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_of_interest)
all_interesting +ylim(0,200)



## what if instead I do the proportions of the communities?


proportional_16s <- as.data.frame(prop.table(as.matrix(normalized_16s_OTUs),  margin=1))
long_prop_16s_otus <- add_date_to_OTU_table_rtrn_long(proportional_16s,
                                                 Jericho_data)
proportional_gp23 <- as.data.frame(prop.table(as.matrix(normalized_gp23_OTUs),  margin=1))
long_prop_gp23_otus <- add_date_to_OTU_table_rtrn_long(proportional_gp23,
                                                      Jericho_data)

long_prop_gp23_otus$OTUid <- gsub("OTU_", "gp23__", long_prop_gp23_otus$OTUid)
long_prop_16s_otus$OTUid <- gsub("OTU_", "S16__", long_prop_16s_otus$OTUid)

## subset by those found in network plot

gp23_prop_of_interest <- subset(long_prop_gp23_otus, OTUid %in% names_of_interest)
S16_prop_of_interest <- subset(long_prop_16s_otus, OTUid %in% names_of_interest)


gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

test_network <- plot_individual_otus_over_time(gp23_prop_of_interest)
test_network
test_network+ylim(0,0.05)
test_network <- plot_individual_otus_over_time(S16_prop_of_interest)
test_network

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_prop_of_interest)
all_interesting +ylim(0,0.2)
all_interesting +ylim(0,0.05)


## try to add in colouring by the type

gp23_prop_of_interest$type <- "gp23"
S16_prop_of_interest$type <- "16s"

gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

plot_individual_otus_over_time_colour_by_type <- function (long_otus) {
 p <- ggplot(long_otus, aes(x=Date, y=value, group=OTUid))+ 
  geom_line(aes(color=type))+
  scale_x_date(breaks = date_breaks("month"),
               labels = date_format("%b"),
               limits = c(as.Date("2010-06-01"),
                          as.Date("2011-07-30")))+
  theme_JAG_presentation()
}


gp23_type <- plot_individual_otus_over_time_colour_by_type(gp23_prop_of_interest)
gp23_type+ylim(0,0.05)


all_interesting <- plot_individual_otus_over_time_colour_by_type(gp23_and_16s_prop_of_interest)
all_interesting
all_interesting +ylim(0,0.2)
all_interesting +ylim(0,0.05)


node_colour = brewer.pal(n = 8, "Dark2")[seq_len(nlevels(as.factor(gp23_and_16s_prop_of_interest$OTUid)))]

node_colour <- colorRampPalette(brewer.pal(n = 12, "Set3"))(nlevels(as.factor(gp23_and_16s_prop_of_interest$OTUid)))


node_colour <- colorRampPalette(brewer.pal(n = 8, "Dark2"))(nlevels(as.factor(gp23_and_16s_prop_of_interest$OTUid)))

names(node_colour) <- levels(as.factor(gp23_and_16s_prop_of_interest$OTUid))


all_interesting <- plot_individual_otus_over_time_colour_by_type(gp23_and_16s_prop_of_interest)
all_interesting
all_interesting +ylim(0,0.2)
all_interesting +ylim(0,0.05)

plot_individual_otus_over_time_colour_by_otu <- function (long_otus, node_colour) {
 p <- ggplot(long_otus, aes(x=Date, y=value, group=OTUid, colour=OTUid))+ 
  geom_line()+
  scale_color_manual(values = node_colour)+
  scale_x_date(breaks = date_breaks("month"),
               labels = date_format("%b"),
               limits = c(as.Date("2010-06-01"),
                          as.Date("2011-07-30")))+
  theme_JAG_presentation()
}

p <- plot_individual_otus_over_time_colour_by_otu(gp23_and_16s_prop_of_interest, node_colour)
p
p +ylim(0,0.2)



### Do another subsection ====

## read in table for analysis
nodes_of_interest <- read.csv("../results/cytoscape_network_group_of_gp23_and_16s_nodes_S16_416_and_first_neighbours.csv")

names_of_interest <- nodes_of_interest$name


long_gp23_otus$OTUid <- gsub("OTU_", "gp23__", long_gp23_otus$OTUid)
long_16s_otus$OTUid <- gsub("OTU_", "S16__", long_16s_otus$OTUid)

## subset by those found in network plot

gp23_of_interest <- subset(long_gp23_otus, OTUid %in% names_of_interest)
S16_of_interest <- subset(long_16s_otus, OTUid %in% names_of_interest)

#ind_gp23_otus_over_time <- plot_individual_otus_over_time(long_gp23_otus)
#ind_gp23_otus_over_time

gp23_and_16s_of_interest <- rbind(gp23_of_interest, S16_of_interest)

test_network <- plot_individual_otus_over_time(gp23_of_interest)
test_network
test_network +ylim(0,50)

test_network <- plot_individual_otus_over_time(S16_of_interest)
test_network

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_of_interest)
all_interesting +ylim(0,200)



## what if instead I do the proportions of the communities?


proportional_16s <- as.data.frame(prop.table(as.matrix(normalized_16s_OTUs),  margin=1))
long_prop_16s_otus <- add_date_to_OTU_table_rtrn_long(proportional_16s,
                                                      Jericho_data)
proportional_gp23 <- as.data.frame(prop.table(as.matrix(normalized_gp23_OTUs),  margin=1))
long_prop_gp23_otus <- add_date_to_OTU_table_rtrn_long(proportional_gp23,
                                                       Jericho_data)

long_prop_gp23_otus$OTUid <- gsub("OTU_", "gp23__", long_prop_gp23_otus$OTUid)
long_prop_16s_otus$OTUid <- gsub("OTU_", "S16__", long_prop_16s_otus$OTUid)

## subset by those found in network plot

gp23_prop_of_interest <- subset(long_prop_gp23_otus, OTUid %in% names_of_interest)
S16_prop_of_interest <- subset(long_prop_16s_otus, OTUid %in% names_of_interest)


gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

test_network_gp23 <- plot_individual_otus_over_time(gp23_prop_of_interest)
test_network
test_network+ylim(0,0.05)
test_network_16s <- plot_individual_otus_over_time(S16_prop_of_interest)
test_network

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_prop_of_interest)
all_interesting +ylim(0,0.2)
all_interesting +ylim(0,0.05)


## try to add in colouring by the type

gp23_prop_of_interest$type <- "gp23"
S16_prop_of_interest$type <- "16s"

gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

plot_individual_otus_over_time_colour_by_type <- function (long_otus) {
 p <- ggplot(long_otus, aes(x=Date, y=value, group=OTUid))+ 
  geom_line(aes(color=type))+
  scale_x_date(breaks = date_breaks("month"),
               labels = date_format("%b"),
               limits = c(as.Date("2010-06-01"),
                          as.Date("2011-07-30")))+
  theme_JAG_presentation()
}


gp23_type <- plot_individual_otus_over_time_colour_by_type(gp23_prop_of_interest)
gp23_type+ylim(0,0.05)


all_interesting <- plot_individual_otus_over_time_colour_by_type(gp23_and_16s_prop_of_interest)
all_interesting
all_interesting +ylim(0,0.2)

pdf("../figures/cyctoscape_network_OTUs_pulled_from_16s_and_gp23_over_time_from_S16_416.pdf", width = 15, height = 11, onefile = FALSE )
grid.arrange(test_network_16s,
             test_network_gp23,
             ncol=1
 )
dev.off()

##### Do it one more time with another set.. #####

## read in table for analysis
nodes_of_interest <- read.csv("../results/cytoscape_network_group_of_gp23_and_16s_nodes_S16_46_and_first_neighbours.csv")

names_of_interest <- nodes_of_interest$name


long_gp23_otus$OTUid <- gsub("OTU_", "gp23__", long_gp23_otus$OTUid)
long_16s_otus$OTUid <- gsub("OTU_", "S16__", long_16s_otus$OTUid)

## subset by those found in network plot

gp23_of_interest <- subset(long_gp23_otus, OTUid %in% names_of_interest)
S16_of_interest <- subset(long_16s_otus, OTUid %in% names_of_interest)

#ind_gp23_otus_over_time <- plot_individual_otus_over_time(long_gp23_otus)
#ind_gp23_otus_over_time

gp23_and_16s_of_interest <- rbind(gp23_of_interest, S16_of_interest)

test_network <- plot_individual_otus_over_time(gp23_of_interest)
test_network
test_network +ylim(0,50)

test_network <- plot_individual_otus_over_time(S16_of_interest)
test_network

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_of_interest)
all_interesting +ylim(0,200)



## what if instead I do the proportions of the communities?


proportional_16s <- as.data.frame(prop.table(as.matrix(normalized_16s_OTUs),  margin=1))
long_prop_16s_otus <- add_date_to_OTU_table_rtrn_long(proportional_16s,
                                                      Jericho_data)
proportional_gp23 <- as.data.frame(prop.table(as.matrix(normalized_gp23_OTUs),  margin=1))
long_prop_gp23_otus <- add_date_to_OTU_table_rtrn_long(proportional_gp23,
                                                       Jericho_data)

long_prop_gp23_otus$OTUid <- gsub("OTU_", "gp23__", long_prop_gp23_otus$OTUid)
long_prop_16s_otus$OTUid <- gsub("OTU_", "S16__", long_prop_16s_otus$OTUid)

## subset by those found in network plot

gp23_prop_of_interest <- subset(long_prop_gp23_otus, OTUid %in% names_of_interest)
S16_prop_of_interest <- subset(long_prop_16s_otus, OTUid %in% names_of_interest)


gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

test_network_gp23 <- plot_individual_otus_over_time(gp23_prop_of_interest)
test_network
test_network+ylim(0,0.05)
test_network_16s <- plot_individual_otus_over_time(S16_prop_of_interest)
test_network

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_prop_of_interest)
all_interesting +ylim(0,0.2)
all_interesting +ylim(0,0.05)


## try to add in colouring by the type

gp23_prop_of_interest$type <- "gp23"
S16_prop_of_interest$type <- "16s"

gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

plot_individual_otus_over_time_colour_by_type <- function (long_otus) {
 p <- ggplot(long_otus, aes(x=Date, y=value, group=OTUid))+ 
  geom_line(aes(color=type))+
  scale_x_date(breaks = date_breaks("month"),
               labels = date_format("%b"),
               limits = c(as.Date("2010-06-01"),
                          as.Date("2011-07-30")))+
  theme_JAG_presentation()
}


gp23_type <- plot_individual_otus_over_time_colour_by_type(gp23_prop_of_interest)
gp23_type+ylim(0,0.05)


all_interesting <- plot_individual_otus_over_time_colour_by_type(gp23_and_16s_prop_of_interest)
all_interesting
all_interesting +ylim(0,0.2)

pdf("../figures/cyctoscape_network_OTUs_pulled_from_16s_and_gp23_over_time_from_S16_46.pdf", width = 15, height = 11, onefile = FALSE )
grid.arrange(test_network_16s,
             test_network_gp23,
             ncol=1
)
dev.off()










## what if I added in the top 2 OTUs of each over time and see what they are doing? ####
## from long_prop_gp23_otus and long_prop_16s_otus

gp23_sorted_by_mean_overall_abundance <- long_prop_gp23_otus %>%
group_by(OTUid) %>%
 dplyr::summarise(m=median(value)) %>%
 dplyr::arrange(desc(m))

## let's take the first two
names_first_two_gp23 <- gp23_sorted_by_mean_overall_abundance$OTUid[1:2]

## add them to the other set
gp23_from_network_and_top_2 <- c(as.character(names_of_interest),names_first_two_gp23)

## do the same with the 16s
S16_sorted_by_mean_overall_abundance <- long_prop_16s_otus %>%
 group_by(OTUid) %>%
 dplyr::summarise(m=median(value)) %>%
 dplyr::arrange(desc(m))

## let's take the first two
names_first_two_16s <- S16_sorted_by_mean_overall_abundance$OTUid[1:2]

## add them to the other set
S16_from_network_and_top_2 <- c(as.character(names_of_interest),names_first_two_16s)


## subset by those found in network plot

gp23_prop_of_interest <- subset(long_prop_gp23_otus, OTUid %in% names_of_interest)
gp23_top_2 <- subset(long_prop_gp23_otus, OTUid %in% names_first_two_gp23)
S16_prop_of_interest <- subset(long_prop_16s_otus, OTUid %in% names_of_interest)
S16_top_2 <- subset(long_prop_16s_otus, OTUid %in% names_first_two_16s)


gp23_and_16s_prop_of_interest_with_top_2 <- rbind(gp23_prop_of_interest_with_top2, S16_prop_of_interest_with_top2)

test_network <- plot_individual_otus_over_time(gp23_prop_of_interest_with_top2)
test_network
test_network+ylim(0,0.05)
test_network <- plot_individual_otus_over_time(S16_prop_of_interest_with_top2)
test_network

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_prop_of_interest_with_top_2)
all_interesting +ylim(0,0.2)
all_interesting +ylim(0,0.05)

## try to add in colouring by the type

gp23_prop_of_interest$type <- "gp23 from network"
gp23_top_2$type <- "gp23 top 2"
S16_prop_of_interest$type <- "16s"
S16_top_2$type <- "16s top 2"

gp23_and_16s_prop_of_interest_with_top2 <- rbind(gp23_prop_of_interest, S16_prop_of_interest, gp23_top_2, S16_top_2)

plot_individual_otus_over_time_colour_by_type <- function (long_otus) {
 p <- ggplot(long_otus, aes(x=Date, y=value, group=OTUid))+ 
  geom_line(aes(color=type), size=1.5)+
  scale_x_date(breaks = date_breaks("month"),
               labels = date_format("%b"),
               limits = c(as.Date("2010-06-01"),
                          as.Date("2011-07-30")))+
  theme_JAG_presentation()
}


all_interesting <- plot_individual_otus_over_time_colour_by_type(gp23_and_16s_prop_of_interest_with_top2)
all_interesting + theme(panel.grid.major.x = element_line(colour="grey20"),
                        panel.grid.major.y=element_blank())


pdf("../figures/cytoscape_network_OTUs_pulled_from_16s_and_gp23_over_time_with_top2.pdf", width = 15, height = 11, onefile = FALSE )
all_interesting + theme(panel.grid.major.x = element_line(colour="grey20"),
                        panel.grid.major.y=element_blank())
dev.off()

## do with another too... ====
## needs to be a node table!!
nodes_of_interest <- read.csv("../results/cytoscape_network_group_of_gp23_and_16s_nodes_S16_26_and_first_neighbours.csv")

names_of_interest <- nodes_of_interest$name


long_gp23_otus$OTUid <- gsub("OTU_", "gp23__", long_gp23_otus$OTUid)
long_16s_otus$OTUid <- gsub("OTU_", "S16__", long_16s_otus$OTUid)

## subset by those found in network plot

gp23_of_interest <- subset(long_gp23_otus, OTUid %in% names_of_interest)
S16_of_interest <- subset(long_16s_otus, OTUid %in% names_of_interest)

#ind_gp23_otus_over_time <- plot_individual_otus_over_time(long_gp23_otus)
#ind_gp23_otus_over_time

gp23_and_16s_of_interest <- rbind(gp23_of_interest, S16_of_interest)

test_network <- plot_individual_otus_over_time(gp23_of_interest)
test_network
test_network +ylim(0,50)

test_network <- plot_individual_otus_over_time(S16_of_interest)
test_network

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_of_interest)
all_interesting +ylim(0,200)



## what if instead I do the proportions of the communities?


proportional_16s <- as.data.frame(prop.table(as.matrix(normalized_16s_OTUs),  margin=1))
long_prop_16s_otus <- add_date_to_OTU_table_rtrn_long(proportional_16s,
                                                      Jericho_data)
proportional_gp23 <- as.data.frame(prop.table(as.matrix(normalized_gp23_OTUs),  margin=1))
long_prop_gp23_otus <- add_date_to_OTU_table_rtrn_long(proportional_gp23,
                                                       Jericho_data)

long_prop_gp23_otus$OTUid <- gsub("OTU_", "gp23__", long_prop_gp23_otus$OTUid)
long_prop_16s_otus$OTUid <- gsub("OTU_", "S16__", long_prop_16s_otus$OTUid)

## subset by those found in network plot

gp23_prop_of_interest <- subset(long_prop_gp23_otus, OTUid %in% names_of_interest)
S16_prop_of_interest <- subset(long_prop_16s_otus, OTUid %in% names_of_interest)


gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

test_network <- plot_individual_otus_over_time(gp23_prop_of_interest)
test_network
test_network+ylim(0,0.05)
test_network <- plot_individual_otus_over_time(S16_prop_of_interest)
test_network

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_prop_of_interest)
all_interesting +ylim(0,0.2)
all_interesting +ylim(0,0.05)


## try to add in colouring by the type

gp23_prop_of_interest$type <- "phage"
S16_prop_of_interest$type <- "bacte"

gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

plot_individual_otus_over_time_colour_by_type <- function (long_otus) {
 p <- ggplot(long_otus, aes(x=Date, y=value, group=OTUid))+ 
  geom_line(aes(color=type))+
  scale_x_date(breaks = date_breaks("month"),
               labels = date_format("%b"),
               limits = c(as.Date("2010-06-01"),
                          as.Date("2011-07-30")))+
  theme_JAG_presentation()
}


gp23_type <- plot_individual_otus_over_time_colour_by_type(gp23_prop_of_interest)
gp23_type
S16_type <- plot_individual_otus_over_time_colour_by_type(S16_prop_of_interest)
S16_type

pdf("../figures/cytoscape_network_OTUs_pulled_from_16s_and_gp23_over_time_from_S16_26.pdf", width = 15, height = 11, onefile = FALSE )

grid.arrange(S16_type + theme(panel.grid.major.x = element_line(colour="grey20"),
                           panel.grid.major.y=element_blank()),
             gp23_type + theme(panel.grid.major.x = element_line(colour="grey20"),
                               panel.grid.major.y=element_blank()),
             ncol=1)
dev.off()



all_interesting <- plot_individual_otus_over_time_colour_by_type(gp23_and_16s_prop_of_interest)
all_interesting
all_interesting +ylim(0,0.2)


## what if I added in the top 2 OTUs of each over time and see what they are doing? ####
## from long_prop_gp23_otus and long_prop_16s_otus

gp23_sorted_by_mean_overall_abundance <- long_prop_gp23_otus %>%
 group_by(OTUid) %>%
 dplyr::summarise(m=median(value)) %>%
 dplyr::arrange(desc(m))

## let's take the first two
names_first_two_gp23 <- gp23_sorted_by_mean_overall_abundance$OTUid[1:2]

## add them to the other set
gp23_from_network_and_top_2 <- c(as.character(names_of_interest),names_first_two_gp23)

## do the same with the 16s
S16_sorted_by_mean_overall_abundance <- long_prop_16s_otus %>%
 group_by(OTUid) %>%
 dplyr::summarise(m=median(value)) %>%
 dplyr::arrange(desc(m))

## let's take the first two
names_first_two_16s <- S16_sorted_by_mean_overall_abundance$OTUid[1:2]

## add them to the other set
S16_from_network_and_top_2 <- c(as.character(names_of_interest),names_first_two_16s)


## subset by those found in network plot

gp23_prop_of_interest <- subset(long_prop_gp23_otus, OTUid %in% names_of_interest)
gp23_top_2 <- subset(long_prop_gp23_otus, OTUid %in% names_first_two_gp23)
S16_prop_of_interest <- subset(long_prop_16s_otus, OTUid %in% names_of_interest)
S16_top_2 <- subset(long_prop_16s_otus, OTUid %in% names_first_two_16s)


gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

test_network <- plot_individual_otus_over_time(gp23_prop_of_interest)
test_network
test_network+ylim(0,0.05)
test_network <- plot_individual_otus_over_time(S16_prop_of_interest)
test_network

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_prop_of_interest)
all_interesting +ylim(0,0.2)
all_interesting +ylim(0,0.05)

## try to add in colouring by the type

gp23_prop_of_interest$type <- "phage from network"
gp23_top_2$type <- "phage top 2"
S16_prop_of_interest$type <- "bacteria from network"
S16_top_2$type <- "bacteria top 2"


gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

gp23_and_16s_prop_of_interest_with_top2 <- rbind(gp23_prop_of_interest, S16_prop_of_interest, gp23_top_2, S16_top_2)
}


all_interesting <- plot_individual_otus_over_time_colour_by_type(gp23_and_16s_prop_of_interest)
all_interesting

all_interesting_with_top2 <- plot_individual_otus_over_time_colour_by_type(gp23_and_16s_prop_of_interest_with_top2)
all_interesting_with_top2 + theme(panel.grid.major.x = element_line(colour="grey20"),
                        panel.grid.major.y=element_blank())


pdf("../figures/cytoscape_network_OTUs_pulled_from_16s_and_gp23_over_time_with_cyano.pdf", width = 15, height = 11, onefile = FALSE )
all_interesting + theme(panel.grid.major.x = element_line(colour="grey20"),
                        panel.grid.major.y=element_blank()) +ylim(0,0.2)
dev.off()


pdf("../figures/OTUs_pulled_from_16s_and_gp23_over_time_with_cyano_with_top2.pdf", width = 15, height = 11, onefile = FALSE )
all_interesting_with_top2 + theme(panel.grid.major.x = element_line(colour="grey20"),
                        panel.grid.major.y=element_blank())
dev.off()


## plot all the correlations over time as example for presentation =====

all_strong_correlations <- read.csv("../results/significant_strong_correlationmerged_all.csv", row.names=1)


## read in table for analysis


names_of_interest <- union(all_strong_correlations$Var1, all_strong_correlations$Var2)


long_gp23_otus$OTUid <- gsub("OTU_", "gp23__", long_gp23_otus$OTUid)
long_16s_otus$OTUid <- gsub("OTU_", "S16__", long_16s_otus$OTUid)

## subset by those found in network plot

gp23_of_interest <- subset(long_gp23_otus, OTUid %in% names_of_interest)
S16_of_interest <- subset(long_16s_otus, OTUid %in% names_of_interest)

#ind_gp23_otus_over_time <- plot_individual_otus_over_time(long_gp23_otus)
#ind_gp23_otus_over_time
gp23_of_interest$type <- "phage with strong correlations"
S16_of_interest$type <- "bacteria with strong correlations"

gp23_and_16s_of_interest <- rbind(gp23_of_interest, S16_of_interest)

test_network_gp23 <- plot_individual_otus_over_time(gp23_of_interest)
test_network_gp23
test_network +ylim(0,500)

test_network_16 <- plot_individual_otus_over_time(S16_of_interest)
test_network_16

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_of_interest)
all_interesting +ylim(0,500)



all_interesting <- plot_individual_otus_over_time_colour_by_type(gp23_and_16s_of_interest)
all_interesting
all_interesting +ylim(0,500)

pdf("../figures/All_strong_correlations_16s_and_gp23_over_time.pdf", width = 15, height = 11, onefile = FALSE )

all_interesting + theme(panel.grid.major.x = element_line(colour="grey20"),
                                  panel.grid.major.y=element_blank())+ylim(0,500)
dev.off()

## DO by proportion


gp23_prop_of_interest <- subset(long_prop_gp23_otus, OTUid %in% names_of_interest)

S16_prop_of_interest <- subset(long_prop_16s_otus, OTUid %in% names_of_interest)

## try to add in colouring by the type

gp23_prop_of_interest$type <- "phage with strong correlations"

S16_prop_of_interest$type <- "bacteria with strong correlations"



gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

test_network_gp23 <- plot_individual_otus_over_time(gp23_prop_of_interest)
test_network_gp23
test_network+ylim(0,0.05)
test_network_16s <- plot_individual_otus_over_time(S16_prop_of_interest)
test_network_16s

all_interesting <- plot_individual_otus_over_time(gp23_and_16s_prop_of_interest)
all_interesting +ylim(0,0.2)
all_interesting +ylim(0,0.05)

all_interesting <- plot_individual_otus_over_time_colour_by_type(gp23_and_16s_prop_of_interest)
all_interesting
all_interesting +ylim(0, 0.3)


pdf("../figures/All_strong_correlations_proportional_16s_and_gp23_over_time.pdf", width = 15, height = 11, onefile = FALSE )

grid.arrange(test_network_gp23+ theme(panel.grid.major.x = element_line(colour="grey20"),
                                      panel.grid.major.y=element_blank()), test_network_16s+ theme(panel.grid.major.x = element_line(colour="grey20"),
                                                                                                   panel.grid.major.y=element_blank()), ncol=1)

dev.off()

## for presenation pull out 2 oTUs that are correlated. 

names <- c("S16__44", "gp23__29")
gp23_prop_of_interest <- subset(long_prop_gp23_otus, OTUid %in% names)

S16_prop_of_interest <- subset(long_prop_16s_otus, OTUid %in% names)


gp23_prop_of_interest$type <- "phage with strong correlations"

S16_prop_of_interest$type <- "bacte with strong correlations"

gp23_and_16s_prop_of_interest <- rbind(gp23_prop_of_interest, S16_prop_of_interest)

gp23_network_type <- plot_individual_otus_over_time_colour_by_type(gp23_prop_of_interest)
S16_network_type <- plot_individual_otus_over_time_colour_by_type(S16_prop_of_interest)

all_interesting <- plot_individual_otus_over_time_colour_by_type(gp23_prop_of_interest)
all_interesting+ geom_line(aes(color=type, alpha=1))

pdf("../figures/Strong_correlations_proportional_16s_and_gp23_over_time_2_OTUs.pdf", width = 15, height = 11, onefile = FALSE )
grid.arrange(S16_network_type +
              geom_line(aes(color=type, alpha=1))+ theme(panel.grid.major.x = element_line(colour="grey20"),                                                        panel.grid.major.y=element_blank()),
             gp23_network_type + theme(panel.grid.major.x = element_line(colour="grey20"),panel.grid.major.y=element_blank()),
             ncol=1)

dev.off()


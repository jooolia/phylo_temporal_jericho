library(ggplot2)
library(plyr)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(scales)
library(vegan)


### FIX HETEROTROPHs

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)
taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)


##########################
## Try to separate by het and phototrophic 18s
## like fig 2 in de Vargas et al 2015

phytoplankton_fam <- c("Diatomea",
                       "Dinophyceae",
                       "Mamiellophyceae",
                       "Chrysophyceae",
                       "Chlorophyceae",
                       "Euglenozoa",
                       "Dictyochophyceae",
                       "Chlorarachniophyta",
                       "Raphidophyceae")
phytoplankton_class <- c("Rhodophyceae",
                         "Chloroplastida")
phytoplantkon_phylum <- c("Haptophyta",
                          "Cryptophyceae")
phytoplankton_order <- c("Dinoflagellata",
                         "Ochrophyta") ## not sure about this, are all dinos phytoplankton?

## algae but maybe not phytoplankton,
algae_class <- c("Chloroplastida")

phytoplankton_tax <- taxonomy_18s %>% 
 filter(Family %in% phytoplankton_fam | Class %in% phytoplankton_class | Phylum %in% phytoplantkon_phylum| Order %in% phytoplankton_order)

colnames(normalized_18s_OTUs) <- gsub(".size.*.", "", colnames(normalized_18s_OTUs))

phytoplankton_otu_table <- droplevels(subset(normalized_18s_OTUs, select=colnames(normalized_18s_OTUs) %in% phytoplankton_tax$otu_number)) 

write.table(phytoplankton_otu_table,"../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv",
            sep="\t",
            col.names=NA)
## masts we don't know...

## phaedoraea, nasselaria, spumellaria, collodaria
protozooplankton_order <- c("Ciliophora")
protozooplankton_family <- c("Acanthoecida",
                             "Foraminifera")

## or coudl try approach that heterotrophs are those that are not phytoplankton
## but this is not quite right....

heterotroph_tax <- taxonomy_18s %>% 
 filter(!(otu_number %in% phytoplankton_tax$otu_number))

heterotroph_otu_table <- droplevels(subset(normalized_18s_OTUs, select=colnames(normalized_18s_OTUs) %in% heterotroph_tax$otu_number)) 

write.table(heterotroph_otu_table,"../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv",
            sep="\t",
            col.names=NA)

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


load(file="../../phyla_palette_euks.txt")
load(file="../../class_palette_euks.txt")
load(file="../../order_palette_euks.txt")
load(file="../../family_palette_euks.txt")
load(file="../../phyla_palette_bac.txt")
load(file="../../class_palette_bac.txt")
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

mid_summer_x <- as.Date("2010-09-22")-as.numeric(as.Date("2010-09-22")-as.Date("2010-06-22"))/2
mid_fall_x <- as.Date("2010-12-22")-as.numeric(as.Date("2010-12-22")-as.Date("2010-09-22"))/2
mid_winter_x <- as.Date("2011-03-22")-as.numeric(as.Date("2011-03-22")-as.Date("2010-12-22"))/2
mid_spring_x <- as.Date("2011-06-22")-as.numeric(as.Date("2011-06-22")-as.Date("2011-03-22"))/2

season_text <- annotate("text",x=c(mid_summer_x, mid_fall_x, mid_winter_x, mid_spring_x), y=1.05, label=c("Summer", "Fall", "Winter", "Spring"))

spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=1)
date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")))

### Get relevant VC numbers ====
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
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

long_18s_otus <- add_date_to_OTU_table_rtrn_long(normalized_18s_OTUs,
                                                 Jericho_data)
names(long_18s_otus)[2] <- "otu_number"
## merge together with taxonomy
long_18s_with_taxonomy <- merge(long_18s_otus, taxonomy_18s, by="otu_number")
## make sure only euks
long_18s_with_taxonomy <- droplevels(subset(long_18s_with_taxonomy, Domain == "Eukaryota"))


long_16s_otus <- add_date_to_OTU_table_rtrn_long(normalized_16s_OTUs,
                                                 Jericho_data)

names(long_16s_otus)[2] <- "otu_number"
long_16s_with_taxonomy <- merge(long_16s_otus, taxonomy_16s, by="otu_number")
## want to remove the chloroplasts
long_16s_with_taxonomy <- long_16s_with_taxonomy %>%
 filter(Class != "Chloroplast")
long_16s_with_taxonomy <- droplevels(subset(long_16s_with_taxonomy, Domain == "Bacteria"))


## need to summarise the data by the counts of each order...
## Plotting functions ####

# want relative abundance by data...
summarise_otus_by_taxon <- function (long_otus_with_taxonomy, taxonomic_level) {
 long_otus_summarise <- ddply(long_otus_with_taxonomy,
                                       c(taxonomic_level, "Date"), 
                                       summarise, 
                                       rel_abun=sum(value))
 return(long_otus_summarise)
}


# long_18s_otus_phyla_over_time <- area_all_taxa_over_time(summarise_otus_by_taxon(long_18s_with_taxonomy, "Phylum"), 
#                                                       phyla_palette_euks, 
#                                                       "Phylum", 
#                                                       "18S by phyla")
# long_18s_otus_phyla_over_time
# 
# long_18s_otus_orders_over_time <- area_all_taxa_over_time(summarise_otus_by_taxon(long_18s_with_taxonomy, 
#                                                                                "Order"),
#                                                        order_palette_euks,
#                                                        "Order",
#                                                        "18S by Order")
# long_18s_otus_orders_over_time

### what if we get rid of the unclassified -have to do this by taxonomic level each time
## Do the same for phyla
summarise_otus_by_phyla_remove_unclassified <- function (long_otus_with_taxonomy) {
 long_otus_phyla_no_unclassified <- filter(long_otus_with_taxonomy,Phylum!="unclassified")
 long_otus_by_phyla_no_unclassified_summarise <- ddply(long_otus_phyla_no_unclassified,
                                                       c("Phylum", "Date"),
                                                       summarise,
                                                       rel_abun=sum(value))
 return(long_otus_by_phyla_no_unclassified_summarise)
}
summarise_otus_by_class_remove_unclassified <- function (long_otus_with_taxonomy) {
 long_otus_phyla_no_unclassified <- filter(long_otus_with_taxonomy,Class!="unclassified")
 long_otus_by_phyla_no_unclassified_summarise <- ddply(long_otus_phyla_no_unclassified,
                                                       c("Class", "Date"),
                                                       summarise,
                                                       rel_abun=sum(value))
 return(long_otus_by_phyla_no_unclassified_summarise)
}

summarise_otus_by_order_remove_unclassified <- function (long_otus_with_taxonomy) {
 long_otus_orders_no_unclassified <- filter(long_otus_with_taxonomy,Order!="unclassified")
 long_otus_by_order_no_unclassified_summarise <- ddply(long_otus_orders_no_unclassified ,
                                                       c("Order", "Date"),
                                                       summarise, 
                                                       rel_abun=sum(value))
 return(long_otus_by_order_no_unclassified_summarise)
}

summarise_otus_by_family_remove_unclassified <- function (long_otus_with_taxonomy) {
 long_otus_orders_no_unclassified <- filter(long_otus_with_taxonomy,Family!="unclassified")
 long_otus_by_order_no_unclassified_summarise <- ddply(long_otus_orders_no_unclassified ,
                                                       c("Family", "Date"),
                                                       summarise, 
                                                       rel_abun=sum(value))
 return(long_otus_by_order_no_unclassified_summarise)
}


plot_no_unclassified_over_time_facet_wrap <- function (long_otus_orders_no_unclassified, Order) {
 ggplot(long_otus_orders_no_unclassified, 
        aes(x=Date, y=rel_abun))+
  season_line +
  season_text+
  spring_bloom_line+
  geom_line()+ 
  facet_wrap(~ Order)+
  date_scaling+
  theme_JAG_presentation()+
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
}

# long_18s_otus_orders_no_unclassified_over_time_facet <- plot_no_unclassified_over_time_facet_wrap(summarise_otus_by_order_remove_unclassified(long_otus_orders_no_unclassified), 
#                                                                                                   Order)
# long_18s_otus_orders_no_unclassified_over_time_facet

find_max_value_of_order <- function (long_otus_orders_no_unclassified) {
 long_otus_orders_no_unclassified_max_value  <- ddply(long_otus_orders_no_unclassified , 
                                                      "Order",
                                                      summarise,
                                                      Order_max_value=max(rel_abun))
 return(long_otus_orders_no_unclassified_max_value)
}

#long_18s_otus_orders_max_value <- find_max_value_of_order(long_18s_otus_orders_no_unclassified)

get_top_5_orders <- function (long_otus_orders_max_values) {
 ordered_smallest_to_largest <- long_otus_orders_max_values %>%
  arrange(Order_max_value)
 Five_largest_rows <- tail(ordered_smallest_to_largest, n=5L)
 Top_5_Orders <- droplevels(Five_largest_rows$Order)
 return(Top_5_Orders)
}


high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)


##### 18s prop #####
## barplots of each over time out of 100

## for missing

proportional_18s <- as.data.frame(prop.table(as.matrix(normalized_18s_OTUs),  margin=1))

long_18s_otus_prop <- add_date_to_OTU_table_rtrn_long(proportional_18s,Jericho_data)

names(long_18s_otus_prop)[2] <- "otu_number"
long_18s_otus_prop$otu_number <- gsub(".size.*.", "", long_18s_otus_prop$otu_number)
long_18s_prop_with_taxonomy <- merge(long_18s_otus_prop, taxonomy_18s, by="otu_number")

long_18s_prop_with_taxonomy  <- droplevels(subset(long_18s_prop_with_taxonomy , Domain == "Eukaryota"))

prop_18s_summarised_by_phyla <- summarise_otus_by_taxon(long_18s_prop_with_taxonomy, "Phylum")
prop_18s_summarised_by_class <- summarise_otus_by_taxon(long_18s_prop_with_taxonomy, "Class")
prop_18s_summarised_by_order <- summarise_otus_by_taxon(long_18s_prop_with_taxonomy, "Order")
prop_18s_summarised_by_family <- summarise_otus_by_taxon(long_18s_prop_with_taxonomy, "Family")

barplot_all_taxa_over_time <- function (long_otus_by_order,
                                        palette, 
                                        taxonomic_level,
                                        title="") {
 p <- ggplot(long_otus_by_order, 
             aes_string(x="Date", y="rel_abun",
                        group=taxonomic_level))+ 
  season_line +
  season_text +
  spring_bloom_line+
  geom_bar(stat="identity",
             aes_string(fill=taxonomic_level),
             colour=line_colour)+
  date_scaling+
  scale_fill_manual(values=palette)+
 # ggtitle(title)+
  theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 16))
 return(p)
}

area_all_taxa_over_time <- function (long_otus_by_order, 
                                     palette,
                                     taxonomic_level,
                                     title="") {
 p <- ggplot(long_otus_by_order,
             aes_string(x="Date", y="rel_abun",
                        group=taxonomic_level))+
  season_line +
  season_text +
  spring_bloom_line+
  geom_area(position="stack", 
              aes_string(fill=taxonomic_level))+
  date_scaling+
  scale_fill_manual(values=palette)+
  #ggtitle(title)+
  theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())
}

#### Proportional plots 18s #####
## both area and barplots

proportional_barplot_18s_phyla <- barplot_all_taxa_over_time(prop_18s_summarised_by_phyla,
                                                             phyla_palette_euks, 
                                                             "Phylum",
                                                             "18S by phyla")


proportional_barplot_18s_phyla


proportional_area_18s_phyla <- area_all_taxa_over_time(prop_18s_summarised_by_phyla,
                                                             phyla_palette_euks, 
                                                        "Phylum",
                                                       "18S by phyla")

proportional_barplot_18s_class <- barplot_all_taxa_over_time(prop_18s_summarised_by_class,
                                                             class_palette_euks, 
                                                             "Class",
                                                             "18S by class")
proportional_area_18s_class <-area_all_taxa_over_time(prop_18s_summarised_by_class,
                                                      class_palette_euks,
                                                      "Class",
                                                      "18S by class")
proportional_barplot_18s_order <- barplot_all_taxa_over_time(prop_18s_summarised_by_order,
                                                             order_palette_euks, 
                                                             "Order",
                                                             "18S by order")
proportional_area_18s_order <-area_all_taxa_over_time(prop_18s_summarised_by_order,
                                                      order_palette_euks,
                                                      "Order",
                                                      "18S by order")
proportional_barplot_18s_family <- barplot_all_taxa_over_time(prop_18s_summarised_by_family,
                                                             family_palette_euks, 
                                                             "Family",
                                                             "18S by family")
proportional_area_18s_family <-area_all_taxa_over_time(prop_18s_summarised_by_family,
                                                      family_palette_euks,
                                                      "Family",
                                                      "18S by family")

### remove unclassified and then plot ####
# Do by Phylum

long_18s_prop_with_taxonomy_no_opisthokont <- subset(long_18s_prop_with_taxonomy, Phylum!="Opisthokonta")
prop_18s_summarised_by_phyla_no_unclassified <- summarise_otus_by_phyla_remove_unclassified(long_18s_prop_with_taxonomy_no_opisthokont)
#subset(long_18s_prop_with_taxonomy, Phylum!="Opisthokonta")
proportional_barplot_18s_phyla_no_unclassified <- barplot_all_taxa_over_time(prop_18s_summarised_by_phyla_no_unclassified,
                                                                             phyla_palette_euks,
                                                                             "Phylum",
                                                                             "18S by phyla no unclassified")
## make some without the high res VCx
long_18s_prop_with_taxonomy_no_high_res <- subset(long_18s_prop_with_taxonomy, !(VC_number %in% high_res_vcs))
prop_18s_summarised_by_phyla_no_high_res <- summarise_otus_by_phyla_remove_unclassified(long_18s_prop_with_taxonomy_no_high_res)
proportional_barplot_18s_phyla_no_high_res <- barplot_all_taxa_over_time(prop_18s_summarised_by_phyla_no_high_res ,
                                                                             phyla_palette_euks,
                                                                             "Phylum",
                                                                             "18S by phyla no unclassified no high res")

## only high res
long_18s_prop_with_taxonomy_only_high_res <- subset(long_18s_prop_with_taxonomy, VC_number %in% high_res_vcs)
prop_18s_summarised_by_phyla_only_high_res <- summarise_otus_by_phyla_remove_unclassified(long_18s_prop_with_taxonomy_only_high_res)
proportional_barplot_18s_phyla_only_high_res <- barplot_all_taxa_over_time(prop_18s_summarised_by_phyla_only_high_res ,
                                                                         phyla_palette_euks,
                                                                         "Phylum",
                                                                         "18S by phyla no unclassified only high res")+
 scale_x_date(breaks = date_breaks("week"), #labels = date_format("%b"),
              limits = c(as.Date("2011-01-15"),
                         as.Date("2011-02-15")))


proportional_area_18s_phyla_no_unclassified <-area_all_taxa_over_time(prop_18s_summarised_by_phyla_no_unclassified, 
                                                                      phyla_palette_euks, 
                                                                      "Phylum",
                                                                      "18S by phyla no unclassified")

prop_18s_summarised_by_class_no_unclassified <- summarise_otus_by_class_remove_unclassified(long_18s_prop_with_taxonomy_no_opisthokont)
proportional_barplot_18s_class_no_unclassified <- barplot_all_taxa_over_time(prop_18s_summarised_by_class_no_unclassified,
                                                                             class_palette_euks,
                                                                             "Class",
                                                                             "18S by class no unclassified")
proportional_area_18s_class_no_unclassified <-area_all_taxa_over_time(prop_18s_summarised_by_class_no_unclassified, 
                                                                      class_palette_euks, 
                                                                      "Class",
                                                                      "18S by class no unclassified")
 
prop_18s_summarised_by_class_no_high_res <- summarise_otus_by_class_remove_unclassified(long_18s_prop_with_taxonomy_no_high_res)
proportional_barplot_18s_class_no_high_res <- barplot_all_taxa_over_time(prop_18s_summarised_by_class_no_high_res ,
                                                                         class_palette_euks,
                                                                         "Class",
                                                                         "18S by class no unclassified no high res")

prop_18s_summarised_by_order_no_unclassified <- summarise_otus_by_order_remove_unclassified(long_18s_prop_with_taxonomy_no_opisthokont)
proportional_barplot_18s_order_no_unclassified <- barplot_all_taxa_over_time(prop_18s_summarised_by_order_no_unclassified,
                                                                       order_palette_euks,
                                                                       "Order",
                                                                       "18S by order no unclassified")

prop_18s_summarised_by_order_no_high_res <- summarise_otus_by_order_remove_unclassified(long_18s_prop_with_taxonomy_no_high_res)
proportional_barplot_18s_order_no_high_res <- barplot_all_taxa_over_time(prop_18s_summarised_by_order_no_high_res ,
                                                                         order_palette_euks,
                                                                         "Order",
                                                                         "18S by order no unclassified no high res")
proportional_area_18s_order_no_unclassified <-area_all_taxa_over_time(na.omit(prop_18s_summarised_by_order_no_unclassified),
                                                                      order_palette_euks,
                                                                      "Order",
                                                                      "18S by order no unclassified")

prop_18s_summarised_by_family_no_unclassified <- summarise_otus_by_family_remove_unclassified(long_18s_prop_with_taxonomy_no_opisthokont)
proportional_barplot_18s_family_no_unclassified <- barplot_all_taxa_over_time(prop_18s_summarised_by_family_no_unclassified, family_palette_euks, "Family", "18S by family no unclassified")

prop_18s_summarised_by_family_no_high_res <- summarise_otus_by_family_remove_unclassified(long_18s_prop_with_taxonomy_no_high_res)

proportional_barplot_18s_family_no_high_res <- barplot_all_taxa_over_time(prop_18s_summarised_by_family_no_high_res,
                                                                         family_palette_euks,
                                                                         "Family",
                                                                         "18S by family no unclassified no high res")

proportional_area_18s_family_no_unclassified <-area_all_taxa_over_time(na.omit(prop_18s_summarised_by_family_no_unclassified), family_palette_euks, "Family", "18S by family no unclassified")


### Print proportional barplots 18s ####
pdf(paste0(figures_dir,"18S_proportional_barplots_time_series%03d.pdf"), width = 15, height = 11, onefile = FALSE)
proportional_barplot_18s_phyla
proportional_barplot_18s_phyla_no_unclassified 

proportional_barplot_18s_phyla_no_high_res
proportional_barplot_18s_phyla_only_high_res

proportional_barplot_18s_class
proportional_barplot_18s_class_no_unclassified

proportional_barplot_18s_order
proportional_barplot_18s_order_no_high_res

proportional_barplot_18s_family
proportional_barplot_18s_family_no_high_res 

dev.off()

### Print proportional areas 18s ####
pdf(paste0(figures_dir,"18S_proportional_area_time_series%03d.pdf"), width = 15, height = 11, onefile = FALSE)
proportional_area_18s_phyla
proportional_area_18s_phyla_no_unclassified
proportional_area_18s_class
proportional_area_18s_class_no_unclassified
proportional_area_18s_order
proportional_area_18s_order_no_unclassified
#proportional_area_18s_order_no_unclassified_top_5
proportional_area_18s_family
proportional_area_18s_family_no_unclassified
dev.off()

### 16s data

#### proportional data. ####
##### 16s #####

proportional_16s <- as.data.frame(prop.table(as.matrix(normalized_16s_OTUs),  margin=1))

long_16s_otus_prop <- add_date_to_OTU_table_rtrn_long(proportional_16s,Jericho_data)
#long_16s_otus_prop <- droplevels(subset(long_16s_otus_prop, value > 0))

names(long_16s_otus_prop)[2] <- "otu_number"
long_16s_otus_prop$otu_number <- gsub(".size.*.", "", long_16s_otus_prop$otu_number)
long_16s_prop_with_taxonomy <- merge(long_16s_otus_prop, taxonomy_16s, by="otu_number")


prop_16s_summarised_by_phyla <- droplevels(summarise_otus_by_taxon(long_16s_prop_with_taxonomy, "Phylum"))
prop_16s_summarised_by_class <- droplevels(summarise_otus_by_taxon(long_16s_prop_with_taxonomy, "Class"))
prop_16s_summarised_by_order <- droplevels(summarise_otus_by_taxon(long_16s_prop_with_taxonomy, "Order"))
prop_16s_summarised_by_family <- droplevels(summarise_otus_by_taxon(long_16s_prop_with_taxonomy, "Family"))

proportional_barplot_16s_phyla <- barplot_all_taxa_over_time(prop_16s_summarised_by_phyla,phyla_palette_bac, "Phylum","16S by phyla")
proportional_barplot_16s_phyla
proportional_barplot_16s_class <- barplot_all_taxa_over_time(prop_16s_summarised_by_class,class_palette_bac, "Class", "16S by class")
proportional_barplot_16s_class
proportional_barplot_16s_order <- barplot_all_taxa_over_time(prop_16s_summarised_by_order,order_palette_bac, "Order", "16S by order")
proportional_barplot_16s_order
proportional_barplot_16s_family <- barplot_all_taxa_over_time(prop_16s_summarised_by_family,family_palette_bac, "Family", "16S by family")
proportional_barplot_16s_family
proportional_area_16s_phyla <- area_all_taxa_over_time(prop_16s_summarised_by_phyla, phyla_palette_bac, "Phylum", "16S by phyla")
proportional_area_16s_phyla
proportional_area_16s_class <- area_all_taxa_over_time(prop_16s_summarised_by_class, class_palette_bac, "Class","16S by class")
proportional_area_16s_class
proportional_area_16s_order <- area_all_taxa_over_time(prop_16s_summarised_by_order, order_palette_bac, "Order", "16S by order")
proportional_area_16s_order
proportional_area_16s_family <- area_all_taxa_over_time(prop_16s_summarised_by_family, family_palette_bac, "Family","16S by family")
proportional_area_16s_family

## remove unclassified
prop_16s_summarised_by_phyla_no_unclassified <- summarise_otus_by_phyla_remove_unclassified(long_16s_prop_with_taxonomy)
proportional_barplot_16s_phyla_no_unclassified <- barplot_all_taxa_over_time(prop_16s_summarised_by_phyla_no_unclassified, phyla_palette_bac, "Phylum",  "16S by phyla no unclassified")
proportional_barplot_16s_phyla_no_unclassified 


## make some without the high res VCx
long_16s_prop_with_taxonomy_no_high_res <- subset(long_16s_prop_with_taxonomy, !(VC_number %in% high_res_vcs))
prop_16s_summarised_by_phyla_no_high_res <- summarise_otus_by_phyla_remove_unclassified(long_16s_prop_with_taxonomy_no_high_res)
proportional_barplot_16s_phyla_no_high_res <- barplot_all_taxa_over_time(prop_16s_summarised_by_phyla_no_high_res ,
                                                                         phyla_palette_bac,
                                                                         "Phylum",
                                                                         "16S by phyla no unclassified no high res")

## only high res
long_16s_prop_with_taxonomy_only_high_res <- subset(long_16s_prop_with_taxonomy, VC_number %in% high_res_vcs)
prop_16s_summarised_by_phyla_only_high_res <- summarise_otus_by_phyla_remove_unclassified(long_16s_prop_with_taxonomy_only_high_res)
proportional_barplot_16s_phyla_only_high_res <- barplot_all_taxa_over_time(prop_16s_summarised_by_phyla_only_high_res ,
                                                                           phyla_palette_bac,
                                                                           "Phylum",
                                                                           "16S by phyla no unclassified only high res")+
 scale_x_date(breaks = date_breaks("week"), #labels = date_format("%b"),
              limits = c(as.Date("2011-01-15"),
                         as.Date("2011-02-15")))

proportional_area_16s_phyla_no_unclassified <- area_all_taxa_over_time(prop_16s_summarised_by_phyla_no_unclassified, phyla_palette_bac, "Phylum","16S by phyla no unclassified")
proportional_area_16s_phyla_no_unclassified


prop_16s_summarised_by_class_no_unclassified <- summarise_otus_by_class_remove_unclassified(long_16s_prop_with_taxonomy)
proportional_barplot_16s_class_no_unclassified <- barplot_all_taxa_over_time(prop_16s_summarised_by_class_no_unclassified, class_palette_bac, "Class","16S by class no unclassified")
proportional_barplot_16s_class_no_unclassified 

prop_16s_summarised_by_class_no_high_res <- summarise_otus_by_class_remove_unclassified(long_16s_prop_with_taxonomy_no_high_res)
proportional_barplot_16s_class_no_high_res <- barplot_all_taxa_over_time(prop_16s_summarised_by_class_no_high_res ,
                                                                         class_palette_bac,
                                                                         "Class",
                                                                         "16S by class no unclassified no high res")


proportional_area_16s_class_no_unclassified <- area_all_taxa_over_time(prop_16s_summarised_by_class_no_unclassified, class_palette_bac, "Class","16S by class no unclassified")
proportional_area_16s_class_no_unclassified

prop_16s_summarised_by_order_no_unclassified <- summarise_otus_by_order_remove_unclassified(long_16s_prop_with_taxonomy)
proportional_barplot_16s_order_no_unclassified <- barplot_all_taxa_over_time(prop_16s_summarised_by_order_no_unclassified, order_palette_bac, "Order", "16S by order no unclassified")
proportional_barplot_16s_order_no_unclassified 
proportional_area_16s_order_no_unclassified <- area_all_taxa_over_time(prop_16s_summarised_by_order_no_unclassified, order_palette_bac, "Order", "16S by order no unclassified")
proportional_area_16s_order_no_unclassified


prop_16s_summarised_by_order_no_high_res <- summarise_otus_by_order_remove_unclassified(long_16s_prop_with_taxonomy_no_high_res)
proportional_barplot_16s_order_no_high_res <- barplot_all_taxa_over_time(prop_16s_summarised_by_order_no_high_res ,
                                                                         order_palette_bac,
                                                                         "Order",
                                                                         "16S by order no unclassified no high res")



prop_16s_summarised_by_family_no_unclassified <- summarise_otus_by_family_remove_unclassified(long_16s_prop_with_taxonomy)
proportional_barplot_16s_family_no_unclassified <- barplot_all_taxa_over_time(prop_16s_summarised_by_family_no_unclassified, family_palette_bac, "Family", "16S by family no unclassified")
proportional_barplot_16s_family_no_unclassified 

proportional_area_16s_family_no_unclassified <- area_all_taxa_over_time(prop_16s_summarised_by_family_no_unclassified, family_palette_bac, "Family", "16S by family no unclassified")
proportional_area_16s_family_no_unclassified


prop_16s_summarised_by_family_no_high_res <- summarise_otus_by_family_remove_unclassified(long_16s_prop_with_taxonomy_no_high_res)
proportional_barplot_16s_family_no_high_res <- barplot_all_taxa_over_time(prop_16s_summarised_by_family_no_high_res ,
                                                                         family_palette_bac,
                                                                         "Family",
                                                                         "16S by family no unclassified no high res")



### Print proportional barplots 16s ####
pdf(paste0(figures_dir,"16S_proportional_barplots_time_series%03d.pdf"), width = 15, height = 11, onefile = FALSE)
proportional_barplot_16s_phyla
proportional_barplot_16s_phyla_no_unclassified
proportional_barplot_16s_phyla_no_high_res # using in paper
proportional_barplot_16s_phyla_only_high_res
proportional_barplot_16s_class
proportional_barplot_16s_class_no_high_res
proportional_barplot_16s_order
proportional_barplot_16s_order_no_high_res 
proportional_barplot_16s_family
proportional_barplot_16s_family_no_high_res
dev.off()

### Print proportional areas 16s  ####
pdf(paste0(figures_dir,"16S_proportional_area_time_series%03d.pdf"), width = 15, height = 11, onefile = FALSE)
proportional_area_16s_phyla
proportional_area_16s_phyla_no_unclassified
proportional_area_16s_class
proportional_area_16s_class_no_unclassified
proportional_area_16s_order
proportional_area_16s_order_no_unclassified
proportional_area_16s_family
proportional_area_16s_family_no_unclassified
dev.off()

for (ind_phylum in unique(long_16s_prop_with_taxonomy$Phylum)) {
 #print(ind_phylum)
 phyla_subset <- droplevels(subset(long_16s_prop_with_taxonomy, Phylum == ind_phylum))
pdf(paste0(figures_dir,"16S_time_series_phylum_",ind_phylum,"_by_specific_orders_barplots%03d.pdf"), width = 15, height = 11, onefile = FALSE)
  for (ind_class in unique(phyla_subset$Class)){
   #print(ind_class)
   class_subset <- droplevels(subset(phyla_subset, Class == ind_class))
   prop_sum_order <- summarise_otus_by_taxon(class_subset, "Order")
   #print(prop_sum_order)
   p <- barplot_all_taxa_over_time(prop_sum_order,
                              order_palette_bac, 
                              "Order", 
                              paste("16S by", ind_phylum, "phylum",ind_class, "class separated by order")
   )
   print(p)
  } 
dev.off()
}

## do by area
for (ind_phylum in unique(long_16s_prop_with_taxonomy$Phylum)) {
 #print(ind_phylum)
 phyla_subset <- droplevels(subset(long_16s_prop_with_taxonomy, Phylum == ind_phylum))
 pdf(paste0(figures_dir,"16S_time_series_phylum_",ind_phylum,"_by_specific_orders_areaplots%03d.pdf"), width = 15, height = 11, onefile=FALSE)
 for (ind_class in unique(phyla_subset$Class)){
  #print(ind_class)
  class_subset <- droplevels(subset(phyla_subset, Class == ind_class))
  prop_sum_order <- summarise_otus_by_taxon(class_subset, "Order")
  #print(prop_sum_order)
  p <- area_all_taxa_over_time(prop_sum_order,
                                  order_palette_bac, 
                                  "Order", 
                                  paste("16S by", ind_phylum, "phylum",ind_class, "class separated by order")
  )
  print(p)
 } 
 dev.off()
}

## want to divide up by phyla and then look at the orders within a specific class.

for (ind_phylum in unique(long_18s_prop_with_taxonomy$Phylum)) {
 #print(ind_phylum)
 phyla_subset <- droplevels(subset(long_18s_prop_with_taxonomy, Phylum == ind_phylum))
 pdf(paste0(figures_dir, "18S_time_series_phylum_",ind_phylum,"_by_specific_orders_barplots%03d.pdf"), width = 15, height = 11, onefile = FALSE)
 for (ind_class in unique(phyla_subset$Class)){
  #print(ind_class)
  class_subset <- droplevels(subset(phyla_subset, Class == ind_class))
  prop_sum_order <- summarise_otus_by_taxon(class_subset, "Order")
 # print(prop_sum_order)
  p <- barplot_all_taxa_over_time(prop_sum_order,
                                  order_palette_euks, 
                                  "Order", 
                                  paste("18S by", ind_phylum, "phylum",ind_class, "class separated by order")
  )
  print(p)
 } 
 dev.off()
}

## do by area
## couldn't this replace all the other stuff??
for (ind_phylum in unique(long_18s_prop_with_taxonomy$Phylum)) {
 phyla_subset <- droplevels(subset(long_18s_prop_with_taxonomy, Phylum == ind_phylum))
 pdf(paste0(figures_dir,"18S_time_series_phylum_",ind_phylum,"_by_specific_orders_areaplots%03d.pdf"), width = 15, height = 11, onefile = FALSE)
 for (ind_class in unique(phyla_subset$Class)){
 # print(ind_class)
  class_subset <- droplevels(subset(phyla_subset, Class == ind_class))
  prop_sum_order <- summarise_otus_by_taxon(class_subset, "Order")
  #print(prop_sum_order)
  p <- area_all_taxa_over_time(prop_sum_order,
                               order_palette_euks, 
                               "Order", 
                               paste("18S by", ind_phylum, "phylum",ind_class, "class separated by order")
  )
  print(p)
 } 
 dev.off()
}

### Print proportional areas 16s subsections ####
# pdf("../figures/16S_proportional_area_time_series_subsection%03d.pdf", width = 15, height = 6, onefile = FALSE)
# proportional_area_16s_class_no_unclassified_Proteo_phylum
# proportional_area_16s_class_no_unclassified_Actino_phylum
# proportional_area_16s_class_no_unclassified_Bacteroi_phylum
# proportional_area_16s_order_no_unclassified_Cyano_phylum
# 
# dev.off()


# ### Print proportional areas 16s subsections ####
# pdf(paste("../figures/16S_proportional_area_time_series_subsection_order%03d.pdf", sep=""), width = 15, height = 6, onefile = FALSE)
# proportional_area_16s_order_no_unclassified_Proteo_phylum
# proportional_area_16s_order_no_unclassified_Actino_phylum
# proportional_area_16s_order_no_unclassified_Bacteroi_phylum
# #proportional_area_16s_genus_no_unclassified_Cyano_phylum
# dev.off()

plot_otus_from_specic_tax_group <- function (long_otus_with_taxonomy, title) {
 
 colourCount  <- length(unique(long_otus_with_taxonomy$otu_number))
 getPalette  <-  colorRampPalette(brewer.pal(11, "Spectral"))
 
 ggplot(long_otus_with_taxonomy,
        aes(x=Date, y=value,group=otu_number))+
  season_line +
  season_text +
  spring_bloom_line+
  geom_area(aes(fill=otu_number),
            position="stack")+
  scale_fill_manual(values = getPalette(colourCount))+
  date_scaling+
 # ggtitle(title)+
 theme_JAG_presentation()+
  theme(axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())
}

pdf(paste0(figures_dir,"18S_time_series_by_order_specific_OTUs%03d.pdf"), width = 15, height = 11, onefile = FALSE)
#plot_otus_from_specic_tax_group(droplevels(subset(long_18s_prop_with_taxonomy,  Order == "Chlorophyta")) , "18S only Chlorophyta")
plot_otus_from_specic_tax_group(droplevels(subset(long_18s_prop_with_taxonomy,  Order == "Ochrophyta")), "18S only Ochrophyta")
plot_otus_from_specic_tax_group(droplevels(subset(long_18s_prop_with_taxonomy,  Order == "Dinoflagellata")), "18S only Dinoflagellata")
dev.off()

pdf(paste0(figures_dir,"16S_time_series_by_order_specific_OTUs%03d.pdf"), width = 15, height = 11, onefile = FALSE)
plot_otus_from_specic_tax_group(droplevels(subset(long_16s_prop_with_taxonomy,  Order == "Flavobacteriales")), "16S only Flavobacteriales" )
plot_otus_from_specic_tax_group(droplevels(subset(long_16s_prop_with_taxonomy,  Phylum == "Cyanobacteria")), "16S only Cyanobacteria" )
plot_otus_from_specic_tax_group(droplevels(subset(long_16s_prop_with_taxonomy,  Phylum == "Proteobacteria")), "16S only Proteobacteria" )
plot_otus_from_specic_tax_group(droplevels(subset(long_16s_prop_with_taxonomy,  Phylum == "Bacteroidetes")), "16S only Bacteroidetes" )
dev.off()

## Use summarised data to look at this one. 
plot_otus_from_specific_order <- function (long_otus_with_taxonomy,tax_palette, title) {
 ggplot(long_otus_with_taxonomy, 
        aes(x=Date, y=rel_abun,group=Order))+
  season_line +
  season_text +
  spring_bloom_line+
  geom_area(aes(fill=Order),
            position="stack",
            colour=line_colour)+
  scale_fill_manual(values = tax_palette)+
  date_scaling+
  #ggtitle(title)+
  theme_JAG_presentation()+
  theme(axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())
}

# prop_16s_summarised_by_order_Bac <- droplevels(summarise_otus_by_taxon(droplevels(subset(long_16s_prop_with_taxonomy, Phylum == "Bacteroidetes")), "Order"))
# 
# area_Bacteroidetes <- plot_otus_from_specific_order(prop_16s_summarised_by_order_Bac, order_palette_bac, "16S only Bacteroidetes" )
# 
# prop_16s_summarised_by_order_Pro <- droplevels(summarise_otus_by_taxon(droplevels(subset(long_16s_prop_with_taxonomy, Phylum == "Proteobacteria")),"Order"))
# 
# area_Proteobacteria <- plot_otus_from_specific_order(prop_16s_summarised_by_order_Pro, order_palette_bac, "16S only Proteobacteria" )
# prop_16s_summarised_by_order_Cyan <- droplevels(summarise_otus_by_taxon(droplevels(subset(long_16s_prop_with_taxonomy, Phylum == "Cyanobacteria")), "Order"))
# 
# area_Cyan <- plot_otus_from_specific_order(prop_16s_summarised_by_order_Cyan,
#                                            order_palette_bac,
#                                            "16S only Cyanobacteria" )
# 
# 
# pdf("../figures/16S_time_series_by_order_classes%03d.pdf", width = 15, height = 6, onefile = FALSE)
# area_Bacteroidetes
# area_Proteobacteria
# area_Cyan
# dev.off()

plot_otus_from_specific_class <- function (long_otus_with_taxonomy, tax_palette, title) {
 ggplot(long_otus_with_taxonomy,
        aes(x=Date, y=rel_abun,group=Family))+
  season_line +
  season_text +
  spring_bloom_line+
  geom_area(aes(fill=Family),
            position="stack")+
  scale_fill_manual(values = tax_palette)+
  date_scaling+
 # ggtitle(title)+
 theme_JAG_presentation()+
  theme(axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())
}

# prop_18s_summarised_by_order_Chlor <- droplevels(summarise_otus_by_taxon(droplevels(subset(long_18s_prop_with_taxonomy, Order == "Chlorophyta")),
#                                                                          "Family"))
# area_Chlorophyta <- plot_otus_from_specific_class(prop_18s_summarised_by_order_Chlor,
#                                                   family_palette_euks,
#                                                   "18S only Chlorophyta")
# prop_18s_summarised_by_order_Och <- droplevels(summarise_otus_by_taxon(droplevels(subset(long_18s_prop_with_taxonomy,Order == "Ochrophyta")),"Family"))
# 
# area_Ochrophyta <- plot_otus_from_specific_class(prop_18s_summarised_by_order_Och,family_palette_euks,"18S only Ochrophyta")
# 
# prop_18s_summarised_by_order_Dino <- droplevels(summarise_otus_by_taxon(droplevels(subset(long_18s_prop_with_taxonomy,Order == "Dinoflagellata")), "Family"))
# 
# area_Dinoflagellata <- plot_otus_from_specific_class(prop_18s_summarised_by_order_Dino, 
#                                                      family_palette_euks,
#                                                      "18S only Dinoflagellata")
# 
# 
# pdf("../figures/18S_time_series_by_Family_orders%03d.pdf", width = 15, height = 6, onefile = FALSE)
# area_Chlorophyta
# area_Ochrophyta
# area_Dinoflagellata
# dev.off()
# 

##########################################
## Try to separate by het and phototrophic 18s
## like fig 2 in de Vargas et al 2015

phytoplankton_fam <- c("Diatomea",
                       "Dinophyceae",
                       "Mamiellophyceae",
                       "Chrysophyceae",
                       "Chlorophyceae",
                       "Euglenozoa",
                       "Dictyochophyceae",
                       "Chlorarachniophyta",
                       "Raphidophyceae")
phytoplankton_class <- c("Rhodophyceae",
                         "Chloroplastida")
phytoplantkon_phylum <- c("Haptophyta",
                          "Cryptophyceae")
phytoplankton_order <- c("Dinoflagellata",
                         "Ochrophyta") ## not sure about this, are all dinos phytoplankton?

## algae but maybe not phytoplankton,
algae_class <- c("Chloroplastida")

## masts we don't know...
phytoplankton_tax <- taxonomy_18s %>% 
 filter(Family %in% phytoplankton_fam | Class %in% phytoplankton_class | Phylum %in% phytoplantkon_phylum| Order %in% phytoplankton_order)

phytoplankton_otus <- droplevels(subset(long_18s_prop_with_taxonomy, otu_number %in% phytoplankton_tax$otu_number)) 

prop_18s_summarised_by_family_phytoplankton <- droplevels(summarise_otus_by_taxon(phytoplankton_otus,
                                                                                 "Family"))

prop_18s_phytos_summarised_by_order <- droplevels(summarise_otus_by_taxon(phytoplankton_otus, "Order"))

prop_18s_phytos_summarised_by_class <- droplevels(summarise_otus_by_taxon(phytoplankton_otus, "Class"))
prop_18s_phytos_summarised_by_phylum <- droplevels(summarise_otus_by_taxon(phytoplankton_otus, "Phylum"))

area_all_taxa_over_time(prop_18s_phytos_summarised_by_order,
                             order_palette_euks, 
                             "Order", 
                             "18S byphylumclass separated by order")

area_all_taxa_over_time(prop_18s_summarised_by_family_phytoplankton,
                        family_palette_euks, 
                        "Family", 
                        "18S byphylumclass separated by family")

area_all_taxa_over_time(prop_18s_phytos_summarised_by_class,
                        class_palette_euks, 
                        "Class", 
                        "18S phytos by class")

area_all_taxa_over_time(prop_18s_phytos_summarised_by_phylum,
                        phyla_palette_euks, 
                        "Phylum", 
                        "18S phytos by phylum")

## phaedoraea, nasselaria, spumellaria, collodaria
protozooplankton_order <- c("Ciliophora")
protozooplankton_family <- c("Acanthoecida",
                              "Foraminifera")

#protozooplankton_genus <- c(
                            
protozooplankton_tax <- taxonomy_18s %>% 
 filter(Family %in% protozooplankton_family | Order %in% protozooplankton_order)

protozooplankton_otus <- droplevels(subset(long_18s_prop_with_taxonomy, otu_number %in% protozooplankton_tax$otu_number)) 

prop_18s_summarised_by_family_protozooplankton <- droplevels(summarise_otus_by_taxon(protozooplankton_otus,"Family"))

prop_18s_protozoop_summarised_by_order <- droplevels(summarise_otus_by_taxon(protozooplankton_otus,"Order"))

prop_18s_protozoop_summarised_by_class <- droplevels(summarise_otus_by_taxon(protozooplankton_otus, "Class"))
prop_18s_protozoop_summarised_by_phylum <- droplevels(summarise_otus_by_taxon(protozooplankton_otus, "Phylum"))

area_all_taxa_over_time(prop_18s_protozoop_summarised_by_order,
                        order_palette_euks, 
                        "Order", 
                        "18S protozoop separated by order")

area_all_taxa_over_time(prop_18s_summarised_by_family_protozooplankton,
                        family_palette_euks, 
                        "Family", 
                        "18S protozoop separated by family")

area_all_taxa_over_time(prop_18s_protozoop_summarised_by_class,
                        class_palette_euks, 
                        "Class", 
                        "18S protozoop by class")

area_all_taxa_over_time(prop_18s_protozoop_summarised_by_phylum,
                        phyla_palette_euks, 
                        "Phylum", 
                        "18S protozoop by phylum")

metazooplantkon_order <- c("Metazoa",
                           "Apicomplexa"
                           )

## or coudl try approach that heterotrophs are those that are not phytoplankton

heterotroph_tax <- taxonomy_18s %>% 
 filter(!(otu_number %in% phytoplankton_tax$otu_number))

heterotroph_otus <- droplevels(subset(long_18s_prop_with_taxonomy, otu_number %in% heterotroph_tax$otu_number)) 

prop_18s_summarised_by_family_het <- droplevels(summarise_otus_by_taxon(heterotroph_otus,
                                                                                     "Family"))
prop_18s_het_summarised_by_order <- droplevels(summarise_otus_by_taxon(heterotroph_otus,
                                                                             "Order"))
prop_18s_het_summarised_by_class <- droplevels(summarise_otus_by_taxon(heterotroph_otus,
                                                                             "Class"))
prop_18s_het_summarised_by_phylum <- droplevels(summarise_otus_by_taxon(heterotroph_otus,
                                                                              "Phylum"))

area_all_taxa_over_time(prop_18s_het_summarised_by_order,
                        order_palette_euks, 
                        "Order", 
                        "18S het separated by order")

area_all_taxa_over_time(prop_18s_summarised_by_family_het,
                        family_palette_euks, 
                        "Family", 
                        "18S protozoop separated by family")

area_all_taxa_over_time(prop_18s_het_summarised_by_class,
                        class_palette_euks, 
                        "Class", 
                        "18S protozoop by class")

area_all_taxa_over_time(prop_18s_het_summarised_by_phylum,
                        phyla_palette_euks, 
                        "Phylum", 
                        "18S protozoop by phylum")
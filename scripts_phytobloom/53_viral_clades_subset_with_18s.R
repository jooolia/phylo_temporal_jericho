

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
## virus groups correlate to env params?
library(Hmisc)
library(stringr)
library(cowplot)



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
                          colour="grey90",
                          size=1.5)
spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=1)
date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")))

date_scaling_bloom <-   scale_x_date(limits = c(as.Date("2011-06-20"), as.Date("2011-07-06")))

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


#MPL_group_sums <- read.csv("../results/MPL_group_sums_by_site_prop.csv")
#gp23_group_sums <- read.csv("../results/gp23_group_sums_by_site_prop.csv")

### do it over the whole time

normalized_18s_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_normalized_18S.tsv", row.names="VC_number")
colnames(normalized_18s_OTUs) <- gsub(".size.*.", "", colnames(normalized_18s_OTUs))


taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)

melted_otus <- melt(as.matrix(normalized_18s_OTUs))

melted_otus$Order <- taxonomy_18s$Order[match(melted_otus$Var2, taxonomy_18s$otu_number)]
## remove unclassified
melted_otus <- subset(melted_otus, Order != "unclassified"  )

ord_melted <- melted_otus %>%
 group_by(Var1,Order) %>%
 summarise(total=sum(value))

## Ok need to put it back together again.
## Families on top
ord_sum_otus <- dcast(ord_melted, Var1 ~  Order)
## need to make them rownames again
row.names(ord_sum_otus) <- ord_sum_otus$Var1
ord_sum_otus <- ord_sum_otus[,-1]
S18_ord_sums <- ord_sum_otus 

write.csv(ord_sum_otus, file="../results/normalized_18s_summarized_by_order.csv")

melted_otus <- melt(as.matrix(normalized_18s_OTUs))

 melted_otus$Family <- taxonomy_18s$Family[match(melted_otus$Var2, taxonomy_18s$otu_number)]
# ## remove unclassified
melted_otus <- subset(melted_otus, Family != "unclassified"  )
# 
fam_melted <- melted_otus %>%
 group_by(Var1,Family) %>%
 summarise(total=sum(value))
# 
# ## Ok need to put it back together again.
# ## Families on top
fam_sum_otus <- dcast(fam_melted, Var1 ~  Family)
# ## need to make them rownames again
 row.names(fam_sum_otus) <- fam_sum_otus$Var1
 fam_sum_otus <- fam_sum_otus[,-1]
 S18_fam_sums <- fam_sum_otus
 
 write.csv(fam_sum_otus, file="../results/normalized_18s_summarized_by_family.csv")
 

## or actually want to pull out the OTus that are in the group and also in the Raphidophyte group....

normalized_16s_OTUs <- read.delim("../data/OTU_table_bloom_with_all_time_series_16S_R1.tsv", row.names="VC_number") 
colnames(normalized_16s_OTUs) <- gsub(".size.*.", "", colnames(normalized_16s_OTUs))
proportional_18s <- as.data.frame(prop.table(as.matrix(normalized_18s_OTUs),  margin=1))
proportional_16s <- as.data.frame(prop.table(as.matrix(normalized_16s_OTUs),  margin=1))

taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)

melted_otus <- melt(as.matrix(normalized_16s_OTUs))

melted_otus$Order <- taxonomy_16s$Order[match(melted_otus$Var2, taxonomy_16s$otu_number)]
## remove unclassified
melted_otus <- subset(melted_otus, Order != "unclassified"  )

ord_melted <- melted_otus %>%
 group_by(Var1,Order) %>%
 summarise(total=sum(value))

## Ok need to put it back together again.
## Families on top
ord_sum_otus <- dcast(ord_melted, Var1 ~  Order)
## need to make them rownames again
row.names(ord_sum_otus) <- ord_sum_otus$Var1
ord_sum_otus <- ord_sum_otus[,-1]
S16_ord_sums <- ord_sum_otus 

write.csv(ord_sum_otus, file="../results/normalized_16s_summarized_by_order.csv")

melted_otus <- melt(as.matrix(normalized_16s_OTUs))

melted_otus$Family <- taxonomy_16s$Family[match(melted_otus$Var2, taxonomy_16s$otu_number)]
# ## remove unclassified
melted_otus <- subset(melted_otus, Family != "unclassified"  )
# 
fam_melted <- melted_otus %>%
 group_by(Var1,Family) %>%
 summarise(total=sum(value))
# 
# ## Ok need to put it back together again.
# ## Families on top
fam_sum_otus <- dcast(fam_melted, Var1 ~  Family)
# ## need to make them rownames again
row.names(fam_sum_otus) <- fam_sum_otus$Var1
fam_sum_otus <- fam_sum_otus[,-1]
S16_fam_sums <- fam_sum_otus

write.csv(fam_sum_otus, file="../results/normalized_16s_summarized_by_family.csv")



## try to split up by taxonomy


##########################
## Try to separate by het and phototrophic 18s
## like fig 2 in de Vargas et al 2015

raphidophyte_fam <- "Raphidophyceae"

diatom_fam <- "Diatomea"

ciliate_ord <- "Ciliophora"

dinoflagellate_ord <- "Dinoflagellata"

raphidophyte_tax <- taxonomy_18s %>% 
 filter(Family %in% raphidophyte_fam)

diatom_tax <- taxonomy_18s %>% 
 filter(Family %in% diatom_fam)

ciliate_tax <- taxonomy_18s %>%
 filter(Order %in% ciliate_ord)

dinoflagellate_tax <- taxonomy_18s %>%
 filter(Order %in% dinoflagellate_ord)


colnames(proportional_18s) <- gsub(".size.*.", "", colnames(proportional_18s))

raphidophyte_otu_table <- droplevels(subset(proportional_18s, select=colnames(proportional_18s) %in% raphidophyte_tax$otu_number)) 
write.csv(raphidophyte_otu_table, "../results/normalized_18s_raphidophyte_otus.csv")

diatom_otu_table <- droplevels(subset(proportional_18s, select=colnames(proportional_18s) %in% diatom_tax$otu_number)) 

ciliate_otu_table <- droplevels(subset(proportional_18s, select=colnames(proportional_18s) %in% ciliate_tax$otu_number)) 

dino_otu_table <- droplevels(subset(proportional_18s, select=colnames(proportional_18s) %in% dinoflagellate_tax$otu_number)) 

names(S18_ord_sums)[1] <- "VC"
names(S18_fam_sums)[1] <- "VC"

S16_ord_sums <- read.csv("../results/normalized_16s_summarized_by_order.csv")

names(S16_ord_sums)[1] <- "VC"

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis_all_time_series.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

parameters_to_exclude <- c("season",
                           # "Average_NO3_NO2",
                           "Standard_error_NO3_NO2",
                           "Secchi_disk_disappears",
                           "Secchi_disk_reappears",
                           #"VC_number",
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

summary(params_Jericho)
#params_Jericho <- params_Jericho[,-12]
## get those in the these samples 

## RdRp

# proportional_MPL_groupA <- read.csv("../results/proportional_MPL_groupA.csv", row.names = 1)
# proportional_MPL_groupB <- read.csv("../results/proportional_MPL_groupB.csv", row.names = 1)
# #proportional_MPL_groupC <- read.csv("../results/proportional_MPL_groupC.csv", row.names = 1)
# proportional_MPL_groupD <- read.csv("../results/proportional_MPL_groupD.csv", row.names = 1)
# proportional_MPL_groupE <- read.csv("../results/proportional_MPL_groupE.csv", row.names = 1)
# #proportional_MPL_groupF <- read.csv("../results/proportional_MPL_groupF.csv", row.names = 1)
# proportional_MPL_groupG <- read.csv("../results/proportional_MPL_groupG.csv", row.names = 1)
# proportional_MPL_groupH <- read.csv("../results/proportional_MPL_groupH.csv", row.names = 1)

## first plot the different stuff



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
# long_rdrp_groupA_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupA,
#                                                   Jericho_data)
# long_rdrp_groupB_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupB,
#                                                          Jericho_data)
# long_rdrp_groupE_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupE,
#                                                          Jericho_data)
# long_rdrp_groupG_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupG,
#                                                          Jericho_data)
# long_rdrp_groupH_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupH,
#                                                          Jericho_data)
long_raphidophyte_otus <- add_date_to_OTU_table_rtrn_long(raphidophyte_otu_table,
                                                         Jericho_data)
long_diatom_otus <- add_date_to_OTU_table_rtrn_long(diatom_otu_table,
                                                          Jericho_data)

long_ciliate_otus <- add_date_to_OTU_table_rtrn_long(ciliate_otu_table,
                                                    Jericho_data)

long_dino_otus <- add_date_to_OTU_table_rtrn_long(dino_otu_table,
                                                     Jericho_data)

## Plot single OTUs over time
plot_individual_otus_over_time <- function (long_otus, title, type_OTU) {
 #cols <- colorRampPalette(brewer.pal(8,"Dark2"))(nrow(long_otus))
 colourCount  <- length(unique(long_otus$OTUid))
mypalette  <-  brewer.pal(10, "Spectral")
 if (type_OTU == "virus") {
  specific_palette <- mypalette[1:5] 
 }
 else if (type_OTU == "host"){
  specific_palette <- mypalette[6:10]
 }
 getPalette  <-  colorRampPalette(specific_palette)
 sampling_points <- long_otus$Date
 p <- ggplot(droplevels(long_otus), aes(x=Date, y=value, group=OTUid))+
  season_line +
  spring_bloom_line+
  HAKA_bloom+
  geom_area(aes(fill=OTUid), position="stack", colour=line_colour)+
  date_scaling +
  scale_fill_manual(values = getPalette(colourCount), guide = FALSE)+
  theme_JAG_presentation(base_size=18)+
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  ylab("Relative abundance")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(face="bold", size=20),
        axis.text.y = element_text(size=18),
        axis.line=element_line(size = 0.3))
#  +
#   ggtitle(title)
}

## Plot single OTUs over time
plot_individual_otus_over_time_only_bloom <- function (long_otus, title, type_OTU) {
 #cols <- colorRampPalette(brewer.pal(8,"Dark2"))(nrow(long_otus))
 colourCount  <- length(unique(long_otus$OTUid))
 mypalette  <-  brewer.pal(10, "Spectral")
 if (type_OTU == "virus") {
  specific_palette <- mypalette[1:5] 
 }
 else if (type_OTU == "host"){
  specific_palette <- mypalette[6:10]
 }
 getPalette  <-  colorRampPalette(specific_palette)
 sampling_points <- long_otus$Date
 p <- ggplot(droplevels(long_otus), aes(x=Date, y=value, group=OTUid))+
  HAKA_bloom+
  season_line +
  spring_bloom_line+
  geom_area(aes(fill=OTUid), position="stack", colour=line_colour)+
  date_scaling_bloom +
  scale_fill_manual(values = getPalette(colourCount), guide = FALSE)+
  annotate("point", 
           x = as.Date(sampling_points),
           y = -0.02,
           colour = "grey",
           shape= 17,
           # solid=TRUE,
           size = 5)+
  theme_JAG_presentation(base_size=18)+
  ylab("Relative abundance")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(face="bold", size=20),
        axis.text.y = element_text(size=18),
        axis.line=element_line(size = 0.3))
#  +
#   ggtitle(title)
}

filtered_long_raphidophyte <- long_raphidophyte_otus

filtered_long_diatom <- long_diatom_otus

filtered_long_ciliate <- long_ciliate_otus

filtered_long_dino <- long_dino_otus

# group_A_rdrp <- plot_individual_otus_over_time(long_rdrp_groupA_otus, "RdRp group a", "virus")
# group_B_rdrp <- plot_individual_otus_over_time(long_rdrp_groupB_otus, "rdrp group b", "virus")
# group_E_rdrp <- plot_individual_otus_over_time(long_rdrp_groupE_otus, "rdrp group E", "virus")
# group_G_rdrp <- plot_individual_otus_over_time(long_rdrp_groupG_otus, "rdrp group G", "virus")
# group_H_rdrp <- plot_individual_otus_over_time(long_rdrp_groupH_otus, "rdrp group H", "virus")

raphidophyte_over_time_matching_rdrp<- plot_individual_otus_over_time(filtered_long_raphidophyte , "Raphidophyte", "host")

raphidophyte_over_time_bloom <- plot_individual_otus_over_time_only_bloom(filtered_long_raphidophyte , "Raphidophyte", "host")

diatom_over_time_matching_rdrp<- plot_individual_otus_over_time(filtered_long_diatom , "Diatom", "host")

diatom_over_time_bloom<- plot_individual_otus_over_time_only_bloom(filtered_long_diatom , "Diatom", "host")

ciliate_over_time_matching_rdrp<- plot_individual_otus_over_time(filtered_long_ciliate, "Ciliate", "host")

ciliate_over_time_bloom<- plot_individual_otus_over_time_only_bloom(filtered_long_ciliate, "Ciliate", "host")

dino_over_time_matching_rdrp<- plot_individual_otus_over_time(filtered_long_dino, "Dinoflagellate", "host")

dino_over_time_bloom<- plot_individual_otus_over_time_only_bloom(filtered_long_dino, "Dinoflagellate", "host")

pdf(file="../figures/Raphidophytes_over_all_times.pdf",onefile = FALSE)
raphidophyte_over_time_matching_rdrp
dev.off()

pdf(file="../figures/Raphidophytes_over_bloom.pdf",onefile = FALSE)
raphidophyte_over_time_bloom
dev.off()

pdf(file="../figures/Ciliates_over_all_times.pdf",onefile = FALSE)
ciliate_over_time_matching_rdrp
dev.off()

pdf(file="../figures/Ciliates_over_bloom.pdf",onefile = FALSE)
ciliate_over_time_bloom
dev.off()

pdf(file="../figures/Diatoms_over_all_times.pdf", onefile = FALSE)
diatom_over_time_matching_rdrp
dev.off()

pdf(file="../figures/Diatoms_over_bloom.pdf", onefile = FALSE)
diatom_over_time_bloom
dev.off()

pdf(file="../figures/Dinoflagellates_over_all_times.pdf", onefile = FALSE)
dino_over_time_matching_rdrp
dev.off()

pdf(file="../figures/Dinoflagellates_over_bloom.pdf", onefile = FALSE)
dino_over_time_bloom
dev.off()


pdf(file="../figures/Spefic_taxa_over_bloom_and_all_time.pdf", onefile = FALSE, width = 20, height = 20)
plot_grid(raphidophyte_over_time_bloom+xlab(NULL)+
           theme(axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank()),
          raphidophyte_over_time_matching_rdrp+xlab(NULL)+
           theme(axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank()),
          ciliate_over_time_bloom+xlab(NULL)+
           theme(axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank()),
          ciliate_over_time_matching_rdrp+xlab(NULL)+
           theme(axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank()),
          diatom_over_time_bloom+xlab(NULL)+
           theme(axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank()),
          diatom_over_time_matching_rdrp+xlab(NULL)+
           theme(axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank()),
          dino_over_time_bloom,
          dino_over_time_matching_rdrp+
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()),
          ncol = 2,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          label_size = 24)
dev.off()


# pdf(file="../figures/group_G_with_Diatoms.pdf", width = 20, height = 17,onefile = FALSE)
# grid.arrange(diatom_over_time_matching_rdrp,
#              group_G_rdrp,
#              ncol=1)
# dev.off()


# pdf(file="../figures/group_H_with_Diatoms.pdf", width = 20, height = 17,onefile = FALSE)
# grid.arrange(diatom_over_time_matching_rdrp,
#              group_H_rdrp,
#              ncol=1)
# dev.off()



### what about with gp23 group I and cyanobacteria

## gp23
# 
# proportional_gp23_groupA <- read.csv("../results/proportional_gp23_groupA.csv", row.names = 1)
# proportional_gp23_groupB <- read.csv("../results/proportional_gp23_groupB.csv", row.names = 1)
# proportional_gp23_groupC <- read.csv("../results/proportional_gp23_groupC.csv", row.names = 1)
# proportional_gp23_groupD <- read.csv("../results/proportional_gp23_groupD.csv", row.names = 1)
# proportional_gp23_groupE <- read.csv("../results/proportional_gp23_groupE.csv", row.names = 1)
# proportional_gp23_groupF <- read.csv("../results/proportional_gp23_groupF.csv", row.names = 1)
# proportional_gp23_groupG <- read.csv("../results/proportional_gp23_groupG.csv", row.names = 1)
# proportional_gp23_groupH <- read.csv("../results/proportional_gp23_groupH.csv", row.names = 1)
# proportional_gp23_groupI <- read.csv("../results/proportional_gp23_groupI.csv", row.names = 1)
# #proportional_gp23_groupJ <- read.csv("../results/proportional_gp23_groupJ.csv", row.names = 1)
### here!
cyanobacteria_phylum <- c("Cyanobacteria")
cyanobacteria_class <- c("Cyanobacteria")

#diatom_fam <- c("Diatomea") ## not sure about this, are all dinos phytoplankton?

cyano_tax <- taxonomy_16s %>% 
 filter(Class %in% cyanobacteria_class )

colnames(proportional_16s) <- gsub(".size.*.", "", colnames(proportional_16s))

cyanobacteria_otu_table <- droplevels(subset(proportional_16s, select=colnames(proportional_16s) %in% cyano_tax$otu_number)) 

## try phylum

cyano_tax_phyla <- taxonomy_16s %>% 
 filter(Phylum %in% cyanobacteria_phylum )

cyanobacteria_otu_table_phyla <- droplevels(subset(proportional_16s, select=colnames(proportional_16s) %in% cyano_tax_phyla$otu_number)) 
write.csv(cyanobacteria_otu_table_phyla, "../results/normalized_16s_cyano_otus.csv")
# 
# ## Read in OTU table and get long format with Date
# long_gp23_groupI_otus <- add_date_to_OTU_table_rtrn_long(proportional_gp23_groupI,
#                                                          Jericho_data)
# 
# #long_rdrp_groupB_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupB,
#                                                         Jericho_data)
long_cyano_otus <- add_date_to_OTU_table_rtrn_long(cyanobacteria_otu_table,
                                                          Jericho_data)

#                                                         Jericho_data)
long_cyano_otus_phyla <- add_date_to_OTU_table_rtrn_long(cyanobacteria_otu_table_phyla,
                                                   Jericho_data)

#filtered_long_cyano <- long_cyano_otus %>% 
# filter(Date < as.Date("2011-06-01") )
# 
# group_I_gp23 <- plot_individual_otus_over_time(long_gp23_groupI_otus, "gp23 group I", "virus")
# #group_B_rdrp <- plot_individual_otus_over_time(long_rdrp_groupB_otus, "rdrp group b")
# 
cyano_over_time_matching_gp23 <- plot_individual_otus_over_time(long_cyano_otus , "Class cyano", "host")

cyano_over_time_matching_gp23 <- cyano_over_time_matching_gp23+
 theme(axis.title.x = element_blank(),
       axis.text.x  = element_blank()
 )

cyano_over_time_matching_gp23_phyla <- plot_individual_otus_over_time(long_cyano_otus_phyla,
                                                                      "phyla cyano", "host") +
 theme(axis.title.x = element_blank(),
       axis.text.x  = element_blank()
)

# pdf(file="../figures/group_I_with_Cyanos.pdf", width = 20, height = 17, onefile = FALSE)
# plot_grid(cyano_over_time_matching_gp23,
#           group_I_gp23,
#           align = "v",
#           ncol=1,
#           labels = c("A", "B"),
#           label_size = 20)
# dev.off()
# 
# pdf(file="../figures/group_I_with_Cyanos_phyla.pdf", width = 20, height = 17, onefile = FALSE)
# plot_grid(cyano_over_time_matching_gp23_phyla,
#           group_I_gp23,
#           align = "v",
#           ncol = 1,
#           labels = c("A", "B"),
#           label_size = 20)
# dev.off()


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


MPL_group_sums <- read.csv("../results/MPL_group_sums_by_site_prop.csv")
gp23_group_sums <- read.csv("../results/gp23_group_sums_by_site_prop.csv")

S18_ord_sums <- read.csv("../results/normalized_18s_summarized_by_order.csv")
S18_fam_sums <- read.csv("../results/normalized_18s_summarized_by_family.csv")

## or actually want to pull out the OTus that are in the group and also in the Raphidophyte group....

S16_class_sums <- read.csv("../results/normalized_16s_summarized_by_class.csv")


normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 
proportional_18s <- as.data.frame(prop.table(as.matrix(normalized_18s_OTUs),  margin=1))
proportional_16s <- as.data.frame(prop.table(as.matrix(normalized_16s_OTUs),  margin=1))

taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)
taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)


## try to split up by taxonomy


##########################
## Try to separate by het and phototrophic 18s
## like fig 2 in de Vargas et al 2015

raphidophyte_fam <- c("Raphidophyceae")

diatom_fam <- c("Diatomea") ## not sure about this, are all dinos phytoplankton?

raphidophyte_tax <- taxonomy_18s %>% 
 filter(Family %in% raphidophyte_fam)

diatom_tax <- taxonomy_18s %>% 
 filter(Family %in% diatom_fam)

colnames(proportional_18s) <- gsub(".size.*.", "", colnames(proportional_18s))

raphidophyte_otu_table <- droplevels(subset(proportional_18s, select=colnames(proportional_18s) %in% raphidophyte_tax$otu_number)) 
write.csv(raphidophyte_otu_table, "../results/normalized_18s_raphidophyte_otus.csv")

diatom_otu_table <- droplevels(subset(proportional_18s, select=colnames(proportional_18s) %in% diatom_tax$otu_number)) 


names(S18_ord_sums)[1] <- "VC"
names(S18_fam_sums)[1] <- "VC"

S16_ord_sums <- read.csv("../results/normalized_16s_summarized_by_order.csv")

names(S16_ord_sums)[1] <- "VC"

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
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

proportional_MPL_groupA <- read.csv("../results/proportional_MPL_groupA.csv", row.names = 1)
proportional_MPL_groupB <- read.csv("../results/proportional_MPL_groupB.csv", row.names = 1)
#proportional_MPL_groupC <- read.csv("../results/proportional_MPL_groupC.csv", row.names = 1)
proportional_MPL_groupD <- read.csv("../results/proportional_MPL_groupD.csv", row.names = 1)
proportional_MPL_groupE <- read.csv("../results/proportional_MPL_groupE.csv", row.names = 1)
#proportional_MPL_groupF <- read.csv("../results/proportional_MPL_groupF.csv", row.names = 1)
proportional_MPL_groupG <- read.csv("../results/proportional_MPL_groupG.csv", row.names = 1)
proportional_MPL_groupH <- read.csv("../results/proportional_MPL_groupH.csv", row.names = 1)

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
long_rdrp_groupA_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupA,
                                                  Jericho_data)
long_rdrp_groupB_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupB,
                                                         Jericho_data)
long_rdrp_groupE_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupE,
                                                         Jericho_data)
long_rdrp_groupG_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupG,
                                                         Jericho_data)
long_rdrp_groupH_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupH,
                                                         Jericho_data)
long_raphidophyte_otus <- add_date_to_OTU_table_rtrn_long(raphidophyte_otu_table,
                                                         Jericho_data)
long_diatom_otus <- add_date_to_OTU_table_rtrn_long(diatom_otu_table,
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
 p <- ggplot(droplevels(long_otus), aes(x=Date, y=value, group=OTUid))+
  season_line +
  spring_bloom_line+
  geom_area(aes(fill=OTUid), position="stack", colour=line_colour)+
  date_scaling +
  scale_fill_manual(values = getPalette(colourCount), guide = FALSE)+
  theme_JAG_presentation()+
  ylab("Relative abundance")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
 #+

}

filtered_long_raphidophyte <- long_raphidophyte_otus %>% 
 filter(Date < as.Date("2011-06-01") )

filtered_long_diatom <- long_diatom_otus %>% 
 filter(Date < as.Date("2011-06-01") )

group_A_rdrp <- plot_individual_otus_over_time(long_rdrp_groupA_otus, "RdRp group a", "virus")
group_B_rdrp <- plot_individual_otus_over_time(long_rdrp_groupB_otus, "rdrp group b", "virus")
group_E_rdrp <- plot_individual_otus_over_time(long_rdrp_groupE_otus, "rdrp group E", "virus")
group_G_rdrp <- plot_individual_otus_over_time(long_rdrp_groupG_otus, "rdrp group G", "virus")
group_H_rdrp <- plot_individual_otus_over_time(long_rdrp_groupH_otus, "rdrp group H", "virus")

raphidophyte_over_time_matching_rdrp<- plot_individual_otus_over_time(filtered_long_raphidophyte , "Raphidophyte", "host")
raphidophyte_over_time_matching_rdrp <- raphidophyte_over_time_matching_rdrp+  theme(axis.title.x = element_blank(),
       axis.text.x  = element_blank()
       )
diatom_over_time_matching_rdrp<- plot_individual_otus_over_time(filtered_long_diatom , "Diatom", "host")

pdf(paste0(figures_dir,"group_A_with_Raphidophytes.pdf"), width = 14, height = 14,onefile = FALSE)
plot_grid(raphidophyte_over_time_matching_rdrp,
          group_A_rdrp,
          ncol = 1,
          align = "v",
          labels = c("A", "B"),
          label_size = 20
)
dev.off()


## try to do correlation
## want the sums of these
sums_raphidophyte_otu_table <- rowSums(raphidophyte_otu_table)
sums_proportional_MPL_groupA <- rowSums(proportional_MPL_groupA)
## compare those that are the same
union_vcs_A_and_B <- intersect(names(sums_raphidophyte_otu_table), names(sums_proportional_MPL_groupA))
## need to add dates and sort 
Jericho_A <- as.data.frame(subset(sums_raphidophyte_otu_table, names(sums_raphidophyte_otu_table) %in% union_vcs_A_and_B))
names(Jericho_A) <- "sum"
Jericho_A$Date <- Jericho_data$Date[match(rownames(Jericho_A), as.character(Jericho_data$VC_number))]
Jericho_data$Date[match(names(Jericho_A), as.character(Jericho_data$VC_number))]
# Jericho_A <- as.data.frame(subset(sums_raphidophyte_otu_table, names(sums_raphidophyte_otu_table) %in% union_vcs_A_and_B))
# names(Jericho_A) <- "sum"
# Jericho_A$Date <- Jericho_data$Date[match(rownames(Jericho_A), as.character(Jericho_data$VC_number))]
Jericho_B <- as.data.frame(subset(sums_proportional_MPL_groupA, names(sums_proportional_MPL_groupA) %in% union_vcs_A_and_B))
names(Jericho_B) <- "sum"
Jericho_B$Date <- Jericho_data$Date[match(rownames(Jericho_B), as.character(Jericho_data$VC_number))]

## do with lag
Jericho_A_date_ordered <- Jericho_A[order(Jericho_A$Date),]
Jericho_B_date_ordered <- Jericho_B[order(Jericho_B$Date),]

## straight correlation
cor(Jericho_A_date_ordered$sum, Jericho_B_date_ordered$sum)
## time-lagged correlation
cor(Jericho_A_date_ordered$sum, lag(Jericho_B_date_ordered$sum), use = "pairwise.complete.obs")


pdf(paste0(figures_dir,"group_B_with_Raphidophytes.pdf"), width = 14, height = 14,onefile = FALSE)
grid.arrange(raphidophyte_over_time_matching_rdrp,
             group_B_rdrp,
             ncol=1)
dev.off()

pdf(paste0(figures_dir,"group_E_with_Diatoms.pdf"), width = 14, height = 14,onefile = FALSE)
grid.arrange(diatom_over_time_matching_rdrp,
             group_E_rdrp,
             ncol=1)
dev.off()

pdf(paste0(figures_dir,"group_G_with_Diatoms.pdf"),width = 14, height = 14, onefile = FALSE)
grid.arrange(diatom_over_time_matching_rdrp,
             group_G_rdrp,
             ncol=1)
dev.off()


pdf(paste0(figures_dir,"group_H_with_Diatoms.pdf"),width = 14, height = 14, onefile = FALSE)
grid.arrange(diatom_over_time_matching_rdrp,
             group_H_rdrp,
             ncol=1)
dev.off()





### compare groups to S18 and env to MPL ####
correlation_between_amplicons_and_print_table <- function (group_sums_1,group_sums_1_name, group_sums_2, group_sums_2_name, table_print_name) {
 group_sums_1 <- raphidophyte_otu_table
 group_sums_1_name <- "raph_otus"
 colnames(group_sums_1) <- gsub("OTU",paste0(group_sums_1_name, "_OTU"), colnames(group_sums_1) )
 group_sums_2 <- proportional_MPL_groupA
 colnames(group_sums_2) <- gsub("OTU",paste0(group_sums_2_name, "_OTU"), colnames(group_sums_2) )
 group_sums_1$Date <- Jericho_data$Date[match(rownames(group_sums_1), Jericho_data$VC)]
 ## want to add the S18 matrix to the viral fam one. 
 group_sums_1_with_env <- merge(as.matrix(group_sums_1),as.matrix(params_Jericho), by="Date")
 group_sums_1_with_env$Date <- as.Date(group_sums_1_with_env$Date)
 
 union_vcs <- intersect(rownames(group_sums_2), group_sums_1_with_env$VC)
 
 group_sums_2_jer <- droplevels(subset(group_sums_2,rownames(group_sums_2) %in% union_vcs))
 
 group_sums_1_with_env_jer <- droplevels(subset(group_sums_1_with_env, group_sums_1_with_env$VC%in% union_vcs))
 
 spearman <- rcorr(as.matrix(group_sums_2_jer), as.matrix(group_sums_1_with_env_jer[,-(1)]),type="spearman") 
 correlations <- spearman$r
 p <- spearman$P
 
 melted_cor <- melt(correlations)
 melted_p <- melt(p)
 mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
 melted_mystars <- melt(mystars)
 
 melted_together <- cbind( melted_cor,p_value=melted_p$value, stars=melted_mystars$value)
 #melted_together <- na.omit(melted_together) ## gets rid of the leftover diagonals
# names(melted_together)
 
 melted_together$cor_with_stars <- paste(round(melted_together$value, digits=3), " ", melted_together$stars)
 
 melted_together <- na.omit(melted_together) ## gets rid of the leftover diagonals
 names(melted_together) 
 
 #filter for correlation about 0.6 and p value less than 0.01
 #filtered_data <- filter(melted_together, p_value <= 0.01  & abs(value) > 0.5)
 
 #filtered_data_groups <- melted_together[str_split_fixed(melted_together$Var1, "_", n=2)[,1]=="sum" |str_split_fixed(melted_together$Var1, "_", n=2)[,1]=="sum",]
 
 significant_cor <- dcast(melted_together, Var1 ~Var2, value.var="cor_with_stars")
 write.csv(significant_cor, file=table_print_name)
}

correlation_between_amplicons_and_print_table(raphidophyte_otu_table,"raphidophyte", proportional_MPL_groupA, "mpl_group_A", "../results/RdRp_groupA_with_raphidopytes_cor.csv")


### what about with gp23 group I and cyanobacteria

## gp23

proportional_gp23_groupA <- read.csv("../results/proportional_gp23_groupA.csv", row.names = 1)
proportional_gp23_groupB <- read.csv("../results/proportional_gp23_groupB.csv", row.names = 1)
proportional_gp23_groupC <- read.csv("../results/proportional_gp23_groupC.csv", row.names = 1)
proportional_gp23_groupD <- read.csv("../results/proportional_gp23_groupD.csv", row.names = 1)
proportional_gp23_groupE <- read.csv("../results/proportional_gp23_groupE.csv", row.names = 1)
proportional_gp23_groupF <- read.csv("../results/proportional_gp23_groupF.csv", row.names = 1)
proportional_gp23_groupG <- read.csv("../results/proportional_gp23_groupG.csv", row.names = 1)
proportional_gp23_groupH <- read.csv("../results/proportional_gp23_groupH.csv", row.names = 1)
proportional_gp23_groupI <- read.csv("../results/proportional_gp23_groupI.csv", row.names = 1)
#proportional_gp23_groupJ <- read.csv("../results/proportional_gp23_groupJ.csv", row.names = 1)
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

## Read in OTU table and get long format with Date
long_gp23_groupI_otus <- add_date_to_OTU_table_rtrn_long(proportional_gp23_groupI,
                                                         Jericho_data)

#long_rdrp_groupB_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupB,
#                                                         Jericho_data)
long_cyano_otus <- add_date_to_OTU_table_rtrn_long(cyanobacteria_otu_table,
                                                          Jericho_data)

#                                                         Jericho_data)
long_cyano_otus_phyla <- add_date_to_OTU_table_rtrn_long(cyanobacteria_otu_table_phyla,
                                                   Jericho_data)

#filtered_long_cyano <- long_cyano_otus %>% 
# filter(Date < as.Date("2011-06-01") )

group_I_gp23 <- plot_individual_otus_over_time(long_gp23_groupI_otus, "gp23 group I", "virus")
#group_B_rdrp <- plot_individual_otus_over_time(long_rdrp_groupB_otus, "rdrp group b")

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

pdf(paste0(figures_dir,"group_I_with_Cyanos.pdf"), width = 14, height = 14, onefile = FALSE)
plot_grid(cyano_over_time_matching_gp23,
          group_I_gp23,
          align = "v",
          ncol=1,
          labels = c("A", "B"),
          label_size = 20)
dev.off()



## try to do correlation
## want the sums of these
sums_cyano_otu_table <- rowSums(cyanobacteria_otu_table_phyla)
sums_proportional_gp23_groupI <- rowSums(proportional_gp23_groupI)
## compare those that are the same
union_vcs_A_and_B <- intersect(names(sums_cyano_otu_table), names(sums_proportional_gp23_groupI))
## need to add dates and sort 
Jericho_data$Date[match(names(Jericho_A), as.character(Jericho_data$VC_number))]
Jericho_A <- as.data.frame(subset(sums_cyano_otu_table, names(sums_cyano_otu_table) %in% union_vcs_A_and_B))
names(Jericho_A) <- "sum"
Jericho_A$Date <- Jericho_data$Date[match(rownames(Jericho_A), as.character(Jericho_data$VC_number))]
Jericho_B <- as.data.frame(subset(sums_proportional_gp23_groupI, names(sums_proportional_gp23_groupI) %in% union_vcs_A_and_B))
names(Jericho_B) <- "sum"
Jericho_B$Date <- Jericho_data$Date[match(rownames(Jericho_B), as.character(Jericho_data$VC_number))]

## do with lag
Jericho_A_date_ordered <- Jericho_A[order(Jericho_A$Date),]
Jericho_B_date_ordered <- Jericho_B[order(Jericho_B$Date),]

## straight correlation
cor(Jericho_A_date_ordered$sum, Jericho_B_date_ordered$sum)
## time-lagged correlation
cor(Jericho_A_date_ordered$sum, lag(Jericho_B_date_ordered$sum), use = "pairwise.complete.obs")





pdf(paste0(figures_dir,"group_I_with_Cyanos_phyla.pdf"), width = 14, height = 14, onefile = FALSE)
plot_grid(cyano_over_time_matching_gp23_phyla,
          group_I_gp23,
          align = "v",
          ncol = 1,
          labels = c("A", "B"),
          label_size = 20)
dev.off()

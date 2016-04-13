## Trees showing the OTus that are being oligtypes

## Author: Julia Gustavsen
##Date: 24 June 2013
##Script to re-create the phylogenetic tree with heatmap that I had made with iTOL.
## Will have to line up the heatmap with the phylogenetic tree. I am adding zeros to the OTU table for the isolates and cloned sequences obtained from Genbank.
#Last updated: 6 September 2014


library(ggplot2)
#library(gplots)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(ape)
library(phytools)
library(scales)
library(ggtree)
library(Hmisc)



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

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_MPL_OTUs_top_10 <- read.delim("../data/OTU_table_Jericho_time_series_MPL_normalized_top_10.tsv", row.names=1)
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_18s_top_10 <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_top_10.tsv", row.names = 1)

normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_gp23_OTUs_top_10 <- read.delim("../data/OTU_table_Jericho_time_series_gp23_normalized_top_10.tsv", row.names = 1)
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number")
normalized_16s_top_10 <- read.delim("../data/OTU_table_Jericho_time_series_16s_R1_normalized_top_10.tsv", row.names = 1)


Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)


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

plot_tree_root_and_ladderize <- function (tree_file, root, plot_tree_file_name) {
 tree_amplicon <- read.tree(tree_file)
 tree_amplicon <- root(tree_amplicon, root, resolve.root=TRUE)
 tree_amplicon  <- ladderize(tree_amplicon)
 pdf(plot_tree_file_name,width = 8, height = 11 )
 plot(tree_amplicon,cex=0.3)
 dev.off()
 #exporting tree.
 write.tree(tree_amplicon, "../results/temp_tree.tree")
 #re-import to be able to match up the tip labels with the OTU table. 
 reordered_tree <- read.tree("../results/temp_tree.tree")
 tree_amplicon <- reordered_tree
 return(tree_amplicon)
}

RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                                 "../figures/RdRp_95_miseq_data_with_env_iso_and_ref_Raxml.pdf")

pdf("../figures/RdRp_tree_with_OTUs_for_oligotyping.pdf")
 ggtree(RdRp_95_ref_tree, ladderize=FALSE) + geom_text(aes(label=label), subset=.(RdRp_95_ref_tree$tip.label %in% colnames(normalized_MPL_OTUs_top_10)), hjust=-.2) 
 dev.off()
 
 
gp23_95_ref_tree <- plot_tree_root_and_ladderize("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
                                                    root = "Enterobacteria_phage_T4",
                                                    "../figures/gp23_95_miseq_data_with_env_iso_and_ref_Fasttree.pdf")
 

pdf("../figures/gp23_tree_with_OTUs_for_oligotyping.pdf")
ggtree(gp23_95_ref_tree, ladderize=FALSE) + geom_text(aes(label=label), subset=.(gp23_95_ref_tree$tip.label %in% colnames(normalized_gp23_OTUs_top_10)), hjust=-.2) 
dev.off()


plot_tree_and_ladderize <- function (tree_file, plot_tree_file_name) {
 tree_amplicon <- read.tree(tree_file)
 tree_amplicon  <- ladderize(tree_amplicon)
 #pdf(plot_tree_file_name,width = 8, height = 11 )
 plot(tree_amplicon,cex=0.3)
 #dev.off()
 #exporting tree.
 write.tree(tree_amplicon, "../results/temp_tree.tree")
 #re-import to be able to match up the tip labels with the OTU table. 
 reordered_tree <- read.tree("../results/temp_tree.tree")
 tree_amplicon <- reordered_tree
 return(tree_amplicon)
}
S18_tree <- plot_tree_and_ladderize("../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    "../figures/S18_97_miseq_data_Raxml.pdf")


colnames(normalized_18s_top_10) <- gsub("\\.", "_", colnames(normalized_18s_top_10))
pdf("../figures/18s_tree_with_OTUs_for_oligotyping.pdf")
ggtree(S18_tree, ladderize=FALSE) + geom_text(aes(label=label), subset=.(S18_tree$tip.label %in% colnames(normalized_18s_top_10)), hjust=-.2) 
dev.off()


S16_tree <- plot_tree_and_ladderize("../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    "../figures/S16_97_miseq_data_Raxml.pdf")

colnames(normalized_16s_top_10) <- gsub("\\.", "_", colnames(normalized_16s_top_10))
pdf("../figures/16s_tree_with_OTUs_for_oligotyping.pdf")
ggtree(S16_tree, ladderize=FALSE) + geom_text(aes(label=label), subset=.(S16_tree$tip.label %in% colnames(normalized_16s_top_10)), hjust=-.2) 
dev.off()


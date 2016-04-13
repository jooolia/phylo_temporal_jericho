## Author: Julia Gustavsen
##Date: 24 June 2013
##Script to re-create the phylogenetic tree with heatmap that I had made with iTOL.
## Will have to line up the heatmap with the phylogenetic tree. I am adding zeros to the OTU table for the isolates and cloned sequences obtained from Genbank.
#Last updated: 6 September 2014


library(ggplot2)
#library(gplots)
library(gridExtra)
library(vegan)
library(RColorBrewer)
library(picante)
library(xtable)
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
 figures_dir <- "../figures/"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
  point_colour <- "white"
  heatmap_fill <- "black"
  low_heatmap <- "grey"
  figures_dir <- "../figures_pres/"
 }
}

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number")


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


#### Drawing the phylogenetic tree ####

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

plot_tree <- function(tree, plot_tree_file_name){
 
 pdf(plot_tree_file_name,width = 8, height = 11 )
 
 plot(tree,cex=0.5)
 
 dev.off()
}


add_isolates_to_otu_table <- function (tree, OTU_table) {
 #### Adding tree tips as zeros to the OTU table ####
 #For tree tips not in OTU table I can add a column of zeros to the matrix. 
 #New copy of OTU table
 OTU_table_with_isolates <- OTU_table
 
 #Vector I would add in place for the isolates based on the number of sites in the sample. 
 sample_length <- length(rownames(OTU_table_with_isolates))
 
 #Make vector of 0s based on number of sites. 
 zeros <- rep(0, sample_length)
 
 #Insert vector in table based on number of tips in tree that are not found in the OTU table. 
 for(i in (tree$tip.label)) {
  if (is.element(i, colnames(OTU_table_with_isolates))) {
   print(i)
  }
  else {
   print("no")
   print(i)
   OTU_table_with_isolates <-cbind(OTU_table_with_isolates, zeros)
   colnames(OTU_table_with_isolates)[dim(OTU_table_with_isolates)[2]] <- paste(i)   
  }
 } 
 #### order matrix so that it is ordered like the tree ####
 OTU_table_with_isolates_reordered <- OTU_table_with_isolates[,rev(tree$tip.label)]
 return(OTU_table_with_isolates_reordered)
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



phylo_diversity_plots_and_table <- function (OTU_table_with_isolates,
                                             reordered_tree,
                                             filename_bar_plot_pd,
                                             filename_table_pd) {
 #OTU_table_with_isolates <- MPL_otu_table_with_isolates
 #reordered_tree <- RdRp_95_ref_tree
 RdRp.pd <- pd(OTU_table_with_isolates, reordered_tree, include.root=FALSE)
 RdRp.pd$site <- rownames(RdRp.pd)
 ## Add in dates:
 RdRp.pd$Date <- Jericho_data$Date[match(RdRp.pd$site, Jericho_data$VC_number)]
 
 RdRp.psd <- psd(OTU_table_with_isolates, reordered_tree)
 RdRp.psd$Date <- Jericho_data$Date[match(rownames(RdRp.psd), Jericho_data$VC_number)]
 #dm_comun <- cophenetic(reordered_tree) ## distance matrix of the tree
 #ses.mpd.comun=ses.mpd(OTU_table_with_isolates, dm_comun, runs=100, iterations=100)
 ##some other metrics
 #ses.mntd(OTU_table_with_isolates, reordered_tree, "taxa.labels")
 #mpd(OTU_table_with_isolates, reordered_tree)
 ## phylogenetic beta diversity
 
 OTU_table_with_isolates_subset <- subset(OTU_table_with_isolates, select= (colnames(OTU_table_with_isolates) %in% reordered_tree$tip.label))
 
 phylogenetic_beta <- as.matrix(comdist(OTU_table_with_isolates_subset, cophenetic(reordered_tree)))
 write.csv(phylogenetic_beta , file=paste0(filename_table_pd, "beta_div.csv"))
 
 #ses.mntd(normalized_MPL_OTUs_subset, cophenetic(RdRp_95_only_miseq),null.model="taxa.labels")
 
 # Calculates mean nearest taxon distance for taxa in a community
 tree_otus.mntd <- mntd(OTU_table_with_isolates_subset, cophenetic(reordered_tree))
 names(tree_otus.mntd) <- rownames(OTU_table_with_isolates_subset)
 tree_otus.mntd_melted <- melt(tree_otus.mntd)
 
 tree_otus.mntd_melted$Date <- Jericho_data$Date[match(rownames(tree_otus.mntd_melted), Jericho_data$VC_number)]
 
 ## mean pairwise distance...Calculates mean pairwise distance separating taxa in a community
 tree_otus.mpd <- mpd(OTU_table_with_isolates_subset, cophenetic(reordered_tree))
 names(tree_otus.mpd) <- rownames(OTU_table_with_isolates_subset)
 tree_otus.mpd_melted <- melt(tree_otus.mpd)
 
 tree_otus.mpd_melted$Date <- Jericho_data$Date[match(rownames(tree_otus.mpd_melted), Jericho_data$VC_number)]
 
 
 pdf(filename_bar_plot_pd, width = 11, height = 8, onefile = FALSE)
 pd <- ggplot(RdRp.pd, aes(x=Date, y=PD))+
  #geom_bar(stat="identity", fill="#1efffe")+
  geom_line(colour=line_colour)+
  geom_point(colour="#1efffe")+
  season_line +
  season_text +
  spring_bloom_line+
  date_scaling +
  ggtitle("phylogentic diversity")+
  theme_JAG_presentation()
 print(pd)  
 sr <- ggplot(RdRp.pd, aes(x=Date, y=SR))+
  # geom_bar(stat="identity",fill="#1e20ff")+
  geom_line(colour=line_colour)+
  geom_point(colour="#1e20ff")+
  season_line +
  spring_bloom_line+
  date_scaling +
  ggtitle("species richness")+
  theme_JAG_presentation()
 print(sr)
 mntd <- ggplot(tree_otus.mntd_melted, aes(x=Date, y=value))+
  geom_bar(stat="identity",fill="steelblue3")+
  season_line +
  spring_bloom_line+
  date_scaling +
  ggtitle("mean nearest taxa distance")+
  theme_JAG_presentation()
 print(mntd)
 
 mpd <- ggplot(tree_otus.mpd_melted, aes(x=Date, y=value))+
  geom_bar(stat="identity",fill="turquoise")+
  season_line +
  spring_bloom_line+
  date_scaling +
  ggtitle("mean phylogenetic distance")+
  theme_JAG_presentation()
 print(mpd)
 
 psr <- ggplot(RdRp.psd, aes(x=Date, y=PSR))+
  geom_bar(stat="identity",fill="turquoise4")+
  season_line +
  spring_bloom_line+
  date_scaling +
  ggtitle("phylogentic  species richness")+
  theme_JAG_presentation()
 print(psr)
 
 psc <- ggplot(RdRp.psd, aes(x=Date, y=PSC))+
  geom_bar(stat="identity",fill="aquamarine")+
  season_line +
  spring_bloom_line+
  date_scaling +
  ggtitle("phylogentic species clustering")+
  theme_JAG_presentation()
 print(psc)
 
 dev.off()
 
 RdRp.pd.table <-xtable(RdRp.pd)
 
 print(RdRp.pd.table, 
       floating = FALSE,
       type="html", 
       file=paste0(filename_table_pd, ".html"))
 
 write.csv(RdRp.pd.table, file=paste0(filename_table_pd, ".csv"))
 ## want to return the pd and sr
 return(list(pd,sr))
 
}




### Run on trees!

### MPL ##### 

#### MPL With references ####

RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                                 paste0(figures_dir,"RdRp_95_miseq_data_with_env_iso_and_ref_Raxml.pdf"))



MPL_otu_table_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs)

MPL_table_PD <- phylo_diversity_plots_and_table(MPL_otu_table_with_isolates,
                                                RdRp_95_ref_tree,
                                                paste0(figures_dir,"RdRp_PD_with_OTUs_by_site%03d.pdf"),
                                                "../results/RdRp_phylogenetic_diversity_table")


high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)
normalized_MPL_OTUs_no_high_res <- subset(normalized_MPL_OTUs, !(rownames(normalized_MPL_OTUs) %in% high_res_vcs))

## only high res

normalized_MPL_OTUs_only_high_res <- subset(normalized_MPL_OTUs, rownames(normalized_MPL_OTUs) %in% high_res_vcs)


RdRp_OTU_annotations <- read.csv("../results/RdRp_groups_with_OTUs.csv")


### gp23  ####

#### gp23  ith references ####

gp23_95_ref_tree <- read.tree("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result")

gp23_95_ref_tree <- root( gp23_95_ref_tree, "Enterobacteria_phage_T4", resolve.root=TRUE)

gp23_OTU_annotations <- read.csv("../results/gp23_groups_with_OTUs.csv")


gp23_otu_table_with_isolates <- add_isolates_to_otu_table(gp23_95_ref_tree, normalized_gp23_OTUs)

gp23_table_PD <-  phylo_diversity_plots_and_table(gp23_otu_table_with_isolates,
                                                  gp23_95_ref_tree,
                                                  paste0(figures_dir,"gp23_PD_with_OTUs_by_site.pdf"),
                                                  "../results/gp23_phylogenetic_diversity_table")

normalized_gp23_OTUs_no_high_res <- subset(normalized_gp23_OTUs, !(rownames(normalized_gp23_OTUs) %in% high_res_vcs))

normalized_gp23_OTUs_only_high_res <- subset(normalized_gp23_OTUs, rownames(normalized_gp23_OTUs) %in% high_res_vcs)


### 18s ####

plot_tree_and_ladderize <- function (tree_file, plot_tree_file_name) {
 tree_amplicon <- read.tree(tree_file)
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

S18_tree <- plot_tree_and_ladderize("../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    paste0(figures_dir,"S18_97_miseq_data_Fasttree.pdf"))

S18_tree$tip.label <- gsub("_size_.*$", "",S18_tree$tip.label)
colnames(normalized_18s_OTUs) <- gsub(".size.*.", "",colnames(normalized_18s_OTUs))

taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)



## fix so it matches!!!!

normalized_18s_OTUs_no_high_res <- subset(normalized_18s_OTUs, !(rownames(normalized_18s_OTUs) %in% high_res_vcs))

## only high res

normalized_18s_OTUs_only_high_res <- subset(normalized_18s_OTUs, rownames(normalized_18s_OTUs) %in% high_res_vcs)


S18_table_PD <- phylo_diversity_plots_and_table(normalized_18s_OTUs,
                                                S18_tree, 
                                                paste0(figures_dir,"S18_PD_with_OTUs_by_site%03d.pdf"),
                                                "../results/S18_phylogenetic_diversity_table")


### Top 18s 100 ####

normalized_18s_OTUs_top_100 <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_top_100.tsv",row.names=1)

normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv",                                  row.names=1)


normalized_18s_not_phytos_to_remove <- normalized_18s_OTUs[,!(colnames(normalized_18s_OTUs) %in% colnames(normalized_18s_OTUs_phytos))]
S18_97_ref_tree_phytos <- drop.tip(S18_tree, colnames(normalized_18s_not_phytos_to_remove ))

S18_table_PD_phytos <-  phylo_diversity_plots_and_table(normalized_18s_OTUs_phytos,
                                                        S18_97_ref_tree_phytos,  
                                                        paste0(figures_dir,"S18_PD_phytos_by_site%03d.pdf"),
                                                        "../results/S18_phytos_phylogenetic_diversity_table")

## hetero

normalized_18s_not_hetero_to_remove <- normalized_18s_OTUs[,!(colnames(normalized_18s_OTUs) %in% colnames(normalized_18s_OTUs_hetero))]
S18_97_ref_tree_hetero <- drop.tip(S18_tree, colnames(normalized_18s_not_hetero_to_remove ))


S18_table_PD_hetero <-  phylo_diversity_plots_and_table(normalized_18s_OTUs_hetero,
                                                        S18_97_ref_tree_hetero,  
                                                        paste0(figures_dir,"S18_PD_hetero_by_site%03d.pdf"),
                                                        "../results/S18_hetero_phylogenetic_diversity_table")



#### 16s ####

S16_tree <- plot_tree_and_ladderize("../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    paste0(figures_dir,"S16_97_miseq_data_Fasttree.pdf"))

S16_tree$tip.label <- gsub("_size_.*$", "",S16_tree$tip.label)
colnames(normalized_16s_OTUs) <- gsub(".size.*.", "",colnames(normalized_16s_OTUs))


normalized_16s_OTUs_no_high_res <- subset(normalized_16s_OTUs, !(rownames(normalized_16s_OTUs) %in% high_res_vcs))


## only high res

normalized_16s_OTUs_only_high_res <- subset(normalized_16s_OTUs, rownames(normalized_16s_OTUs) %in% high_res_vcs)



S16_table_PD <- phylo_diversity_plots_and_table(normalized_16s_OTUs,
                                                S16_tree, 
                                                paste0(figures_dir,"S16_PD_with_OTUs_by_site%03d.pdf"),
                                                "../results/S16_phylogenetic_diversity_table")

## plot all together
pdf(paste0(figures_dir,"phylogenetic_diversity_over_time_all_amplicons_no_AVS_line_graphs%03d.pdf"), width = 20, height = 15, onefile = FALSE)
grid.arrange(
 gp23_table_PD[[2]],
 MPL_table_PD[[2]],
 gp23_table_PD[[1]],
 MPL_table_PD[[1]],
 S16_table_PD[[2]],
 S18_table_PD[[2]],
 S16_table_PD[[1]],
 S18_table_PD[[1]],
 S18_table_PD_hetero[[2]],
 S18_table_PD_phytos[[2]],
 S18_table_PD_hetero[[1]],
 S18_table_PD_phytos[[1]],
 ncol=2)
dev.off()
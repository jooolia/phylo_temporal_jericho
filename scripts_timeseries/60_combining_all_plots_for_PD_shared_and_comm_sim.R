#combine the pd, the shared OTUs and the diversity over time


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
library(plyr)
library(dplyr)
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
 heatmap_fill <- "white"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
  point_colour <- "white"
  heatmap_fill <- "black"
 }
}

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number")

normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv", row.names=1)


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
                          size=2.5)
spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=2)
date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")))


#### Drawing the phylogenetic tree ####

plot_tree_root_and_ladderize <- function (tree_file, root, plot_tree_file_name) {
 tree_amplicon <- read.tree(tree_file)
 tree_amplicon <- root(tree_amplicon, root, resolve.root=TRUE)
 tree_amplicon  <- ladderize(tree_amplicon)
 plot(tree_amplicon,cex=0.3)
 #exporting tree.
 write.tree(tree_amplicon, "../results/temp_tree.tree")
 #re-import to be able to match up the tip labels with the OTU table. 
 reordered_tree <- read.tree("../results/temp_tree.tree")
 tree_amplicon <- reordered_tree
 return(tree_amplicon)
}

plot_tree <- function(tree, plot_tree_file_name){
 plot(tree,cex=0.5)
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

make_heatmap_based_on_tree_ordering <- function (OTU_table_with_isolates,
                                                 tree, 
                                                 plot_filename) {
 #### Heatmap with log scaled values ####
 #otu_table$VC_number = rownames(otu_table)
 #long_otus <- melt(otu_table, id="VC_number", variable.name="OTUid")
 ## Add in dates:
 OTU_table_long <- melt(as.matrix(OTU_table_with_isolates))
 OTU_table_long$value[OTU_table_long$value == 0] <- NA
 OTU_table_long$Date <- Jericho_data$Date[match(OTU_table_long$Var1, Jericho_data$VC_number)]
 
 OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(tree$tip.label))
 OTU_heatmap <- ggplot(OTU_table_long, 
                       aes(Date,Var2)) + 
  season_line +
  spring_bloom_line+
  # OTU_heatmap <- ggplot(long_OTUs, aes(x=Date,y=OTUid)) + 
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = heatmap_fill,
                      high = "#40627C",
                      na.value = heatmap_fill ,
                      trans = "log"
  ) +
  theme_JAG_presentation()
 ##fix up formatting
 base_size=9
 
 OTU_heatmap_formatted <- OTU_heatmap + 
  labs(x = "",  y = "") + 
  date_scaling
 #print out
 print(OTU_heatmap_formatted)
 return(OTU_heatmap_formatted)
}

tree_plot_and_then_heatmap_print <- function (tree_file, root, plot_tree_pdf,OTU_table,heatmap_plot_pdf) {
 tree_refs_env <- plot_tree_root_and_ladderize(tree_file, root, plot_tree_pdf)
 OTU_table_with_isolates <- add_isolates_to_otu_table(tree_refs_env,OTU_table )
 make_heatmap_based_on_tree_ordering(OTU_table_with_isolates,
                                     tree_refs_env,
                                     heatmap_plot_pdf)
}

#### Richness across trees ####

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

 OTU_table_with_isolates_subset <- subset(OTU_table_with_isolates, select= (colnames(OTU_table_with_isolates) %in% reordered_tree$tip.label))
 
 phylogenetic_beta <- as.matrix(comdist(OTU_table_with_isolates_subset, cophenetic(reordered_tree)))
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

 pd <- ggplot(RdRp.pd, aes(x=Date, y=PD))+
  season_line +
  spring_bloom_line+
    geom_line(colour=line_colour)+
  geom_point(colour="#1efffe")+
  date_scaling +
  theme_JAG_presentation(base_size=18)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=20),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size=18),
        axis.line=element_line(size = 0.3))+
 # ggtitle("phylogenetic diversity")+
  ylab("Phylogenetic \ndiversity")

 print(pd)  
 sr <- ggplot(RdRp.pd, aes(x=Date, y=SR))+
  # geom_bar(stat="identity",fill="#1e20ff")+
  season_line +
  spring_bloom_line+
    geom_line(colour=line_colour)+
  geom_point(colour="#1e20ff")+
  date_scaling +
  theme_JAG_presentation()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=20),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size=18),
        axis.line=element_line(size = 0.3))+
  ylab("Species \nrichness")
  #ggtitle("species richness")

 print(sr)
 mntd <- ggplot(tree_otus.mntd_melted, aes(x=Date, y=value))+
  season_line +
  spring_bloom_line+
    geom_bar(stat="identity",fill="steelblue3")+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.line=element_line(size = 0.3))+
  date_scaling +
  ggtitle("mean nearest taxa distance")+
  theme_JAG_presentation()
 print(mntd)
 
 mpd <- ggplot(tree_otus.mpd_melted, aes(x=Date, y=value))+
  season_line +
  spring_bloom_line+
    geom_bar(stat="identity",fill="turquoise")+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.line=element_line(size = 0.3))+
  date_scaling +
  ggtitle("mean phylogenetic distance")+
  theme_JAG_presentation()
 print(mpd)
 
 psr <- ggplot(RdRp.psd, aes(x=Date, y=PSR))+
  season_line +
  spring_bloom_line+
    geom_bar(stat="identity",fill="turquoise4")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=20),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size=18),
        axis.line=element_line(size = 0.3))+
  date_scaling +
  ggtitle("phylogenetic  species richness")+
 theme_JAG_presentation(base_size=18)+
  ylab("Species richness")
 print(psr)
 
 psc <- ggplot(RdRp.psd, aes(x=Date, y=PSC))+
  season_line +
  spring_bloom_line+
    geom_bar(stat="identity",fill="aquamarine")+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.line=element_line(size = 0.3))+
  date_scaling +
  ggtitle("phylogenetic species clustering")+
  theme_JAG_presentation()
 print(psc)

 
 RdRp.pd.table <-xtable(RdRp.pd)

 ## want to return the pd and sr
 return(list(pd,sr))
 
}

### Run on trees!

### MPL ##### 

#### MPL With references ####

RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                                 "../figures/RdRp_95_miseq_data_with_env_iso_and_ref_Raxml.pdf")



MPL_otu_table_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs)

MPL_table_PD <- phylo_diversity_plots_and_table(MPL_otu_table_with_isolates,
                                                RdRp_95_ref_tree,
                                                "../figures/RdRp_PD_with_OTUs_by_site%03d.pdf",
                                                "../results/RdRp_phylogenetic_diversity_table")


high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)
normalized_MPL_OTUs_no_high_res <- subset(normalized_MPL_OTUs, !(rownames(normalized_MPL_OTUs) %in% high_res_vcs))

### gp23  ####

#### gp23  ith references ####

gp23_95_ref_tree <- read.tree("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result")

gp23_95_ref_tree <- root( gp23_95_ref_tree, "Enterobacteria_phage_T4", resolve.root=TRUE)

## having problems with memory for this tree...

gp23_otu_table_with_isolates <- add_isolates_to_otu_table(gp23_95_ref_tree, normalized_gp23_OTUs)

gp23_table_PD <-  phylo_diversity_plots_and_table(gp23_otu_table_with_isolates,
                                                  gp23_95_ref_tree,
                                                  "../figures/gp23_PD_with_OTUs_by_site.pdf",
                                                  "../results/gp23_phylogenetic_diversity_table")

normalized_gp23_OTUs_no_high_res <- subset(normalized_gp23_OTUs, !(rownames(normalized_gp23_OTUs) %in% high_res_vcs))

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
                                    "../figures/S18_97_miseq_data_Fasttree.pdf")

S18_tree$tip.label <- gsub("_size_.*$", "",S18_tree$tip.label)
colnames(normalized_18s_OTUs) <- gsub(".size.*.", "",colnames(normalized_18s_OTUs))

taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)

## fix so it matches!!!!
high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)
normalized_18s_OTUs_no_high_res <- subset(normalized_18s_OTUs, !(rownames(normalized_18s_OTUs) %in% high_res_vcs))

## only high res

normalized_18s_OTUs_only_high_res <- subset(normalized_18s_OTUs, rownames(normalized_18s_OTUs) %in% high_res_vcs)


S18_table_PD <- phylo_diversity_plots_and_table(normalized_18s_OTUs,
                                                S18_tree, 
                                                "../figures/S18_PD_with_OTUs_by_site%03d.pdf",
                                                "../results/S18_phylogenetic_diversity_table")


### Top 18s 100 ####

normalized_18s_OTUs_top_100 <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_top_100.tsv",row.names=1)

normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv",                                  row.names=1)
 
normalized_18s_not_phytos_to_remove <- normalized_18s_OTUs[,!(colnames(normalized_18s_OTUs) %in% colnames(normalized_18s_OTUs_phytos))]
S18_97_ref_tree_phytos <- drop.tip(S18_tree, colnames(normalized_18s_not_phytos_to_remove ))


S18_table_PD_phytos <-  phylo_diversity_plots_and_table(normalized_18s_OTUs_phytos,
                                                        S18_97_ref_tree_phytos,  
                                                        "../figures/S18_PD_phytos_by_site%03d.pdf",
                                                        "../results/S18_phytos_phylogenetic_diversity_table")

## hetero

normalized_18s_not_hetero_to_remove <- normalized_18s_OTUs[,!(colnames(normalized_18s_OTUs) %in% colnames(normalized_18s_OTUs_hetero))]
S18_97_ref_tree_hetero <- drop.tip(S18_tree, colnames(normalized_18s_not_hetero_to_remove ))

plot_tree(S18_97_ref_tree_hetero,                                                         "../figures/S18_97_hetero_miseq_data_Fasttree.pdf")

S18_table_PD_hetero <-  phylo_diversity_plots_and_table(normalized_18s_OTUs_hetero,
                                                        S18_97_ref_tree_hetero,  
                                                        "../figures/S18_PD_hetero_by_site%03d.pdf",
                                                        "../results/S18_hetero_phylogenetic_diversity_table")



#### 16s ####

S16_tree <- plot_tree_and_ladderize("../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    "../figures/S16_97_miseq_data_Fasttree.pdf")

S16_tree$tip.label <- gsub("_size_.*$", "",S16_tree$tip.label)
colnames(normalized_16s_OTUs) <- gsub(".size.*.", "",colnames(normalized_16s_OTUs))

taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)


normalized_16s_OTUs_no_high_res <- subset(normalized_16s_OTUs, !(rownames(normalized_16s_OTUs) %in% high_res_vcs))

S16_table_PD <- phylo_diversity_plots_and_table(normalized_16s_OTUs,
                                                S16_tree, 
                                                "../figures/S16_PD_with_OTUs_by_site%03d.pdf",
                                                "../results/S16_phylogenetic_diversity_table")



##
## make distance matrices over time and perform regression


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


## ok maybe easier to change to dates before hand
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

Env_data_for_merging <- Jericho_data [,c("Date","VC_number")]
row.names(Env_data_for_merging) <- Env_data_for_merging$VC_number

## will just merge OTUs with date. GOod for plotting then. 
get_community_similarity_plot_over_time <- function (normalized_otus) {
 new_table <- cbind(normalized_otus, Env_data_for_merging[, "Date"][match(rownames(normalized_otus), rownames(Env_data_for_merging))])
 row.names(new_table) <- new_table[,(dim(new_table)[2])]
 new_table <- new_table[,-(dim(new_table)[2])]
 
 ## what about community dissimilarity between times?
 ## want to look at the diss between following times. Time A to Time B
 OTU_diss <- vegdist(new_table)
 OTU_matrix <- as.matrix(OTU_diss)
 
 sorted_by_date_OTU_matrix <- OTU_matrix[sort(row.names(OTU_matrix)),sort(colnames(OTU_matrix))]
 
 community_over_time <- c("Date", "Community_diss_between_dates")
 
 ## or just a loop through all the dates using the indexing to get at them
 for (sample_date in row.names(sorted_by_date_OTU_matrix)){
  print(sample_date)
  index_time_a <- (which(row.names(sorted_by_date_OTU_matrix)==sample_date))
  ## use the next row in the sorted matrix!
  index_time_b <- index_time_a - 1
  print(index_time_b)
  ## at time 2 how does it compare to the last time
  print(sorted_by_date_OTU_matrix[index_time_a, index_time_b])
  comm_diss_between_dates <- sorted_by_date_OTU_matrix[index_time_a, index_time_b]
  row_add <- c(sample_date, comm_diss_between_dates)
  community_over_time <- rbind(community_over_time, row_add)
 }
 
 ## don't want the first few because the first is an artifact and the other is the header names
 colnames(community_over_time) <- community_over_time[1,]
 
 community_over_time <- community_over_time[-c(1:2),1:2]
 community_over_time <- as.data.frame(community_over_time)
 str(community_over_time)
 community_over_time$Date <- as.Date(community_over_time$Date)
 community_over_time$Community_diss_between_dates  <- as.numeric(as.character(community_over_time$Community_diss_between_dates))
 community_over_time$Community_sim_between_dates <- 1-community_over_time$Community_diss_between_dates
 community_over_time$group <- "all"
 str(community_over_time)
 
 p_all <- ggplot(community_over_time, aes(x=Date, y=Community_sim_between_dates, group=group))+ 
  season_line +
  spring_bloom_line+
  geom_line(colour=line_colour)+
  geom_point(colour="grey")+
  theme_JAG_presentation(base_size=18)+
  theme(axis.title = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.line=element_line(size = 0.3))+
  ggtitle(NULL)+ 
  ylab("Community similarity \n between dates")+
  date_scaling+
  ylim(0,1)
 
 return(p_all)
 
}


MPL_div <- get_community_similarity_plot_over_time(normalized_MPL_OTUs)
AVS_div <- get_community_similarity_plot_over_time(normalized_AVS_OTUs)
gp23_div <- get_community_similarity_plot_over_time(normalized_gp23_OTUs)
S18_div <- get_community_similarity_plot_over_time(normalized_18s_OTUs)

S16_div <- get_community_similarity_plot_over_time(normalized_16s_OTUs)

## next 
## shared OTUs over time

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)
## Summarise 

count_16s_otus <- adply(t(normalized_16s_OTUs), 2,function(x)sum(x>0))
count_18s_otus <-adply(t(normalized_18s_OTUs), 2, function(x)sum(x>0))
count_gp23_otus <-adply(t(normalized_gp23_OTUs), 2, function(x)sum(x>0))
count_MPL_otus <-adply(t(normalized_MPL_OTUs), 2, function(x)sum(x>0))

get_percent_shared_over_time <- function (normalized_otus) {
 #normalized_otus <- normalized_MPL_OTUs
 new_table <- cbind(normalized_otus, Env_data_for_merging[, "Date"][match(rownames(normalized_otus), rownames(Env_data_for_merging))])
 row.names(new_table) <- new_table[,(dim(new_table)[2])]
 new_table <- new_table[,-(dim(new_table)[2])]
 
 shared_16s <- betadiver(new_table, triangular = FALSE)
 #this is what I wanted....how many are shared between pools!
 #so now I need to sum and then average over time....need some sort aaply here....eep.
 
 each_time_16s <- as.matrix(shared_16s$a + shared_16s$b +shared_16s$c)

 percent_shared <- (as.matrix(shared_16s$a)/each_time_16s)*100
 #shared_16s$a #these are the shared OTUs between each pool

 
 ## what about community dissimilarity between times?
 ## want to look at the diss between following times. Time A to Time B
# OTU_diss <- vegdist(new_table)
 OTU_matrix <- as.matrix(percent_shared)
 
 sorted_by_date_OTU_matrix <- OTU_matrix[sort(row.names(OTU_matrix)),sort(colnames(OTU_matrix))]
 
 community_over_time <- c("Date", "Shared_between_dates")
 
 ## or just a loop through all the dates using the indexing to get at them
 for (sample_date in row.names(sorted_by_date_OTU_matrix)){
  print(sample_date)
  index_time_a <- (which(row.names(sorted_by_date_OTU_matrix)==sample_date))
  ## use the next row in the sorted matrix!
  index_time_b <- index_time_a - 1
  print(index_time_b)
  ## at time 2 how does it compare to the last time
  print(sorted_by_date_OTU_matrix[index_time_a, index_time_b])
  comm_percent_between_dates <- sorted_by_date_OTU_matrix[index_time_a, index_time_b]
  row_add <- c(sample_date, comm_percent_between_dates)
  community_over_time <- rbind(community_over_time, row_add)
 }
 
 ## don't want the first few because the first is an artifact and the other is the header names
 colnames(community_over_time) <- community_over_time[1,]
 
 community_over_time <- community_over_time[-c(1:2),1:2]
 community_over_time <- as.data.frame(community_over_time)
 str(community_over_time)
 community_over_time$Date <- as.Date(community_over_time$Date)
 community_over_time$Shared_between_dates  <- as.numeric(as.character(community_over_time$Shared_between_dates))
 community_over_time$group <- "all"
 str(community_over_time)
 
 
 p_all <- ggplot(community_over_time, aes(x=Date, y=Shared_between_dates, group=group))+ 
  season_line +
  spring_bloom_line+
  geom_line(colour=line_colour)+
  geom_point(colour="green")+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=18),
        axis.line=element_line(size = 0.3))+
  ggtitle(NULL)+ 
  date_scaling+ylim(c(0,100))
 return(p_all)
 
}


MPL_shared <- get_percent_shared_over_time(normalized_MPL_OTUs)
AVS_shared <- get_percent_shared_over_time(normalized_AVS_OTUs)
gp23_shared <- get_percent_shared_over_time(normalized_gp23_OTUs)
S18_shared <- get_percent_shared_over_time(normalized_18s_OTUs)
S16_shared <- get_percent_shared_over_time(normalized_16s_OTUs)


## plot all together
pdf("../figures/PD_Shared_OTUS_comm_sim_over_time_all_amplicons_no_AVS_line_graphs%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(
 gp23_table_PD[[2]],
 MPL_table_PD[[2]],
 gp23_table_PD[[1]],
 MPL_table_PD[[1]],
 gp23_div, 
 MPL_div, 
# gp23_shared,
# MPL_shared,
 S16_table_PD[[2]],
 S18_table_PD[[2]],
 S16_table_PD[[1]],
 S18_table_PD[[1]],
 S16_div,
 S18_div,
 #S16_shared,
 #S18_shared, 
 ncol=2)
dev.off()

pdf("../figures/PD_Shared_OTUS_comm_sim_over_time_all_amplicons_no_AVS_line_graphs_labelled%03d.pdf", width = 20, height = 17, onefile=FALSE)
plot_grid(
 gp23_table_PD[[2]],
 MPL_table_PD[[2]],
 gp23_table_PD[[1]],
 MPL_table_PD[[1]],
 gp23_div, 
 MPL_div, 
 # gp23_shared,
 # MPL_shared,
 S16_table_PD[[2]],
 S18_table_PD[[2]],
 S16_table_PD[[1]],
 S18_table_PD[[1]],
 S16_div,
 S18_div,
 #S16_shared,
 #S18_shared, 
 ncol = 2,
 labels = c("A", "B","","","","", "C", "D", "", "" , "", ""),
 label_size=25,
 align = "v")
dev.off()


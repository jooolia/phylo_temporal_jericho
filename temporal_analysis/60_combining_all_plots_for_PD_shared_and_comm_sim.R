## combine the pd, the shared OTUs and the diversity over time

## Author: Julia Gustavsen

library(ggplot2)
library(vegan)
library(picante) ## ape loaded by picante
library(reshape2)
library(scales)
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

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv",
                                  row.names="VC_number")
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv",
                                  row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",
                                   row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number")

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv",
                         row.names=1)
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

plot_tree_root_and_ladderize <- function (tree_file,
                                          root,
                                          plot_tree_file_name) {
  tree_amplicon <- read.tree(tree_file)
  tree_amplicon <- root(tree_amplicon,
                        root,
                        resolve.root=TRUE)
  tree_amplicon  <- ladderize(tree_amplicon)
  plot(tree_amplicon,cex=0.3)
  #exporting tree.
  write.tree(tree_amplicon,
             "../results/temp_tree.tree")
  #re-import to be able to match up the tip labels with the OTU table. 
  reordered_tree <- read.tree("../results/temp_tree.tree")
  return(reordered_tree)
}

add_isolates_to_otu_table <- function (tree, OTU_table) {
  #### Adding tree tips as zeros to the OTU table ####
  ## For tree tips not in OTU table I can add a column of zeros to the matrix. 
  ## New copy of OTU table
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
      OTU_table_with_isolates <-cbind(OTU_table_with_isolates,
                                      zeros)
      colnames(OTU_table_with_isolates)[dim(OTU_table_with_isolates)[2]] <- paste(i)   
    }
  } 
  #### order matrix so that it is ordered like the tree ####
  OTU_table_with_isolates_reordered <- OTU_table_with_isolates[,rev(tree$tip.label)]
  return(OTU_table_with_isolates_reordered)
}

#### Richness across trees ####

phylo_diversity_plots_and_table <- function (OTU_table_with_isolates,
                                             reordered_tree,
                                             filename_bar_plot_pd,
                                             filename_table_pd) {
  #OTU_table_with_isolates <- MPL_otu_table_with_isolates
  #reordered_tree <- RdRp_95_ref_tree
  RdRp.pd <- pd(OTU_table_with_isolates,
                reordered_tree,
                include.root = FALSE)
  RdRp.pd$site <- rownames(RdRp.pd)
  ## Add in dates:
  RdRp.pd$Date <- Jericho_data$Date[match(RdRp.pd$site,
                                          Jericho_data$VC_number)]
  
  OTU_table_with_isolates_subset <- subset(OTU_table_with_isolates,
                                           select = (colnames(OTU_table_with_isolates) %in% reordered_tree$tip.label))
  
  
  pd <- ggplot(RdRp.pd,
               aes(x = Date,
                   y = PD))+
    season_line+
    # spring_bloom_line+
    geom_line(colour = line_colour)+
    geom_point(colour = "#1efffe")+
    date_scaling +
    theme_JAG_presentation(base_size = 18)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold",
                                      size = 20),
          axis.text.x  = element_blank(),
          axis.text.y = element_text(size = 18),
          axis.line=element_line(size = 0.3))+
    ylab("Phylogenetic \ndiversity")
  
  print(pd) 
  
  sr <- ggplot(RdRp.pd,
               aes(x = Date,
                   y = SR))+
    season_line +
    # spring_bloom_line+
    geom_line(colour = line_colour)+
    geom_point(colour = "#1e20ff")+
    date_scaling +
    theme_JAG_presentation()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold",
                                      size = 20),
          axis.text.x  = element_blank(),
          axis.text.y = element_text(size = 18),
          axis.line = element_line(size = 0.3))+
    ylab("Species \nrichness")
  
  print(sr)
  
  ## want to return the pd and sr
  return(list(pd,sr))
}

### Run on trees!

### MPL ##### 

#### MPL With references ####

RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                                 "../figures/RdRp_95_miseq_data_with_env_iso_and_ref_Raxml.pdf")


MPL_otu_table_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree,
                                                         normalized_MPL_OTUs)

MPL_table_PD <- phylo_diversity_plots_and_table(MPL_otu_table_with_isolates,
                                                RdRp_95_ref_tree,
                                                "../figures/RdRp_PD_with_OTUs_by_site%03d.pdf",
                                                "../results/RdRp_phylogenetic_diversity_table")


high_res_vcs <- c(1198,
                  1199,
                  1200,
                  1201,
                  1202)
normalized_MPL_OTUs_no_high_res <- subset(normalized_MPL_OTUs,
                                          !(rownames(normalized_MPL_OTUs) %in% high_res_vcs))

### gp23  ####

#### gp23  ith references ####

gp23_95_ref_tree <- read.tree("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result")

gp23_95_ref_tree <- root( gp23_95_ref_tree,
                          "Enterobacteria_phage_T4",
                          resolve.root = TRUE)

## having problems with memory for this tree...

gp23_otu_table_with_isolates <- add_isolates_to_otu_table(gp23_95_ref_tree,
                                                          normalized_gp23_OTUs)

gp23_table_PD <-  phylo_diversity_plots_and_table(gp23_otu_table_with_isolates,
                                                  gp23_95_ref_tree,
                                                  "../figures/gp23_PD_with_OTUs_by_site.pdf",
                                                  "../results/gp23_phylogenetic_diversity_table")

normalized_gp23_OTUs_no_high_res <- subset(normalized_gp23_OTUs,
                                           !(rownames(normalized_gp23_OTUs) %in% high_res_vcs))

### 18s ####

plot_tree_and_ladderize <- function (tree_file,
                                     plot_tree_file_name) {
  tree_amplicon <- read.tree(tree_file)
  tree_amplicon <- ladderize(tree_amplicon)
  pdf(plot_tree_file_name,
      width = 8,
      height = 11 )
  plot(tree_amplicon,
       cex = 0.3)
  dev.off()
  
  write.tree(tree_amplicon,
             "../results/temp_tree.tree")
  ## re-import to be able to match up the tip labels
  ## with the OTU table. 
  reordered_tree <- read.tree("../results/temp_tree.tree")
  return(reordered_tree)
}

S18_tree <- plot_tree_and_ladderize("../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    "../figures/S18_97_miseq_data_Fasttree.pdf")

S18_tree$tip.label <- gsub("_size_.*$",
                           "",
                           S18_tree$tip.label)

colnames(normalized_18s_OTUs) <- gsub(".size.*.",
                                      "",
                                      colnames(normalized_18s_OTUs))

normalized_18s_OTUs_no_high_res <- subset(normalized_18s_OTUs,
                                          !(rownames(normalized_18s_OTUs) %in% high_res_vcs))

## only high res

normalized_18s_OTUs_only_high_res <- subset(normalized_18s_OTUs,
                                            rownames(normalized_18s_OTUs) %in% high_res_vcs)


S18_table_PD <- phylo_diversity_plots_and_table(normalized_18s_OTUs,
                                                S18_tree, 
                                                "../figures/S18_PD_with_OTUs_by_site%03d.pdf",
                                                "../results/S18_phylogenetic_diversity_table")


#### 16s ####

S16_tree <- plot_tree_and_ladderize("../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    "../figures/S16_97_miseq_data_Fasttree.pdf")

S16_tree$tip.label <- gsub("_size_.*$",
                           "",
                           S16_tree$tip.label)

colnames(normalized_16s_OTUs) <- gsub(".size.*.",
                                      "",
                                      colnames(normalized_16s_OTUs))

normalized_16s_OTUs_no_high_res <- subset(normalized_16s_OTUs,
                                          !(rownames(normalized_16s_OTUs) %in% high_res_vcs))

S16_table_PD <- phylo_diversity_plots_and_table(normalized_16s_OTUs,
                                                S16_tree, 
                                                "../figures/S16_PD_with_OTUs_by_site%03d.pdf",
                                                "../results/S16_phylogenetic_diversity_table")


## make distance matrices over time and perform regression

Env_data_for_merging <- Jericho_data [,c("Date","VC_number")]
row.names(Env_data_for_merging) <- Env_data_for_merging$VC_number

## will just merge OTUs with date. GOod for plotting then. 
get_community_similarity_plot_over_time <- function (normalized_otus) {
  new_table <- cbind(normalized_otus,
                     Env_data_for_merging[, "Date"][match(rownames(normalized_otus),
                                                          rownames(Env_data_for_merging))])
  
  row.names(new_table) <- new_table[,(dim(new_table)[2])]
  new_table <- new_table[,-(dim(new_table)[2])]
  
  ## want to look at the diss between following times. Time A to Time B
  OTU_diss <- vegdist(new_table)
  OTU_matrix <- as.matrix(OTU_diss)
  
  sorted_by_date_OTU_matrix <- OTU_matrix[sort(row.names(OTU_matrix)),
                                          sort(colnames(OTU_matrix))]
  
  community_over_time <- c("Date",
                           "Community_diss_between_dates")
  
  ## or just a loop through all the dates using the indexing to get at them
  for (sample_date in row.names(sorted_by_date_OTU_matrix)){
    print(sample_date)
    index_time_a <- (which(row.names(sorted_by_date_OTU_matrix) == sample_date))
    ## use the next row in the sorted matrix!
    index_time_b <- index_time_a - 1
    print(index_time_b)
    ## at time 2 how does it compare to the last time
    print(sorted_by_date_OTU_matrix[index_time_a,
                                    index_time_b])
    
    comm_diss_between_dates <- sorted_by_date_OTU_matrix[index_time_a,
                                                         index_time_b]
    row_add <- c(sample_date,
                 comm_diss_between_dates)
    
    community_over_time <- rbind(community_over_time,
                                 row_add)
  }
  
  ## don't want the first few because the first is an artifact and the other is the header names
  colnames(community_over_time) <- community_over_time[1,]
  
  community_over_time <- community_over_time[-c(1:2),
                                             1:2]
  
  community_over_time <- as.data.frame(community_over_time)
  community_over_time$Date <- as.Date(community_over_time$Date)
  community_over_time$Community_diss_between_dates  <- as.numeric(as.character(community_over_time$Community_diss_between_dates))
  community_over_time$Community_sim_between_dates <- 1 - community_over_time$Community_diss_between_dates
  community_over_time$group <- "all"
  str(community_over_time)
  
  p_all <- ggplot(community_over_time,
                  aes(x = Date, 
                      y = Community_sim_between_dates,
                      group = group))+ 
    season_line +
    # spring_bloom_line+
    geom_line(colour = line_colour)+
    geom_point(colour = "grey")+
    theme_JAG_presentation(base_size = 18)+
    theme(axis.title = element_text(face = "bold",
                                    size = 20),
          axis.text.x  = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.line=element_line(size = 0.3))+
    ggtitle(NULL)+ 
    ylab("Community similarity \n between dates\n")+
    xlab("\nDate")+
    date_scaling+
    ylim(0,1)
  
  return(p_all)
  
}


MPL_div <- get_community_similarity_plot_over_time(normalized_MPL_OTUs)
gp23_div <- get_community_similarity_plot_over_time(normalized_gp23_OTUs)
S18_div <- get_community_similarity_plot_over_time(normalized_18s_OTUs)
S16_div <- get_community_similarity_plot_over_time(normalized_16s_OTUs)


pdf("../figures/PD_Shared_OTUS_comm_sim_over_time_all_amplicons_line_graphs_labelled%03d.pdf",
    width = 20,
    height = 16)
print(
plot_grid(
  gp23_table_PD[[2]],
  MPL_table_PD[[2]],
  gp23_table_PD[[1]],
  MPL_table_PD[[1]],
  gp23_div+xlab(""), 
  MPL_div+xlab(""), 
  S16_table_PD[[2]],
  S18_table_PD[[2]],
  S16_table_PD[[1]],
  S18_table_PD[[1]],
  S16_div,
  S18_div,
  ncol = 2,
  labels = c("A", "B","","","","", "C", "D", "", "" , "", ""),
  label_size=25,
  align = "v")
)
dev.off()

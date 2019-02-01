## Author: Julia Gustavsen
##Date: 24 June 2013
##Script to re-create the phylogenetic tree with bubbleplot that I had made with iTOL.
## Will have to line up the bubbleplot with the phylogenetic tree. I am adding zeros to the OTU table for the isolates and cloned sequences obtained from Genbank.
#Last updated: 6 September 2014


library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(ape)
library(phytools)
library(scales)
library(ggtree)
library(Hmisc)
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
  bubbleplot_fill <- "white"
  low_bubbleplot <- "yellow"
  figures_dir <- "../figures/"
} else {
  print("Cool you passed a nice theme file to this script")
  source(inputFile)
  if (inputFile == "../../JAG_black_presentation.R"){
    path_colour <- "white"
    line_colour <- "white"
    point_colour <- "white"
    bubbleplot_fill <- "black"
    low_bubbleplot <- "grey"
    figures_dir <- "../figures_pres/"
  }
}

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number")


## first plot the different stuff
high_res_vcs <- c(1198,
                  1199,
                  1200,
                  1201,
                  1202)
Jericho_data <- subset(Jericho_data,
                       !(Jericho_data$VC_number %in% high_res_vcs))

Jericho_no_high_res_vcs <- Jericho_data$VC_number

missing_MPL_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_MPL_OTUs))]
missing_MPL_dates <- Jericho_data$Date[match(missing_MPL_samples,
                                             Jericho_data$VC_number)]

missing_gp23_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_gp23_OTUs))]
missing_gp23_dates <- Jericho_data$Date[match(missing_gp23_samples,
                                              Jericho_data$VC_number)]

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

season_text <- annotate("text",
                        x=c(mid_summer_x,
                            mid_fall_x, 
                            mid_winter_x,
                            mid_spring_x),
                        y=1.05,
                        label=c("Summer",
                                "Fall",
                                "Winter", 
                                "Spring"))


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

bubbleplot_by_tree_ordering <- function (OTU_table_with_isolates,
                                                 tree, 
                                                 plot_filename) {
  #### bubbleplot with log scaled values ####
  #otu_table$VC_number = rownames(otu_table)
  #long_otus <- melt(otu_table, id="VC_number", variable.name="OTUid")
  ## Add in dates:
  OTU_table_long <- melt(as.matrix(OTU_table_with_isolates))
  OTU_table_long$value[OTU_table_long$value == 0] <- NA
  OTU_table_long$Date <- Jericho_data$Date[match(OTU_table_long$Var1, Jericho_data$VC_number)]
  
  OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(tree$tip.label))
  str(OTU_table_long$value)
  ## want it to be something like the annotations or taxonomic annotations. 
  OTU_bubbleplot_formatted <- ggplot(OTU_table_long, 
                                  aes(Date,Var2)) + 
    season_line +
    geom_point(aes(size = value)) +
    # scale_fill_gradient(name = "Relative\n abundance",
    #                     low = "lightblue",
    #                     high = "#40627C",
    #                     na.value = bubbleplot_fill ,
    #                     trans = "log",
    #                     breaks=c(min(OTU_table_long$value, na.rm = TRUE), 
    #                              sqrt(max(OTU_table_long$value, na.rm= TRUE)),
    #                              max(OTU_table_long$value, na.rm = TRUE)),
    #                     labels=c(min(OTU_table_long$value, na.rm = TRUE),
    #                              round(sqrt(max(OTU_table_long$value, na.rm= TRUE)),digits = -1),
    #                              round(max(OTU_table_long$value, na.rm = TRUE), digits = -2)
    #                     )) +
    labs(y = "") + 
    date_scaling +
    theme_JAG_presentation(base_size = 18)
  
  pdf(plot_filename, onefile=FALSE, width = 20, height = 11)
  print(OTU_bubbleplot_formatted)
  dev.off()
  return(OTU_bubbleplot_formatted)
}

bubbleplot_by_tree_ordering_high_res_view <- function (OTU_table_with_isolates,
                                                               tree, 
                                                               plot_filename) {
  OTU_table_long <- melt(as.matrix(OTU_table_with_isolates))
  OTU_table_long$value[OTU_table_long$value == 0] <- NA
  OTU_table_long$Date <- Jericho_data$Date[match(OTU_table_long$Var1, Jericho_data$VC_number)]
  OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(tree$tip.label))
  
  ## want it to be something like the annotations or taxonomic annotations. 
  OTU_bubbleplot <- ggplot(OTU_table_long, 
                        aes(Date,Var2)) + 
    season_line +
    spring_bloom_line+
    geom_point(aes(size = value)) +
    # scale_fill_gradient(name = "Relative\n abundance",
    #                     low = "lightblue",
    #                     high = "#40627C",
    #                     na.value = bubbleplot_fill ,
    #                     trans = "log",
    #                     breaks=c(min(OTU_table_long$value, na.rm = TRUE), 
    #                              sqrt(max(OTU_table_long$value, na.rm= TRUE)),
    #                              max(OTU_table_long$value, na.rm = TRUE)),
    #                     labels=c(min(OTU_table_long$value, na.rm = TRUE),
    #                              round(sqrt(max(OTU_table_long$value, na.rm= TRUE)),digits = -1),
    #                              round(max(OTU_table_long$value, na.rm = TRUE), digits = -2)
    #                     ))+
    theme_JAG_presentation()
  ##fix up formatting
  base_size=9
  
  OTU_bubbleplot_formatted <- OTU_bubbleplot + 
    labs(x = "",  y = "") + scale_x_date(breaks = date_breaks("week"),
                                         limits = c(as.Date("2011-01-15"),
                                                    as.Date("2011-02-15")))
  #print out
  
  pdf(plot_filename, onefile=FALSE, width = 20, height = 11)
  print(OTU_bubbleplot_formatted)
  dev.off()
  return(OTU_bubbleplot_formatted)
}


tree_plot_and_then_bubbleplot_print <- function (tree_file, root, plot_tree_pdf,OTU_table,bubbleplot_plot_pdf) {
  tree_refs_env <- plot_tree_root_and_ladderize(tree_file, root, plot_tree_pdf)
  OTU_table_with_isolates <- add_isolates_to_otu_table(tree_refs_env,OTU_table )
  bubbleplot_by_tree_ordering(OTU_table_with_isolates,
                                      tree_refs_env,
                                      bubbleplot_plot_pdf)
}


make_label_for_bubbleplot_based_on_tree_ordering_with_annotations <- function (OTU_table_with_isolates,
                                                                            tree, 
                                                                            plot_filename, 
                                                                            annotations,
                                                                            number_of_groups_in_annotation,
                                                                            colour_palette) {
  
  #Insert vector in table based on number of tips in tree that are not found in the OTU table. 
  for(i in colnames((OTU_table_with_isolates))) {
    if (is.element(i, annotations$values)) {
      print(i)
    }
    else {
      print("no")
      print(i)
      new_row <- data.frame(X="X", group="groupNull",values= i)
      annotations <-rbind(annotations,new_row)
    }
  } 
  #### order matrix so that it is ordered like the tree ####
  annotations_reordered <- annotations[order(match(tree$tip.label,annotations$values)),]
  annotations_reordered$values <- factor(annotations_reordered$values, levels = factor(tree$tip.label))
  
  make_hcb <- function(data, var, name = NULL, fillScale = NULL, ...) {
    ## from Suann Holmes
    ## http://statweb.stanford.edu/~susan/papers/Pregnancy/PNAS_Vaginal_Analysis.html
    #data <- annotations_reordered
    hcb <- ggplot(data=data, aes_string(x=1, y="values", fill="group")) + 
      geom_raster()+
      xlab(NULL) + ylab(NULL) +
      theme(axis.title=element_blank(), axis.ticks=element_blank()) +
      theme(axis.text.x=element_blank())+
      theme(axis.text.y=element_blank())
    #+   annotate("text", x=1, y=data$values, label = data$group)
    
    if(!is.null(fillScale)) hcb <- hcb + fillScale
    return(hcb)
  }
  
  #CSTColors <- brewer.pal(number_of_groups_in_annotation,"Paired")# Length 6 for consistency with pre-revision CST+ coloration
  names(colour_palette) <- levels(annotations$group)
  CSTColorScale <- scale_colour_manual(name = "group",
                                       values = colour_palette)
  CSTFillScale <- scale_fill_manual(name = "group",
                                    values = colour_palette)
  
  hcb <- make_hcb(annotations_reordered,
                  "group",
                  name="group",
                  fillScale = CSTFillScale)
  
  pdf(plot_filename,
      width = 20,
      height = 11)
  print(hcb)
  
  dev.off()
  return(hcb)
}


### Run on trees!

### MPL ##### 

#### MPL With references ####

tree_plot_and_then_bubbleplot_print("../results/RAxML_bipartitions.RdRptree",
                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                 plot_tree_pdf =paste0(figures_dir,"RdRp_95_miseq_data_with_env_iso_and_ref_Raxml.pdf"),
                                 OTU_table = normalized_MPL_OTUs,
                                 paste0(figures_dir,"RdRp_tree_bubbleplot_ggplot.pdf"))


RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                                 paste0(figures_dir,"RdRp_95_miseq_data_with_env_iso_and_ref_Raxml.pdf"))

MPL_otu_table_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs)

high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)
normalized_MPL_OTUs_no_high_res <- subset(normalized_MPL_OTUs, !(rownames(normalized_MPL_OTUs) %in% high_res_vcs))
normalized_MPL_OTUs_no_high_res_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs_no_high_res)


testing_MPL_bubbleplot <- bubbleplot_by_tree_ordering(normalized_MPL_OTUs_no_high_res_with_isolates,
                                                           RdRp_95_ref_tree, 
                                                           paste0(figures_dir,"RdRp_miseq_data_bubbleplot_no_high_res.pdf"))


## only high res

normalized_MPL_OTUs_only_high_res <- subset(normalized_MPL_OTUs, rownames(normalized_MPL_OTUs) %in% high_res_vcs)

normalized_MPL_OTUs_only_high_res_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs_only_high_res)


normalized_MPL_OTUs_no_high_res_with_isolates <- as.data.frame(prop.table(normalized_MPL_OTUs_no_high_res_with_isolates)*100)

## no high res
add_to_this_MPL_bubbleplot <- bubbleplot_by_tree_ordering(add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs_no_high_res_with_isolates),
                                                               RdRp_95_ref_tree, 
                                                               paste0(figures_dir,"RdRp_miseq_data_bubbleplot_no_high_res_test.pdf"))

RdRp_OTU_annotations <- read.csv("../results/RdRp_groups_with_OTUs.csv")

## want to order like
mypalette<-brewer.pal(7,"Dark2")
MPL_colour <-  c(mypalette[1],
                 "grey",
                 mypalette[2:7],
                 "white")



p1<- add_to_this_MPL_bubbleplot+  theme(axis.title.y=element_blank(), axis.ticks=element_blank()) +
  theme(axis.text.y=element_blank())+
  annotate("text",
           x = missing_MPL_dates,
           y = 0.02,
           label="x",
           colour = "grey",
           size=7)+guides(size=guide_legend(title="Relative\nAbundance (%)"))


p2 <- make_label_for_bubbleplot_based_on_tree_ordering_with_annotations(normalized_MPL_OTUs_no_high_res_with_isolates,
                                                                     RdRp_95_ref_tree, 
                                                                     paste0(figures_dir,"RdRp_miseq_data_bubbleplot_no_high_res_label.pdf"), 
                                                                     RdRp_OTU_annotations,
                                                                     9,
                                                                     MPL_colour) +theme(legend.position="left")

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxHeight)
gp2$heights[2:5] <- as.list(maxHeight)

pdf(paste0(figures_dir,"MPL_bubbleplot_with_label.pdf"), onefile=FALSE, width = 20, height = 11)
grid.arrange( gp2,gp1, ncol=2,widths=c(1/10,9/10))
dev.off()


### gp23  ####

#### gp23  ith references ####

#  gp23_95_ref_tree <- read.tree("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result")
#  
#  gp23_95_ref_tree <- root( gp23_95_ref_tree, "Enterobacteria_phage_T4", resolve.root=TRUE)
#  
gp23_95_ref_tree <- plot_tree_root_and_ladderize("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
                                                 root = "Enterobacteria_phage_T4",
                                                 paste0(figures_dir,"gp23_95_miseq_data_with_env_iso_and_ref_Raxml.pdf"))


gp23_OTU_annotations <- read.csv("../results/gp23_groups_with_OTUs.csv")


## having problems with memory for this tree...

tree_plot_and_then_bubbleplot_print("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
                                 root = "Enterobacteria_phage_T4",
                                 plot_tree_pdf =paste0(figures_dir,"gp23_95_miseq_data_with_env_iso_and_ref_Raxml.pdf"),
                                 OTU_table = normalized_gp23_OTUs,
                                 paste0(figures_dir,"gp23_tree_bubbleplot_ggplot.pdf"))

gp23_otu_table_with_isolates <- add_isolates_to_otu_table(gp23_95_ref_tree, normalized_gp23_OTUs)

normalized_gp23_OTUs_no_high_res <- subset(gp23_otu_table_with_isolates, !(rownames(gp23_otu_table_with_isolates) %in% high_res_vcs))

normalized_gp23_OTUs_no_high_res <- as.data.frame(prop.table(normalized_gp23_OTUs_no_high_res)*100)


testing_gp23_bubbleplot <- bubbleplot_by_tree_ordering(normalized_gp23_OTUs_no_high_res,
                                                            gp23_95_ref_tree, 
                                                            paste0(figures_dir,"gp23_miseq_data_bubbleplot_no_high_res.pdf"))




mypalette<-brewer.pal(9,"Set3")
gp23_colour <- c("grey", mypalette, "white") 


p1<- testing_gp23_bubbleplot+  theme(axis.title.y=element_blank(), axis.ticks=element_blank()) +
  theme(axis.text.y=element_blank())+
  annotate("text",
           x = missing_gp23_dates,
           y = 0.02,
           label="x",
           colour = "grey",
           size=7)+guides(size=guide_legend(title="Relative\nAbundance (%)"))


p2 <- make_label_for_bubbleplot_based_on_tree_ordering_with_annotations(normalized_gp23_OTUs_no_high_res,
                                                                     gp23_95_ref_tree, 
                                                                     paste0(figures_dir,"gp23_miseq_data_bubbleplot_no_high_res_label.pdf"), 
                                                                     gp23_OTU_annotations,10,
                                                                     gp23_colour) +theme(legend.position="left")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxHeight)
gp2$heights[2:5] <- as.list(maxHeight)

pdf(paste0(figures_dir,"gp23_bubbleplot_with_label.pdf"), onefile=FALSE, width = 20, height = 11)
grid.arrange( gp2,gp1, ncol=2,widths=c(1/10,9/10))
dev.off()



normalized_gp23_OTUs_only_high_res <- subset(normalized_gp23_OTUs, rownames(normalized_gp23_OTUs) %in% high_res_vcs)


bubbleplot_by_tree_ordering_high_res_view(normalized_gp23_OTUs_only_high_res,
                                                  gp23_95_ref_tree, 
                                                  paste0(figures_dir,"gp23_95_miseq_data_bubbleplot_only_high_res.pdf"))




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

bubbleplot_by_tree_ordering_with_tax <- function (OTU_table_with_isolates,
                                                          tree, 
                                                          plot_filename, 
                                                          taxonomy) {
  #OTU_table_with_isolates <- normalized_16s_OTUs
  #taxonomy <- taxonomy_16s
  #tree <- S16_tree
  
  OTU_table_long <- melt(as.matrix(OTU_table_with_isolates))
  OTU_table_long$value[OTU_table_long$value == 0] <- NA
  OTU_table_long$Date <- Jericho_data$Date[match(OTU_table_long$Var1, Jericho_data$VC_number)]
  OTU_table_long$Phylum <- taxonomy$Phylum[match(OTU_table_long$Var2, taxonomy$otu_number)]
  
  OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(tree$tip.label))
  #OTU_table_long <- within(OTU_table_long, Phylum <- reorder(OTU_table_long, Var2))
  
  tree_tips_in_table <- subset(tree$tip.label, tree$tip.label %in% OTU_table_long$Var2)
  unique_and_in_tree <- subset(unique(OTU_table_long$Var2),unique(OTU_table_long$Var2)%in% tree_tips_in_table )
  
  OTU_table_redone <- unique_and_in_tree[order(unique_and_in_tree, tree_tips_in_table)]
  
  OTU_table_Phylum <- taxonomy$Phylum[match(OTU_table_redone, taxonomy$otu_number)]
  
  #OTU_table_Phylum <- reorder(OTU_table_Phylum, OTU_table_redone)
  
  
  print(head(OTU_table_long))
  
  #reordered_phyla <- reorder(OTU_table_long$Phylum, OTU_table_long$Var2)
  
  OTU_bubbleplot <- ggplot(OTU_table_long, 
                        aes(Date,Var2)) +
    season_line +
    spring_bloom_line+
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = bubbleplot_fill,
                        high = "#40627C",
                        na.value = bubbleplot_fill ,
                        trans = "log"
    ) +
    theme_JAG_presentation()+ 
    labs(x = "",  y = "") + 
    date_scaling + scale_y_discrete(labels=rev(OTU_table_Phylum)) 
  
  
  #new_bubbleplot <- OTU_bubbleplot 
  pdf(plot_filename, width = 20, height = 11)
  print(OTU_bubbleplot)
  dev.off()
  return(OTU_bubbleplot)
  
}


label_bar_for_bubbleplot <- function (taxonomy, tree) {
  
  data_for_bar <- as.data.frame(taxonomy$Phylum[match(tree$tip.label, taxonomy$otu_number)])
  names(data_for_bar)[1] <- "phylum"
  data_for_bar$phylum <- factor(data_for_bar$phylum, levels=unique(as.character(data_for_bar$phylum)))
  
  return(ggplot(data_for_bar, aes(x=1,y=phylum, fill=phylum))+geom_bar(stat="identity")+
           xlab(NULL) + ylab(NULL) +
           theme(axis.title=element_blank(), axis.ticks=element_blank()) +
           theme(axis.text.x=element_blank())+
           theme(axis.text.y=element_blank())
  )
} 



S18_tree <- plot_tree_and_ladderize("../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    paste0(figures_dir,"S18_97_miseq_data_Fasttree.pdf"))

S18_tree$tip.label <- gsub("_size_.*$", "",S18_tree$tip.label)
colnames(normalized_18s_OTUs) <- gsub(".size.*.", "",colnames(normalized_18s_OTUs))

bubbleplot_by_tree_ordering(normalized_18s_OTUs,
                                    S18_tree, 
                                    paste0(figures_dir,"S18_97_miseq_data_bubbleplot.pdf"))


testing_18s_bubbleplot <- bubbleplot_by_tree_ordering(normalized_18s_OTUs,
                                                           S18_tree, 
                                                           paste0(figures_dir,"S18_97_miseq_data_bubbleplot.pdf"))

taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)

pdf(paste0(figures_dir,"S18_97_miseq_data_bubbleplot_bar_label.pdf"), width=10, height=10)
label_bar_for_bubbleplot(taxonomy_18s,
                      S18_tree)
dev.off()



testing_18s_bubbleplot_with_phyla <- bubbleplot_by_tree_ordering_with_tax(normalized_18s_OTUs,
                                                                               S18_tree, 
                                                                               paste0(figures_dir,"S18_97_miseq_data_bubbleplot_phyla_names.pdf"),
                                                                               taxonomy_18s)

## fix so it matches!!!!

normalized_18s_OTUs_no_high_res <- subset(normalized_18s_OTUs, !(rownames(normalized_18s_OTUs) %in% high_res_vcs))
testing_18s_bubbleplot <- bubbleplot_by_tree_ordering(normalized_18s_OTUs_no_high_res,
                                                           S18_tree, 
                                                           paste0(figures_dir,"S18_97_miseq_data_bubbleplot_no_high_res.pdf"))

bubbleplot_by_tree_ordering_with_tax(normalized_18s_OTUs_no_high_res,
                                             S18_tree, 
                                             paste0(figures_dir,"S18_97_miseq_data_bubbleplot_no_high_res_phyla_names.pdf"),
                                             taxonomy_18s)


p1<- testing_18s_bubbleplot +  theme(axis.title.y=element_blank(), axis.ticks=element_blank()) +
  theme(axis.text.y=element_blank())

p2 <-label_bar_for_bubbleplot(taxonomy_18s,
                           S18_tree) +theme(legend.position="left")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxHeight)
gp2$heights[2:5] <- as.list(maxHeight)

grid.arrange( gp2,gp1, ncol=2,widths=c(1/10,9/10))

pdf(paste0(figures_dir,"18s_bubbleplot_with_label.pdf"), onefile=FALSE, width = 20, height = 11)
#grid.arrange( gp2,gp1, ncol=2,widths=c(2/10,8/10))
p1

dev.off()



#### 16s ####

S16_tree <- plot_tree_and_ladderize("../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    paste0(figures_dir,"S16_97_miseq_data_Fasttree.pdf"))

S16_tree$tip.label <- gsub("_size_.*$", "",S16_tree$tip.label)
colnames(normalized_16s_OTUs) <- gsub(".size.*.", "",colnames(normalized_16s_OTUs))


bubbleplot_by_tree_ordering(normalized_16s_OTUs,
                                    S16_tree, 
                                    paste0(figures_dir,"S16_97_miseq_data_bubbleplot.pdf"))

taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)

testing_16s_bubbleplot <- bubbleplot_by_tree_ordering_with_tax(normalized_16s_OTUs,
                                                                    S16_tree, 
                                                                    paste0(figures_dir,"S16_97_miseq_data_bubbleplot_with_tax.pdf"),
                                                                    taxonomy_16s)


normalized_16s_OTUs_no_high_res <- subset(normalized_16s_OTUs, !(rownames(normalized_16s_OTUs) %in% high_res_vcs))
testing_16s_bubbleplot_no_high_res <- bubbleplot_by_tree_ordering(normalized_16s_OTUs_no_high_res,
                                                                       S16_tree, 
                                                                       paste0(figures_dir,"S16_97_miseq_data_bubbleplot_no_high_res.pdf"))



# 
# testing_16s_bubbleplot_no_high_res <- bubbleplot_by_tree_ordering_with_tax(normalized_16s_OTUs_no_high_res,
#                                                                                 S16_tree, 
#                                                                                 paste0(figures_dir,"S16_97_miseq_data_bubbleplot_no_high_res.pdf"),
#                                                                                 taxonomy_16s)
# 
# 
# p1<- testing_16s_bubbleplot_no_high_res+  theme(axis.title.y=element_blank(), axis.ticks=element_blank()) + theme(axis.text.y=element_blank())
# 
# p2 <-label_bar_for_bubbleplot(taxonomy_16s,
#                            S16_tree) +theme(legend.position="left")
# gp1<-ggplotGrob(p1)
# gp2<-ggplotGrob(p2)  
# 
# maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
# gp1$heights[2:5] <- as.list(maxHeight)
# gp2$heights[2:5] <- as.list(maxHeight)
# 
# grid.arrange( gp2,gp1, ncol=2,widths=c(1/10,9/10))
# 
# pdf(paste0(figures_dir,"16s_bubbleplot_with_label.pdf"), onefile=FALSE, width = 20, height = 11)
# #grid.arrange( gp2,gp1, ncol=2,widths=c(2/10,8/10))
# p1
# dev.off()
# 

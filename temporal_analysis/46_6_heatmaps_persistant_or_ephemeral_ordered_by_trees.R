## Author: Julia Gustavsen
##Date: 24 June 2013
##Script to re-create the phylogenetic tree with heatmap that I had made with iTOL.
## Will have to line up the heatmap with the phylogenetic tree. I am adding zeros to the OTU table for the isolates and cloned sequences obtained from Genbank.
#Last updated: 6 September 2014
# 
# TO_DO
# Indicate missing samples


library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(ape)
library(phytools)
library(scales)
library(ggtree)
library(Hmisc)
library(plyr)
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
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
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
spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=1)
date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")))


high_res_vcs <- c(1198,
                  1199,
                  1200,
                  1201,
                  1202)


## want Jericho_dates that are not in MPL table and high res
Jericho_no_high_res_vcs <- Jericho_data$VC_number[!(Jericho_data$VC_number %in% high_res_vcs)]

## want to annotate missing data days. 
## want to annotate those days using star
missing_MPL_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_MPL_OTUs))]
missing_MPL_dates <- Jericho_data$Date[match(missing_MPL_samples,
                                             Jericho_data$VC_number)]

missing_gp23_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_gp23_OTUs))]
missing_gp23_dates <- Jericho_data$Date[match(missing_gp23_samples,
                                              Jericho_data$VC_number)]


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

make_heatmap_based_on_tree_ordering <- function (OTU_table_with_isolates,
                                                 tree, 
                                                 plot_filename) {
 OTU_table_long <- melt(as.matrix(OTU_table_with_isolates))
 OTU_table_long$value[OTU_table_long$value == 0] <- NA
 OTU_table_long$Date <- Jericho_data$Date[match(OTU_table_long$Var1, Jericho_data$VC_number)]
 
 OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(tree$tip.label))
 
 ## want it to be something like the annotations or taxonomic annotations. 
 OTU_heatmap_formatted <- ggplot(OTU_table_long, 
                                 aes(Date,Var2)) + 

  # spring_bloom_line+
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(name = "Relative\n abundance\n",
                      low = "lightblue",
                      high = "#40627C",
                      na.value = heatmap_fill ,
                      trans = "log",
                      breaks=c(min(OTU_table_long$value, na.rm = TRUE), 
                               sqrt(max(OTU_table_long$value, na.rm= TRUE)),
                               max(OTU_table_long$value, na.rm = TRUE)),
                      labels=c(min(OTU_table_long$value, na.rm = TRUE),
                               round(sqrt(max(OTU_table_long$value, na.rm= TRUE)),digits = -1),
                               round(max(OTU_table_long$value, na.rm = TRUE), digits = -2)
                      )) +
  labs(x = "",  y = "") + 
   season_line +
  date_scaling+ 
 theme_JAG_presentation()
 ##fix up formatting
 
 
 pdf(plot_filename, onefile=FALSE, width = 20, height = 11)
 print(OTU_heatmap_formatted)
 dev.off()
 return(OTU_heatmap_formatted)
}


make_heatmap_with_persistant_and_ephemeral_otus <- function (OTU_table_with_isolates_pers,
                                                             OTU_table_with_isolates_eph,
                                                 tree, 
                                                 plot_filename) {
 #### Heatmap with log scaled values ####
 #otu_table$VC_number = rownames(otu_table)
 #long_otus <- melt(otu_table, id="VC_number", variable.name="OTUid")
 ## Add in dates:
 OTU_table_long_pers <- melt(as.matrix(OTU_table_with_isolates_pers))
 OTU_table_long_pers$value[OTU_table_long_pers$value == 0] <- NA
 OTU_table_long_pers$Date <- Jericho_data$Date[match(OTU_table_long_pers$Var1, Jericho_data$VC_number)]
 OTU_table_long_pers$Var2 <- factor(OTU_table_long_pers$Var2, levels = factor(tree$tip.label))
 OTU_table_long_pers$group <- rep("pers", dim(OTU_table_long_pers)[1])
 
 OTU_table_long_eph <- melt(as.matrix(OTU_table_with_isolates_eph))
 OTU_table_long_eph$value[OTU_table_long_eph$value == 0] <- NA
 OTU_table_long_eph$Date <- Jericho_data$Date[match(OTU_table_long_eph$Var1, Jericho_data$VC_number)]
 OTU_table_long_eph$Var2 <- factor(OTU_table_long_eph$Var2, levels = factor(tree$tip.label))
 OTU_table_long_eph$group <- rep("eph", dim(OTU_table_long_eph)[1])

 OTU_table_pers_and_eph <- rbind(OTU_table_long_pers, OTU_table_long_eph)
 print(OTU_table_pers_and_eph)
 
 ## want it to be something like the annotations or taxonomic annotations. 
 OTU_heatmap <- ggplot(OTU_table_pers_and_eph, 
                       aes(Date,Var2)) + 
  # spring_bloom_line+
  geom_tile(fill = "white") +
  geom_tile(data = na.omit(OTU_table_pers_and_eph),
            aes(fill = group),
            size = 4) +
  scale_fill_discrete(labels = c("Ephemeral",
                                 "Persistent"))+
  theme_JAG_presentation()
 ##fix up formatting

 OTU_heatmap_formatted <- OTU_heatmap + 
  labs(y = "",
       fill="",
       x="\nDate") + 
  date_scaling+
   season_line +
  theme(legend.text = element_text(size = 20))
 #print out
 
 pdf(plot_filename, onefile=FALSE, width = 20, height = 11)
 print(OTU_heatmap_formatted)
 dev.off()
 return(OTU_heatmap_formatted)
}


make_heatmap_with_persistant_medium_and_ephemeral_otus <- function (OTU_table_with_isolates_pers,
                                                                    OTU_table_with_isolates_med,
                                                             OTU_table_with_isolates_eph,
                                                             tree, 
                                                             plot_filename) {
 #### Heatmap with log scaled values ####
 #otu_table$VC_number = rownames(otu_table)
 #long_otus <- melt(otu_table, id="VC_number", variable.name="OTUid")
 ## Add in dates:
 OTU_table_long_pers <- melt(as.matrix(OTU_table_with_isolates_pers))
 OTU_table_long_pers$value[OTU_table_long_pers$value == 0] <- NA
 OTU_table_long_pers$Date <- Jericho_data$Date[match(OTU_table_long_pers$Var1, Jericho_data$VC_number)]
 OTU_table_long_pers$Var2 <- factor(OTU_table_long_pers$Var2, levels = factor(tree$tip.label))
 OTU_table_long_pers$group <- rep("pers", dim(OTU_table_long_pers)[1])
 
 OTU_table_long_med <- melt(as.matrix(OTU_table_with_isolates_med))
 OTU_table_long_med$value[OTU_table_long_med$value == 0] <- NA
 OTU_table_long_med$Date <- Jericho_data$Date[match(OTU_table_long_med$Var1, Jericho_data$VC_number)]
 OTU_table_long_med$Var2 <- factor(OTU_table_long_med$Var2, levels = factor(tree$tip.label))
 OTU_table_long_med$group <- rep("med", dim(OTU_table_long_med)[1])
 print(OTU_table_long_med)
 OTU_table_long_eph <- melt(as.matrix(OTU_table_with_isolates_eph))
 OTU_table_long_eph$value[OTU_table_long_eph$value == 0] <- NA
 OTU_table_long_eph$Date <- Jericho_data$Date[match(OTU_table_long_eph$Var1, Jericho_data$VC_number)]
 OTU_table_long_eph$Var2 <- factor(OTU_table_long_eph$Var2, levels = factor(tree$tip.label))
 OTU_table_long_eph$group <- rep("eph", dim(OTU_table_long_eph)[1])
 
 OTU_table_pers_and_eph <- rbind(OTU_table_long_pers, OTU_table_long_med,OTU_table_long_eph)
 #print(OTU_table_pers_and_eph)
 
 ## want it to be something like the annotations or taxonomic annotations. 
 OTU_heatmap <- ggplot(OTU_table_pers_and_eph, 
                       aes(Date,Var2)) + 
  # spring_bloom_line+
  geom_tile(fill = "white") +
  geom_tile(data = na.omit(OTU_table_pers_and_eph),
            aes(fill = group),
            size=4) +
  #   scale_fill_gradient(low = "lightblue",
  #                       high = "#40627C",
  #                       na.value = heatmap_fill ,
  #                       trans = "log"
  #   ) +
  #   geom_tile(aes(fill = group)) +
  #   scale_fill_manual(values = c("pers" = "red","eph" = "pink"))
  #   +
  
  
  theme_JAG_presentation()
 ##fix up formatting
 base_size=9
 
 OTU_heatmap_formatted <- OTU_heatmap + 
  labs(x = "",
       y = "")+ 
  date_scaling+
   season_line +
  theme(legend.text = element_text(size = 20))
 #print out
 
 pdf(plot_filename, onefile=FALSE, width = 20, height = 11)
 print(OTU_heatmap_formatted)
 dev.off()
 return(OTU_heatmap_formatted)
}



make_heatmap_based_on_tree_ordering_high_res_view <- function (OTU_table_with_isolates,
                                                               tree, 
                                                               plot_filename) {
 OTU_table_long <- melt(as.matrix(OTU_table_with_isolates))
 OTU_table_long$value[OTU_table_long$value == 0] <- NA
 OTU_table_long$Date <- Jericho_data$Date[match(OTU_table_long$Var1, Jericho_data$VC_number)]
 OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(tree$tip.label))
 
 ## want it to be something like the annotations or taxonomic annotations. 
 OTU_heatmap <- ggplot(OTU_table_long, 
                       aes(Date,Var2)) + 
  # spring_bloom_line+
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(name="Relative\nabundance",
                      low = "lightblue",
                      high = "#40627C",
                      na.value = heatmap_fill ,
                      trans = "log"
  ) +
   season_line +
  theme_JAG_presentation()
 ##fix up formatting
 base_size=9
 
 OTU_heatmap_formatted <- OTU_heatmap + 
  labs(x = "",  y = "") + 
  scale_x_date(breaks = date_breaks("week"),
                                       limits = c(as.Date("2011-01-15"),
                                                  as.Date("2011-02-15")))+
  theme(legend.text = element_text(size = 20))
 #print out
 
 pdf(plot_filename, onefile=FALSE, width = 20, height = 11)
 print(OTU_heatmap_formatted)
 dev.off()
 return(OTU_heatmap_formatted)
}


tree_plot_and_then_heatmap_print <- function (tree_file,
                                              root,
                                              plot_tree_pdf,
                                              OTU_table,
                                              heatmap_plot_pdf) {
 tree_refs_env <- plot_tree_root_and_ladderize(tree_file, root, plot_tree_pdf)
 OTU_table_with_isolates <- add_isolates_to_otu_table(tree_refs_env,
                                                      OTU_table )
 make_heatmap_based_on_tree_ordering(OTU_table_with_isolates,
                                     tree_refs_env,
                                     heatmap_plot_pdf)
}

make_label_for_heatmap_based_on_tree_ordering_with_annotations <- function (OTU_table_with_isolates,
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
   new_row <- data.frame(X = "X",
                         group = "",
                         values = i)
   annotations <-rbind(annotations,
                       new_row)
   #colnames(OTU_table_with_isolates)[dim(OTU_table_with_isolates)[2]] <- paste(i)   
  }
 } 
 #### order matrix so that it is ordered like the tree ####
 annotations_reordered <- annotations[order(match(tree$tip.label,annotations$values)),]
 annotations_reordered$values <- factor(annotations_reordered$values,
                                        levels = factor(tree$tip.label))
 
 make_hcb <- function(data, var, name = NULL, fillScale = NULL, ...) {
  ## from Suann Holmes
  ## http://statweb.stanford.edu/~susan/papers/Pregnancy/PNAS_Vaginal_Analysis.html
  #data <- annotations_reordered
  hcb <- ggplot(data = data,
                aes_string(x = 1,
                           y = "values",
                           fill = "group")) + 
   geom_raster()+
   xlab(NULL)+
   ylab(NULL) +
   theme(axis.title=element_blank(), 
         axis.ticks=element_blank(),
   axis.text.x=element_blank(),
   axis.text.y=element_blank(),
   axis.line = element_blank(),
   legend.text = element_text(size = 18)) 
  #+
  #    theme(axis.text.y=element_text(size=8, face="bold")) +
  #    theme(plot.margin=unit(c(0,0,0,0),"lines"), 
  #          axis.ticks.margin = unit(0,"null"), ...) +
  #    guides(fill=F)
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

 return(hcb)
}


### Run on trees!

### MPL ##### 

#### MPL With references ####

RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                                 "../figures/RdRp_95_miseq_data_with_env_iso_and_ref_Raxml.pdf")

MPL_otu_table_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs)

get_otu_table_highlighting_persistant_otus <- function (otu_table_with_isolates) {
 ### so be persistant OTU needs to be there 90% of time and ephemral less than 20% of time. For now show on same heatmap so can compare. 
  otu_table_with_isolates_for_persistant <-  otu_table_with_isolates # by OTU-in column
  min_percent_sites <- 0.90  # 90% of times for persistant OTUs
  number_sites <- dim(otu_table_with_isolates_for_persistant)[1] 
  min_number_sites <- min_percent_sites * number_sites
  ## get columns where OTUs are in more than 3 sites
  
  count_otus <- adply(otu_table_with_isolates_for_persistant, 2,function(x)sum(x>0))
  persistant_otus <- count_otus[count_otus$V1 > min_number_sites,]
  
  ## or do I just want to put all the other ones as 0?
  otu_table_with_isolates_for_persistant[which(!(colnames(otu_table_with_isolates_for_persistant) %in% persistant_otus$X1))] <- 0
  return(otu_table_with_isolates_for_persistant)
}

get_otu_table_highlighting_ephemeral_otus <- function (otu_table_with_isolates) {
 ### so be persistant OTU needs to be there 90% of time and ephemral less than 20% of time. For now show on same heatmap so can compare. 
 otu_table_with_isolates_for_ephemeral <-  otu_table_with_isolates # by OTU-in column
 min_percent_sites <- 0.20  # 90% of times for persistant OTUs
 number_sites <- dim(otu_table_with_isolates_for_ephemeral)[1] 
 min_number_sites <- min_percent_sites * number_sites
 ## get columns where OTUs are in more than 3 sites
 
 count_otus <- adply(otu_table_with_isolates_for_ephemeral, 2,function(x)sum(x>0))
 ephemeral_otus <- count_otus[count_otus$V1 < min_number_sites,]
 
 ## or do I just want to put all the other ones as 0?
 otu_table_with_isolates_for_ephemeral[which(!(colnames(otu_table_with_isolates_for_ephemeral) %in% ephemeral_otus$X1))] <- 0
 return(otu_table_with_isolates_for_ephemeral)
}

get_otu_table_highlighting_medium_otus <- function (otu_table_with_isolates) {
 ### so be persistant OTU needs to be there 90% of time and ephemral less than 20% of time. For now show on same heatmap so can compare. 
 otu_table_with_isolates_for_med <-  otu_table_with_isolates # by OTU-in column
 min_percent_sites <- 0.20  # 90% of times for persistant OTUs
 number_sites <- dim(otu_table_with_isolates_for_med)[1] 
 min_number_sites <- min_percent_sites * number_sites
 ## get columns where OTUs are in more than 3 sites
 
 max_percent_sites <- 0.90  # 90% of times for persistant OTUs
 max_number_sites <- max_percent_sites * number_sites
 
 count_otus <- adply(otu_table_with_isolates_for_med, 2,function(x)sum(x>0))
 med_otus <- count_otus[(count_otus$V1 > min_number_sites) & (count_otus$V1 < max_number_sites),]
 print(med_otus)
 ## or do I just want to put all the other ones as 0?
 otu_table_with_isolates_for_med[which(!(colnames(otu_table_with_isolates_for_med) %in% med_otus$X1))] <- 0
 return(otu_table_with_isolates_for_med)
}


MPL_otu_table_with_isolates_for_persistant <- get_otu_table_highlighting_persistant_otus(MPL_otu_table_with_isolates)
MPL_otu_table_with_isolates_for_ephemeral <- get_otu_table_highlighting_ephemeral_otus(MPL_otu_table_with_isolates)

make_heatmap_based_on_tree_ordering(MPL_otu_table_with_isolates_for_persistant,
                                                           RdRp_95_ref_tree, 
                                                           "../figures/RdRp_miseq_data_heatmap_only_persistant.pdf")
make_heatmap_based_on_tree_ordering(MPL_otu_table_with_isolates_for_ephemeral,
                                    RdRp_95_ref_tree, 
                                    "../figures/RdRp_miseq_data_heatmap_only_ephemeral.pdf")


high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)
normalized_MPL_OTUs_no_high_res <- subset(normalized_MPL_OTUs, !(rownames(normalized_MPL_OTUs) %in% high_res_vcs))
normalized_MPL_OTUs_no_high_res_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs_no_high_res)

make_heatmap_based_on_tree_ordering(normalized_MPL_OTUs_no_high_res_with_isolates,
                                                           RdRp_95_ref_tree, 
                                                           "../figures/RdRp_miseq_data_heatmap_no_high_res.pdf")

MPL_otu_table_with_isolates_for_persistant_no_high_res <- get_otu_table_highlighting_persistant_otus(normalized_MPL_OTUs_no_high_res_with_isolates)
MPL_otu_table_with_isolates_for_ephemeral_no_high_res <- get_otu_table_highlighting_ephemeral_otus(normalized_MPL_OTUs_no_high_res_with_isolates)

MPL_otu_table_with_isolates_for_med_no_high_res <- get_otu_table_highlighting_medium_otus(normalized_MPL_OTUs_no_high_res_with_isolates)

make_heatmap_based_on_tree_ordering(MPL_otu_table_with_isolates_for_persistant_no_high_res,
                                    RdRp_95_ref_tree, 
                                    "../figures/RdRp_miseq_data_heatmap_only_persistant_no_high_res.pdf")
make_heatmap_based_on_tree_ordering(MPL_otu_table_with_isolates_for_ephemeral_no_high_res,
                                    RdRp_95_ref_tree, 
                                    "../figures/RdRp_miseq_data_heatmap_only_ephemeral_no_high_res.pdf")

MPL_pers_and_eph_heatmap <- make_heatmap_with_persistant_and_ephemeral_otus(MPL_otu_table_with_isolates_for_persistant_no_high_res,
                                                MPL_otu_table_with_isolates_for_ephemeral_no_high_res,
                                                RdRp_95_ref_tree, 
                                                "../figures/RdRp_miseq_data_heatmap_persistant_and_ephemeral_no_high_res.pdf")

Jericho_no_high_res_vcs <- Jericho_data$VC_number[!(Jericho_data$VC_number %in% high_res_vcs)]

missing_MPL_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_MPL_OTUs_no_high_res_with_isolates))]
missing_MPL_dates <- Jericho_data$Date[match(missing_MPL_samples,
                                             Jericho_data$VC_number)]

MPL_pers_and_eph_heatmap <- MPL_pers_and_eph_heatmap+ annotate("text",
                                                               x = missing_MPL_dates,
                                                               y = 0.00,
                                                               label="*",
                                                               size = 12)


make_heatmap_with_persistant_medium_and_ephemeral_otus(MPL_otu_table_with_isolates_for_persistant_no_high_res,
                                                       MPL_otu_table_with_isolates_for_med_no_high_res,
                                                       MPL_otu_table_with_isolates_for_ephemeral_no_high_res,
                                                       RdRp_95_ref_tree, 
                                                       "../figures/RdRp_miseq_data_heatmap_persistant_medium_and_ephemeral_no_high_res.pdf")

RdRp_OTU_annotations <- read.csv("../results/RdRp_groups_with_OTUs.csv")

## want to order like
mypalette<-brewer.pal(7,"Dark2")
MPL_colour <- c(mypalette[1],
                "grey",
                mypalette[2:7],
                "white")

p1 <- MPL_pers_and_eph_heatmap +  
  theme(axis.title.y = element_blank(),
        axis.ticks = element_blank())+
 theme(axis.text.y = element_blank())

p2 <- make_label_for_heatmap_based_on_tree_ordering_with_annotations(normalized_MPL_OTUs_no_high_res_with_isolates,
                                                                     RdRp_95_ref_tree, 
                                                                     "../figures/RdRp_miseq_data_heatmap_no_high_res_label.pdf", 
                                                                     RdRp_OTU_annotations,
                                                                     9,
                                                                     MPL_colour) +theme(legend.position="left")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxHeight)
gp2$heights[2:5] <- as.list(maxHeight)

pdf("../figures/MPL_heatmap_with_pers_and_eph_with_label.pdf", onefile=FALSE, width = 20, height = 11)
grid.arrange( gp2,gp1, ncol=2,widths=c(4/20,16/20))
dev.off()

MPL_overall_heatmap <- make_heatmap_based_on_tree_ordering(normalized_MPL_OTUs_no_high_res_with_isolates,
                                    RdRp_95_ref_tree, 
                                    "../figures/RdRp_miseq_data_heatmap_no_high_res.pdf")

p1<- MPL_overall_heatmap + 
  xlab("Date")+
  theme(axis.title.y=element_blank(), axis.ticks=element_blank()) +
  theme(axis.text.y=element_blank())+
  annotate("text",
           x = missing_MPL_dates,
           y = 0.00,
           label="*",
           size=10)

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxHeight)
gp2$heights[2:5] <- as.list(maxHeight)

## need to change the gp23 as well. 
pdf("../figures/MPL_heatmap_with_label.pdf", onefile=FALSE, width = 9, height = 4)
grid.arrange(gp2,
             gp1,
             ncol = 2,
             widths = c(6/40,34/40))
dev.off()

### gp23  ####

gp23_95_ref_tree <- plot_tree_root_and_ladderize("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
                                                 root = "Enterobacteria_phage_T4",
                                                 "../figures/gp23_95_miseq_data_with_env_iso_and_ref_Raxml.pdf")

## having problems with memory for this tree...

tree_plot_and_then_heatmap_print("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
                                 root = "Enterobacteria_phage_T4",
                                 plot_tree_pdf ="../figures/gp23_95_miseq_data_with_env_iso_and_ref_Raxml.pdf",
                                 OTU_table = normalized_gp23_OTUs,
                                 "../figures/gp23_tree_heatmap_ggplot.pdf")

gp23_otu_table_with_isolates <- add_isolates_to_otu_table(gp23_95_ref_tree, normalized_gp23_OTUs)

gp23_otu_table_with_isolates_for_persistant <- get_otu_table_highlighting_persistant_otus(gp23_otu_table_with_isolates)
gp23_otu_table_with_isolates_for_ephemeral <- get_otu_table_highlighting_ephemeral_otus(gp23_otu_table_with_isolates)

make_heatmap_based_on_tree_ordering(gp23_otu_table_with_isolates_for_persistant,
                                    gp23_95_ref_tree, 
                                    "../figures/gp23_miseq_data_heatmap_only_persistant.pdf")
make_heatmap_based_on_tree_ordering(gp23_otu_table_with_isolates_for_ephemeral,
                                    gp23_95_ref_tree, 
                                    "../figures/gp23_miseq_data_heatmap_only_ephemeral.pdf")


normalized_gp23_OTUs_no_high_res_with_isolates <- subset(gp23_otu_table_with_isolates, !(rownames(gp23_otu_table_with_isolates) %in% high_res_vcs))
make_heatmap_based_on_tree_ordering(normalized_gp23_OTUs_no_high_res_with_isolates,
                                                            gp23_95_ref_tree, 
                                                            "../figures/gp23_miseq_data_heatmap_no_high_res.pdf")

## need to make this with the label...






gp23_otu_table_with_isolates_for_persistant_no_high_res <- get_otu_table_highlighting_persistant_otus(normalized_gp23_OTUs_no_high_res_with_isolates)
gp23_otu_table_with_isolates_for_ephemeral_no_high_res <- get_otu_table_highlighting_ephemeral_otus(normalized_gp23_OTUs_no_high_res_with_isolates)

make_heatmap_based_on_tree_ordering(gp23_otu_table_with_isolates_for_persistant_no_high_res,
                                    gp23_95_ref_tree, 
                                    "../figures/gp23_miseq_data_heatmap_only_persistant_no_high_res.pdf")
make_heatmap_based_on_tree_ordering(gp23_otu_table_with_isolates_for_ephemeral_no_high_res,
                                    gp23_95_ref_tree, 
                                    "../figures/gp23_miseq_data_heatmap_only_ephemeral_no_high_res.pdf")

gp23_pers_and_eph_heatmap <- make_heatmap_with_persistant_and_ephemeral_otus(gp23_otu_table_with_isolates_for_persistant_no_high_res,
                                                gp23_otu_table_with_isolates_for_ephemeral_no_high_res,
                                                gp23_95_ref_tree, 
                                                "../figures/gp23_miseq_data_heatmap_persistant_and_ephemeral_no_high_res.pdf")


missing_gp23_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_gp23_OTUs_no_high_res_with_isolates))]

missing_gp23_dates <- Jericho_data$Date[match(missing_gp23_samples,
                                             Jericho_data$VC_number)]

gp23_pers_and_eph_heatmap <- gp23_pers_and_eph_heatmap+ annotate("text", x = missing_gp23_dates, y = 0.00, label="*", size=12)

gp23_otu_table_with_isolates_for_med_no_high_res <- get_otu_table_highlighting_medium_otus(normalized_gp23_OTUs_no_high_res_with_isolates)

make_heatmap_with_persistant_medium_and_ephemeral_otus(gp23_otu_table_with_isolates_for_persistant_no_high_res,
                                                       gp23_otu_table_with_isolates_for_med_no_high_res,
                                                       gp23_otu_table_with_isolates_for_ephemeral_no_high_res,
                                                       gp23_95_ref_tree, 
                                                       "../figures/gp23_miseq_data_heatmap_persistant_medium_and_ephemeral_no_high_res.pdf")



gp23_OTU_annotations <- read.csv("../results/gp23_groups_with_OTUs.csv")

mypalette<-brewer.pal(9,"Set3")
gp23_colour <- c("grey", mypalette[1:8], "white") 

p1<- gp23_pers_and_eph_heatmap+  theme(axis.title.y=element_blank(), axis.ticks=element_blank()) +
 theme(axis.text.y=element_blank())

p2 <- make_label_for_heatmap_based_on_tree_ordering_with_annotations(normalized_gp23_OTUs_no_high_res_with_isolates,
                                                                     gp23_95_ref_tree, 
                                                                     "../figures/gp23_miseq_data_heatmap_no_high_res_label.pdf", 
                                                                     gp23_OTU_annotations,10,
                                                                     gp23_colour)+ theme(legend.position = "left")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxHeight)
gp2$heights[2:5] <- as.list(maxHeight)

pdf("../figures/gp23_heatmap_with_pers_and_eph_with_label.pdf", onefile=FALSE, width = 20, height = 11)
grid.arrange( gp2,gp1, ncol=2,widths=c(3/20,17/20))
dev.off()

### all gp23

gp23_heatmap_all <- make_heatmap_based_on_tree_ordering(gp23_otu_table_with_isolates,
                                    gp23_95_ref_tree, 
                                    "../figures/gp23_miseq_data_heatmap.pdf")

gp23_heatmap_all_no_high_res <- make_heatmap_based_on_tree_ordering(normalized_gp23_OTUs_no_high_res_with_isolates,
                                    gp23_95_ref_tree, 
                                    "../figures/gp23_miseq_data_heatmap_no_high_res.pdf")

p1 <- gp23_heatmap_all_no_high_res+
  xlab("Date")+
  theme(axis.title.y = element_blank(),
        axis.ticks = element_blank())+
  theme(axis.text.y = element_blank())+
  annotate("text",
           x = missing_gp23_dates,
           y = 0.00,
           label="*",
           size=10)

p2 <- make_label_for_heatmap_based_on_tree_ordering_with_annotations(normalized_gp23_OTUs_no_high_res_with_isolates,
                                                                     gp23_95_ref_tree, 
                                                                     "../figures/gp23_miseq_data_heatmap_no_high_res_label.pdf", 
                                                                     gp23_OTU_annotations,10,
                                                                     gp23_colour) +theme(legend.position="left")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxHeight)
gp2$heights[2:5] <- as.list(maxHeight)

pdf("../figures/gp23_heatmap_with_label.pdf",
    onefile = FALSE,
    width = 9,
    height = 4)
grid.arrange(gp2,
              gp1,
              ncol = 2,
              widths = c(6/40,34/40))
dev.off()


normalized_gp23_OTUs_only_high_res <- subset(normalized_gp23_OTUs, rownames(normalized_gp23_OTUs) %in% high_res_vcs)


make_heatmap_based_on_tree_ordering_high_res_view(normalized_gp23_OTUs_only_high_res,
                                                  gp23_95_ref_tree, 
                                                  "../figures/gp23_95_miseq_data_heatmap_only_high_res.pdf")

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

make_heatmap_based_on_tree_ordering_with_tax <- function (OTU_table_with_isolates,
                                                          tree, 
                                                          plot_filename, 
                                                          taxonomy) {
 
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
 
 print(head(OTU_table_long))
 
 #reordered_phyla <- reorder(OTU_table_long$Phylum, OTU_table_long$Var2)
 
 OTU_heatmap <- ggplot(OTU_table_long, 
                       aes(Date,Var2)) +
  season_line +
  #spring_bloom_line+
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = heatmap_fill,
                      high = "#40627C",
                      na.value = heatmap_fill ,
                      trans = "log"
  ) +
  theme_JAG_presentation()+ 
  labs(x = "",  y = "") + 
  date_scaling + scale_y_discrete(labels=rev(OTU_table_Phylum)) 
 
 
 #new_heatmap <- OTU_heatmap 
 pdf(plot_filename, width = 20, height = 11)
 print(OTU_heatmap)
 dev.off()
 return(OTU_heatmap)
 
}


label_bar_for_heatmap <- function (taxonomy, tree) {
 
 data_for_bar <- as.data.frame(taxonomy$Phylum[match(tree$tip.label, taxonomy$otu_number)])
 names(data_for_bar)[1] <- "phylum"
 data_for_bar$phylum <- factor(data_for_bar$phylum, levels=unique(as.character(data_for_bar$phylum)))
 
 return(ggplot(data_for_bar, aes(x=1,y=phylum, fill=phylum))+geom_bar(stat="identity"))
} 



S18_tree <- plot_tree_and_ladderize("../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    "../figures/S18_97_miseq_data_Fasttree.pdf")

S18_tree$tip.label <- gsub("_size_.*$", "",S18_tree$tip.label)
colnames(normalized_18s_OTUs) <- gsub(".size.*.", "",colnames(normalized_18s_OTUs))

S18_otu_table_for_persistant <- get_otu_table_highlighting_persistant_otus(normalized_18s_OTUs)
S18_otu_table_for_ephemeral <- get_otu_table_highlighting_ephemeral_otus(normalized_18s_OTUs)

make_heatmap_based_on_tree_ordering(S18_otu_table_for_persistant,
                                    S18_tree, 
                                    "../figures/S18_miseq_data_heatmap_only_persistant.pdf")
make_heatmap_based_on_tree_ordering(S18_otu_table_for_ephemeral,
                                    S18_tree, 
                                    "../figures/S18_miseq_data_heatmap_only_ephemeral.pdf")


make_heatmap_based_on_tree_ordering(normalized_18s_OTUs,
                                                           S18_tree, 
                                                           "../figures/S18_97_miseq_data_heatmap.pdf")

taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)

pdf("../figures/S18_97_miseq_data_heatmap_bar_label.pdf", width=10, height=10)
label_bar_for_heatmap(taxonomy_18s,
                      S18_tree)
dev.off()

testing_18s_heatmap_with_phyla <- make_heatmap_based_on_tree_ordering_with_tax(normalized_18s_OTUs,
                                                                               S18_tree, 
                                                                               "../figures/S18_97_miseq_data_heatmap_phyla_names.pdf",
                                                                               taxonomy_18s)

## fix so it matches!!!!

normalized_18s_OTUs_no_high_res <- subset(normalized_18s_OTUs, !(rownames(normalized_18s_OTUs) %in% high_res_vcs))
testing_18s_heatmap <- make_heatmap_based_on_tree_ordering(normalized_18s_OTUs_no_high_res,
                                                           S18_tree, 
                                                           "../figures/S18_97_miseq_data_heatmap_no_high_res.pdf")


S18_otu_table_for_persistant_no_high_res <- get_otu_table_highlighting_persistant_otus(normalized_18s_OTUs_no_high_res )
S18_otu_table_for_ephemeral_no_high_res  <- get_otu_table_highlighting_ephemeral_otus(normalized_18s_OTUs_no_high_res )

make_heatmap_based_on_tree_ordering(S18_otu_table_for_persistant_no_high_res,
                                    S18_tree, 
                                    "../figures/S18_miseq_data_heatmap_only_persistant_no_high_res.pdf")
make_heatmap_based_on_tree_ordering(S18_otu_table_for_ephemeral_no_high_res,
                                    S18_tree, 
                                    "../figures/S18_miseq_data_heatmap_only_ephemeral_no_high_res.pdf")


make_heatmap_with_persistant_and_ephemeral_otus(S18_otu_table_for_persistant_no_high_res,
                                                S18_otu_table_for_ephemeral_no_high_res,
                                                S18_tree, 
                                                "../figures/S18_miseq_data_heatmap_persistant_and_ephemeral_no_high_res.pdf")

make_heatmap_based_on_tree_ordering_with_tax(normalized_18s_OTUs_no_high_res,
                                             S18_tree, 
                                             "../figures/S18_97_miseq_data_heatmap_no_high_res_phyla_names.pdf",
                                             taxonomy_18s)


S18_otu_table_for_med_no_high_res <- get_otu_table_highlighting_medium_otus(normalized_18s_OTUs_no_high_res)

## only high res

normalized_18s_OTUs_only_high_res <- subset(normalized_18s_OTUs, rownames(normalized_18s_OTUs) %in% high_res_vcs)


make_heatmap_based_on_tree_ordering_high_res_view(normalized_18s_OTUs_only_high_res,
                                                  S18_tree, 
                                                  "../figures/S18_97_miseq_data_heatmap_only_high_res.pdf")


make_heatmap_based_on_tree_ordering_with_tax(normalized_18s_OTUs_only_high_res,
                                             S18_tree, 
                                             "../figures/S18_97_miseq_data_heatmap_only_high_res_phyla_names.pdf",
                                             taxonomy_18s) 


#### 16s ####

S16_tree <- plot_tree_and_ladderize("../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                    "../figures/S16_97_miseq_data_Fasttree.pdf")

S16_tree$tip.label <- gsub("_size_.*$", "",S16_tree$tip.label)
colnames(normalized_16s_OTUs) <- gsub(".size.*.", "",colnames(normalized_16s_OTUs))


make_heatmap_based_on_tree_ordering(normalized_16s_OTUs,
                                    S16_tree, 
                                    "../figures/S16_97_miseq_data_heatmap.pdf")


S16_otu_table_for_persistant <- get_otu_table_highlighting_persistant_otus(normalized_16s_OTUs)
S16_otu_table_for_ephemeral <- get_otu_table_highlighting_ephemeral_otus(normalized_16s_OTUs)

make_heatmap_based_on_tree_ordering(S16_otu_table_for_persistant,
                                    S16_tree, 
                                    "../figures/S16_miseq_data_heatmap_only_persistant.pdf")
make_heatmap_based_on_tree_ordering(S16_otu_table_for_ephemeral,
                                    S16_tree, 
                                    "../figures/S16_miseq_data_heatmap_only_ephemeral.pdf")


taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)

testing_16s_heatmap <- make_heatmap_based_on_tree_ordering_with_tax(normalized_16s_OTUs,
                                                                    S16_tree, 
                                                                    "../figures/S16_97_miseq_data_heatmap_with_tax.pdf",
                                                                    taxonomy_16s)


normalized_16s_OTUs_no_high_res <- subset(normalized_16s_OTUs, !(rownames(normalized_16s_OTUs) %in% high_res_vcs))
testing_16s_heatmap_no_high_res <- make_heatmap_based_on_tree_ordering(normalized_16s_OTUs_no_high_res,
                                                                       S16_tree, 
                                                                       "../figures/S16_97_miseq_data_heatmap_no_high_res.pdf")

S16_otu_table_for_persistant_no_high_res <- get_otu_table_highlighting_persistant_otus(normalized_16s_OTUs_no_high_res )
S16_otu_table_for_ephemeral_no_high_res  <- get_otu_table_highlighting_ephemeral_otus(normalized_16s_OTUs_no_high_res )

make_heatmap_based_on_tree_ordering(S16_otu_table_for_persistant_no_high_res,
                                    S16_tree, 
                                    "../figures/S16_miseq_data_heatmap_only_persistant_no_high_res.pdf")
make_heatmap_based_on_tree_ordering(S16_otu_table_for_ephemeral_no_high_res,
                                    S16_tree, 
                                    "../figures/S16_miseq_data_heatmap_only_ephemeral_no_high_res.pdf")

make_heatmap_with_persistant_and_ephemeral_otus(S16_otu_table_for_persistant_no_high_res,
                                                S16_otu_table_for_ephemeral_no_high_res,
                                                S16_tree, 
                                                "../figures/S16_miseq_data_heatmap_persistant_and_ephemeral_no_high_res.pdf")

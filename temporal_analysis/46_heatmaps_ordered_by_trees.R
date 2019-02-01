## Author: Julia Gustavsen
##Date: 24 June 2013
##Script to re-create the phylogenetic tree with heatmap that I had made with iTOL.
## Will have to line up the heatmap with the phylogenetic tree. I am adding zeros to the OTU table for the isolates and cloned sequences obtained from Genbank.
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
str(OTU_table_long$value)
## want it to be something like the annotations or taxonomic annotations. 
OTU_heatmap_formatted <- ggplot(OTU_table_long, 
                                aes(Date,Var2)) + 
 season_line +
 # annotate("text",x=c(mid_summer_x, mid_fall_x, mid_winter_x, mid_spring_x), y=tree$tip.label[1], label=c("Summer", "Fall", "Winter", "Spring")) +
 spring_bloom_line+
 geom_tile(aes(fill = value)) +
 scale_fill_gradient(name = "Relative\n abundance",
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
 labs(y = "") + 
 date_scaling +
theme_JAG_presentation(base_size = 18)

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
  season_line +
  spring_bloom_line+
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(name = "Relative\n abundance",
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
                      ))+
  theme_JAG_presentation()
 ##fix up formatting
 base_size=9
 
 OTU_heatmap_formatted <- OTU_heatmap + 
  labs(x = "",  y = "") + scale_x_date(breaks = date_breaks("week"),
                                       limits = c(as.Date("2011-01-15"),
                                                  as.Date("2011-02-15")))
 #print out
 
 pdf(plot_filename, onefile=FALSE, width = 20, height = 11)
 print(OTU_heatmap_formatted)
 dev.off()
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
plot_otus_by_date_on_tree <- function (OTU_table_with_isolates,
                                       reordered_tree,
                                       plot_filename) {
 test_otu_table <- rev(OTU_table_with_isolates) 
 row.names(test_otu_table)<- Jericho_data$Date[match(row.names(test_otu_table),
                                                     Jericho_data$VC_number)]
 nrow(test_otu_table)
 pdf(plot_filename, 
     width = 80, 
     height = 40,
     onefile = FALSE)
 
 par(mfrow = c(3,
               (nrow(test_otu_table)/3)))
 for (i in sort(rownames(test_otu_table))) {
  
  plot(reordered_tree,
       show.tip.label = FALSE,
       main = i)
  tiplabels(tip = which(test_otu_table[i, ] > 0), 
            pch = 19,
            cex = 1.5, 
            col ="#1E90FF")
  #legend("left" , i, bty = "n")
 }
 
 dev.off()
 
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


make_heatmap_facet <- function (OTU_table_with_isolates,
                                tree,
                                group_annotations,
                                plot_filename) {
 
 #OTU_table_with_isolates <- normalized_gp23_OTUs_no_high_res
 #group_annotations <- gp23_OTU_annotations
 #tree <- gp23_95_ref_tree
 #   for(i in colnames((OTU_table_with_isolates))) {
 #    if (is.element(i, annotations$values)) {
 #     print(i)
 #    }
 #    else {
 #     print("no")
 #     print(i)
 #     new_row <- data.frame(X="X", group="groupNull",values= i)
 #     annotations <-rbind(annotations,new_row)
 #    }
 #   } 
 #   #### order matrix so that it is ordered like the tree ####
 #   annotations_reordered <- annotations[order(match(tree$tip.label,annotations$values)),]
 #   annotations_reordered$values <- factor(annotations_reordered$values, levels = factor(tree$tip.label))
 
 
 
 OTU_table_long <- melt(as.matrix(OTU_table_with_isolates))
 OTU_table_long$group <- group_annotations$group[match(OTU_table_long$Var2,
                                                       group_annotations$values)]
 OTU_table_long$value[OTU_table_long$value == 0] <- NA
 OTU_table_long$Date <- Jericho_data$Date[match(OTU_table_long$Var1, Jericho_data$VC_number)]
 OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(tree$tip.label))
 
 ## want it to be something like the annotations or taxonomic annotations. 
 OTU_heatmap <- ggplot(OTU_table_long, 
                       aes(Date,Var2)) + 
  facet_grid(group~.,
             space="free",
             scales="free") + 
  #season_line +
  #spring_bloom_line+
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(name = "Relative\n abundance",
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
                      ))+
  theme_JAG_presentation()
 ##fix up formatting
 base_size=9
 
 OTU_heatmap_formatted <- OTU_heatmap 
 #+ 
 #    labs(x = "",  y = "") + scale_x_date(breaks = date_breaks("week"),
 #                                         limits = c(as.Date("2011-01-15"),
 #                                                    as.Date("2011-02-15")))
 #print out
 
 pdf(plot_filename, onefile=FALSE, width = 20, height = 11)
 print(OTU_heatmap_formatted)
 dev.off()
 return(OTU_heatmap_formatted)
}


### Run on trees!

### MPL ##### 

#### MPL With references ####

tree_plot_and_then_heatmap_print("../results/RAxML_bipartitions.RdRptree",
                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                 plot_tree_pdf =paste0(figures_dir,"RdRp_95_miseq_data_with_env_iso_and_ref_Raxml.pdf"),
                                 OTU_table = normalized_MPL_OTUs,
                                 paste0(figures_dir,"RdRp_tree_heatmap_ggplot.pdf"))


RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                  root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                                  paste0(figures_dir,"RdRp_95_miseq_data_with_env_iso_and_ref_Raxml.pdf"))

MPL_otu_table_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs)

plot_otus_by_date_on_tree(MPL_otu_table_with_isolates,
                          RdRp_95_ref_tree,
                          paste0(figures_dir,"RdRp_95_miseq_data_with_env_iso_and_ref_Raxml_with_otus_many_trees%03d.pdf"))

high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)
normalized_MPL_OTUs_no_high_res <- subset(normalized_MPL_OTUs, !(rownames(normalized_MPL_OTUs) %in% high_res_vcs))
normalized_MPL_OTUs_no_high_res_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs_no_high_res)


testing_MPL_heatmap <- make_heatmap_based_on_tree_ordering(normalized_MPL_OTUs_no_high_res_with_isolates,
                                                           RdRp_95_ref_tree, 
                                                           paste0(figures_dir,"RdRp_miseq_data_heatmap_no_high_res.pdf"))


## only high res

normalized_MPL_OTUs_only_high_res <- subset(normalized_MPL_OTUs, rownames(normalized_MPL_OTUs) %in% high_res_vcs)

normalized_MPL_OTUs_only_high_res_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs_only_high_res)


## no high res
add_to_this_MPL_heatmap <- make_heatmap_based_on_tree_ordering(add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs_no_high_res_with_isolates),
                                                           RdRp_95_ref_tree, 
                                                           paste0(figures_dir,"RdRp_miseq_data_heatmap_no_high_res_test.pdf"))

RdRp_OTU_annotations <- read.csv("../results/RdRp_groups_with_OTUs.csv")

## want to order like
mypalette<-brewer.pal(7,"Dark2")
MPL_colour <-  c(mypalette[1],
                 "grey",
                 mypalette[2:7],
                 "white")



p1<- add_to_this_MPL_heatmap+  theme(axis.title.y=element_blank(), axis.ticks=element_blank()) +
   theme(axis.text.y=element_blank())

p2 <- make_label_for_heatmap_based_on_tree_ordering_with_annotations(normalized_MPL_OTUs_no_high_res_with_isolates,
                                                                      RdRp_95_ref_tree, 
                                                                     paste0(figures_dir,"RdRp_miseq_data_heatmap_no_high_res_label.pdf"), 
                                                                      RdRp_OTU_annotations,
                                                                     9,
                                                                     MPL_colour) +theme(legend.position="left")

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxHeight)
gp2$heights[2:5] <- as.list(maxHeight)

pdf(paste0(figures_dir,"MPL_heatmap_with_label.pdf"), onefile=FALSE, width = 20, height = 11)
grid.arrange( gp2,gp1, ncol=2,widths=c(1/10,9/10))
dev.off()



make_heatmap_based_on_tree_ordering_high_res_view(normalized_MPL_OTUs_only_high_res_with_isolates,
                                                           RdRp_95_ref_tree, 
                                                  paste0(figures_dir,"MPL_95_miseq_data_heatmap_only_high_res.pdf"))


# ggtree(RdRp_95_only_miseq, ladderize=FALSE) + theme_tree2()
# ggtree(RdRp_95_only_miseq)+geom_point(aes(shape=isTip, color=isTip), size=3)
# ggtree(RdRp_95_only_miseq, layout="fan") + geom_text(aes(label=label, angle=angle), size=3, color="purple", vjust=-0.3)
# gzoom(RdRp_95_only_miseq, grep("OTU_135", RdRp_95_only_miseq$tip.label))
# #gzoom(RdRp_95_only_miseq, "OTU_1")
# ggtree(RdRp_95_only_miseq) + geom_text(aes(label=node))
# 
# ggtree(RdRp_95_only_miseq) + geom_text(aes(label=node)) 
# #p <- ggtree(RdRp_95_only_miseq) + geom_tiplab()
# #annotation_clade(p, node=816, "selected clade", offset.text=2) %>% 
# #annotation_clade(node=1223, "selected clade", offset.text=2)
# 
# #ggtree(RdRp_95_only_miseq) %>% ggtree:hilight(node=816, fill="steelblue", alpha=.6) 
# 
# #cp <- ggtree(RdRp_95_only_miseq) %>% collapse(node=1223)
# #cp + geom_point(subset=.(node == 816), size=5, shape=23, fill="steelblue")
# tree <- groupClade(RdRp_95_only_miseq, node=816)
# ggtree(tree, aes(color=group, linetype=group))

gheatmap_jag <- function (p, OTU_table_with_isolates,tree, offset = 1, width = 10, low = "#A0B0BE", high = "#40627C", 
          color = "white", colnames = TRUE, font.size = 4) 
{
#  p <- p
#  OTU_table_with_isolates <- normalized_MPL_OTUs
#  tree <- RdRp_95_only_miseq
#  #width=1
#  offset=0
#  color=line_colour
 width <- width * (p$data$x %>% range %>% diff)/ncol(OTU_table_with_isolates)
 isTip <- x <- y <- variable <- value <- from <- to <- NULL
 df = p$data
 df = df[df$isTip, ]
 start <- max(df$x) + offset
 dd <- as.data.frame(t(OTU_table_with_isolates[,df$label[order(df$y)] ]))
  dd$y <- sort(df$y)
  dd$lab <- rownames(dd)
  dd <- melt(dd, id = c("lab", "y"))
  #dd$variable <- reorder(dd$variable, sort(as.numeric(as.character(dd$variable))), max)
  #dd$variable<- factor(dd$variable, levels = sort(as.numeric(as.character(dd$variable))))
  dd$variable <- factor( dd$variable, levels=sort(as.numeric(as.character(levels(dd$variable)))))
  dd$Date <- Jericho_data$Date[match(dd$variable, Jericho_data$VC_number)]
  
  ## need to sort by VC number or use dates...
  
  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from = dd$Date, to = V2)
  mapping <- unique(mapping)
  dd$x <- V2
#  OTU_table_long <- melt(as.matrix(OTU_table_with_isolates))
#  OTU_table_long$value[OTU_table_long$value == 0] <- NA
#  OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(tree$tip.label))

 # p2 <- p + geom_tile(data = dd, aes(as.factor(variable), lab, fill = value), color = color)
  p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), color = color)
  
  if (is(dd$value, "numeric")) {
   p2 <- p2 + scale_fill_gradient(low = low, high = high, 
                                  na.value = line_colour,trans = "log")
  }
  else {
   p2 <- p2 + scale_fill_discrete(na.value = line_colour)
  }
  if (colnames) {
   p2 <- p2 + geom_text(data = mapping, aes(x = to, label = from), 
                        y = 0,angle=90,  size = font.size)
  }
  p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())

 p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())
 p2 <- p2 + guides(fill = guide_legend(override.aes = list(colour = NULL)))
 attr(p2, "mapping") <- mapping
 return(p2)
}

# p <- ggtree(RdRp_95_only_miseq, colour=line_colour)+theme_JAG_presentation()
# p <- p + geom_tiplab(size=3)
# gheatmap_jag(p, normalized_MPL_OTUs, RdRp_95_only_miseq)

p <- ggtree(RdRp_95_ref_tree, colour=line_colour,ladderize = FALSE)+theme_JAG_presentation()
#p <- p + geom_tiplab(size=2, align=TRUE, linetype = 0)

pdf(paste0(figures_dir,"RdRp_95_heatmap_with_tree_with_isolates.pdf"), width = 50, height = 40, onefile = FALSE)
print(gheatmap_jag(p, MPL_otu_table_with_isolates, RdRp_95_ref_tree, color=line_colour))
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
 
  tree_plot_and_then_heatmap_print("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
                                   root = "Enterobacteria_phage_T4",
                                   plot_tree_pdf =paste0(figures_dir,"gp23_95_miseq_data_with_env_iso_and_ref_Raxml.pdf"),
                                   OTU_table = normalized_gp23_OTUs,
                                   paste0(figures_dir,"gp23_tree_heatmap_ggplot.pdf"))
 # 
#   gp23_95_ref_tree <- plot_tree_root_and_ladderize("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
#                                                    root = "Enterobacteria_phage_T4",
#                                                    "../figures/gp23_95_miseq_data_with_env_iso_and_ref_Fasttree.pdf")
#  
 
 gp23_otu_table_with_isolates <- add_isolates_to_otu_table(gp23_95_ref_tree, normalized_gp23_OTUs)
 
 #plot_otus_by_date_on_tree(gp23_otu_table_with_isolates,
  #                         gp23_95_ref_tree,
   #                        "../figures/gp23_95_miseq_data_with_env_iso_and_ref_Fasttree_with_otus_many_trees.pdf")
 

 normalized_gp23_OTUs_no_high_res <- subset(normalized_gp23_OTUs, !(rownames(normalized_gp23_OTUs) %in% high_res_vcs))
 testing_gp23_heatmap <- make_heatmap_based_on_tree_ordering(normalized_gp23_OTUs_no_high_res,
                                                            gp23_95_ref_tree, 
                                                            paste0(figures_dir,"gp23_miseq_data_heatmap_no_high_res.pdf"))
 
#### test out facetting the heatmap
 
 

 
 testing_gp23_heatmap_facet <- make_heatmap_facet(normalized_gp23_OTUs_no_high_res,
                                            gp23_95_ref_tree,
                                            gp23_OTU_annotations,
                                            "test.pdf")
 
  

 
 
 mypalette<-brewer.pal(9,"Set3")
 gp23_colour <- c("grey", mypalette, "white") 
 
 
 p1<- testing_gp23_heatmap+  theme(axis.title.y=element_blank(), axis.ticks=element_blank()) +
  theme(axis.text.y=element_blank())
 
 p2 <- make_label_for_heatmap_based_on_tree_ordering_with_annotations(normalized_gp23_OTUs_no_high_res,
                                                                      gp23_95_ref_tree, 
                                                                      paste0(figures_dir,"gp23_miseq_data_heatmap_no_high_res_label.pdf"), 
                                                                      gp23_OTU_annotations,10,
                                                                      gp23_colour) +theme(legend.position="left")
 gp1<-ggplotGrob(p1)
 gp2<-ggplotGrob(p2)  
 
 maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
 gp1$heights[2:5] <- as.list(maxHeight)
 gp2$heights[2:5] <- as.list(maxHeight)
 
 pdf(paste0(figures_dir,"gp23_heatmap_with_label.pdf"), onefile=FALSE, width = 20, height = 11)
 grid.arrange( gp2,gp1, ncol=2,widths=c(1/10,9/10))
 dev.off()
 
 
 
 normalized_gp23_OTUs_only_high_res <- subset(normalized_gp23_OTUs, rownames(normalized_gp23_OTUs) %in% high_res_vcs)
 
 
 make_heatmap_based_on_tree_ordering_high_res_view(normalized_gp23_OTUs_only_high_res,
                                                            gp23_95_ref_tree, 
                                                   paste0(figures_dir,"gp23_95_miseq_data_heatmap_only_high_res.pdf"))
 
 
 
 p <- ggtree(gp23_95_ref_tree, ladderize = FALSE,colour=line_colour)+theme_JAG_presentation()
 p <- p + geom_tiplab(size=3, align=TRUE, linetype = 0)
 
 pdf(paste0(figures_dir,"gp23_95_heatmap_with_tree_with_isolates.pdf"), width = 50, height = 40, onefile = FALSE)
 print(gheatmap_jag(p, gp23_otu_table_with_isolates, gp23_95_ref_tree, width = 20, color=line_colour))
 dev.off()
 
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
  
  OTU_heatmap <- ggplot(OTU_table_long, 
                        aes(Date,Var2)) +
   season_line +
   spring_bloom_line+
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
 
make_heatmap_based_on_tree_ordering(normalized_18s_OTUs,
                                                  S18_tree, 
                                    paste0(figures_dir,"S18_97_miseq_data_heatmap.pdf"))


testing_18s_heatmap <- make_heatmap_based_on_tree_ordering(normalized_18s_OTUs,
                                    S18_tree, 
                                    paste0(figures_dir,"S18_97_miseq_data_heatmap.pdf"))

taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)

pdf(paste0(figures_dir,"S18_97_miseq_data_heatmap_bar_label.pdf"), width=10, height=10)
label_bar_for_heatmap(taxonomy_18s,
                      S18_tree)
dev.off()



testing_18s_heatmap_with_phyla <- make_heatmap_based_on_tree_ordering_with_tax(normalized_18s_OTUs,
                                                           S18_tree, 
                                                           paste0(figures_dir,"S18_97_miseq_data_heatmap_phyla_names.pdf"),
                                                           taxonomy_18s)

## fix so it matches!!!!

normalized_18s_OTUs_no_high_res <- subset(normalized_18s_OTUs, !(rownames(normalized_18s_OTUs) %in% high_res_vcs))
testing_18s_heatmap <- make_heatmap_based_on_tree_ordering(normalized_18s_OTUs_no_high_res,
                                                           S18_tree, 
                                                           paste0(figures_dir,"S18_97_miseq_data_heatmap_no_high_res.pdf"))

make_heatmap_based_on_tree_ordering_with_tax(normalized_18s_OTUs_no_high_res,
                                    S18_tree, 
                                    paste0(figures_dir,"S18_97_miseq_data_heatmap_no_high_res_phyla_names.pdf"),
                                    taxonomy_18s)


p1<- testing_18s_heatmap +  theme(axis.title.y=element_blank(), axis.ticks=element_blank()) +
 theme(axis.text.y=element_blank())

p2 <-label_bar_for_heatmap(taxonomy_18s,
                           S18_tree) +theme(legend.position="left")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxHeight)
gp2$heights[2:5] <- as.list(maxHeight)

grid.arrange( gp2,gp1, ncol=2,widths=c(1/10,9/10))

pdf(paste0(figures_dir,"18s_heatmap_with_label.pdf"), onefile=FALSE, width = 20, height = 11)
#grid.arrange( gp2,gp1, ncol=2,widths=c(2/10,8/10))
p1

dev.off()




# ## only high res
# 
# normalized_18s_OTUs_only_high_res <- subset(normalized_18s_OTUs, rownames(normalized_18s_OTUs) %in% high_res_vcs)
# 
# 
# make_heatmap_based_on_tree_ordering_high_res_view(normalized_18s_OTUs_only_high_res,
#                                                            S18_tree, 
#                                                            "../figures/S18_97_miseq_data_heatmap_only_high_res.pdf")
# 
# 
# pdf("../figures/S18_97_miseq_data_heatmap_only_high_res_phyla_names.pdf")
# make_heatmap_based_on_tree_ordering_with_tax(normalized_18s_OTUs_only_high_res,
#                                              S18_tree, 
#                                              "../figures/S18_97_miseq_data_heatmap_only_high_res_phyla_names.pdf",
#                                              taxonomy_18s) +
#  scale_x_date(breaks = date_breaks("week"),
#               limits = c(as.Date("2011-01-15"),
#                          as.Date("2011-02-15")))
# dev.off()



plot_otus_by_date_on_tree(normalized_18s_OTUs,
                           S18_tree, 
                          paste0(figures_dir,"S18_97_miseq_data_Fasttree_with_otus_many_trees%03d.pdf"))
 

### Top 18s 100 ####
 
 normalized_18s_OTUs_top_100 <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_top_100.tsv",row.names=1)
 
# normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
# normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv",                                  row.names=1)
#  
 
#  ## want only OTUs in normalized_mpls 100 and the refs..
#  ## so want to get the names of OTUs not in the 100
#  colnames(normalized_18s_OTUs_top_100) <- gsub(".size.*.", "",colnames(normalized_18s_OTUs_top_100))
#  normalized_18s_to_remove <- normalized_18s_OTUs[,!(colnames(normalized_18s_OTUs) %in% colnames(normalized_18s_OTUs_top_100))]
# S18_97_ref_tree_top_100 <- drop.tip(S18_tree, colnames(normalized_18s_to_remove))
#  
#  plot(S18_97_ref_tree_top_100)
#  
#  plot_tree(S18_97_ref_tree_top_100,                                                         "../figures/S18_97_top_100_miseq_data_Fasttree.pdf")
#  
#  make_heatmap_based_on_tree_ordering(normalized_18s_OTUs_top_100,
#                                      S18_97_ref_tree_top_100, 
#                                      "../figures/S18_tree_top_100_heatmap_ggplot.pdf")
#  
#  plot_otus_by_date_on_tree(normalized_18s_OTUs_top_100,
#                            S18_97_ref_tree_top_100, 
#                            "../figures/S18_97_miseq_data_top_100_Fasttree_with_otus_many_trees%03d.pdf")
#  
#  ## phytos
# #  
#  normalized_18s_not_phytos_to_remove <- normalized_18s_OTUs[,!(colnames(normalized_18s_OTUs) %in% colnames(normalized_18s_OTUs_phytos))]
#  S18_97_ref_tree_phytos <- drop.tip(S18_tree, colnames(normalized_18s_not_phytos_to_remove ))
#  
#  plot(S18_97_ref_tree_phytos)
#  
#  plot_tree(S18_97_ref_tree_phytos,                                                         "../figures/S18_97_phytos_miseq_data_Fasttree.pdf")
#  
#  make_heatmap_based_on_tree_ordering(normalized_18s_OTUs_phytos,
#                                      S18_97_ref_tree_phytos, 
#                                      "../figures/S18_tree_phytos_heatmap_ggplot.pdf")
#  
#  plot_otus_by_date_on_tree(normalized_18s_OTUs_phytos,
#                            S18_97_ref_tree_phytos, 
#                            "../figures/S18_97_phytos_Fasttree_with_otus_many_trees%03d.pdf")
#  
# 
#  ## hetero
#  
#  normalized_18s_not_hetero_to_remove <- normalized_18s_OTUs[,!(colnames(normalized_18s_OTUs) %in% colnames(normalized_18s_OTUs_hetero))]
#  S18_97_ref_tree_hetero <- drop.tip(S18_tree, colnames(normalized_18s_not_hetero_to_remove ))
#  
#  plot(S18_97_ref_tree_hetero)
#  
#  plot_tree(S18_97_ref_tree_hetero,                                                         "../figures/S18_97_hetero_miseq_data_Fasttree.pdf")
#  
#  make_heatmap_based_on_tree_ordering(normalized_18s_OTUs_hetero,
#                                      S18_97_ref_tree_hetero, 
#                                      "../figures/S18_tree_hetero_heatmap_ggplot.pdf")
#  
#  plot_otus_by_date_on_tree(normalized_18s_OTUs_hetero,
#                            S18_97_ref_tree_hetero, 
#                            "../figures/S18_97_hetero_Fasttree_with_otus_many_trees%03d.pdf")
#  

 
 
#### 16s ####
 
 S16_tree <- plot_tree_and_ladderize("../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree",
                                     paste0(figures_dir,"S16_97_miseq_data_Fasttree.pdf"))
 
 S16_tree$tip.label <- gsub("_size_.*$", "",S16_tree$tip.label)
 colnames(normalized_16s_OTUs) <- gsub(".size.*.", "",colnames(normalized_16s_OTUs))
 

  make_heatmap_based_on_tree_ordering(normalized_16s_OTUs,
                                     S16_tree, 
                                     paste0(figures_dir,"S16_97_miseq_data_heatmap.pdf"))
 
  taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)
  
 testing_16s_heatmap <- make_heatmap_based_on_tree_ordering_with_tax(normalized_16s_OTUs,
                                                            S16_tree, 
                                                            paste0(figures_dir,"S16_97_miseq_data_heatmap_with_tax.pdf"),
                                                            taxonomy_16s)
 

 normalized_16s_OTUs_no_high_res <- subset(normalized_16s_OTUs, !(rownames(normalized_16s_OTUs) %in% high_res_vcs))
 testing_16s_heatmap_no_high_res <- make_heatmap_based_on_tree_ordering(normalized_16s_OTUs_no_high_res,
                                                            S16_tree, 
                                                            paste0(figures_dir,"S16_97_miseq_data_heatmap_no_high_res.pdf"))
 
 


 testing_16s_heatmap_no_high_res <- make_heatmap_based_on_tree_ordering_with_tax(normalized_16s_OTUs_no_high_res,
                                                                        S16_tree, 
                                                                        paste0(figures_dir,"S16_97_miseq_data_heatmap_no_high_res.pdf"),
                                                                        taxonomy_16s)
 pdf(paste0(figures_dir,"S16_97_miseq_data_heatmap_no_high_res_phyla_names_label.pdf"), width = 10, height = 11)
label_bar_for_heatmap(taxonomy_16s, S16_tree)
dev.off()
 

p1<- testing_16s_heatmap_no_high_res+  theme(axis.title.y=element_blank(), axis.ticks=element_blank()) + theme(axis.text.y=element_blank())

p2 <-label_bar_for_heatmap(taxonomy_16s,
                           S16_tree) +theme(legend.position="left")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxHeight = grid::unit.pmax(gp1$heights[2:5], gp2$heights[2:5])
gp1$heights[2:5] <- as.list(maxHeight)
gp2$heights[2:5] <- as.list(maxHeight)

grid.arrange( gp2,gp1, ncol=2,widths=c(1/10,9/10))

pdf(paste0(figures_dir,"16s_heatmap_with_label.pdf"), onefile=FALSE, width = 20, height = 11)
#grid.arrange( gp2,gp1, ncol=2,widths=c(2/10,8/10))
p1
dev.off()


## only high res

# normalized_16s_OTUs_only_high_res <- subset(normalized_16s_OTUs, rownames(normalized_16s_OTUs) %in% high_res_vcs)
# 
# 
# make_heatmap_based_on_tree_ordering_high_res_view(normalized_16s_OTUs_only_high_res,
#                                                            S16_tree, 
#                                                            "../figures/S16_97_miseq_data_heatmap_only_high_res.pdf")
# 


 
 plot_otus_by_date_on_tree(normalized_16s_OTUs,
                           S16_tree, 
                           paste0(figures_dir,"S16_97_miseq_data_Fasttree_with_otus_many_trees%03d.pdf"))
 


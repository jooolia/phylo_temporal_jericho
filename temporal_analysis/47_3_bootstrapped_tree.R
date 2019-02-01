## big tree bootstrap

## need to find an efficient way to do this:
library(ggplot2)
library(ape)
library(scales)
library(ggtree)
library(geiger)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)

args <- commandArgs(TRUE)
inputFile <- args[1]
## test to see if input file is given, so I can decide whether to use this argument or the orginal one. 
## test out the black manuscript file

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
                          size=2)

mid_summer_x <- as.Date("2010-09-22") - as.numeric(as.Date("2010-09-22") - as.Date("2010-06-22"))/2
mid_fall_x <- as.Date("2010-12-22") - as.numeric(as.Date("2010-12-22") - as.Date("2010-09-22"))/2
mid_winter_x <- as.Date("2011-03-22") - as.numeric(as.Date("2011-03-22") -  as.Date("2010-12-22"))/2
mid_spring_x <- as.Date("2011-06-22") - as.numeric(as.Date("2011-06-22") - as.Date("2011-03-22"))/2

season_text <- annotate("text",
                        x = c(mid_summer_x,
                              mid_fall_x,
                              mid_winter_x,
                              mid_spring_x),
                        y = 1.05,
                        label = c("Summer",
                                  "Fall",
                                  "Winter",
                                  "Spring"))

spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour = "green",
                                size = 2)
date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")))



#http://www.r-phylo.org/wiki/HowTo/DataTreeManipulation#Is_there_a_shorthand_way_to_refer_to_a_specific_list_of_taxa_.28for_example.2C_all_members_of_a_particular_group.29.3F

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv",
                                  row.names="VC_number")
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv",
                                  row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",
                                   row.names="VC_number")

normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv",
                                  row.names="VC_number")

### maybe change these to prop abundance???

proportional_MPL <- as.data.frame(prop.table(as.matrix(normalized_MPL_OTUs),
                                             margin=1))
proportional_gp23 <- as.data.frame(prop.table(as.matrix(normalized_gp23_OTUs),
                                              margin=1))

### Get relevant VC numbers ====
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv",
                         row.names=1)
## Reformat date
Jericho_data$Date <- as.Date(Jericho_data$Date)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv",
                             nrows=61)

high_res_vcs <- c(1198,
                  1199,
                  1200,
                  1201,
                  1202)


## want Jericho_dates that are not in MPL table and high res
Jericho_no_high_res_vcs <- Jericho_data$VC_number[!(Jericho_data$VC_number %in% high_res_vcs)]

## want to annotate missing data days. 
## want to annotate those days using star
missing_MPL_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(proportional_MPL))]
missing_MPL_dates <- Jericho_data$Date[match(missing_MPL_samples,
                                             Jericho_data$VC_number)]

missing_gp23_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(proportional_gp23))]
missing_gp23_dates <- Jericho_data$Date[match(missing_gp23_samples,
                                              Jericho_data$VC_number)]

plot_tree_root_and_ladderize <- function (tree_file,
                                          root) {
  tree_amplicon <- read.tree(tree_file)
  tree_amplicon <- root(tree_amplicon,
                        root,
                        resolve.root = TRUE)
  tree_amplicon  <- ladderize(tree_amplicon)
  plot(tree_amplicon,
       cex = 0.3)
  # exporting tree.
  write.tree(tree_amplicon, "../results/temp_tree.tree")
  # re-import to be able to match up the tip labels with the OTU table. 
  reordered_tree <- read.tree("../results/temp_tree.tree")
  return(reordered_tree)
}


## full tree with bootstrap
mypalette<-brewer.pal(7,
                      "Dark2")

RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762")


RdRp_95_ref_tree_with_bootstrap <- apeBoot(RdRp_95_ref_tree,
                                           RdRp_95_ref_tree$node.label)

RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap[as.numeric(as.character(RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap)) <50] <- NA
RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap[as.character(RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap) == "Root"] <- NA
RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap <- as.numeric(as.character(RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap))


tree2 <- groupClade(RdRp_95_ref_tree_with_bootstrap,
                    node = c(1768,
                             1764,
                             1760,
                             1701,
                             1613,
                             1700,
                             1450,
                             1017))

bootstrap_tree <- ggtree(tree2,
                         aes(color = group),
                         ladderize = FALSE,
                         colour = line_colour)+
  scale_color_manual(values = c("black",
                                "grey",
                                mypalette[1:7]),
                     guide = FALSE
  )+
  scale_y_continuous(labels = NULL)+
  scale_x_continuous(labels = NULL)+
  theme_JAG_presentation()+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  geom_treescale(x = 0.05,
                 y = 2,
                 offset = 2.5)+
  geom_point(aes(fill = bootstrap,
                 alpha = bootstrap),
             size = 3,
             na.rm = TRUE,
             shape = 21,
             colour = "white")+
  scale_fill_continuous(low = 'grey',
                        high = 'black',
                        guide = "legend")+
  scale_alpha_continuous(range = c(0.7,
                                   1),
                         na.value = 0)+
  theme(legend.position = "right")+
  theme(plot.margin = unit(c(1,5,1,1),
                           "cm"))


pdf(paste0(figures_dir,"MPL_tree_with_bootstrap.pdf"),
    width = 12,
    height = 17,
    onefile = FALSE)

bootstrap_tree

dev.off()

## add on
bootstrap_with_points <- bootstrap_tree+ 
  geom_point2(aes(subset = (!grepl("OTU|Culley|0\\.[:digit:]+|^[0-9]{2}$",
                                   p$data$label) & isTip)),
              color = "grey30",
              shape = 18,
              size = 2.2)+
  geom_point2(aes(subset = ( grepl("Culley",
                                   p$data$label) & isTip)),
              color = "green",
              size = 2,
              shape = 24)+
  geom_point2(aes(subset = ( grepl("OTU", p$data$label) & isTip)),
              shape = 23,
              fill = "steelblue",
              colour = NA,
              size = 1 )

bootstrap_with_points_and_clade_labels <- bootstrap_with_points+
  geom_cladelabel(node = 1768,
                  "A",
                  angle = 0,
                  offset.text = 0.12,
                  barsize = 4, 
                  fontsize = 6,
                  offset = 2.0,
                  color = "grey") + 
  geom_cladelabel(node = 1764,
                  "B",
                  angle = 0,
                  offset.text = 0.12,
                  barsize = 4,
                  fontsize = 6,
                  offset = 1.9242,
                  color = mypalette[1]) + 
  geom_cladelabel(node = 1760,
                  "C",
                  angle = 0,
                  color = mypalette[2],
                  barsize = 4,
                  fontsize = 6,
                  offset = 1.546,
                  offset.text = -0.3) + 
  geom_cladelabel(node = 1701,
                  "D",
                  angle = 0,
                  barsize = 4,
                  fontsize = 6,
                  offset = 0.325,
                  offset.text = -0.3,
                  color = mypalette[3]) + 
  geom_cladelabel(node = 1613,
                  "E",
                  angle = 0,
                  barsize = 4,
                  fontsize = 6,
                  offset = 1.482,
                  offset.text = 0.12,
                  color = mypalette[4]) + 
  geom_cladelabel(node = 1700,
                  "F",
                  angle = 0,
                  offset.text = 0.12,
                  barsize = 4,
                  fontsize = 6,
                  offset = 1.355,
                  color = mypalette[5]) + 
  geom_cladelabel(node = 1450,
                  "G",
                  angle = 0,
                  offset.text = 0.12,
                  barsize = 4,
                  fontsize = 6,
                  offset = 0.78,
                  color = mypalette[6]) +
  geom_cladelabel(node = 1017,
                  "H",
                  angle = 0,
                  offset.text = 0.12,
                  barsize = 4,
                  fontsize = 6,
                  offset = 1.065,
                  color = mypalette[7])

names_tips <- c("Isolates",
                "Culley et al 2003 and\n Culley and Steward, 2007",
                "OTUs")
cols_tips <- c("grey30",
               "green",
               "steelblue")

length_cols <-c(rep(1,
                    length(cols_tips)))  

df_cols <- data.frame(names_tips,
                      cols_tips,length_cols)

colouring_bars <- as.character(df_cols$cols_tips)
names(colouring_bars) <- df_cols$names_tips

labeling_tree <- ggplot(df_cols,
                        aes(x = names_tips,
                            y = length_cols,
                            group = names_tips))+
  geom_bar(stat = "identity",
           aes(fill = names_tips))+
  scale_fill_manual(values = colouring_bars,
                    name = "")



## get only legend

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

## can't get legend to work
length_cols <-c(rep(1,
                    length(cols_tips)))  

df_cols <- data.frame(names_tips,
                      cols_tips,length_cols)

colouring_bars <- as.character(df_cols$cols_tips)
names(colouring_bars) <- df_cols$names_tips

labeling_tree <- ggplot(df_cols,
                        aes(x = names_tips,
                            y = length_cols,
                            group = names_tips))+
  geom_bar(stat = "identity",
           aes(fill = names_tips))+
  scale_fill_manual(values = colouring_bars,
                    name = "")

legend <- g_legend(labeling_tree)


pdf(paste0(figures_dir,"MPL_tree_with_bootstrap_and_points.pdf"),
    width = 12,
    height = 17,
    onefile = FALSE)

grid.arrange(bootstrap_with_points_and_clade_labels,
             legend,
             ncol = 2,
             widths = c(19/20,1/20))
dev.off()


################# gp23 

mypalette<-brewer.pal(9,
                      "Set3")
colours_for_groups <- colorRampPalette(mypalette)(15)

gp23_95_ref_tree <- plot_tree_root_and_ladderize("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
                                                 root = "Enterobacteria_phage_T4")

###### legend stuff
names_tips <- c("Isolates",
                "Chow and Fuhrman, 2012",
                "Filee et al, 2005",
                "Jia et al, 2007",
                "Lopez-Bueno et al, 2009",
                "Bellas and Anesio, 2013",
                "Knapik and Prentice, 2011",
                "Liu et al, 2012",
                "Butina et al, 2013",
                "Comeau et al, 2010",
                "Sandaa and Kristiansen, 2007",
                "Mabizela and Litthauer, 2009",
                "OTUs from this study")
cols_tips <- c("grey30",
               "green",
               "brown",
               "lightblue",
               "darkgreen",
               "purple",
               "turquoise",
               "lightgrey",
               "beige",
               "darkred",
               "pink",
               "darkslateblue",
               "steelblue")


#####

length_cols <-c(rep(1,
                    length(cols_tips)))  

df_cols <- data.frame(names_tips,
                      cols_tips,length_cols)

colouring_bars <- as.character(df_cols$cols_tips)
names(colouring_bars) <- df_cols$names_tips

labeling_tree <- ggplot(df_cols,
                        aes(x = names_tips,
                            y = length_cols,
                            group = names_tips))+
  geom_bar(stat = "identity",
           aes(fill = names_tips))+
  scale_fill_manual(values = colouring_bars,
                    name = "")

## get only legend

length_cols <-c(rep(1,
                    length(cols_tips)))  

df_cols <- data.frame(names_tips,
                      cols_tips,length_cols)

colouring_bars <- as.character(df_cols$cols_tips)
names(colouring_bars) <- df_cols$names_tips

labeling_tree <- ggplot(df_cols,
                        aes(x = names_tips,
                            y = length_cols,
                            group = names_tips))+
  geom_bar(stat = "identity",
           aes(fill = names_tips))+
  scale_fill_manual(values = colouring_bars,
                    name = "")

legend <- g_legend(labeling_tree)

#####



pdf(paste0(figures_dir,"gp23_95_env_OTUs_legend.pdf"), width = 30, height = 40, onefile = FALSE)
grid.arrange(legend)
dev.off()

colour_tips_of_trees <- function (tree,
                                  node_to_extract) {
  new_subtree <- extract.clade(tree,
                               node = node_to_extract)
  #print(new_subtree$tip.labels)
  
  new_tree_bootstrap <- gp23_95_ref_tree$node.label[match(new_subtree$tip.label,
                                                          gp23_95_ref_tree$tip.label)][-length(new_subtree$tip.label)]
  
  new_subtree <- apeBoot(new_subtree,
                         new_tree_bootstrap)
  
  new_subtree@bootstrap$bootstrap[as.numeric(as.character(new_subtree@bootstrap$bootstrap)) <50] <- NA
  new_subtree@bootstrap$bootstrap[as.character(new_subtree@bootstrap$bootstrap) == "Root"] <- NA
  new_subtree@bootstrap$bootstrap <- as.numeric(as.character(new_subtree@bootstrap$bootstrap))
  
  p <- ggtree(new_subtree, ladderize = FALSE,colour=line_colour) +
    theme_JAG_presentation()
  
  colour_stuff <- p+
    geom_tiplab(aes(subset=(!grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                   p$data$label) & isTip)),
                color="grey30",
                size=2.2) +
    geom_point2(aes(subset=( grepl("Chow",
                                   p$data$label) & isTip)),
                color="green",
                size=2,
                shape = 24)+
    geom_point2(aes(subset=( grepl("Filee",
                                   p$data$label) & isTip)),
                color="brown",
                size=2,
                shape = 24)+
    geom_point2(aes(subset=( grepl("Jia",
                                   p$data$label) & isTip)),
                color="lightblue",
                size=2,
                shape = 24)+ 
    geom_point2(aes(subset=( grepl("Lopez-Bueno|LopezBueno",
                                   p$data$label) & isTip)),
                color="darkgreen",
                size=2,
                shape = 24)+
    geom_point2(aes(subset=( grepl("Bellas",
                                   p$data$label) & isTip)),
                color="purple",
                size=2,
                shape = 24)+
    geom_point2(aes(subset=( grepl("Knapik",
                                   p$data$label) & isTip)),
                color="turquoise",
                size=2,
                shape = 24)+
    geom_point2(aes(subset=( grepl("Liu",
                                   p$data$label) & isTip)),
                color="lightgrey",
                size=2,
                shape = 24)+
    geom_point2(aes(subset=( grepl("Butina",
                                   p$data$label) & isTip)),
                color="beige",
                size=2,
                shape = 24)+ 
    geom_point2(aes(subset=( grepl("Comeau",
                                   p$data$label) & isTip)),
                color="darkred",
                size=2,
                shape = 24)+ 
    geom_point2(aes(subset=( grepl("Sandaa",
                                   p$data$label) & isTip)),
                color="pink",
                size=2,
                shape = 24)+
    geom_point2(aes(subset=( grepl("Mabiz",
                                   p$data$label) & isTip)),
                color="darkslateblue",
                size=2,
                shape = 24)+
    geom_point2(aes(subset=( grepl("OTU", p$data$label) & isTip)), shape=23, fill="steelblue", size=1 )+
    geom_tiplab(aes(subset=(grepl("OTU", p$data$label) & isTip)),
                color="black",
                size=1.8)+
    geom_point(aes(colour=bootstrap,alpha=bootstrap),
               size=3,na.rm=TRUE)+
    scale_colour_continuous(low='grey', high='black')+
    scale_alpha_continuous(range = c(0.7, 1), na.value=0)
  return(colour_stuff)
}


gp23_95_ref_tree_with_bootstrap <- apeBoot(gp23_95_ref_tree,
                                           gp23_95_ref_tree$node.label)

gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap[as.numeric(as.character(gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap)) <50] <- NA
gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap[as.character(gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap) == "Root"] <- NA
gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap <- as.numeric(as.character(gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap))


tree2 <- groupClade(gp23_95_ref_tree_with_bootstrap,
                    node=c(4472,
                           4363,
                           4252,
                           4217,
                           4147,
                           3905,
                           3295,
                           2955,
                           2295))

bootstrap_tree <- ggtree(tree2,
                         aes(color=group),
                         ladderize = FALSE,
                         colour=line_colour)+
  scale_color_manual(values=c("black","grey",
                              mypalette[1:8]),
                     guide = FALSE
  )+
  scale_y_continuous(labels=NULL)+
  scale_x_continuous(labels=NULL)+
  theme_JAG_presentation()+
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank())+
  geom_treescale(x=0.05,
                 y=2,
                 offset=2.5)+
  geom_point(aes(fill=bootstrap,
                 alpha=bootstrap),
             size=3,
             na.rm=TRUE,
             shape = 21,
             colour="white")+
  scale_fill_continuous(low='grey',
                        high='black',
                        guide = "legend")+
  scale_alpha_continuous(range = c(0.7, 1),
                         na.value=0)+
  theme(legend.position="right")+
  theme(plot.margin = unit(c(1,5,1,1),
                           "cm"))


pdf(paste0(figures_dir,"gp23_tree_with_bootstrap.pdf"), width = 12, height = 17, onefile = FALSE)
bootstrap_tree
dev.off()


## add on
bootstrap_with_points <- bootstrap_tree+ 
  geom_point2(aes(subset=(!grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                 p$data$label) & isTip)),
              color="grey30",
              shape=18,
              size=2.2) +
  geom_point2(aes(subset=( grepl("Chow",
                                 p$data$label) & isTip)),
              color="green",
              shape=23,
              size=2)+
  geom_point2(aes(subset=( grepl("Filee",
                                 p$data$label) & isTip)),
              color="brown",
              shape=23,
              size=2)+
  geom_point2(aes(subset=( grepl("Jia",
                                 p$data$label) & isTip)),
              color="lightblue",
              shape=23,
              size=2)+ 
  geom_point2(aes(subset=( grepl("Lopez-Bueno|LopezBueno",
                                 p$data$label) & isTip)),
              color="darkgreen",
              shape=23,
              size=2)+
  geom_point2(aes(subset=( grepl("Bellas",
                                 p$data$label) & isTip)),
              color="purple",
              shape=23,
              size=2)+
  geom_point2(aes(subset=( grepl("Knapik",
                                 p$data$label) & isTip)),
              color="turquoise",
              shape=23,
              size=2)+
  geom_point2(aes(subset=( grepl("Liu",
                                 p$data$label) & isTip)),
              color="lightgrey",
              shape=23,
              size=2)+
  geom_point2(aes(subset=( grepl("Butina",
                                 p$data$label) & isTip)),
              color="beige",
              shape=23,
              size=2)+ 
  geom_point2(aes(subset=( grepl("Comeau",
                                 p$data$label) & isTip)),
              color="darkred",
              shape=23,
              size=2)+ 
  geom_point2(aes(subset=( grepl("Sandaa",
                                 p$data$label) & isTip)),
              color="pink",
              shape=23,
              size=2)+
  geom_point2(aes(subset=( grepl("Mabiz",
                                 p$data$label) & isTip)),
              color="darkslateblue",
              shape=23,
              size=2)+
  geom_point2(aes(subset=( grepl("OTU",
                                 p$data$label) & isTip)),
              shape=23,
              fill="steelblue",
              colour=NA,
              size=1 )



bootstrap_with_points_and_clade_labels <-  bootstrap_with_points+geom_cladelabel(node=4472,
                                                                                 "A",
                                                                                 offset = 2.415,
                                                                                 barsize=2,
                                                                                 fontsize=6,
                                                                                 angle=0,
                                                                                 colour="grey") + 
  geom_cladelabel(node=4363,
                  "B",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  #offset.text=-0.1,
                  offset = 1.5,
                  color= mypalette[1]) +
  geom_cladelabel( node=4252,
                   "C",
                   barsize=2,
                   fontsize=6,
                   angle=0,
                   #offset.text=0.12,
                   offset = 1.425,
                   color= mypalette[2]) + 
  geom_cladelabel( node=4217,
                   "D",
                   barsize=2,
                   fontsize=6,
                   angle=0,
                   #offset.text=-0.1, 
                   color= mypalette[3]) + 
  geom_cladelabel( node=4147,
                   "E",  
                   barsize=2,
                   fontsize=6,
                   angle=0,
                   #offset.text=0.12,
                   offset = 1.4,
                   color= mypalette[4]) + 
  geom_cladelabel( node=3905,
                   "F",
                   barsize=2,
                   fontsize=6,
                   angle=0,
                   offset=1.165,
                   # offset.text=-0.1,
                   color= mypalette[5]) +
  geom_cladelabel( node=3295,
                   "G",
                   barsize=2,
                   fontsize=6,
                   angle=0,
                   offset=1.3,
                   # offset.text=0.12,
                   color= mypalette[6]) + 
  geom_cladelabel(node=2955,
                  "H", 
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=0.62,
                  color= mypalette[7]) + 
  geom_cladelabel(node=2295, 
                  "I",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=0.89,
                  # offset.text=0.12,
                  color= mypalette[8]) 


pdf(paste0(figures_dir,
           "gp23_tree_with_bootstrap_and_points.pdf"),
    width = 12,
    height = 17,
    onefile = FALSE)
grid.arrange(bootstrap_with_points_and_clade_labels,
             legend,
             ncol=2,
             widths=c(19/20,
                      1/20))
dev.off()


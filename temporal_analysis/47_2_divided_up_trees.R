## divided up trees. 

## need to find an efficient way to do this:
library(ggplot2)
library(ape)
library(scales)
library(ggtree)
library(geiger)
library(plyr)
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

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv",
                                  row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",
                                   row.names="VC_number")

proportional_MPL <- as.data.frame(prop.table(as.matrix(normalized_MPL_OTUs),
                                             margin=1))
proportional_gp23 <- as.data.frame(prop.table(as.matrix(normalized_gp23_OTUs),
                                              margin=1))

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

colour_tips_of_trees <- function (tree,
                                  node_to_extract) {
  new_subtree <- extract.clade(tree,
                               node = node_to_extract)
  
  new_tree_bootstrap <- RdRp_95_ref_tree$node.label[match(new_subtree$tip.label,
                                                          RdRp_95_ref_tree$tip.label)][-length(new_subtree$tip.label)]
  
  new_subtree <- apeBoot(new_subtree, 
                         new_tree_bootstrap)
  
  
  new_subtree@bootstrap$bootstrap[as.numeric(as.character(new_subtree@bootstrap$bootstrap)) <50] <- NA
  new_subtree@bootstrap$bootstrap[as.character(new_subtree@bootstrap$bootstrap) == "Root"] <- NA
  new_subtree@bootstrap$bootstrap <- as.numeric(as.character(new_subtree@bootstrap$bootstrap))
  
  
  p <- ggtree(new_subtree,
              ladderize = FALSE,
              colour = line_colour)
  
  colour_stuff <- p+
    geom_tiplab(aes(subset = (!grepl("OTU|Culley|0\\.[:digit:]+|^[0-9]{2}$",
                                     p$data$label) & isTip)),
                color = "grey30",
                size = 2.2) +
    geom_point2(aes(subset = ( grepl("Culley",
                                     p$data$label) & isTip)),
                color = "green",
                size = 2,
                shape = 24)+
    geom_point2(aes(subset = ( grepl("OTU",
                                     p$data$label) & isTip)),
                shape = 23,
                fill = "steelblue",
                size = 1 )+
    geom_tiplab(aes(subset = (grepl("OTU",
                                    p$data$label) & isTip)),
                color = "black",
                size = 1.8)+
    geom_point(aes(colour = bootstrap,
                   alpha = bootstrap),
               size = 3,
               na.rm = TRUE)+
    scale_colour_continuous(low = 'grey',
                            high = 'black')+
    scale_alpha_continuous(range = c(0.7, 1),
                           na.value = 0)
  
  return(colour_stuff)
}

## get only legend

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

names_tips <- c("Isolates",
                "Culley_2003_and_2007",
                "OTUs")
cols_tips <- c("grey30",
               "green",
               "steelblue")


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

pdf(paste0(figures_dir,"RdRp_95_env_OTUs_legend.pdf"), 
    width = 8,
    height = 10,
    onefile = FALSE)

grid.arrange(legend)

dev.off()


RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762")


tree <- groupClade(RdRp_95_ref_tree,
                   node = c(1768,
                            1764,
                            1760,
                            1701,
                            1613,
                            1700,
                            1450,
                            1017)
)

ggtree(tree, ladderize = FALSE) +
  geom_text2(aes(subset = !isTip,
                 label = node),
             hjust = -.3) +
  geom_tiplab()

mypalette<-brewer.pal(7,
                      "Dark2")

p2a <- ggtree(tree,
              aes(color = group),
              ladderize = FALSE)+
  scale_color_manual(values = c(line_colour,
                                "grey",
                                mypalette[1:7]),
                     guide = FALSE)+
  scale_y_continuous(labels = NULL)+
  scale_x_continuous(labels = NULL)+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) 

# ## top part
tips_to_drop <- tips(RdRp_95_ref_tree,
                     node = 1006)
top_tree_part <- drop.tip(RdRp_95_ref_tree,
                          tips_to_drop)

## bootstrap values for the top part of the tree. have to exclude the last because it is no longer a node
top_tree_bootstrap <- RdRp_95_ref_tree$node.label[match(top_tree_part$tip.label,
                                                        RdRp_95_ref_tree$tip.label)][-length(top_tree_part$tip.label)]

top_tree_part <- apeBoot(top_tree_part,
                         top_tree_bootstrap)

top_tree_part@bootstrap$bootstrap[as.numeric(as.character(top_tree_part@bootstrap$bootstrap)) <50] <- NA
top_tree_part@bootstrap$bootstrap[as.character(top_tree_part@bootstrap$bootstrap) == "Root"] <- NA
top_tree_part@bootstrap$bootstrap <- as.numeric(as.character(top_tree_part@bootstrap$bootstrap))



top_tree <- ggtree(top_tree_part,
                   ladderize = FALSE,
                   colour = line_colour)+
  #geom_tiplab(size = 3)+
  geom_treescale(x = 0.05,
                 y = 2,
                 offset = 2.5)+
  geom_tiplab(aes(subset = (!grepl("OTU|Culley|0\\.[:digit:]+|^[0-9]{2}$",
                                   p$data$label) & isTip)),
              color = "grey30",
              size = 2.2) +
  geom_point2(aes(subset = (grepl("Culley",
                                  p$data$label) & isTip)),
              color = "green",
              size = 2)+
  geom_point2(aes(subset = (grepl("OTU",
                                  p$data$label) & isTip)),
              shape = 23,
              fill = "steelblue",
              size = 1)+
  geom_point(aes(colour = bootstrap,
                 alpha = bootstrap),
             size = 3,
             na.rm = TRUE)+
  scale_colour_continuous(low = 'grey',
                          high = 'black')+
  scale_alpha_continuous(range = c(0.7, 1),
                         na.value = 0)+
  theme(plot.margin = unit(c(1,
                             5,
                             1,
                             1),
                           "cm"))

tree_clade_A <- colour_tips_of_trees(tree,
                                     1768)+
  geom_treescale()

small_tree_A_highlighted <- p2a+
  geom_cladelabel(node = 1768,
                  "A",
                  offset = 1.0,
                  barsize = 1,
                  fontsize = 4,
                  angle = 0,
                  colour = "grey")+
  theme(plot.margin = unit(c(1,1,1,1),
                           "cm"))


tree_clade_B <- colour_tips_of_trees(tree,
                                     1764)+
  geom_treescale()

small_tree_B_highlighted <- p2a+
  geom_cladelabel(node = 1764,
                  "B",
                  offset = 1.0,
                  barsize = 1,
                  fontsize = 4,
                  angle = 0,
                  color = mypalette[1])+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

tree_clade_C <- colour_tips_of_trees(tree,
                                     1760)+
  geom_treescale()

small_tree_C_highlighted <- p2a+
  geom_cladelabel(node = 1760,
                  "C",
                  offset = 1.0,
                  barsize = 1,
                  fontsize = 4,
                  angle = 0,
                  color = mypalette[2])+
  theme(plot.margin = unit(c(1,1,1,1),
                           "cm"))

tree_clade_D <- colour_tips_of_trees(tree, 
                                     1701)+
  geom_treescale()

small_tree_D_highlighted <- p2a+
  geom_cladelabel(node = 1701,
                  "D",
                  offset = 0.3,
                  barsize = 1,
                  fontsize = 4,
                  angle = 0,
                  color = mypalette[3])+
  theme(plot.margin = unit(c(1,1,1,1),
                           "cm"))


tree_clade_E <- colour_tips_of_trees(tree,
                                     1613)+
  geom_treescale()

small_tree_E_highlighted <- p2a+
  geom_cladelabel(node = 1613,
                  "E",
                  offset = 1.0,
                  barsize = 1,
                  fontsize = 4,
                  angle = 0,
                  color = mypalette[4])+
  theme(plot.margin = unit(c(1,1,1,1),
                           "cm"))


tree_clade_F <- colour_tips_of_trees(tree,
                                     1700)+
  geom_treescale()

small_tree_F_highlighted <- p2a+
  geom_cladelabel(node = 1700,
                  "F",
                  offset = 1.0,
                  barsize = 1,
                  fontsize = 4,
                  angle = 0,
                  color = mypalette[5])+
  theme(plot.margin = unit(c(1,1,1,1),
                           "cm"))

tree_clade_G <- colour_tips_of_trees(tree,
                                     1450)+
  geom_treescale()

small_tree_G_highlighted <- p2a+
  geom_cladelabel(node = 1450,
                  "G",
                  offset = 0.5,
                  barsize = 1,
                  fontsize = 4,
                  angle = 0,
                  color = mypalette[6])+
  theme(plot.margin = unit(c(1,1,1,1),
                           "cm"))

tree_clade_H <- colour_tips_of_trees(tree,
                                     1017)+
  geom_treescale()

small_tree_H_highlighted <- p2a+
  geom_cladelabel(node = 1017,
                  "H",
                  offset = 0.8,
                  barsize = 1,
                  fontsize = 4,
                  angle = 0,
                  color = mypalette[7])+
  theme(plot.margin = unit(c(1,1,1,1),
                           "cm"))



# 
# ## make one with all of the stuff on it for the manuscript
# ## tips getting cut off...in ggtree
pdf(paste0(figures_dir,"MPL_zoomed_in_tree%0d.pdf"),
    width = 12,
    height = 17,
    onefile = FALSE)
plot_grid(p2a,
          top_tree,
          small_tree_A_highlighted,
          tree_clade_A,
          ncol=2,
          nrow=2,
          rel_heights = c(1.5,1),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_B_highlighted,
          tree_clade_B,
          small_tree_C_highlighted,
          tree_clade_C,
          ncol = 2,
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_D_highlighted,
          tree_clade_D,
          small_tree_E_highlighted,
          tree_clade_E,
          ncol = 2,
          nrow = 2,
          rel_heights = c(1,1.5),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_F_highlighted,
          tree_clade_F,
          small_tree_G_highlighted,
          tree_clade_G,
          ncol = 2,
          rel_heights = c(1,2),
          rel_widths = c(1/8,7/8),
          label_size = 20)


plot_grid(small_tree_H_highlighted,
          tree_clade_H,
          ncol = 2,
          rel_widths = c(1/8,7/8),
          label_size = 20)

dev.off()


########################## gp23 #############

mypalette<-brewer.pal(9,
                      "Set3")
colours_for_groups <- colorRampPalette(mypalette)(15)

gp23_95_ref_tree <- plot_tree_root_and_ladderize("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
                                                 root = "Enterobacteria_phage_T4")

tree <- groupClade(gp23_95_ref_tree, node=c(4472,4363, 4252,4217,4147,3905, 3295,2955,2295))

p2a <- ggtree(tree,
              aes(color=group),
              ladderize = FALSE)+
  scale_color_manual(values=c(line_colour,"grey",
                              mypalette[1:8]),
                     guide = FALSE
                     #,labels=c("Ref","A","B","C","D","E","F","G","H", "I")
  )+
  scale_y_continuous(labels=NULL)+
  scale_x_continuous(labels=NULL)+
  theme_JAG_presentation()+
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        axis.line = element_blank()
  ) 


colour_tips_of_trees <- function (tree, node_to_extract) {
  new_subtree <- extract.clade(tree, node = node_to_extract)
  #print(new_subtree$tip.labels)
  
  new_tree_bootstrap <- gp23_95_ref_tree$node.label[match(new_subtree$tip.label,gp23_95_ref_tree$tip.label)][-length(new_subtree$tip.label)]
  
  new_subtree <- apeBoot(new_subtree, new_tree_bootstrap)
  
  new_subtree@bootstrap$bootstrap[as.numeric(as.character(new_subtree@bootstrap$bootstrap)) <50] <- NA
  new_subtree@bootstrap$bootstrap[as.character(new_subtree@bootstrap$bootstrap) == "Root"] <- NA
  new_subtree@bootstrap$bootstrap <- as.numeric(as.character(new_subtree@bootstrap$bootstrap))
  
  p <- ggtree(new_subtree, ladderize = FALSE,colour=line_colour) +
    theme_JAG_presentation()
  
  colour_stuff <- p+ geom_tiplab(aes(subset=(!grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz", p$data$label) & isTip)),
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



ggtree(tree, ladderize = FALSE) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

## top part
tips_to_drop <- tips(tree, node=2288) 
top_tree_part <- drop.tip(tree, tips_to_drop)


## bootstrap values for the top part of the tree. have to exclude the last because it is no longer a node
top_tree_bootstrap <- gp23_95_ref_tree$node.label[match(top_tree_part$tip.label,gp23_95_ref_tree$tip.label)][-length(top_tree_part$tip.label)]

top_tree_part <- apeBoot(top_tree_part, top_tree_bootstrap)

top_tree_part@bootstrap$bootstrap[as.numeric(as.character(top_tree_part@bootstrap$bootstrap)) <50] <- NA
top_tree_part@bootstrap$bootstrap[as.character(top_tree_part@bootstrap$bootstrap) == "Root"] <- NA
top_tree_part@bootstrap$bootstrap <- as.numeric(as.character(top_tree_part@bootstrap$bootstrap))




# grep("OTU", top_tree_part$tip.label, invert = TRUE, value = TRUE)
# 
# groupOTU(p2a, c(top_tree_part$tip.label))
# p2a+ geom_text("test")
#  geom_segment2(mapping = top_tree_part$tip.label, data = p2a)

top_tree <- ggtree(top_tree_part, ladderize = FALSE,colour=line_colour)+
  #geom_tiplab(size=3)+
  geom_treescale(x=0.05,y=2, offset=2.5)+geom_tiplab(aes(subset=(!grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz", p$data$label) & isTip)),
                                                     color="grey30",
                                                     size=2.2) +
  geom_point2(aes(subset=( grepl("Chow",
                                 p$data$label) & isTip)),
              color="green",
              size=2)+
  geom_point2(aes(subset=( grepl("Filee",
                                 p$data$label) & isTip)),
              color="brown",
              size=2)+
  geom_point2(aes(subset=( grepl("Jia",
                                 p$data$label) & isTip)),
              color="lightblue",
              size=2)+ 
  geom_point2(aes(subset=( grepl("Lopez-Bueno|LopezBueno",
                                 p$data$label) & isTip)),
              color="darkgreen",
              size=2)+
  geom_point2(aes(subset=( grepl("Bellas",
                                 p$data$label) & isTip)),
              color="purple",
              size=2)+
  geom_point2(aes(subset=( grepl("Knapik",
                                 p$data$label) & isTip)),
              color="turquoise",
              size=2)+
  geom_point2(aes(subset=( grepl("Liu",
                                 p$data$label) & isTip)),
              color="lightgrey",
              size=2)+
  geom_point2(aes(subset=( grepl("Butina",
                                 p$data$label) & isTip)),
              color="beige",
              size=2)+ 
  geom_point2(aes(subset=( grepl("Comeau",
                                 p$data$label) & isTip)),
              color="darkred",
              size=2)+ 
  geom_point2(aes(subset=( grepl("Sandaa",
                                 p$data$label) & isTip)),
              color="pink",
              size=2)+
  geom_point2(aes(subset=( grepl("Mabiz",
                                 p$data$label) & isTip)),
              color="darkslateblue",
              size=2)+
  geom_point2(aes(subset=( grepl("OTU", p$data$label) & isTip)), shape=23, fill="steelblue", size=1 )+
  geom_tiplab(aes(subset=(grepl("OTU", p$data$label) & isTip)),
              color="black",
              size=1.8)+
  geom_point(aes(colour=bootstrap,alpha=bootstrap),
             size=3,na.rm=TRUE)+
  scale_colour_continuous(low='grey', high='black')+
  scale_alpha_continuous(range = c(0.7, 1), na.value=0)+
  theme_JAG_presentation()+
  theme(plot.margin = unit(c(1,5,1,1), "cm"))

#par(mar=c(5,4,4,16))
tree_clade_A <- colour_tips_of_trees(tree, 4472)+
  geom_treescale()

small_tree_A_highlighted <- p2a+
  geom_cladelabel(node=4472,
                  "A",
                  offset = 1.8,
                  barsize=1,
                  fontsize=4,
                  angle=0,
                  colour="grey"
  )

tree_clade_B <- colour_tips_of_trees(tree,4363)+
  geom_treescale(x=0.05,y=2, offset=2.5)

small_tree_B_highlighted <- p2a +
  geom_cladelabel(node=4363, 
                  "B",
                  barsize=1,
                  fontsize=4,
                  angle=0,
                  #offset.text=-0.1,
                  offset = 0.8,
                  color= mypalette[1]
  )+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

tree_clade_C <- colour_tips_of_trees(tree,4252)+
  geom_treescale(x=0.05,y=2, offset=2.5)+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

small_tree_C_highlighted <- p2a +
  geom_cladelabel( node=4252, "C",
                   barsize=1,
                   fontsize=4,
                   angle=0,
                   #offset.text=0.12,
                   offset = 0.8,
                   color= mypalette[2])+
  theme(plot.margin = unit(c(1,1.2,1,1), "cm"))

tree_clade_D <- colour_tips_of_trees(tree,4217)+
  geom_treescale(x=2.2,y=2, offset=1.5)+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

small_tree_D_highlighted <- p2a +
  geom_cladelabel( node=4217, "D", 
                   barsize=1,
                   fontsize=4,
                   angle=0,
                   #offset.text=-0.1, 
                   color= mypalette[3])+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

tree_clade_E <- colour_tips_of_trees(tree,4147)+
  geom_treescale(x=0.8,y=2, offset=1.5)

small_tree_E_highlighted <- p2a +
  geom_cladelabel( node=4147, 
                   "E",  
                   barsize=1,
                   fontsize=4,
                   angle=0,
                   #offset.text=0.12,
                   offset = 0.7,
                   color= mypalette[4])

tree_clade_F <- colour_tips_of_trees(tree,3905)+
  geom_treescale(offset=5)

small_tree_F_highlighted <- p2a +
  geom_cladelabel( node=3905, "F",
                   barsize=1,
                   fontsize=4,
                   angle=0,
                   offset=0.6,
                   # offset.text=-0.1,
                   color= mypalette[5])

## needs to be split

clade_g <- extract.clade(tree, 3295)
ggtree(clade_g, ladderize = FALSE) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

tips_to_drop <- tips(clade_g, node=671) 
clade_g_b <- drop.tip(clade_g, tips_to_drop)
plot(clade_g_b)
plot(clade_g)
# clade_g_a <- extract.clade(clade_g, 669)
# plot(clade_g_a)

## bootstrap values for the top part of the tree. have to exclude the last because it is no longer a node
clade_g_b_bootstrap <- gp23_95_ref_tree$node.label[match(clade_g_b$tip.label,gp23_95_ref_tree$tip.label)][-length(clade_g_b$tip.label)]

clade_g_b <- apeBoot(clade_g_b, clade_g_b_bootstrap)

clade_g_b@bootstrap$bootstrap[as.numeric(as.character(clade_g_b@bootstrap$bootstrap)) <50] <- NA
clade_g_b@bootstrap$bootstrap[as.character(clade_g_b@bootstrap$bootstrap) == "Root"] <- NA
clade_g_b@bootstrap$bootstrap <- as.numeric(as.character(clade_g_b@bootstrap$bootstrap))


tree_clade_G_a <- colour_tips_of_trees(clade_g, 669)+
  geom_treescale(offset=15)

small_tree_G_highlighted <- p2a +
  geom_cladelabel( node=3295, "G",
                   barsize=1,
                   fontsize=4,
                   angle=0,
                   offset=0.6,
                   # offset.text=0.12,
                   color= mypalette[6])

tree_clade_G_b <- ggtree(clade_g_b, ladderize = FALSE,colour=line_colour)+
  #geom_tiplab(size=3)+
  geom_treescale(x=0.05,y=2, offset=2.5)+geom_tiplab(aes(subset=(!grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz", p$data$label) & isTip)),
                                                     color="grey30",
                                                     size=2.2) +
  geom_point2(aes(subset=( grepl("Chow",
                                 p$data$label) & isTip)),
              color="green",
              size=2)+
  geom_point2(aes(subset=( grepl("Filee",
                                 p$data$label) & isTip)),
              color="brown",
              size=2)+
  geom_point2(aes(subset=( grepl("Jia",
                                 p$data$label) & isTip)),
              color="lightblue",
              size=2)+ 
  geom_point2(aes(subset=( grepl("Lopez-Bueno|LopezBueno",
                                 p$data$label) & isTip)),
              color="darkgreen",
              size=2)+
  geom_point2(aes(subset=( grepl("Bellas",
                                 p$data$label) & isTip)),
              color="purple",
              size=2)+
  geom_point2(aes(subset=( grepl("Knapik",
                                 p$data$label) & isTip)),
              color="turquoise",
              size=2)+
  geom_point2(aes(subset=( grepl("Liu",
                                 p$data$label) & isTip)),
              color="lightgrey",
              size=2)+
  geom_point2(aes(subset=( grepl("Butina",
                                 p$data$label) & isTip)),
              color="beige",
              size=2)+ 
  geom_point2(aes(subset=( grepl("Comeau",
                                 p$data$label) & isTip)),
              color="darkred",
              size=2)+ 
  geom_point2(aes(subset=( grepl("Sandaa",
                                 p$data$label) & isTip)),
              color="pink",
              size=2)+
  geom_point2(aes(subset=( grepl("Mabiz",
                                 p$data$label) & isTip)),
              color="darkslateblue",
              size=2)+
  geom_point2(aes(subset=( grepl("OTU", p$data$label) & isTip)),
              shape=23,
              fill="steelblue",
              size=1 )+ 
  geom_tiplab(aes(subset=(grepl("OTU", p$data$label) & isTip)),
              color="black",
              size=1.8)+
  geom_point(aes(colour=bootstrap,alpha=bootstrap),
             size=3,na.rm=TRUE)+
  scale_colour_continuous(low='grey', high='black')+
  scale_alpha_continuous(range = c(0.7, 1), na.value=0)


tree_clade_H <- colour_tips_of_trees(tree,2955)+
  geom_treescale(x=0.1,y=2, offset=7)

small_tree_H_highlighted <- p2a +
  geom_cladelabel(node=2955,
                  "H", 
                  barsize=1,
                  fontsize=4,
                  angle=0,
                  offset=0.2,
                  color= mypalette[7])

## needs to be split

clade_i <- extract.clade(tree, 2295)
ggtree(clade_i, ladderize = FALSE,colour=line_colour) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

tips_to_drop <- tips(clade_i, node=671) 
clade_i_b <- drop.tip(clade_i, tips_to_drop)


## bootstrap values for the top part of the tree. have to exclude the last because it is no longer a node
clade_i_b_bootstrap <- gp23_95_ref_tree$node.label[match(clade_i_b$tip.label,gp23_95_ref_tree$tip.label)][-length(clade_i_b$tip.label)]

clade_i_b <- apeBoot(clade_i_b, clade_i_b_bootstrap)

clade_i_b@bootstrap$bootstrap[as.numeric(as.character(clade_i_b@bootstrap$bootstrap)) <50] <- NA
clade_i_b@bootstrap$bootstrap[as.character(clade_i_b@bootstrap$bootstrap) == "Root"] <- NA
clade_i_b@bootstrap$bootstrap <- as.numeric(as.character(clade_i_b@bootstrap$bootstrap))




tree_clade_I_a <- colour_tips_of_trees(clade_i,671)+
  geom_treescale(x=0.1,y=2, offset=15)

tree_clade_I_b <- ggtree(clade_i_b, ladderize = FALSE, colour=line_colour)+
  #geom_tiplab(size=3)+
  geom_treescale(x=0.05,y=2, offset=2.5)+
  geom_tiplab(aes(subset=(!grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz", p$data$label) & isTip)),
              color="grey30",
              size=2.2) +
  geom_point2(aes(subset=( grepl("Chow",
                                 p$data$label) & isTip)),
              color="green",
              size=2)+
  geom_point2(aes(subset=( grepl("Filee",
                                 p$data$label) & isTip)),
              color="brown",
              size=2)+
  geom_point2(aes(subset=( grepl("Jia",
                                 p$data$label) & isTip)),
              color="lightblue",
              size=2)+ 
  geom_point2(aes(subset=( grepl("Lopez-Bueno|LopezBueno",
                                 p$data$label) & isTip)),
              color="darkgreen",
              size=2)+
  geom_point2(aes(subset=( grepl("Bellas",
                                 p$data$label) & isTip)),
              color="purple",
              size=2)+
  geom_point2(aes(subset=( grepl("Knapik",
                                 p$data$label) & isTip)),
              color="turquoise",
              size=2)+
  geom_point2(aes(subset=( grepl("Liu",
                                 p$data$label) & isTip)),
              color="lightgrey",
              size=2)+
  geom_point2(aes(subset=( grepl("Butina",
                                 p$data$label) & isTip)),
              color="beige",
              size=2)+ 
  geom_point2(aes(subset=( grepl("Comeau",
                                 p$data$label) & isTip)),
              color="darkred",
              size=2)+ 
  geom_point2(aes(subset=( grepl("Sandaa",
                                 p$data$label) & isTip)),
              color="pink",
              size=2)+
  geom_point2(aes(subset=( grepl("Mabiz",
                                 p$data$label) & isTip)),
              color="darkslateblue",
              size=2)+
  geom_point2(aes(subset=( grepl("OTU", p$data$label) & isTip)),
              shape=23,
              fill="steelblue",
              size=1 )+
  geom_tiplab(aes(subset=(grepl("OTU", p$data$label) & isTip)),
              color="black",
              size=1.8)+
  geom_point(aes(colour=bootstrap,alpha=bootstrap),
             size=3,na.rm=TRUE)+
  scale_colour_continuous(low='grey', high='black')+
  scale_alpha_continuous(range = c(0.7, 1), na.value=0)


small_tree_I_highlighted <- p2a +
  geom_cladelabel( node=2295, "I",
                   barsize=1,
                   fontsize=4,
                   angle=0,
                   offset=0.2,
                   # offset.text=0.12,
                   color= mypalette[8])


## make one with all of the stuff on it for the manuscript
## tips getting cut off...in ggtree
pdf(paste0(figures_dir,"gp23_zoomed_in_tree%0d.pdf"), width = 12, height = 17, onefile = FALSE)
plot_grid(p2a,
          top_tree,
          small_tree_A_highlighted, 
          tree_clade_A,
          ncol=2,
          nrow=2,
          rel_heights=c(2,1),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_B_highlighted,
          tree_clade_B,
          ncol=2,
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_C_highlighted,
          tree_clade_C,
          ncol=2,
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_D_highlighted,
          tree_clade_D,
          small_tree_E_highlighted,
          tree_clade_E,
          ncol=2,
          nrow = 2,
          rel_heights=c(1,1.5),
          rel_widths = c(1/8,7/8),
          #scale=c(0.8,0.9,0.9), #because keep missing text
          label_size = 20)

plot_grid(small_tree_F_highlighted,
          tree_clade_F,
          ncol=2,
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_G_highlighted,
          tree_clade_G_a,
          ncol=2,
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(NULL,
          tree_clade_G_b,
          ncol=2,
          rel_widths = c(1/8,7/8),
          label_size = 20)


plot_grid(small_tree_H_highlighted,
          tree_clade_H,
          ncol=2,
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_I_highlighted,
          tree_clade_I_a,
          ncol=2,
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(NULL,
          tree_clade_I_b,
          ncol=2,
          rel_widths = c(1/8,7/8),
          label_size = 20)
dev.off()


## need to find an efficient way to do this:
library(ggplot2)
library(ape)
library(phytools)
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

# adding in the seasons and the spring bloom for the ggplots
season_line <- geom_vline(xintercept = as.numeric(c(as.Date("2010-03-22"), 
                                                    as.Date("2010-06-22"), 
                                                    as.Date("2010-09-22"),
                                                    as.Date("2010-12-22"),
                                                    as.Date("2011-03-22"),
                                                    as.Date("2011-06-22"))),
                          colour="grey",
                          size=2)

mid_summer_x <- as.Date("2010-09-22")-as.numeric(as.Date("2010-09-22")-as.Date("2010-06-22"))/2
mid_fall_x <- as.Date("2010-12-22")-as.numeric(as.Date("2010-12-22")-as.Date("2010-09-22"))/2
mid_winter_x <- as.Date("2011-03-22")-as.numeric(as.Date("2011-03-22")-as.Date("2010-12-22"))/2
mid_spring_x <- as.Date("2011-06-22")-as.numeric(as.Date("2011-06-22")-as.Date("2011-03-22"))/2

season_text <- annotate("text",x=c(mid_summer_x, mid_fall_x, mid_winter_x, mid_spring_x), y=1.05, label=c("Summer", "Fall", "Winter", "Spring"))

spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=2)
date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")))



#http://www.r-phylo.org/wiki/HowTo/DataTreeManipulation#Is_there_a_shorthand_way_to_refer_to_a_specific_list_of_taxa_.28for_example.2C_all_members_of_a_particular_group.29.3F

#The package geiger offers the function node.leaves, which facilitates subgroup designation if you know the basal node of a group. The basic syntax is:
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number")

### maybe change these to prop abundance???

proportional_MPL <- as.data.frame(prop.table(as.matrix(normalized_MPL_OTUs),  margin=1))
proportional_gp23 <- as.data.frame(prop.table(as.matrix(normalized_gp23_OTUs),  margin=1))

### Get relevant VC numbers ====
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
## Reformat date
Jericho_data$Date <- as.Date(Jericho_data$Date)

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)

high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)
## want to annotate missing data days. 

#missing dates in MPL

## don't want to care about the high res_vcs
## want Jericho_dates that are not in MPL table and high res
Jericho_no_high_res_vcs <- Jericho_data$VC_number[!(Jericho_data$VC_number %in% high_res_vcs)]

missing_MPL_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(proportional_MPL))]
missing_MPL_dates <- Jericho_data$Date[match(missing_MPL_samples,
                                             Jericho_data$VC_number)]

missing_gp23_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(proportional_gp23))]
missing_gp23_dates <- Jericho_data$Date[match(missing_gp23_samples,
                                              Jericho_data$VC_number)]
## want to annotate those days using star


plot_tree_root_and_ladderize <- function (tree_file, root) {
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


AVS_90_ref_tree <- plot_tree_root_and_ladderize("../results/AVS_90_miseq_data_with_env_iso_and_ref_align_clustal_edited.tree",
                                                root = "Ectocarpus_siliculosus_virus_1_gi_13177282")
reference_tips <- AVS_90_ref_tree$tip.label[!(AVS_90_ref_tree$tip.label %in% colnames(normalized_AVS_OTUs))]
AVS_90_only_miseq <- drop.tip(AVS_90_ref_tree, reference_tips)
plot(AVS_90_only_miseq)
ggtree(AVS_90_only_miseq) + geom_text(aes(label=node))
## so could use these to look at the OTU table
groupA <- tips(AVS_90_only_miseq, 14)
groupB <- tips(AVS_90_only_miseq, 13)
groupC <- tips(AVS_90_only_miseq, 15)
groupD <- tips(AVS_90_only_miseq, 8)

## would like to summarize somehow....
sum_of_group_otus <- function (normalized_OTUs, group, group_string) {
 
 test_group <- subset(normalized_OTUs, select=colnames(normalized_OTUs) %in% group)
 test_row <- as.data.frame(rowSums(test_group))
 names(test_row)[1] <- paste0("sum_of_",group_string )
 test <- adply(test_group, 1, sum)
 names(test)[3] <- paste0("sum_of_",group_string )
 return( test_row)
}

## first rename the column names in the OTU table to something better
barplot_of_top_groups_over_time <- function (all_groups_merged, 
                                             Jericho_data,
                                             title,
                                             colours_custom) {
 melt_all_together <- melt(all_groups_merged)
 ## Add in dates:
 melt_all_together$Date <- Jericho_data$Date[match(melt_all_together$VC, Jericho_data$VC_number)]
 
 group_barplot <- ggplot(melt_all_together, aes(x=Date,y=value, group=variable)) + 
  season_line +
  season_text +
  spring_bloom_line+
  geom_bar(aes(fill=variable), stat="identity", position="stack", colour=line_colour)+
  scale_fill_manual(values=colours_custom)+
  theme_JAG_presentation(base_size=18)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  date_scaling +
  ylab("Relative abundance")+
  guides(fill=guide_legend(title="group"))
 #+
 #guide_legend(title="group from tree")+
 #ggtitle(title)
 
 return(group_barplot)
}

### MPL ####

RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762")

pu <- ggtree(RdRp_95_ref_tree, ladderize = FALSE) 
ggtree(RdRp_95_ref_tree, ladderize = FALSE) +
 geom_text(aes(label=node), hjust=0.3)+
 geom_tiplab()
## not working right now. 
# +
#  geom_label(aes(label=node.label),
#            color="red",
#            size=10,
#            vjust=-0.3) # this does bootstrap

## so could use these to look at the OTU table
groupA <- tips(RdRp_95_ref_tree, 1768) ## includes HAKA
groupB <- tips(RdRp_95_ref_tree, 1764) ## env and Aster glacialis
groupC <- tips(RdRp_95_ref_tree, 1760) ## could I add to the previous one...
groupD <- tips(RdRp_95_ref_tree, 1701) ##CUlley seqs, Marine from SF
groupE <- tips(RdRp_95_ref_tree,1613) ## Diatom viruses
groupF <- tips(RdRp_95_ref_tree,1700) ## mini.....
groupG <- tips(RdRp_95_ref_tree,1450) ## low bootstrap...
groupH <- tips(RdRp_95_ref_tree,1017)

groupA_data_frame <- data.frame("group"="A",
                                values=groupA)
groupB_data_frame <- data.frame("group"="B",
                                values=groupB)
groupC_data_frame <- data.frame("group"="C",
                                values=groupC)
groupD_data_frame <- data.frame("group"="D",
                                values=groupD)
groupE_data_frame <- data.frame("group"="E",
                                values=groupE)
groupF_data_frame <- data.frame("group"="F",
                                values=groupF)
groupG_data_frame <- data.frame("group"="G",
                                values=groupG)
groupH_data_frame <- data.frame("group"="H",
                                values=groupH)

group_data_frame <- rbind(groupA_data_frame, groupB_data_frame)
group_data_frame <- rbind(group_data_frame, groupC_data_frame)
group_data_frame <- rbind(group_data_frame, groupD_data_frame)
group_data_frame <- rbind(group_data_frame, groupE_data_frame)
group_data_frame <- rbind(group_data_frame, groupF_data_frame)
group_data_frame <- rbind(group_data_frame, groupG_data_frame)
group_data_frame <- rbind(group_data_frame, groupH_data_frame)

## make table for annotation for use in script 46
write.csv(group_data_frame, "../results/RdRp_groups_with_OTUs.csv") 

## missing OTUs 119 -node splits from 1014 and 377 from slip 1015

mypalette<-brewer.pal(7,"Dark2")
#colours_for_groups <- colorRampPalette(mypalette)(10)

p <- ggtree(RdRp_95_ref_tree, ladderize=FALSE)+
 geom_text(aes(label=node))
p <- ggtree(RdRp_95_ref_tree,
            ladderize = FALSE,
            colour=line_colour)+
 scale_y_continuous(labels=NULL)+
 scale_x_continuous(labels=NULL)+
 theme_JAG_presentation()+
 theme(axis.title=element_blank(),
       axis.ticks=element_blank(),
       panel.grid = element_blank()) 


pdf(paste0(figures_dir,"RdRp_95_ref_free_plain.pdf"),
    width = 30,
    height = 15,
    onefile = FALSE)
p
dev.off()

p2 <- p+ geom_cladelabel(node=1768,
                         "A",
                         angle=0,
                         offset.text=0.12,
                         barsize=4, 
                         fontsize=6,
                         offset=2,
                         color="grey") + 
 geom_cladelabel(node=1764,
                 "B",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.9242,
                 color=mypalette[1]) + 
 geom_cladelabel(node=1760,
                 "C",
                 angle=0,
                 color=mypalette[2],
                 barsize=4,
                 fontsize=6,
                 offset=1.546,
                 offset.text=-0.3) + 
 geom_cladelabel(node=1701,
                 "D",
                 angle=0,
                 barsize=4,
                 fontsize=6,
                 offset=0.325,
                 offset.text=-0.3,
                 color=mypalette[3]) + 
 geom_cladelabel(node=1613,
                 "E",
                 angle=0,
                 barsize=4,
                 fontsize=6,
                 offset=1.482,
                 offset.text=0.12,
                 color=mypalette[4]) + 
 geom_cladelabel(node=1700,
                 "F",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.355,
                 color=mypalette[5]) + 
 geom_cladelabel(node=1450,
                 "G",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=0.78,
                 color=mypalette[6]) +
 geom_cladelabel(node=1017,
                 "H",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.065,
                 color=mypalette[7])

pdf(paste0(figures_dir,"RdRp_95_ref_free_annotated_with_groups.pdf"),
    width = 30,
    height = 15,
    onefile = FALSE)

p2

dev.off()



tree <- groupClade(RdRp_95_ref_tree, node=c(1768,1764,1760,1701,1613,1700,1450,1017))
p2a <- ggtree(tree,
              aes(color=group),ladderize = FALSE)+
 scale_color_manual(values=c(line_colour,"grey",
                             mypalette[1:7]),
                    guide = FALSE
 )+
 scale_y_continuous(labels=NULL)+
 scale_x_continuous(labels=NULL)+
 theme(axis.title=element_blank(),
       axis.ticks=element_blank(),
       panel.grid = element_blank()
 ) 

p3a <-  p2a+ geom_cladelabel(
 node=1768,
 "A",
 angle=0,
 offset.text=0.12,
 barsize=4, 
 fontsize=6,
 offset=2,
 color="grey") + 
 geom_cladelabel(node=1764,
                 "B",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.9242,
                 color=mypalette[1]) + 
 geom_cladelabel(node=1760,
                 "C",
                 angle=0,
                 color=mypalette[2],
                 barsize=4,
                 fontsize=6,
                 offset=1.546,
                 offset.text=-0.3) + 
 geom_cladelabel(node=1701,
                 "D",
                 angle=0,
                 barsize=4,
                 fontsize=6,
                 offset=0.325,
                 offset.text=-0.3,
                 color=mypalette[3]) + 
 geom_cladelabel(node=1613,
                 "E",
                 angle=0,
                 barsize=4,
                 fontsize=6,
                 offset=1.482,
                 offset.text=0.12,
                 color=mypalette[4]) + 
 geom_cladelabel(node=1700,
                 "F",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.355,
                 color=mypalette[5]) + 
 geom_cladelabel(node=1450,
                 "G",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=0.78,
                 color=mypalette[6]) +
 geom_cladelabel(node=1017,
                 "H",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.065,
                 color=mypalette[7])
RdRP_annotated_tree <- p3a
RdRP_annotated_tree
#gheatmap_jag(p2, normalized_MPL_OTUs, RdRp_95_ref_tree, color=line_colour)


p2+geom_point2(aes(subset=isTip,
                   label=label),
               size=4,
               colour="steelblue")



p3 <- p2+geom_point2(aes(subset=(!grepl("OTU|Culley|0\\.[:digit:]+",
                                        RdRP_annotated_tree$data$label)& isTip),
                         label=label),
                     size=4,
                     # hjust=2,
                     colour="steelblue") +
 geom_point2(aes(subset=(grepl("OTU", RdRP_annotated_tree$data$label) & isTip)),
             color="#FDAC4F",
             shape=8,
             size=2) +
 geom_point2(aes(subset=(grepl("Culley", RdRP_annotated_tree$data$label) & isTip)),
             color="purple",
             shape=8,
             size=2)


pdf(paste0(figures_dir,"RdRp_95_ref_free_annotated_with_groups_and_env_ref_sequences.pdf"), width = 30, height = 15, onefile = FALSE)
print(p3)
dev.off()



p4 <- p2+ geom_text2(aes(subset=( !grepl("OTU|Culley|0\\.[:digit:]+",
                                         RdRp_95_ref_tree$tip.label)  & isTip),
                         label=label),
                     color="grey30",
                     size=2.5,
                     hjust = -.1) +
 geom_point2(aes(subset=(grepl("OTU", RdRP_annotated_tree$data$label) & isTip)),
             color="#FDAC4F",
             shape=8,
             size=2) +
 geom_point2(aes(subset=(grepl("Culley", RdRP_annotated_tree$data$label) & isTip)),
             color="purple",
             shape=8,
             size=2)
p4

pdf(paste0(figures_dir,"RdRp_95_ref_free_text_of_refs_groups_and_env_sequences.pdf"), width = 30, height = 40, onefile = FALSE)
print(p4)
dev.off()

(p + geom_text(aes(label=node)) +geom_tiplab())

## need to collapse some clades. Do annotations after I guess
collapsed_tree <- p %>% collapse(node=1025) %>% collapse(node=1462) %>% collapse(node=1615) %>% collapse(node=1415) %>% 
 collapse(node=1608) %>% 
 collapse(node=1447) %>% 
 collapse(node=1407) %>% 
 collapse(node=1384) %>% 
 collapse(node=1603) %>% 
 collapse(node=1755) %>% 
 collapse(node=1596) %>% 
 collapse(node=1566) %>% 
 collapse(node=1442) %>% 
 collapse(node=1399) %>% 
 collapse(node=1381) %>% 
 collapse(node=1412) %>% 
 collapse(node=1606) %>% 
 collapse(node=1716) %>% 
 # collapse(node=823) %>% 
 collapse(node=1410) %>% 
 collapse(node=1753) %>% 
 collapse(node=1565) %>% 
 collapse(node=1767) %>% 
 collapse(node=1563)


collapsed_tree_with_nodes <- collapsed_tree +
 geom_point2(aes(subset=(node == 1025)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1462)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1615)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1415)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1608)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1447)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1407)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1384)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1603)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1755)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1596)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1566)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1442)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1399)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1381)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1412)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1606)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1716)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1410)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1753)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1565)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1767)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 1563)), size=2, shape=23, fill="steelblue")

collapsed_with_names <-  collapsed_tree_with_nodes +
 geom_tiplab(aes(subset=(node=1563)),
             size=2,
             label=paste0("N=",
                          length(getDescendants(RdRp_95_ref_tree,
                                                1563))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1025)), size=2, label=paste0("N=", length(getDescendants(RdRp_95_ref_tree, 1025))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1462)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1462))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1615)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1615))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1415)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1415))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1608)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1608))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1447)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1447))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1407)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1407))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1384)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1384))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1603)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1603))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1755)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1755))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1596)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1596))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1566)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1566))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1442)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1442))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1399)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1399))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1381)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1381))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1412)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1412))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1606)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1606))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1716)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1716))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1410)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1410))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1753)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1753))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1565)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1565))," OTUs"))+
 geom_tiplab(aes(subset=(node == 1767)), size=2, label=paste0("N=",                           length(getDescendants(RdRp_95_ref_tree, 1767))," OTUs"))+
 geom_tiplab(aes(subset=( !grepl("OTU|Culley|0\\.[:digit:]+",
                                 RdRp_95_ref_tree$tip.label)  & isTip),
                 label=label),
             color="grey30",
             size=2.3
 )+
 geom_point2(aes(subset=(grepl("Culley", RdRP_annotated_tree$data$label) & isTip)),
             color="purple",
             shape=8,
             size=1)+
 geom_point2(aes(subset=(grepl("OTU", RdRP_annotated_tree$data$label) & isTip)),
             shape=23, fill="steelblue",
             size=1)

collapsed_with_names 

pdf(paste0(figures_dir,"RdRp_95_ref_collapsed_OTUs.pdf"), width = 10, height = 30, onefile = FALSE)
print(collapsed_with_names )
dev.off()

## add annotations: not working!
test <- collapsed_with_names +  geom_cladelabel(
 node=1768,
 "A",
 angle=0,
 offset.text=0.12,
 barsize=4, 
 fontsize=6,
 offset=0,
 color="grey") + 
 geom_cladelabel(node=1764,
                 "B",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.9242,
                 color=mypalette[1]) + 
 geom_cladelabel(node=1760,
                 "C",
                 angle=0,
                 color=mypalette[2],
                 barsize=4,
                 fontsize=6,
                 offset=1.546,
                 offset.text=-0.3) + 
 geom_cladelabel(node=1701,
                 "D",
                 angle=0,
                 barsize=4,
                 fontsize=6,
                 offset=0.325,
                 offset.text=-0.3,
                 color=mypalette[3]) + 
 geom_cladelabel(node=1613,
                 "E",
                 angle=0,
                 barsize=4,
                 fontsize=6,
                 offset=1.482,
                 offset.text=0.12,
                 color=mypalette[4]) + 
 geom_cladelabel(node=1700,
                 "F",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.355,
                 color=mypalette[5]) + 
 geom_cladelabel(node=1450,
                 "G",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=0.78,
                 color=mypalette[6]) +
 geom_cladelabel(node=1017,
                 "H",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.065,
                 color=mypalette[7])


p5 <- p2+ geom_text2(aes(subset=( !grepl("OTU|Culley|0\\.[:digit:]+",
                                         RdRp_95_ref_tree$tip.label)  & isTip),
                         label=label),
                     color="grey30",
                     size=2.5,
                     hjust = -.1)  +
 geom_point2(aes(subset=(grepl("OTU", RdRP_annotated_tree$data$label) & isTip)),
             color="#FDAC4F",
             shape=8,
             size=2)  +
 geom_point2(aes(subset=(grepl("Culley", RdRP_annotated_tree$data$label) & isTip)),
             color="purple",
             shape=8,
             size=2)
p5

pdf(paste0(figures_dir,"RdRp_95_ref_free_text_of_env_seqs_groups_and_refs_sequences.pdf"), width = 30, height = 40, onefile = FALSE)
print(p5)
dev.off()



names_tips <- c("Isolates",
                "Culley_2003_and_2007",
                "OTUs")
cols_tips <- c("grey30",
               "green",
               "steelblue")


## can't get legend to work
length_cols <-c(rep(1,length(cols_tips)))  

df_cols <- data.frame(names_tips, cols_tips,length_cols)

colouring_bars <- as.character(df_cols$cols_tips)
names(colouring_bars) <- df_cols$names_tips

labeling_tree <- ggplot(df_cols, aes(x=names_tips, y=length_cols,group=names_tips))+
 geom_bar(stat="identity", aes(fill=names_tips))+
 scale_fill_manual(values=colouring_bars, name="")

## get only legend
library(gridExtra)
g_legend<-function(a.gplot){
 tmp <- ggplot_gtable(ggplot_build(a.gplot))
 leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
 legend <- tmp$grobs[[leg]]
 legend
}

legend <- g_legend(labeling_tree)

pdf(paste0(figures_dir,"RdRp_95_env_OTUs_legend.pdf"), width = 8, height = 10, onefile = FALSE)
grid.arrange(legend)
dev.off()


pdf(paste0(figures_dir,"RdRp_95_ref_collapsed_OTUs.pdf"), width = 30, height = 40, onefile = FALSE)
grid.arrange(legend, collapsed_with_names, 
             ncol=2, nrow=1, widths=c(1/6,5/6))
dev.off()

 
tree_tips_coloured_and_clades_annotated <- p2+ 
 geom_tiplab(aes(subset=(!grepl("OTU|Culley|0\\.[:digit:]+|^[0-9]{2}$",
                                collapsed_tree_with_nodes$data$label) & isTip)),
             color="grey30",
             size=2.2) +
 geom_point2(aes(subset=( grepl("Culley",
                                p$data$label) & isTip)),
             color="green",
             size=2)+
 geom_point2(aes(subset=( grepl("OTU", p$data$label) & isTip)),
             shape=23,
             fill="steelblue", 
             size=1 )


colour_tips_of_trees <- function (tree, node_to_extract) {
 new_subtree <- extract.clade(tree, node = node_to_extract)
 #print(new_subtree$tip.labels)
 
 new_tree_bootstrap <- RdRp_95_ref_tree$node.label[match(new_subtree$tip.label,RdRp_95_ref_tree$tip.label)][-length(new_subtree$tip.label)]
 
 new_subtree <- apeBoot(new_subtree, new_tree_bootstrap)
 
 
 new_subtree@bootstrap$bootstrap[as.numeric(as.character(new_subtree@bootstrap$bootstrap)) <50] <- NA
 new_subtree@bootstrap$bootstrap[as.character(new_subtree@bootstrap$bootstrap) == "Root"] <- NA
 new_subtree@bootstrap$bootstrap <- as.numeric(as.character(new_subtree@bootstrap$bootstrap))
 
 
 p <- ggtree(new_subtree, ladderize = FALSE, colour=line_colour)
 
 colour_stuff <- p+
  geom_tiplab(aes(subset=(!grepl("OTU|Culley|0\\.[:digit:]+|^[0-9]{2}$", p$data$label) & isTip)),
              color="grey30",
              size=2.2) +
  geom_point2(aes(subset=( grepl("Culley",
                                 p$data$label) & isTip)),
              color="green",
              size=2,
              shape = 24)+
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
 
 return(colour_stuff)
}

# 
# ## view clade doesn't provide all the information
# 
ggtree(tree, ladderize = FALSE) +
 geom_text2(aes(subset=!isTip, label=node),
            hjust=-.3) +
 geom_tiplab()

# ## top part
tips_to_drop <- tips(RdRp_95_ref_tree, node=1006)
top_tree_part <- drop.tip(RdRp_95_ref_tree, tips_to_drop)

grep("OTU", top_tree_part$tip.label, invert = TRUE, value = TRUE)


## bootstrap values for the top part of the tree. have to exclude the last because it is no longer a node
top_tree_bootstrap <- RdRp_95_ref_tree$node.label[match(top_tree_part$tip.label,RdRp_95_ref_tree$tip.label)][-length(top_tree_part$tip.label)]

top_tree_part <- apeBoot(top_tree_part, top_tree_bootstrap)

top_tree_part@bootstrap$bootstrap[as.numeric(as.character(top_tree_part@bootstrap$bootstrap)) <50] <- NA
top_tree_part@bootstrap$bootstrap[as.character(top_tree_part@bootstrap$bootstrap) == "Root"] <- NA
top_tree_part@bootstrap$bootstrap <- as.numeric(as.character(top_tree_part@bootstrap$bootstrap))



top_tree <- ggtree(top_tree_part, ladderize = FALSE,colour=line_colour)+
 #geom_tiplab(size=3)+
 geom_treescale(x=0.05,y=2, offset=2.5)+geom_tiplab(aes(subset=(!grepl("OTU|Culley|0\\.[:digit:]+|^[0-9]{2}$", p$data$label) & isTip)),
                                                    color="grey30",
                                                    size=2.2) +
 geom_point2(aes(subset=( grepl("Culley",
                                p$data$label) & isTip)),
             color="green",
             size=2)+
 geom_point2(aes(subset=( grepl("OTU", p$data$label) & isTip)), shape=23, fill="steelblue", size=1 )+
 geom_point(aes(colour=bootstrap,alpha=bootstrap), size=3,na.rm=TRUE)+
 scale_colour_continuous(low='grey', high='black')+
 scale_alpha_continuous(range = c(0.7, 1), na.value=0)+
 theme(plot.margin = unit(c(1,5,1,1), "cm"))


#par(mar=c(5,4,4,16))
tree_clade_A <- colour_tips_of_trees(tree, 1768)+
 geom_treescale()
# 
small_tree_A_highlighted <- p2a+
 geom_cladelabel(node=1768,
                 "A",
                 offset = 1.0,
                 barsize=1,
                 fontsize=4,
                 angle=0,
                 colour="grey")+
 theme(plot.margin = unit(c(1,1,1,1), "cm"))


tree_clade_B <- colour_tips_of_trees(tree, 1764)+
 geom_treescale()

small_tree_B_highlighted <- p2a+
 geom_cladelabel(node=1764,
                 "B",
                 offset = 1.0,
                 barsize=1,
                 fontsize=4,
                 angle=0,
                 color=mypalette[1])+
 theme(plot.margin = unit(c(1,1,1,1), "cm"))


tree_clade_C <- colour_tips_of_trees(tree, 1760)+
 geom_treescale()

small_tree_C_highlighted <- p2a+
 geom_cladelabel(node=1760,
                 "C",
                 offset = 1.0,
                 barsize=1,
                 fontsize=4,
                 angle=0,
                 color=mypalette[2])+
 theme(plot.margin = unit(c(1,1,1,1), "cm"))

tree_clade_D <- colour_tips_of_trees(tree, 1701)+
 geom_treescale()

small_tree_D_highlighted <- p2a+
 geom_cladelabel(node=1701,
                 "D",
                 offset = 0.3,
                 barsize=1,
                 fontsize=4,
                 angle=0,
                 color=mypalette[3])+
 theme(plot.margin = unit(c(1,1,1,1), "cm"))


tree_clade_E <- colour_tips_of_trees(tree, 1613)+
 geom_treescale()

small_tree_E_highlighted <- p2a+
 geom_cladelabel(node=1613,
                 "E",
                 offset = 1.0,
                 barsize=1,
                 fontsize=4,
                 angle=0,
                 color=mypalette[4])+
 theme(plot.margin = unit(c(1,1,1,1), "cm"))


tree_clade_F <- colour_tips_of_trees(tree, 1700)+
 geom_treescale()

small_tree_F_highlighted <- p2a+
 geom_cladelabel(node=1700,
                 "F",
                 offset = 1.0,
                 barsize=1,
                 fontsize=4,
                 angle=0,
                 color=mypalette[5])+
 theme(plot.margin = unit(c(1,1,1,1), "cm"))

tree_clade_G <- colour_tips_of_trees(tree, 1450)+
 geom_treescale()

small_tree_G_highlighted <- p2a+
 geom_cladelabel(node=1450,
                 "G",
                 offset = 0.5,
                 barsize=1,
                 fontsize=4,
                 angle=0,
                 color=mypalette[6])+
 theme(plot.margin = unit(c(1,1,1,1), "cm"))

tree_clade_H <- colour_tips_of_trees(tree, 1017)+
 geom_treescale()

small_tree_H_highlighted <- p2a+
 geom_cladelabel(node=1017,
                 "H",
                 offset = 0.8,
                 barsize=1,
                 fontsize=4,
                 angle=0,
                 color=mypalette[7])+
 theme(plot.margin = unit(c(1,1,1,1), "cm"))



# 
# ## make one with all of the stuff on it for the manuscript
# ## tips getting cut off...in ggtree
pdf(paste0(figures_dir,"MPL_zoomed_in_tree%0d.pdf"), width = 12, height = 17, onefile = FALSE)
plot_grid(p2a,
          top_tree,
          small_tree_A_highlighted,
          tree_clade_A,
          ncol=2,
          nrow=2,
          #labels = c("A", "B", "C", "D"),
          rel_heights=c(1.5,1),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_B_highlighted,
          tree_clade_B,
          small_tree_C_highlighted,
          tree_clade_C,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_D_highlighted,
          tree_clade_D,
          small_tree_E_highlighted,
          tree_clade_E,
          ncol=2,
          nrow = 2,
          #labels = c("A", "B", "C", "D"),
          rel_heights=c(1,1.5),
          rel_widths = c(1/8,7/8),
          #scale=c(0.8,0.9,0.9), #because keep missing text
          label_size = 20)

plot_grid(small_tree_F_highlighted,
          tree_clade_F,
          small_tree_G_highlighted,
          tree_clade_G,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_heights=c(1,2),
          rel_widths = c(1/8,7/8),
          label_size = 20)


plot_grid(small_tree_H_highlighted,
          tree_clade_H,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_widths = c(1/8,7/8),
          label_size = 20)

dev.off()


## full tree with bootstrap

RdRp_95_ref_tree_with_bootstrap <- apeBoot(RdRp_95_ref_tree, RdRp_95_ref_tree$node.label)

RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap[as.numeric(as.character(RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap)) <50] <- NA
RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap[as.character(RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap) == "Root"] <- NA
RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap <- as.numeric(as.character(RdRp_95_ref_tree_with_bootstrap@bootstrap$bootstrap))


tree2 <- groupClade(RdRp_95_ref_tree_with_bootstrap, node=c(1768,1764,1760,1701,1613,1700,1450,1017))

bootstrap_tree <- ggtree(tree2,
       aes(color=group),ladderize = FALSE,colour=line_colour)+
 scale_color_manual(values=c("black","grey",
                             mypalette[1:7]),
                    guide = FALSE
 )+
 scale_y_continuous(labels=NULL)+
 scale_x_continuous(labels=NULL)+
 theme_JAG_presentation()+
 theme(axis.title=element_blank(),
       axis.ticks=element_blank(),
       panel.grid = element_blank())+
 geom_treescale(x=0.05,y=2, offset=2.5)+
 geom_point(aes(fill=bootstrap,alpha=bootstrap),
            size=3,
            na.rm=TRUE,
            shape = 21,
            colour="white")+
 scale_fill_continuous(low='grey', high='black', guide = "legend")+
 scale_alpha_continuous(range = c(0.7, 1), na.value=0)+
theme(legend.position="right")+
 theme(plot.margin = unit(c(1,5,1,1), "cm"))

 
pdf(paste0(figures_dir,"MPL_tree_with_bootstrap.pdf"), width = 12, height = 17, onefile = FALSE)
bootstrap_tree
dev.off()

## add on
bootstrap_with_points <- bootstrap_tree+ 
 geom_point2(aes(subset=(!grepl("OTU|Culley|0\\.[:digit:]+|^[0-9]{2}$", p$data$label) & isTip)),
             color="grey30",
             shape=18,
             size=2.2)+
 geom_point2(aes(subset=( grepl("Culley",
                                p$data$label) & isTip)),
             color="green",
             size=2,
             shape = 24)+
 geom_point2(aes(subset=( grepl("OTU", p$data$label) & isTip)),
             shape=23,
             fill="steelblue",
             colour=NA,
             size=1 )

bootstrap_with_points_and_clade_labels <- bootstrap_with_points+geom_cladelabel(
 node=1768,
 "A",
 angle=0,
 offset.text=0.12,
 barsize=4, 
 fontsize=6,
 offset=2.0,
 color="grey") + 
 geom_cladelabel(node=1764,
                 "B",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.9242,
                 color=mypalette[1]) + 
 geom_cladelabel(node=1760,
                 "C",
                 angle=0,
                 color=mypalette[2],
                 barsize=4,
                 fontsize=6,
                 offset=1.546,
                 offset.text=-0.3) + 
 geom_cladelabel(node=1701,
                 "D",
                 angle=0,
                 barsize=4,
                 fontsize=6,
                 offset=0.325,
                 offset.text=-0.3,
                 color=mypalette[3]) + 
 geom_cladelabel(node=1613,
                 "E",
                 angle=0,
                 barsize=4,
                 fontsize=6,
                 offset=1.482,
                 offset.text=0.12,
                 color=mypalette[4]) + 
 geom_cladelabel(node=1700,
                 "F",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.355,
                 color=mypalette[5]) + 
 geom_cladelabel(node=1450,
                 "G",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=0.78,
                 color=mypalette[6]) +
 geom_cladelabel(node=1017,
                 "H",
                 angle=0,
                 offset.text=0.12,
                 barsize=4,
                 fontsize=6,
                 offset=1.065,
                 color=mypalette[7])

pdf(paste0(figures_dir,"MPL_tree_with_bootstrap_and_points.pdf"), width = 12, height = 17, onefile = FALSE)
grid.arrange(bootstrap_with_points_and_clade_labels, legend,
             ncol=2,
             widths=c(19/20,1/20))
dev.off()


## would like to summarize somehow....
sum_of_group_otus <- function (normalized_OTUs, group, group_string) {
 group_only_otu <- subset(group,grepl("OTU", group) )
 test_group <- subset(normalized_OTUs, select=colnames(normalized_OTUs) %in% group_only_otu)
 test_row <- as.data.frame(rowSums(test_group))
 names(test_row)[1] <- paste0("sum_of_",group_string )
 test <- adply(test_group, 1, sum)
 names(test)[3] <- paste0("sum_of_",group_string )
 return( test_row)
}


### propotional reads

sum_otus_groupA <- sum_of_group_otus(proportional_MPL, groupA, "A")
sum_otus_groupB <- sum_of_group_otus(proportional_MPL, groupB, "B")
#sum_otus_groupC <- sum_of_group_otus(proportional_MPL, groupC, "C")
sum_otus_groupD <- sum_of_group_otus(proportional_MPL, groupD, "D")
sum_otus_groupE <- sum_of_group_otus(proportional_MPL, groupE, "E")
# sum_otus_groupF <- sum_of_group_otus(proportional_MPL, groupF, "F")
sum_otus_groupG <- sum_of_group_otus(proportional_MPL, groupG, "G")
sum_otus_groupH <- sum_of_group_otus(proportional_MPL, groupH, "H")
#sum_otus_groupI <- sum_of_group_otus(proportional_MPL, groupI, "I")
# sum_otus_groupJ <- sum_of_group_otus(proportional_MPL, groupJ, "J")
#sum_otus_groupK <- sum_of_group_otus(proportional_MPL, groupK, "K")

### want to get the sums of all these groups and then melt and make a plot
l = list(sum_otus_groupA,sum_otus_groupB,
         #sum_otus_groupC,
         sum_otus_groupD, sum_otus_groupE, 
         #sum_otus_groupF, 
         sum_otus_groupG, 
         sum_otus_groupH
         #, sum_otus_groupI
         #,
         # sum_otus_groupJ,
         #sum_otus_groupK
)
all_together_MPL <- Reduce(merge, lapply(l, function(x) data.frame(x, VC = row.names(x))))
write.csv(all_together_MPL, file="../results/MPL_group_sums_by_site_prop.csv")

names(all_together_MPL) <- gsub("sum_of_", "", names(all_together_MPL))

mpl_bar_prop <- barplot_of_top_groups_over_time(all_together_MPL, 
                                                Jericho_data,
                                                "testing out MPL proportional",
                                                c("grey", mypalette[c(-2,-5)]))




mpl_bar_prop_without_high_res <- barplot_of_top_groups_over_time(subset(all_together_MPL,!(VC %in% high_res_vcs)), 
                                                                 Jericho_data,
                                                                 "MPL proportional",
                                                                 c("grey", mypalette[c(-2,-5)]))

mpl_bar_prop_without_high_res <- mpl_bar_prop_without_high_res + annotate("text", x = missing_MPL_dates, y = 0.00, label="*", size=12)

mpl_bar_prop_only_high_res <- barplot_of_top_groups_over_time(all_together_MPL, 
                                                              Jericho_data,
                                                              "MPL proportional",
                                                              c("grey", mypalette[c(-2,-5)]))+ 
 scale_x_date(breaks = date_breaks("week"), 
              #   labels = date_format("%b"),
              limits = c(as.Date("2011-01-15"),
                         as.Date("2011-02-15")))

pdf(paste0(figures_dir,"RdRp_95_ref_groups_over_time_prop.pdf"), width = 30, height = 15, onefile = FALSE)
print(mpl_bar_prop)
dev.off()

### Adding in clustering labels. 

MPL_abundance_clusters <- read.csv("../results/normalized_MPL_otus_clusters.csv")
names(MPL_abundance_clusters) <- c("Date", "cluster")
MPL_abundance_clusters$Date <- as.Date(MPL_abundance_clusters$Date)

p1<-mpl_bar_prop
p2<-ggplot(MPL_abundance_clusters,aes(Date,y=1,fill=factor(cluster)))+geom_tile()+
 scale_y_continuous(expand=c(0,0))+
 date_scaling +
 theme_JAG_presentation(base_size=18)+
 theme(axis.title=element_blank(),
       axis.ticks=element_blank(),
       axis.text=element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()
       # ,       legend.position="none"
 )

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)

## clustering labels 
pdf(paste0(figures_dir,"RdRp_95_ref_groups_over_time_with_label_prop.pdf"), width = 30, height = 15, onefile = FALSE)
grid.arrange(gp1, gp2, ncol=1,heights=c(9/10,1/10))
dev.off()


pdf(paste0(figures_dir,"RdRp_95_ref_groups_over_time_no_high_res_prop.pdf"), width = 30, height = 15, onefile = FALSE)
print(mpl_bar_prop_without_high_res )
dev.off()

MPL_abundance_clusters_no_high_res <- read.csv("../results/normalized_MPL_otus_clusters_no_high_res.csv")

names(MPL_abundance_clusters_no_high_res) <- c("Date", "cluster")
MPL_abundance_clusters_no_high_res$Date <- as.Date(MPL_abundance_clusters_no_high_res$Date)

p1<-mpl_bar_prop_without_high_res
p2<-ggplot(MPL_abundance_clusters_no_high_res,aes(Date,y=1,fill=factor(cluster)))+geom_tile()+
 scale_y_continuous(expand=c(0,0))+
 date_scaling +
 theme_JAG_presentation(base_size=18)+
 theme(axis.title=element_blank(),
       axis.ticks=element_blank(),
       axis.text=element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()
       # ,       legend.position="none"
 )

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)

pdf(paste0(figures_dir,"RdRp_95_ref_groups_over_time_with_label_no_high_res_prop.pdf"), width = 30, height = 15, onefile = FALSE)
grid.arrange(gp1, gp2, ncol=1,heights=c(9/10,1/10))
dev.off()

## add in env clustering label
MPL_env_clusters_no_high_res <- read.csv("../results/normalized_MPL_env_clusters_no_high_res.csv")

names(MPL_env_clusters_no_high_res) <- c("Date", "cluster")
MPL_env_clusters_no_high_res$Date <- as.Date(MPL_env_clusters_no_high_res$Date)
p3<-ggplot(MPL_env_clusters_no_high_res,aes(Date,y=1,fill=factor(cluster)))+geom_tile()+
 scale_y_continuous(expand=c(0,0))+  date_scaling +
 theme_JAG_presentation(base_size=18)+
 theme(axis.title=element_blank(),
       axis.ticks=element_blank(),
       axis.text=element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()
       # ,       legend.position="none"
 )

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  
gp3<-ggplotGrob(p3)  

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5], gp3$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
gp3$widths[2:5] <- as.list(maxWidth)

pdf(paste0(figures_dir,"RdRp_95_ref_groups_over_time_with_label_abun_and_env_no_high_res_prop.pdf"), width = 30, height = 15, onefile = FALSE)
grid.arrange(gp1, gp2,gp3, ncol=1,heights=c(8/10,1/10,1/10))
dev.off()




pdf(paste0(figures_dir,"RdRp_95_ref_groups_over_time_only_high_res_prop.pdf"), width = 30, height = 15, onefile = FALSE)
print(mpl_bar_prop_only_high_res )
dev.off()


## otus tables from the groups 

save_group_OTUs_to_file <- function (proportional_table, group, table_file_name) {
 proportional_table_group <- subset(proportional_table, select=colnames(proportional_table) %in% group)
 write.csv(proportional_table_group, file=paste0("../results/", table_file_name))
}

save_group_OTUs_to_file(proportional_MPL, groupA, "proportional_MPL_groupA.csv") 
save_group_OTUs_to_file(proportional_MPL, groupB, "proportional_MPL_groupB.csv") 
save_group_OTUs_to_file(proportional_MPL, groupC, "proportional_MPL_groupC.csv") 
save_group_OTUs_to_file(proportional_MPL, groupD, "proportional_MPL_groupD.csv") 
save_group_OTUs_to_file(proportional_MPL, groupE, "proportional_MPL_groupE.csv") 
save_group_OTUs_to_file(proportional_MPL, groupF, "proportional_MPL_groupF.csv") 
save_group_OTUs_to_file(proportional_MPL, groupG, "proportional_MPL_groupG.csv") 
save_group_OTUs_to_file(proportional_MPL, groupH, "proportional_MPL_groupH.csv") 


### gp23 ####
mypalette<-brewer.pal(9,"Set3")
colours_for_groups <- colorRampPalette(mypalette)(15)

gp23_95_ref_tree <- plot_tree_root_and_ladderize("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
                                                 root = "Enterobacteria_phage_T4")
#gp23_95_ref_tree$node.label <- NULL



ggtree(gp23_95_ref_tree, ladderize = FALSE) +
 geom_text(aes(label=node), hjust=0.3)+
 geom_tiplab() +
 geom_text2(aes(subset=!isTip,
                label=label),
            color="red",
            size=10,
            vjust=-0.3) 

## so could use these to look at the OTU table
groupA <- tips(gp23_95_ref_tree, 4472)
groupB <- tips(gp23_95_ref_tree, 4363)
groupC <- tips(gp23_95_ref_tree, 4252)
groupD <- tips(gp23_95_ref_tree, 4217)
groupE <- tips(gp23_95_ref_tree,4147)
groupF <- tips(gp23_95_ref_tree,3905)
groupG <- tips(gp23_95_ref_tree,3295)
groupH <- tips(gp23_95_ref_tree,2955)
groupI <- tips(gp23_95_ref_tree,2295)
#  groupJ <- tips(gp23_95_ref_tree,3876) 
#  groupK <- tips(gp23_95_ref_tree,3911)
#  groupL <- tips(gp23_95_ref_tree,4095) 
#  groupM <- tips(gp23_95_ref_tree,4097)
#  groupN <- tips(gp23_95_ref_tree,4392) 
#  groupO <- tips(gp23_95_ref_tree,4409)
#  groupP <- tips(gp23_95_ref_tree,4537) 



groupA_data_frame <- data.frame("group"="A",
                                values=groupA)
groupB_data_frame <- data.frame("group"="B",
                                values=groupB)
groupC_data_frame <- data.frame("group"="C",
                                values=groupC)
groupD_data_frame <- data.frame("group"="D",
                                values=groupD)
groupE_data_frame <- data.frame("group"="E",
                                values=groupE)
groupF_data_frame <- data.frame("group"="F",
                                values=groupF)
groupG_data_frame <- data.frame("group"="G",
                                values=groupG)
groupH_data_frame <- data.frame("group"="H",
                                values=groupH)
groupI_data_frame <- data.frame("group"="I",
                                values=groupI)

group_data_frame <- rbind(groupA_data_frame, groupB_data_frame)
group_data_frame <- rbind(group_data_frame, groupC_data_frame)
group_data_frame <- rbind(group_data_frame, groupD_data_frame)
group_data_frame <- rbind(group_data_frame, groupE_data_frame)
group_data_frame <- rbind(group_data_frame, groupF_data_frame)
group_data_frame <- rbind(group_data_frame, groupG_data_frame)
group_data_frame <- rbind(group_data_frame, groupH_data_frame)
group_data_frame <- rbind(group_data_frame, groupI_data_frame)


## make table for annotation for use in script 46
write.csv(group_data_frame, "../results/gp23_groups_with_OTUs.csv") 

p <- ggtree(gp23_95_ref_tree, ladderize = FALSE)+ geom_text(aes(label=node)) 
p <- ggtree(gp23_95_ref_tree,
            #layout="circular",
            colour=line_colour,
            ladderize = FALSE)+
 theme_JAG_presentation()+
 scale_y_continuous(labels=NULL)+
 scale_x_continuous(labels=NULL)+
 theme(axis.title=element_blank(), axis.ticks=element_blank(), panel.grid = element_blank()) 

pdf(paste0(figures_dir,"gp23_ref_free_plain.pdf"), width = 30, height = 15, onefile = FALSE)
p
dev.off()


p2 <- p+ geom_cladelabel(node=4472,
                         "A",
                         offset = 2.415,
                         barsize=2,
                         fontsize=6,
                         angle=0,
                         colour="grey") + 
 geom_cladelabel(node=4363, "B", barsize=2,
                 fontsize=6,
                 angle=0,
                 #offset.text=-0.1,
                 offset = 1.5,
                 color= mypalette[1]) +
 geom_cladelabel( node=4252, "C",                   barsize=2,
                  fontsize=6,
                  angle=0,
                  #offset.text=0.12,
                  offset = 1.425,
                  color= mypalette[2]) + 
 geom_cladelabel( node=4217, "D", barsize=2,
                  fontsize=6,
                  angle=0,
                  #offset.text=-0.1, 
                  color= mypalette[3]) + 
 geom_cladelabel( node=4147, "E",  
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  #offset.text=0.12,
                  offset = 1.4,
                  color= mypalette[4]) + 
 geom_cladelabel( node=3905, "F",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=1.165,
                  # offset.text=-0.1,
                  color= mypalette[5]) +
 geom_cladelabel( node=3295, "G",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=1.3,
                  # offset.text=0.12,
                  color= mypalette[6]) + 
 geom_cladelabel(node=2955,
                 "H", barsize=2,
                 fontsize=6,
                 angle=0,
                 offset=0.62,
                 color= mypalette[7]) + 
 geom_cladelabel( node=2295, "I",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=0.89,
                  # offset.text=0.12,
                  color= mypalette[8]) 
#%>% 
#   geom_cladelabel( node=3876, "J", color= mypalette[9]) %>% 
#   geom_cladelabel( node=3911, "K", color=colours_for_groups[10]) %>% 
#  geom_cladelabel( node=4095, "L", color=colours_for_groups[11]) %>% 
#  geom_cladelabel( node=4097, "M", color=colours_for_groups[12]) %>% 
#  geom_cladelabel( node=4392, "N", color=colours_for_groups[13]) %>% 
#  geom_cladelabel( node=4409, "O", color=colours_for_groups[14]) %>% 
#  geom_cladelabel( node=4537, "P", color=colours_for_groups[15])
#   
pdf(paste0(figures_dir,"gp23_95_ref_free_annotated_with_groups.pdf"), width = 30, height = 15, onefile = FALSE)
p2
dev.off()


tree <- groupClade(gp23_95_ref_tree, node=c(4472,4363, 4252,4217,4147,3905, 3295,2955,2295))
p2a <- ggtree(tree,
              aes(color=group),ladderize = FALSE,colour=line_colour)+
 scale_color_manual(values=c("black","grey",
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

p3a <- p2a+geom_cladelabel(node=4472,
                           "A",
                           offset = 2.415,
                           barsize=2,
                           fontsize=6,
                           angle=0,
                           colour="grey") + 
 geom_cladelabel(node=4363, "B", barsize=2,
                 fontsize=6,
                 angle=0,
                 #offset.text=-0.1,
                 offset = 1.5,
                 color= mypalette[1]) +
 geom_cladelabel( node=4252, "C",                   barsize=2,
                  fontsize=6,
                  angle=0,
                  #offset.text=0.12,
                  offset = 1.425,
                  color= mypalette[2]) + 
 geom_cladelabel( node=4217, "D", barsize=2,
                  fontsize=6,
                  angle=0,
                  #offset.text=-0.1, 
                  color= mypalette[3]) + 
 geom_cladelabel( node=4147, "E",  
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  #offset.text=0.12,
                  offset = 1.4,
                  color= mypalette[4]) + 
 geom_cladelabel( node=3905, "F",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=1.165,
                  # offset.text=-0.1,
                  color= mypalette[5]) +
 geom_cladelabel( node=3295, "G",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=1.3,
                  # offset.text=0.12,
                  color= mypalette[6]) + 
 geom_cladelabel(node=2955,
                 "H", barsize=2,
                 fontsize=6,
                 angle=0,
                 offset=0.62,
                 color= mypalette[7]) + 
 geom_cladelabel( node=2295, "I",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=0.89,
                  # offset.text=0.12,
                  color= mypalette[8]) 
gp23_ref_tree_annotated <- p3a

p3 <- p2+ geom_point2(aes(subset=( !grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz", p2$data$label) & isTip),
                          label=label),
                      color="blue",
                      size=4) +
 geom_point2(aes(subset=( grepl("OTU", p2$data$label) & isTip)), color="#FDAC4F", shape=8, size=2 ) +
 geom_point2(aes(subset=( grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz", p2$data$label) & isTip)), color="purple", shape=8, size=2 )

pdf(paste0(figures_dir,"gp23_95_ref_free_annotated_with_groups_and_env_ref_sequences.pdf"), width = 30, height = 15, onefile = FALSE)
print(p3)
dev.off()


p4 <-p2+ geom_text2(aes(subset=(!grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                       p2$data$label) & isTip),
                        label=label),
                    color="grey30",
                    size=2.5,
                    hjust = -.1) +
 geom_point2(aes(subset=( grepl("OTU", p2$data$label) & isTip)), color="#FDAC4F", shape=8, size=2 ) +
 geom_point2(aes(subset=( grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz", p2$data$label) & isTip)), color="purple", shape=8, size=2 )
p4

pdf(paste0(figures_dir,"gp23_95_ref_free_text_of_refs_groups_and_env_sequences.pdf"), width = 30, height = 40, onefile = FALSE)
print(p4)
dev.off()

p5 <- p2+ geom_text2(aes(subset=(!grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                        p2$data$label) & isTip),
                         label=label),
                     color="grey30",
                     size=2.5,
                     hjust = -.1) +
 geom_point2(aes(subset=( grepl("OTU", p2$data$label) & isTip)), color="#FDAC4F", shape=8, size=0.5 )+
 geom_text2(aes(subset=( grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz", p2$data$label) & isTip), label=label), color="darkslateblue",  size=2.5, hjust = -0.1 )

p5

pdf(paste0(figures_dir,"gp23_95_ref_free_text_of_env_seqs_groups_and_refs_sequences.pdf"), width = 30, height = 40, onefile = FALSE)
print(p5)
dev.off()


(p + geom_text(aes(label=node)) +geom_tiplab())

## need to collapse some clades. Do annotations after I guess
collapsed_tree <-
 p   %>% collapse(node=2857) %>%
 collapse(node=2569) %>%
 collapse(node=2402) %>%
 collapse(node=3042) %>% 
 collapse(node=4371) %>% 
 collapse(node=3691) %>% 
 collapse(node=3428) %>% 
 collapse(node=4449) %>% 
 collapse(node=4270) %>% 
 collapse(node=4153) %>% 
 collapse(node=3769) %>% 
 collapse(node=3497) %>% 
 collapse(node=2525) %>% 
 collapse(node=4230) %>% 
 collapse(node=4186) %>% 
 collapse(node=4113) %>% 
 collapse(node=3386) %>% 
 collapse(node=4458) %>% 
 collapse(node=4447) %>% 
 collapse(node=4313)



collapsed_tree_with_nodes <- collapsed_tree +
 geom_point2(aes(subset=(node == 2857)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 2569)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 2402)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 3042)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 4371)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 3691)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 3428)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 4449)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 4270)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 4153)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 3769)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 3497)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 2525)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 4230)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 4186)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 4113)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 3386)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 4458)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 4447)), size=2, shape=23, fill="steelblue")+
 geom_point2(aes(subset=(node == 4313)), size=2, shape=23, fill="steelblue")


names_tips <- c("Isolates",
                "Chow",
                "Filee",
                "Jia",
                "Lopez-Bueno",
                "Bellas",
                "Knapik",
                "Liu",
                "Butina",
                "Comeau",
                "Sandaa",
                "Mabiz",
                "OTUs")
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

collapsed_with_names <- collapsed_tree_with_nodes +
 geom_tiplab(aes(subset=(node == 2857)), size=2,label=paste0("N=",length(getDescendants( gp23_95_ref_tree,2857))," OTUs"))+
 geom_tiplab(aes(subset=(node == 2569)), size=2,label=paste0("N=", length(getDescendants( gp23_95_ref_tree,2569))," OTUs"))+
 geom_tiplab(aes(subset=(node == 2402)), size=2,label=paste0("N=",length(getDescendants( gp23_95_ref_tree, 2402))," OTUs"))+
 geom_tiplab(aes(subset=(node == 3042)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 3042))," OTUs"))+
 geom_tiplab(aes(subset=(node == 4371)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 4371))," OTUs"))+
 geom_tiplab(aes(subset=(node == 3691)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 3691))," OTUs"))+
 geom_tiplab(aes(subset=(node == 3428)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 3428))," OTUs"))+
 geom_tiplab(aes(subset=(node == 4449)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 4449))," OTUs"))+
 geom_tiplab(aes(subset=(node == 4270)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 4270))," OTUs"))+
 geom_tiplab(aes(subset=(node == 4153)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 4153))," OTUs"))+
 geom_tiplab(aes(subset=(node == 3769)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 3769))," OTUs"))+
 geom_tiplab(aes(subset=(node == 3497)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 3497))," OTUs"))+
 geom_tiplab(aes(subset=(node == 2525)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 2525))," OTUs"))+
 geom_tiplab(aes(subset=(node == 4230)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 4230))," OTUs"))+
 geom_tiplab(aes(subset=(node == 4186)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 1563))," OTUs"))+
 geom_tiplab(aes(subset=(node == 4113)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 1563))," OTUs"))+
 geom_tiplab(aes(subset=(node == 3386)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 1563))," OTUs"))+
 geom_tiplab(aes(subset=(node == 4458)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 4458))," OTUs"))+
 geom_tiplab(aes(subset=(node == 4447)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 4447))," OTUs"))+
 geom_tiplab(aes(subset=(node == 4313)), size=2, label=paste0("N=", length(getDescendants( gp23_95_ref_tree, 4313))," OTUs"))+
 geom_tiplab(aes(subset=(!grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz|^[0-9]{2}$",
                                collapsed_tree_with_nodes$data$label) & isTip)),
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
 geom_point2(aes(subset=( grepl("OTU", p$data$label) & isTip)), shape=23, fill="steelblue", size=1 )+ theme(legend.background = element_rect())


## can't get legend to work
length_cols <-c(rep(1,length(cols_tips)))  

df_cols <- data.frame(names_tips, cols_tips,length_cols)

colouring_bars <- as.character(df_cols$cols_tips)
names(colouring_bars) <- df_cols$names_tips

labeling_tree <- ggplot(df_cols, aes(x=names_tips, y=length_cols,group=names_tips))+
 geom_bar(stat="identity", aes(fill=names_tips))+
 scale_fill_manual(values=colouring_bars, name="Reference")

## get only legend
library(gridExtra)
g_legend<-function(a.gplot){
 tmp <- ggplot_gtable(ggplot_build(a.gplot))
 leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
 legend <- tmp$grobs[[leg]]
 legend
}

legend <- g_legend(labeling_tree)

pdf(paste0(figures_dir,"gp23_95_env_OTUs_legend.pdf"), width = 30, height = 40, onefile = FALSE)
grid.arrange(legend)
dev.off()


pdf(paste0(figures_dir,"gp23_95_ref_collapsed_OTUs.pdf"), width = 30, height = 40, onefile = FALSE)
grid.arrange(legend, collapsed_with_names, 
             ncol=2, nrow=1, widths=c(1/6,5/6))
dev.off()


## not collapsed

tree_tips_coloured_and_clades_annotated <- p2+ geom_tiplab(aes(subset=(!grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz|^[0-9]{2}$",
                                                                              collapsed_tree_with_nodes$data$label) & isTip)),
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
 geom_point2(aes(subset=( grepl("OTU", p$data$label) & isTip)), shape=23, fill="steelblue", size=1 )


pdf(paste0(figures_dir,"gp23_95_ref_OTUs_coloured_env_OTUs.pdf"), width = 17, height = 25, onefile = FALSE)
grid.arrange(legend, tree_tips_coloured_and_clades_annotated, 
             ncol=2, nrow=1, widths=c(1/6,5/6))
dev.off()


## ok what should actuallyt is subset the tre by groups ####


## would be good to hold all the nodes somewhere. 
## so should rearrange all these plots somewhere. 

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
          #labels = c("A", "B", "C", "D"),
          rel_heights=c(2,1),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_B_highlighted,
          tree_clade_B,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_C_highlighted,
          tree_clade_C,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_D_highlighted,
          tree_clade_D,
          small_tree_E_highlighted,
          tree_clade_E,
          ncol=2,
          nrow = 2,
          #labels = c("A", "B", "C", "D"),
          rel_heights=c(1,1.5),
          rel_widths = c(1/8,7/8),
          #scale=c(0.8,0.9,0.9), #because keep missing text
          label_size = 20)

plot_grid(small_tree_F_highlighted,
          tree_clade_F,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_G_highlighted,
          tree_clade_G_a,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(NULL,
          tree_clade_G_b,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_widths = c(1/8,7/8),
          label_size = 20)


plot_grid(small_tree_H_highlighted,
          tree_clade_H,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(small_tree_I_highlighted,
          tree_clade_I_a,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_widths = c(1/8,7/8),
          label_size = 20)

plot_grid(NULL,
          tree_clade_I_b,
          ncol=2,
          #labels = c("A", "B", "C", "D"),
          rel_widths = c(1/8,7/8),
          label_size = 20)
dev.off()



gp23_95_ref_tree_with_bootstrap <- apeBoot(gp23_95_ref_tree, gp23_95_ref_tree$node.label)

gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap[as.numeric(as.character(gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap)) <50] <- NA
gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap[as.character(gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap) == "Root"] <- NA
gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap <- as.numeric(as.character(gp23_95_ref_tree_with_bootstrap@bootstrap$bootstrap))


tree2 <- groupClade(gp23_95_ref_tree_with_bootstrap, node=c(4472,4363, 4252,4217,4147,3905, 3295,2955,2295))

bootstrap_tree <- ggtree(tree2,
                         aes(color=group),ladderize = FALSE,colour=line_colour)+
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
 geom_treescale(x=0.05,y=2, offset=2.5)+
 geom_point(aes(fill=bootstrap,alpha=bootstrap),
            size=3,
            na.rm=TRUE,
            shape = 21,
            colour="white")+
 scale_fill_continuous(low='grey', high='black', guide = "legend")+
 scale_alpha_continuous(range = c(0.7, 1), na.value=0)+
 theme(legend.position="right")+
 theme(plot.margin = unit(c(1,5,1,1), "cm"))


pdf(paste0(figures_dir,"gp23_tree_with_bootstrap.pdf"), width = 12, height = 17, onefile = FALSE)
bootstrap_tree
dev.off()


## add on
bootstrap_with_points <- bootstrap_tree+ 
 geom_point2(aes(subset=(!grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz", p$data$label) & isTip)),
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
 geom_point2(aes(subset=( grepl("OTU", p$data$label) & isTip)),
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
 geom_cladelabel(node=4363, "B", barsize=2,
                 fontsize=6,
                 angle=0,
                 #offset.text=-0.1,
                 offset = 1.5,
                 color= mypalette[1]) +
 geom_cladelabel( node=4252, "C",                   barsize=2,
                  fontsize=6,
                  angle=0,
                  #offset.text=0.12,
                  offset = 1.425,
                  color= mypalette[2]) + 
 geom_cladelabel( node=4217, "D", barsize=2,
                  fontsize=6,
                  angle=0,
                  #offset.text=-0.1, 
                  color= mypalette[3]) + 
 geom_cladelabel( node=4147, "E",  
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  #offset.text=0.12,
                  offset = 1.4,
                  color= mypalette[4]) + 
 geom_cladelabel( node=3905, "F",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=1.165,
                  # offset.text=-0.1,
                  color= mypalette[5]) +
 geom_cladelabel( node=3295, "G",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=1.3,
                  # offset.text=0.12,
                  color= mypalette[6]) + 
 geom_cladelabel(node=2955,
                 "H", barsize=2,
                 fontsize=6,
                 angle=0,
                 offset=0.62,
                 color= mypalette[7]) + 
 geom_cladelabel( node=2295, "I",
                  barsize=2,
                  fontsize=6,
                  angle=0,
                  offset=0.89,
                  # offset.text=0.12,
                  color= mypalette[8]) 


pdf(paste0(figures_dir,"gp23_tree_with_bootstrap_and_points.pdf"), width = 12, height = 17, onefile = FALSE)
grid.arrange(bootstrap_with_points_and_clade_labels, legend,
             ncol=2,
             widths=c(19/20,1/20))
dev.off()

#+ scale_fill_manual(name="Bar",values=cols)
### Proportional

sum_otus_groupA <- sum_of_group_otus(proportional_gp23, groupA, "A")
sum_otus_groupB <- sum_of_group_otus(proportional_gp23, groupB, "B")
sum_otus_groupC <- sum_of_group_otus(proportional_gp23, groupC, "C")
sum_otus_groupD <- sum_of_group_otus(proportional_gp23, groupD, "D")
sum_otus_groupE <- sum_of_group_otus(proportional_gp23, groupE, "E")
sum_otus_groupF <- sum_of_group_otus(proportional_gp23, groupF, "F")
sum_otus_groupG <- sum_of_group_otus(proportional_gp23, groupG, "G")
sum_otus_groupH <- sum_of_group_otus(proportional_gp23, groupH, "H")
sum_otus_groupI <- sum_of_group_otus(proportional_gp23, groupI, "I")
#  #sum_otus_groupJ <- sum_of_group_otus(proportional_gp23, groupJ, "J") ## only OTU OTU_2428
#  sum_otus_groupK <- sum_of_group_otus(proportional_gp23, groupK, "K")
#  sum_otus_groupL <- sum_of_group_otus(proportional_gp23, groupL, "L")
#  sum_otus_groupM <- sum_of_group_otus(proportional_gp23, groupM, "M")
#  sum_otus_groupN <- sum_of_group_otus(proportional_gp23, groupN, "N")
#  sum_otus_groupO <- sum_of_group_otus(proportional_gp23, groupO, "O")
#  sum_otus_groupP <- sum_of_group_otus(proportional_gp23, groupP, "P")
#  
### want to get the sums of all these groups and then melt and make a plot
l = list(sum_otus_groupA,sum_otus_groupB,sum_otus_groupC, sum_otus_groupD, sum_otus_groupE, sum_otus_groupF, sum_otus_groupG, sum_otus_groupH, sum_otus_groupI
         #,
         #           #sum_otus_groupJ, 
         #           sum_otus_groupK,
         #           sum_otus_groupL,
         #           sum_otus_groupM,
         #           sum_otus_groupN,
         #           sum_otus_groupO,
         #           sum_otus_groupP
)
all_together_gp23 <- Reduce(merge, lapply(l, function(x) data.frame(x, VC = row.names(x))))

write.csv(all_together_gp23, file="../results/gp23_group_sums_by_site_prop.csv")

names(all_together_gp23) <- gsub("sum_of_", "", names(all_together_gp23))

gp23_bar_prop <- barplot_of_top_groups_over_time(all_together_gp23, 
                                                 Jericho_data,
                                                 "testing out gp23 proportional",
                                                 c("grey", mypalette))

gp23_bar_prop_without_high_res <- barplot_of_top_groups_over_time(subset(all_together_gp23,!(VC %in% high_res_vcs)), 
                                                                  Jericho_data,
                                                                  "testing out gp23 proportional",
                                                                  c("grey", mypalette))

gp23_bar_prop_without_high_res <- gp23_bar_prop_without_high_res + annotate("text", x = missing_gp23_dates, y = 0.00, label="*", size=12)

gp23_bar_prop_only_high_res <- barplot_of_top_groups_over_time(all_together_gp23, 
                                                               Jericho_data,
                                                               "testing out gp23 proportional",
                                                               c("grey", mypalette))+
 scale_x_date(breaks = date_breaks("week"), 
              #   labels = date_format("%b"),
              limits = c(as.Date("2011-01-15"),
                         as.Date("2011-02-15")))


pdf(paste0(figures_dir,"gp23_95_ref_groups_over_time_prop.pdf"), width = 30, height = 15, onefile = FALSE)
print(gp23_bar_prop)
dev.off()


gp23_abundance_clusters <- read.csv("../results/normalized_gp23_otus_clusters.csv")
names(gp23_abundance_clusters) <- c("Date", "cluster")
gp23_abundance_clusters$Date <- as.Date(gp23_abundance_clusters$Date)

p1<-gp23_bar_prop
p2<-ggplot(gp23_abundance_clusters,aes(Date,y=1,fill=factor(cluster)))+geom_tile()+
 scale_y_continuous(expand=c(0,0))+  date_scaling + 
 theme_JAG_presentation(base_size=18)+
 theme(axis.title=element_blank(),
       axis.ticks=element_blank(),
       axis.text=element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()
       # ,       legend.position="none"
 )

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)


pdf(paste0(figures_dir,"gp23_95_ref_groups_over_time_with_label_prop.pdf"), width = 30, height = 15, onefile = FALSE)
grid.arrange(gp1, gp2, ncol=1,heights=c(9/10,1/10))
dev.off()


pdf(paste0(figures_dir,"gp23_95_ref_groups_over_time_no_high_res_prop.pdf"), width = 30, height = 15, onefile = FALSE)
print(gp23_bar_prop_without_high_res )
dev.off()

gp23_abundance_clusters_no_high_res <- read.csv("../results/normalized_gp23_otus_clusters_no_high_res.csv")

names(gp23_abundance_clusters_no_high_res) <- c("Date", "cluster")
gp23_abundance_clusters_no_high_res$Date <- as.Date(gp23_abundance_clusters_no_high_res$Date)

p1<-gp23_bar_prop_without_high_res
p2<-ggplot(gp23_abundance_clusters_no_high_res,aes(Date,y=1,fill=factor(cluster)))+geom_tile()+
 scale_y_continuous(expand=c(0,0))+  date_scaling + 
 theme_JAG_presentation(base_size=18)+
 theme(axis.title=element_blank(),
       axis.ticks=element_blank(),
       axis.text=element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()
       # ,       legend.position="none"
 )

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)

pdf(paste0(figures_dir,"gp23_95_ref_groups_over_time_with_label_no_high_res_prop.pdf"), width = 30, height = 15, onefile = FALSE)
grid.arrange(gp1, gp2, ncol=1,heights=c(9/10,1/10))
dev.off()


## add in env clustering label
gp23_env_clusters_no_high_res <- read.csv("../results/normalized_gp23_env_clusters_no_high_res.csv")

names(gp23_env_clusters_no_high_res) <- c("Date", "cluster")
gp23_env_clusters_no_high_res$Date <- as.Date(gp23_env_clusters_no_high_res$Date)
p3<-ggplot(gp23_env_clusters_no_high_res,aes(Date,y=1,fill=factor(cluster)))+geom_tile()+
 scale_y_continuous(expand=c(0,0))+  date_scaling + 
 theme_JAG_presentation(base_size=18)+
 theme(axis.title=element_blank(),
       axis.ticks=element_blank(),
       axis.text=element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()
       # ,       legend.position="none"
 )

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  
gp3<-ggplotGrob(p3)  

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5], gp3$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
gp3$widths[2:5] <- as.list(maxWidth)

pdf(paste0(figures_dir,"gp23_95_ref_groups_over_time_with_label_abun_and_env_no_high_res_prop.pdf"), width = 30, height = 15, onefile = FALSE)
grid.arrange(gp1, gp2,gp3, ncol=1,heights=c(8/10,1/10,1/10))
dev.off()





pdf(paste0(figures_dir,"gp23_95_ref_groups_over_time_only_high_res_prop.pdf"), width = 30, height = 15, onefile = FALSE)
print(gp23_bar_prop_only_high_res )
dev.off()

save_group_OTUs_to_file(proportional_gp23, groupA, "proportional_gp23_groupA.csv") 
save_group_OTUs_to_file(proportional_gp23, groupB, "proportional_gp23_groupB.csv")
save_group_OTUs_to_file(proportional_gp23, groupC, "proportional_gp23_groupC.csv")
save_group_OTUs_to_file(proportional_gp23, groupD, "proportional_gp23_groupD.csv")
save_group_OTUs_to_file(proportional_gp23, groupE, "proportional_gp23_groupE.csv")
save_group_OTUs_to_file(proportional_gp23, groupF, "proportional_gp23_groupF.csv")
save_group_OTUs_to_file(proportional_gp23, groupG, "proportional_gp23_groupG.csv")
save_group_OTUs_to_file(proportional_gp23, groupH, "proportional_gp23_groupH.csv")
save_group_OTUs_to_file(proportional_gp23, groupI, "proportional_gp23_groupI.csv")


## make one with all of the stuff on it for the manuscript
pdf(paste0(figures_dir,"viral_groups_over_time_no_high_res_prop.pdf"), width = 20, height = 17.5, onefile = FALSE)
plot_grid(RdRP_annotated_tree,
          mpl_bar_prop_without_high_res,
          gp23_ref_tree_annotated,
          gp23_bar_prop_without_high_res,
          ncol=2,
          labels = c("A", "B", "C", "D"),
          label_size = 20)
dev.off()




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

## first rename the column names in the OTU table to something better
barplot_of_top_groups_over_time <- function (all_groups_merged, 
                                             Jericho_data,
                                             title,
                                             colours_custom) {
  melt_all_together <- melt(all_groups_merged)
  ## Add in dates:
  melt_all_together$Date <- Jericho_data$Date[match(melt_all_together$VC,
                                                    Jericho_data$VC_number)]
  
  group_barplot <- ggplot(arrange(melt_all_together,
                                  desc(variable)),## this orders the barplot legend in the same way as the barplots
                          aes(x = Date,
                              y = value,
                              group = variable)
  ) + 
    season_line+
    season_text+
   # spring_bloom_line+
    geom_bar(aes(fill = variable),
             stat = "identity",
             position = "stack",
             colour = line_colour)+
    scale_fill_manual(values = colours_custom)+
    theme_JAG_presentation(base_size = 18)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    date_scaling+
    ylab("Relative abundance\n ")+
    xlab("\nDate")+
    guides(fill = guide_legend(title = "group"))
  
  return(group_barplot)
}

### MPL ####

RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762")

ggtree(RdRp_95_ref_tree,
       ladderize = FALSE) +
  geom_text(aes(label=node), hjust=0.3)+
  geom_tiplab()


## so could use these to look at the OTU table
groupA <- tips(RdRp_95_ref_tree,
               1768) ## includes HAKA
groupB <- tips(RdRp_95_ref_tree,
               1764) ## env and Aster glacialis
groupC <- tips(RdRp_95_ref_tree,
               1760) ## could I add to the previous one...
groupD <- tips(RdRp_95_ref_tree,
               1701) ##CUlley seqs, Marine from SF
groupE <- tips(RdRp_95_ref_tree,
               1613) ## Diatom viruses
groupF <- tips(RdRp_95_ref_tree,
               1700) ## mini.....
groupG <- tips(RdRp_95_ref_tree,
               1450) ## low bootstrap...
groupH <- tips(RdRp_95_ref_tree,
               1017)
## leftover tips are I
##all_tips_except_I <- c(groupA, groupB, groupC, groupD, groupE, groupF, groupG, groupH)
##RdRp_95_ref_tree$tip.label[!RdRp_95_ref_tree$tip.label %in% all_tips_except_I]


groupA_data_frame <- data.frame("group" = "A",
                                values = groupA)
groupB_data_frame <- data.frame("group" = "B",
                                values = groupB)
groupC_data_frame <- data.frame("group" = "C",
                                values = groupC)
groupD_data_frame <- data.frame("group" = "D",
                                values = groupD)
groupE_data_frame <- data.frame("group"="E",
                                values = groupE)
groupF_data_frame <- data.frame("group" = "F",
                                values = groupF)
groupG_data_frame <- data.frame("group" = "G",
                                values = groupG)
groupH_data_frame <- data.frame("group" = "H",
                                values = groupH)

group_data_frame <- rbind(groupA_data_frame,
                          groupB_data_frame)
group_data_frame <- rbind(group_data_frame,
                          groupC_data_frame)
group_data_frame <- rbind(group_data_frame,
                          groupD_data_frame)
group_data_frame <- rbind(group_data_frame,
                          groupE_data_frame)
group_data_frame <- rbind(group_data_frame,
                          groupF_data_frame)
group_data_frame <- rbind(group_data_frame,
                          groupG_data_frame)
group_data_frame <- rbind(group_data_frame,
                          groupH_data_frame)

## make table for annotation for use in script 46
write.csv(group_data_frame,
          "../results/RdRp_groups_with_OTUs.csv") 

## missing OTUs 119 -node splits from 1014 and 377 from slip 1015

mypalette<-brewer.pal(7,
                      "Dark2")

p <- ggtree(RdRp_95_ref_tree,
            ladderize = FALSE)+
  geom_text(aes(label = node)) # just to look at labelled nodes

p <- ggtree(RdRp_95_ref_tree,
            ladderize = FALSE,
            colour = line_colour)+
  scale_y_continuous(labels = NULL)+
  scale_x_continuous(labels = NULL)+
  theme_JAG_presentation()+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) 

pdf(paste0(figures_dir,"RdRp_95_ref_free_plain.pdf"),
    width = 30,
    height = 15,
    onefile = FALSE)
p
dev.off()

p2 <- p+ geom_cladelabel(node = 1768,
                         "A",
                         angle = 0,
                         offset.text = 0.12,
                         barsize = 4, 
                         fontsize = 6,
                         offset = 2,
                         color = mypalette[1]) + 
  geom_cladelabel(node = 1764,
                  "B",
                  angle = 0,
                  offset.text = 0.12,
                  barsize = 4,
                  fontsize = 6,
                  offset = 1.9242,
                  color = "grey") + 
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

pdf(paste0(figures_dir,
           "RdRp_95_ref_free_annotated_with_groups.pdf"),
    width = 30,
    height = 15,
    onefile = FALSE)

p2

dev.off()



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
p2a <- ggtree(tree,
              aes(color = group),
              ladderize = FALSE)+
  scale_color_manual(values = c(line_colour,
                                mypalette[1],
                                "grey",
                                mypalette[2:7]),
                     guide = FALSE)+
  scale_y_continuous(labels = NULL)+
  scale_x_continuous(labels = NULL)+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) 

p3a <-  p2a+ geom_cladelabel(node = 1768,
                             "A",
                             angle = 0,
                             offset.text = 0.12,
                             barsize = 4, 
                             fontsize = 6,
                             offset = 2,
                             color = mypalette[1]) + 
  geom_cladelabel(node = 1764,
                  "B",
                  angle = 0,
                  offset.text = 0.12,
                  barsize = 4,
                  fontsize = 6,
                  offset = 1.9242,
                  color = "grey") + 
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
RdRP_annotated_tree <- p3a
RdRP_annotated_tree

p3 <- p2+
  geom_point2(aes(subset = (!grepl("OTU|Culley|0\\.[:digit:]+",
                                   RdRP_annotated_tree$data$label)& isTip),
                  label = label),
              size = 4,
              # hjust = 2,
              colour = "steelblue") +
  geom_point2(aes(subset = (grepl("OTU",
                                  RdRP_annotated_tree$data$label) & isTip)),
              color = "#FDAC4F",
              shape = 8,
              size = 2) +
  geom_point2(aes(subset = (grepl("Culley",
                                  RdRP_annotated_tree$data$label) & isTip)),
              color = "purple",
              shape = 8,
              size = 2)

pdf(paste0(figures_dir,"RdRp_95_ref_free_annotated_with_groups_and_env_ref_sequences.pdf"),
    width = 30,
    height = 15,
    onefile = FALSE)

print(p3)

dev.off()

p4 <- p2+ geom_text2(aes(subset = ( !grepl("OTU|Culley|0\\.[:digit:]+",
                                           RdRp_95_ref_tree$tip.label)  & isTip),
                         label = label),
                     color = "grey30",
                     size = 2.5,
                     hjust = -.1) +
  geom_point2(aes(subset = (grepl("OTU",
                                  RdRP_annotated_tree$data$label) & isTip)),
              color = "#FDAC4F",
              shape = 8,
              size = 2) +
  geom_point2(aes(subset = (grepl("Culley",
                                  RdRP_annotated_tree$data$label) & isTip)),
              color = "purple",
              shape = 8,
              size = 2)
p4

pdf(paste0(figures_dir,
           "RdRp_95_ref_free_text_of_refs_groups_and_env_sequences.pdf"),
    width = 30,
    height = 40,
    onefile = FALSE)

print(p4)

dev.off()

(p + geom_text(aes(label = node)) +geom_tiplab())

p5 <- p2+ geom_text2(aes(subset = ( !grepl("OTU|Culley|0\\.[:digit:]+",
                                           RdRp_95_ref_tree$tip.label)  & isTip),
                         label = label),
                     color = "grey30",
                     size = 2.5,
                     hjust = -.1)  +
  geom_point2(aes(subset = (grepl("OTU",
                                  RdRP_annotated_tree$data$label) & isTip)),
              color = "#FDAC4F",
              shape = 8,
              size = 2)  +
  geom_point2(aes(subset = (grepl("Culley",
                                  RdRP_annotated_tree$data$label) & isTip)),
              color = "purple",
              shape = 8,
              size = 2)
p5

pdf(paste0(figures_dir,
           "RdRp_95_ref_free_text_of_env_seqs_groups_and_refs_sequences.pdf"),
    width = 30,
    height = 40,
    onefile = FALSE)

print(p5)

dev.off()


p6 <- p3a+ geom_text2(aes(subset = ( !grepl("OTU|Culley|0\\.[:digit:]+",
                                           RdRp_95_ref_tree$tip.label)  & isTip),
                         label = label),
                     color = "grey30",
                     size = 2.5,
                     hjust = -.1)  +
  geom_point2(aes(subset = (grepl("Culley",
                                  RdRP_annotated_tree$data$label) & isTip)),
              color = "purple",
              shape = 8,
              size = 2)
p6

pdf(paste0(figures_dir,
           "RdRp_95_ref_free_text_and_refs_sequences.pdf"),
    width = 15,
    height = 20,
    onefile = FALSE)

print(p6)

dev.off()

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

## get only legend

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
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

## would like to summarize somehow....
sum_of_group_otus <- function (normalized_OTUs,
                               group,
                               group_string) {
  group_only_otu <- subset(group,
                           grepl("OTU",
                                 group) )
  test_group <- subset(normalized_OTUs,
                       select = colnames(normalized_OTUs) %in% group_only_otu)
  test_row <- as.data.frame(rowSums(test_group))
  names(test_row)[1] <- paste0("sum_of_",
                               group_string)
  test <- adply(test_group,
                1,
                sum)
  names(test)[3] <- paste0("sum_of_",
                           group_string )
  return(test_row)
}


### propotional reads

sum_otus_groupA <- sum_of_group_otus(proportional_MPL,
                                     groupA,
                                     "A")
sum_otus_groupB <- sum_of_group_otus(proportional_MPL,
                                     groupB,
                                     "B")
#sum_otus_groupC <- sum_of_group_otus(proportional_MPL,
#groupC,
#"C")
sum_otus_groupD <- sum_of_group_otus(proportional_MPL,
                                     groupD, 
                                     "D")
sum_otus_groupE <- sum_of_group_otus(proportional_MPL,
                                     groupE,
                                     "E")
# sum_otus_groupF <- sum_of_group_otus(proportional_MPL, groupF, "F")
sum_otus_groupG <- sum_of_group_otus(proportional_MPL,
                                     groupG,
                                     "G")
sum_otus_groupH <- sum_of_group_otus(proportional_MPL,
                                     groupH,
                                     "H")
#sum_otus_groupI <- sum_of_group_otus(proportional_MPL, groupI, "I")
# sum_otus_groupJ <- sum_of_group_otus(proportional_MPL, groupJ, "J")
#sum_otus_groupK <- sum_of_group_otus(proportional_MPL, groupK, "K")

### want to get the sums of all these groups and then melt and make a plot
l = list(sum_otus_groupA,
         sum_otus_groupB,
         #sum_otus_groupC,
         sum_otus_groupD, 
         sum_otus_groupE, 
         #sum_otus_groupF, 
         sum_otus_groupG, 
         sum_otus_groupH
         #, sum_otus_groupI
         #,
         # sum_otus_groupJ,
         #sum_otus_groupK
)
all_together_MPL <- Reduce(merge,
                           lapply(l, function(x) data.frame(x, VC = row.names(x))))

write.csv(all_together_MPL,
          file="../results/MPL_group_sums_by_site_prop.csv")

names(all_together_MPL) <- gsub("sum_of_",
                                "",
                                names(all_together_MPL))

mpl_bar_prop_without_high_res <- barplot_of_top_groups_over_time(subset(all_together_MPL,!(VC %in% high_res_vcs)), 
                                                                 Jericho_data,
                                                                 "MPL proportional",
                                                                 c(mypalette[1],"grey",
                                                                   mypalette[c(-1,-2,-5)]))

mpl_bar_prop_without_high_res <- mpl_bar_prop_without_high_res+
  annotate("text",
           x = missing_MPL_dates,
           y = 0.00,
           label="*",
           size=14)

pdf(paste0(figures_dir,
           "RdRp_95_ref_groups_over_time_no_high_res_prop.pdf"),
    width = 30,
    height = 15,
    onefile = FALSE)

print(mpl_bar_prop_without_high_res)

dev.off()

## otus tables from the groups 

save_group_OTUs_to_file <- function (proportional_table,
                                     group,
                                     table_file_name) {
  proportional_table_group <- subset(proportional_table,
                                     select = colnames(proportional_table) %in% group)
  write.csv(proportional_table_group,
            file = paste0("../results/",
                          table_file_name))
}

save_group_OTUs_to_file(proportional_MPL,
                        groupA,
                        "proportional_MPL_groupA.csv") 
save_group_OTUs_to_file(proportional_MPL,
                        groupB,
                        "proportional_MPL_groupB.csv") 
save_group_OTUs_to_file(proportional_MPL,
                        groupC,
                        "proportional_MPL_groupC.csv") 
save_group_OTUs_to_file(proportional_MPL,
                        groupD,
                        "proportional_MPL_groupD.csv") 
save_group_OTUs_to_file(proportional_MPL,
                        groupE,
                        "proportional_MPL_groupE.csv") 
save_group_OTUs_to_file(proportional_MPL,
                        groupF,
                        "proportional_MPL_groupF.csv") 
save_group_OTUs_to_file(proportional_MPL,
                        groupG,
                        "proportional_MPL_groupG.csv") 
save_group_OTUs_to_file(proportional_MPL,
                        groupH, 
                        "proportional_MPL_groupH.csv") 


### gp23 ####
mypalette<-brewer.pal(9,
                      "Set3")
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

pdf(paste0(figures_dir,"gp23_ref_free_plain.pdf"),
    width = 30,
    height = 15,
    onefile = FALSE)
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
   
pdf(paste0(figures_dir,"gp23_95_ref_free_annotated_with_groups.pdf"), width = 30, height = 15, onefile = FALSE)
p2
dev.off()


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
 # theme_JAG_presentation()+
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

p3 <- p2+
  geom_point2(aes(subset=( !grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz", p2$data$label) & isTip),
                          label=label),
                      color="blue",
                      size=4) +
  geom_point2(aes(subset=( grepl("OTU", p2$data$label) & isTip)),
              color="#FDAC4F",
              shape=8,
              size=2 ) +
  geom_point2(aes(subset=( grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                 p2$data$label) & isTip)),
              color="purple",
              shape=8,
              size=2)

pdf(paste0(figures_dir,
           "gp23_95_ref_free_annotated_with_groups_and_env_ref_sequences.pdf"),
    width = 30,
    height = 15,
    onefile = FALSE)

print(p3)

dev.off()


p4 <-p2+ geom_text2(aes(subset=(!grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                       p2$data$label) & isTip),
                        label=label),
                    color="grey30",
                    size=2.5,
                    hjust = -.1) +
  geom_point2(aes(subset=( grepl("OTU", p2$data$label) & isTip)),
              color="#FDAC4F",
              shape=8,
              size=2 ) +
  geom_point2(aes(subset=( grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                 p2$data$label) & isTip)),
              color="purple",
              shape=8,
              size=2 )
p4

pdf(paste0(figures_dir,"gp23_95_ref_free_text_of_refs_groups_and_env_sequences.pdf"),
    width = 30,
    height = 40,
    onefile = FALSE)
print(p4)
dev.off()

p4_5 <-p3a+ geom_text2(aes(subset=(!grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                       p2$data$label) & isTip),
                        label=label),
                    color="grey30",
                    size=2.5,
                    hjust = -.1) +
  geom_point2(aes(subset=( grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                 p2$data$label) & isTip)),
              color="purple",
              shape=8,
              size=2 )
p4_5

pdf(paste0(figures_dir,"gp23_95_ref_free_text_of_refs_groups_and_env_sequences.pdf"),
    width = 30,
    height = 40,
    onefile = FALSE)
print(p4_5)
dev.off()

p5 <- p2+
  geom_text2(aes(subset=(!grepl("OTU|Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                        p2$data$label) & isTip),
                         label=label),
                     color="grey30",
                     size=2.5,
                     hjust = -.1) +
  geom_point2(aes(subset=( grepl("OTU", p2$data$label) & isTip)),
              color="#FDAC4F",
              shape=8,
              size=0.5)+
  geom_text2(aes(subset=( grepl("Chow|Filee|Jia|Lopez-Bueno|LopezBueno|Bellas|Knapik|Liu|Butina|Comeau|Sandaa|Mabiz",
                                p2$data$label) & isTip),
                 label=label),
             color="darkslateblue",
             size=2.5,
             hjust = -0.1)

p5

pdf(paste0(figures_dir,"gp23_95_ref_free_text_of_env_seqs_groups_and_refs_sequences.pdf"),
    width = 30,
    height = 40,
    onefile = FALSE)
print(p5)
dev.off()


(p + geom_text(aes(label=node)) +geom_tiplab())



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


## would be good to hold all the nodes somewhere. 
## so should rearrange all these plots somewhere. 


#+ scale_fill_manual(name="Bar",values=cols)
### Proportional

sum_otus_groupA <- sum_of_group_otus(proportional_gp23,
                                     groupA,
                                     "A")
sum_otus_groupB <- sum_of_group_otus(proportional_gp23,
                                     groupB,
                                     "B")
sum_otus_groupC <- sum_of_group_otus(proportional_gp23,
                                     groupC,
                                     "C")
sum_otus_groupD <- sum_of_group_otus(proportional_gp23,
                                     groupD,
                                     "D")
sum_otus_groupE <- sum_of_group_otus(proportional_gp23,
                                     groupE,
                                     "E")
sum_otus_groupF <- sum_of_group_otus(proportional_gp23,
                                     groupF,
                                     "F")
sum_otus_groupG <- sum_of_group_otus(proportional_gp23,
                                     groupG,
                                     "G")
sum_otus_groupH <- sum_of_group_otus(proportional_gp23,
                                     groupH,
                                     "H")
sum_otus_groupI <- sum_of_group_otus(proportional_gp23,
                                     groupI,
                                     "I")
#  #sum_otus_groupJ <- sum_of_group_otus(proportional_gp23, groupJ, "J") ## only OTU OTU_2428
#  sum_otus_groupK <- sum_of_group_otus(proportional_gp23, groupK, "K")
#  sum_otus_groupL <- sum_of_group_otus(proportional_gp23, groupL, "L")
#  sum_otus_groupM <- sum_of_group_otus(proportional_gp23, groupM, "M")
#  sum_otus_groupN <- sum_of_group_otus(proportional_gp23, groupN, "N")
#  sum_otus_groupO <- sum_of_group_otus(proportional_gp23, groupO, "O")
#  sum_otus_groupP <- sum_of_group_otus(proportional_gp23, groupP, "P")
#  
### want to get the sums of all these groups and then melt and make a plot
l = list(sum_otus_groupA,
         sum_otus_groupB,
         sum_otus_groupC,
         sum_otus_groupD,
         sum_otus_groupE,
         sum_otus_groupF,
         sum_otus_groupG,
         sum_otus_groupH,
         sum_otus_groupI
         #,
         #           #sum_otus_groupJ, 
         #           sum_otus_groupK,
         #           sum_otus_groupL,
         #           sum_otus_groupM,
         #           sum_otus_groupN,
         #           sum_otus_groupO,
         #           sum_otus_groupP
)
all_together_gp23 <- Reduce(merge,
                            lapply(l, function(x) data.frame(x, VC = row.names(x))))

write.csv(all_together_gp23,
          file="../results/gp23_group_sums_by_site_prop.csv")

names(all_together_gp23) <- gsub("sum_of_",
                                 "",
                                 names(all_together_gp23))

gp23_bar_prop <- barplot_of_top_groups_over_time(all_together_gp23, 
                                                 Jericho_data,
                                                 "testing out gp23 proportional",
                                                 c("grey",
                                                   mypalette))

gp23_bar_prop_without_high_res <- barplot_of_top_groups_over_time(subset(all_together_gp23,!(VC %in% high_res_vcs)), 
                                                                  Jericho_data,
                                                                  "testing out gp23 proportional",
                                                                  c("grey",
                                                                    mypalette))

gp23_bar_prop_without_high_res <- gp23_bar_prop_without_high_res +
  annotate("text",
           x = missing_gp23_dates,
           y = 0.00,
           label="*",
           size=14)

gp23_bar_prop_only_high_res <- barplot_of_top_groups_over_time(all_together_gp23, 
                                                               Jericho_data,
                                                               "testing out gp23 proportional",
                                                               c("grey",
                                                                 mypalette))+
  scale_x_date(breaks = date_breaks("week"), 
               #   labels = date_format("%b"),
               limits = c(as.Date("2011-01-15"),
                          as.Date("2011-02-15")))


pdf(paste0(figures_dir,
           "gp23_95_ref_groups_over_time_prop.pdf"),
    width = 30,
    height = 15,
    onefile = FALSE)
print(gp23_bar_prop)
dev.off()


pdf(paste0(figures_dir,"gp23_95_ref_groups_over_time_no_high_res_prop.pdf"), width = 30, height = 15, onefile = FALSE)
print(gp23_bar_prop_without_high_res)
dev.off()

pdf(paste0(figures_dir,
           "gp23_95_ref_groups_over_time_only_high_res_prop.pdf"),
    width = 30,
    height = 15,
    onefile = FALSE)

print(gp23_bar_prop_only_high_res)

dev.off()

save_group_OTUs_to_file(proportional_gp23,
                        groupA,
                        "proportional_gp23_groupA.csv") 
save_group_OTUs_to_file(proportional_gp23,
                        groupB, 
                        "proportional_gp23_groupB.csv")
save_group_OTUs_to_file(proportional_gp23,
                        groupC, 
                        "proportional_gp23_groupC.csv")
save_group_OTUs_to_file(proportional_gp23,
                        groupD,
                        "proportional_gp23_groupD.csv")
save_group_OTUs_to_file(proportional_gp23,
                        groupE,
                        "proportional_gp23_groupE.csv")
save_group_OTUs_to_file(proportional_gp23,
                        groupF,
                        "proportional_gp23_groupF.csv")
save_group_OTUs_to_file(proportional_gp23,
                        groupG,
                        "proportional_gp23_groupG.csv")
save_group_OTUs_to_file(proportional_gp23,
                        groupH,
                        "proportional_gp23_groupH.csv")
save_group_OTUs_to_file(proportional_gp23,
                        groupI,
                        "proportional_gp23_groupI.csv")

dev.off()
## make one with all of the stuff on it for the manuscript
pdf(paste0(figures_dir,
           "viral_groups_over_time_no_high_res_prop.pdf"),
    width = 20,
    height = 17.5, 
    onefile = FALSE)

plot_grid(RdRP_annotated_tree,
          mpl_bar_prop_without_high_res,
          gp23_ref_tree_annotated,
          gp23_bar_prop_without_high_res,
          ncol=2,
          labels = c("A", "B", "C", "D"),
          label_size = 20)

dev.off()
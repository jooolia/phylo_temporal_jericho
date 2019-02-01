## plots ordered by groups to replace the heatmaps

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
library(dplyr)
library(tidyr)

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

## margin 1 is by row
normalized_MPL_OTUs <- as.data.frame(prop.table(as.matrix(normalized_MPL_OTUs), 1))*100

#normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")

normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")

normalized_gp23_OTUs <- as.data.frame(prop.table(as.matrix(normalized_gp23_OTUs), 1))*100

#normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number")


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
  pdf(plot_tree_file_name,
      width = 6
      # ,
      # height = 11
      )
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
  pdf(plot_tree_file_name
      ,width = 6
      # , height = 11 
      )
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

### MPL ##### 

#### MPL With references ####

RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RAxML_bipartitions.RdRptree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                                 paste0(figures_dir,"RdRp_95_miseq_data_with_env_iso_and_ref_Raxml.pdf"))

MPL_otu_table_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs)

high_res_vcs <- c(1198, 1199, 1200, 1201, 1202)
normalized_MPL_OTUs_no_high_res <- subset(normalized_MPL_OTUs, !(rownames(normalized_MPL_OTUs) %in% high_res_vcs))
normalized_MPL_OTUs_no_high_res_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs_no_high_res)

## only high res

normalized_MPL_OTUs_only_high_res <- subset(normalized_MPL_OTUs, rownames(normalized_MPL_OTUs) %in% high_res_vcs)

normalized_MPL_OTUs_only_high_res_with_isolates <- add_isolates_to_otu_table(RdRp_95_ref_tree, normalized_MPL_OTUs_only_high_res)




## 4 differenct ways:
#1 - richness of each group
#2 - relative abundance of each group
#3- richness of each group coloured by relative abundance
#4- relative abundance with points coloured by richness

## need to add in the groups...

## where is I? there is no I in the RdRp!!!!!
## it is all of those that are not annotated, is that dangerous??
RdRp_OTU_annotations <- read.csv("../results/RdRp_groups_with_OTUs.csv")

mypalette<-brewer.pal(7,"Dark2")
MPL_colour <-  c(mypalette[1],
                 "grey",
                 mypalette[2:7],
                 "white")


OTU_table_with_isolates <- add_isolates_to_otu_table(
  RdRp_95_ref_tree,
  normalized_MPL_OTUs_no_high_res_with_isolates)
tree <- RdRp_95_ref_tree
  OTU_table_long <- melt(as.matrix(OTU_table_with_isolates))
  OTU_table_long$value[OTU_table_long$value == 0] <- NA
  OTU_table_long$Date <- Jericho_data$Date[match(OTU_table_long$Var1, Jericho_data$VC_number)]

  OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(tree$tip.label))
  ## need to add in "I" groups
  OTU_table_long <- merge(OTU_table_long, RdRp_OTU_annotations, by.x = "Var2", by.y = "values", all.x = TRUE)
  OTU_table_long$Date <- as.Date(OTU_table_long$Date, origin = "1970-01-01")  
  str(OTU_table_long$value)



    
OTU_table_summarised <- OTU_table_long %>% 
  group_by(Date, group) %>%
  filter(value > 0) %>% 
  summarise(richness = n(), sum_rel_abun=sum(value, na.rm = TRUE)) %>% 
  complete(Date,group, fill = list(richness = 0, sum_rel_abun=0)) ## this keeps the 0 observation groups
    
## exclude groups C and F
OTU_table_summarised <- filter(OTU_table_summarised, group != c("C", "F"))

Jericho_no_high_res_vcs <- Jericho_data$VC_number[!(Jericho_data$VC_number %in% high_res_vcs)]

## want to annotate missing data days. 
## want to annotate those days using star
missing_MPL_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_MPL_OTUs))]
missing_MPL_dates <- Jericho_data$Date[match(missing_MPL_samples,
                                             Jericho_data$VC_number)]

missing_gp23_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_gp23_OTUs))]
missing_gp23_dates <- Jericho_data$Date[match(missing_gp23_samples,
                                              Jericho_data$VC_number)]



OTU_table_summarised <- OTU_table_summarised[!is.na(OTU_table_summarised$group),]

## for each missing dates need to have for each group
## inelegant solution, btu it wokrs for now. 
missing_MPL_dates_df <- c()
for (group in unique(OTU_table_summarised$group)){
  print(group)
  for (date in missing_MPL_dates){
    # cat(date)
    df <- data.frame(Date = as.Date(date, origin = "1970-01-01"),
                     group = group,
                     richness = NA,
                     sum_rel_abun = NA)
    missing_MPL_dates_df <- rbind(missing_MPL_dates_df, df)
  }
}

OTU_table_summarised <- bind_rows(OTU_table_summarised, missing_MPL_dates_df)



## plot 1 richness
  
  # ggplot(OTU_table_summarised , 
  #        aes(x=Date,y=richness)) + 
  #   geom_point()+
  #   facet_grid(group~., scales = "free_y")+
  #   date_scaling +
  #   theme_JAG_presentation(base_size = 18)+
  #   ## colour of facet titles?
  #   theme(strip.background =element_rect(fill=MPL_colour))

  ## plot 2 richness with relative abundance
  # 
  # ggplot(OTU_table_summarised , 
  #        aes(x=Date,y=richness)) + 
  #   geom_point(aes(colour=sum_rel_abun))+
  #   facet_grid(group~., scales = "free_y")+
  #   date_scaling +
  #   theme_JAG_presentation(base_size = 18)+
  #   ## colour of facet titles?
  #   theme(strip.background =element_rect(fill=MPL_colour))

  ## plot 3 relative abundance maybe this one?
  ## need the relative abundace to be % 

#ann_text <- data.frame(Date = missing_MPL_dates,sum_rel_abun = 0, lab = "*", group = "H")
pdf(paste0(figures_dir,
               "MPL_richness_over_time_by_group.pdf"),
    width = 5)
  ggplot(OTU_table_summarised , 
         aes(x = Date,
             y = richness)) + 
    season_line+
    geom_point()+
    annotate("text",
             x = missing_MPL_dates,
             y = -0.1,
             label="x",
             size=3,
             colour = "grey")+
    facet_grid(group~., scales = "free_y")+
    date_scaling +
    #theme_JAG_presentation(base_size = 15)+
    theme_JAG_presentation(base_size = 7)+
    theme(strip.background = element_rect(fill = "white"))+
   # theme_bw()+
   # expand_limits(y = 0)+

    ylab("Number of OTUs observed")
  dev.off()
  
  

blank_data <-   OTU_table_summarised %>% 
    group_by(group) %>% 
    summarise(richness = max(richness, na.rm = TRUE)*1.1) %>% 
      mutate(Date =  as.Date("2010-06-22", origin = "1970-01-01") )
  

  pdf(paste0(figures_dir,
             "MPL_richness_over_time_by_group_with_rel_abun.pdf"),
      width = 6,
      height = 4)  
  ggplot(OTU_table_summarised , 
         aes(x=Date,y=richness))+
    geom_hline(yintercept=0, colour = "gray")+
    season_line+
    geom_line()+
    geom_point(aes(colour = sum_rel_abun), size = 2)+
    annotate("text",
             x = missing_MPL_dates,
             y = 0,
             label="x",
             size=3,
             colour = "grey")+
    geom_blank(data = blank_data, aes(x = as.Date(Date, origin = "1970-01-01") , y = richness))+ 
    facet_grid(group~., scales = "free_y")+
    date_scaling +
    scale_colour_viridis_c(name = "Relative\nAbundance\n(%)")+
    theme_JAG_presentation(base_size = 7)+
    theme(strip.background = element_rect(fill = "white"))+
   # theme_bw()+
    ylab("Number of OTUs observed")
dev.off()  



### gp23  ####


gp23_95_ref_tree <- plot_tree_root_and_ladderize("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result",
                                                 root = "Enterobacteria_phage_T4",
                                                 paste0(figures_dir,"gp23_95_miseq_data_with_env_iso_and_ref_Raxml.pdf"))


gp23_OTU_annotations <- read.csv("../results/gp23_groups_with_OTUs.csv")



gp23_otu_table_with_isolates <- add_isolates_to_otu_table(gp23_95_ref_tree, normalized_gp23_OTUs)


normalized_gp23_OTUs_no_high_res <- subset(normalized_gp23_OTUs, !(rownames(normalized_gp23_OTUs) %in% high_res_vcs))


mypalette<-brewer.pal(9,"Set3")
gp23_colour <- c("grey", mypalette, "white") 



OTU_table_with_isolates <- add_isolates_to_otu_table(gp23_95_ref_tree, normalized_gp23_OTUs_no_high_res)
tree <- gp23_95_ref_tree
OTU_table_long <- melt(as.matrix(OTU_table_with_isolates))
OTU_table_long$value[OTU_table_long$value == 0] <- NA
OTU_table_long$Date <- Jericho_data$Date[match(OTU_table_long$Var1, Jericho_data$VC_number)]

OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(tree$tip.label))
## need to add in "I" groups
OTU_table_long <- merge(OTU_table_long, gp23_OTU_annotations, by.x = "Var2", by.y = "values", all.x = TRUE)
OTU_table_long$Date <- as.Date(OTU_table_long$Date, origin = "1970-01-01")  
str(OTU_table_long$value)

## exclude NA group


OTU_table_summarised <- OTU_table_long %>% 
  group_by(Date, group) %>%
  filter(value > 0) %>% 
  summarise(richness = n(), sum_rel_abun=sum(value, na.rm = TRUE)) %>% 
  complete(Date,group, fill = list(richness = 0, sum_rel_abun=0)) ## this keeps the 0 observation groups

OTU_table_summarised <- OTU_table_summarised[!is.na(OTU_table_summarised$group),]


## for each missing dates need to have for each group
## inelegant solution, btu it wokrs for now. 
missing_gp23_dates_df <- c()
 for (group in unique(OTU_table_summarised$group)){
  print(group)
for (date in missing_gp23_dates){
 # cat(date)
  df <- data.frame(Date = as.Date(date, origin = "1970-01-01"),
             group = group,
             richness = NA,
             sum_rel_abun = NA)
  missing_gp23_dates_df <- rbind(missing_gp23_dates_df, df)
}
}

OTU_table_summarised <- bind_rows(OTU_table_summarised, missing_gp23_dates_df)


pdf(paste0(figures_dir,
           "gp23_richness_over_time_by_group.pdf"),
    width = 6)
ggplot(OTU_table_summarised , 
       aes(x=Date,y=richness)) + 
  season_line+
  geom_hline(yintercept=0)+
  geom_point()+
  geom_line()+
  annotate("text",
           x = missing_gp23_dates,
           y = -0.1,
           label="x",
           size=3,
           colour = "grey")+
  facet_grid(group~., scales = "free_y")+
  date_scaling +
  #theme_JAG_presentation(base_size = 15)+
  theme_bw()+
  # expand_limits(y = 0)+
  ylab("Number of OTUs observed")
dev.off()

blank_data <-   OTU_table_summarised %>% 
  group_by(group) %>% 
  summarise(richness = max(richness, na.rm = TRUE)*1.1) %>% 
  mutate(Date =  as.Date("2010-06-22", origin = "1970-01-01") )


pdf(paste0(figures_dir,
           "gp23_richness_over_time_by_group_with_rel_abun.pdf"),
    width = 6,
    height = 6)
ggplot(OTU_table_summarised , 
       aes(x=Date,y=richness)) + 
  season_line+
  geom_hline(yintercept=0, colour = "gray")+
  geom_line()+
  geom_point(aes(colour = sum_rel_abun), size = 2)+
  annotate("text",
           x = missing_gp23_dates,
           y = 0,
           label="x",
           size=3,
           colour = "grey")+
  geom_blank(data = blank_data, aes(x = as.Date(Date, origin = "1970-01-01") , y = richness))+ 
  facet_grid(group~., scales = "free_y")+
  date_scaling +
  scale_colour_viridis_c(name = "Relative\nAbundance\n(%)", expand = expand_scale(mult = c(0.1, 2)))+
  theme_JAG_presentation(base_size = 7)+
  theme(strip.background=element_rect(fill="white"))+
  #theme_linedraw()+
  ylab("Number of OTUs observed")
dev.off()

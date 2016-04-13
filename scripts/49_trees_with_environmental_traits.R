## trees with traits

library(ape)
library(reshape2)
library(ggplot2)
library(phylobase)
library(adephylo)
library(RColorBrewer)
library(gridExtra)
library(plyr)
library(tabplot)

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number")

if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
 path_colour <- "black"
 line_colour <- "black"
 point_colour <- "black"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
  point_colour <- "white"
 }
}



## so what I want to do is generate data from the Jericho env to visualize under what conditions an OTU was seen. 
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)


make_env_plots_based_on_tree <- function (amplicon_tree,
                                          OTU_table,
                                          amplicon_name) {
  
#   convert_otu_table_to_env_mean <- function (OTU_table, env_param, env_col_name ) {
#     ## add in the isolates
#    
#    OTU_table_with_isolates <- OTU_table
#   #  OTU_table_with_isolates <- normalized_MPL_OTUs
#   #  OTU_table <- normalized_MPL_OTUs
#   #  env_param <- Jericho_data$Temperature_YSI
#   #  env_col_name <- "average_temp"
#   #  
#     sample_length <- length(rownames(OTU_table_with_isolates))
#     
#     #Make vector of 0s based on number of sites. 
#     zeros <- rep(0, sample_length)
#     
#     #Insert vector in table based on number of tips in tree that are not found in the OTU table. 
#     OTU_table_with_isolates  <- lapply(amplicon_tree$tip.label,function(i){
#      if (is.element(i, colnames(OTU_table_with_isolates))) {
#       zero_df <- c()
#       #break
#      }
#      else {
#       print("no")
#       print(i)
#       zero_df <- data.frame(zeros) 
#       colnames(zero_df) <- paste(i)
#       return(zero_df) 
#      }
#      })
#     ## remove the empty from the lists
#     OTU_table_with_isolates_non_null <- OTU_table_with_isolates[!sapply(OTU_table_with_isolates,is.null)]
#     OTU_table_with_isolates1 <- do.call(cbind,OTU_table_with_isolates_non_null)
#     OTU_table_with_OTUs_and_isolates <- cbind(OTU_table, OTU_table_with_isolates1)
#     
#     OTU_table_with_isolates_with_env <-   as.data.frame(apply(OTU_table_with_OTUs_and_isolates,
#                                                               c(1,2),
#                                                               function(x){
#                                                                ## replace each value above 0 with the
#                                                                ## value of the environmental param
#      replace(x,
#              x > 0,
#              env_param[match(rownames(OTU_table_with_OTUs_and_isolates)[x], Jericho_data$VC)])
#     }) )
#   
#   
#     ## to get means to work on available samples
#     ## so instead of 0 make sure it is NA. Better for plotting. 
#     OTU_table_with_isolates_with_env[OTU_table_with_isolates_with_env  == 0] <- NA
#     av_temp_OTU <- as.data.frame(colMeans(OTU_table_with_isolates_with_env , na.rm=TRUE))
#   
#     #### order matrix so that it is ordered like the tree ####
#     colnames(av_temp_OTU) <- env_col_name
#     av_temp_OTU$otu_name <- rownames(av_temp_OTU) # don't lose rownames if I have 2 rows
#     av_temp_OTU_reordered <- av_temp_OTU[rev(amplicon_tree$tip.label),]
#     return(av_temp_OTU_reordered)
#   }
#   
#   
#   
#   MPL_otus_av_temp <- convert_otu_table_to_env_mean(OTU_table, Jericho_data$Temperature_YSI, "average_temp")
#   MPL_otus_av_sali <- convert_otu_table_to_env_mean(OTU_table, Jericho_data$Salinity_ppt_YSI , "average_sal")
#   MPL_otus_av_chla <- convert_otu_table_to_env_mean(OTU_table, Jericho_data$Average_chl_a, "average_chla" )
#   MPL_otus_av_VA <- convert_otu_table_to_env_mean(OTU_table, Jericho_data$Average_viral_abundance, "average_VA" )
#   MPL_otus_av_BA <- convert_otu_table_to_env_mean(OTU_table, Jericho_data$Average_bacterial_abundance, "average_BA")
#   MPL_otus_av_PO4 <- convert_otu_table_to_env_mean(OTU_table, Jericho_data$Average_PO4, "average_PO4")
#   MPL_otus_av_NO <- convert_otu_table_to_env_mean(OTU_table, Jericho_data$Average_NO3_NO2, "average_NO")
#   MPL_otus_av_SI <- convert_otu_table_to_env_mean(OTU_table, Jericho_data$Average_SiO2, "average_SI")
#   
#   ### make into function and do with many others
#   
#   all_env_ordered_by_tree <- cbind(MPL_otus_av_temp, MPL_otus_av_sali, MPL_otus_av_chla, MPL_otus_av_VA, MPL_otus_av_BA, MPL_otus_av_PO4, MPL_otus_av_NO, MPL_otus_av_SI)
#   
#   all_env_ordered_by_tree <- all_env_ordered_by_tree[, c("average_temp",
#                                                          "average_sal",
#                                                          "average_chla",
#                                                          "average_VA",
#                                                          "average_BA",
#                                                          "average_PO4",
#                                                          "average_NO",
#                                                          "average_SI")]
#   all_env_scaled_by_tree <- scale(all_env_ordered_by_tree)
#   
#   env_otu_melt <- melt(all_env_scaled_by_tree)
#   
#   ### Heatmap with log scaled values ####
#   
#   env_otu_melt$Var1 <- factor(env_otu_melt$Var1, levels = factor(amplicon_tree$tip.label))
#   
#   #OTU_table_long$Var2 <- factor(OTU_table_long$Var2, levels = factor(sort(OTU_table_long$Var2)))
#   
#    OTU_heatmap <- ggplot(env_otu_melt, 
#                        aes(Var2,Var1)) + 
#     geom_tile(aes(fill = value),
#             colour = "white") +
#    scale_fill_gradient(low = "red",high = "green",na.value = "white") +
#    theme_bw()
#   
#   ##fix up formatting
#   base_size=9
#   
#   OTU_heatmap_formatted <- OTU_heatmap + 
#    theme_bw(base_size = base_size) + 
#    labs(x = "",  y = "") + 
#    theme(axis.ticks = element_blank(),
#          axis.text.y = element_text(size = base_size / 1.2),
#          axis.text.x = element_text(size = base_size * 1.2,
#                                     angle = 330, 
#                                     hjust = 0,
#                                     colour = "grey20"))
#   
#   pdf(paste0("../figures/",amplicon_name,"_tree_with_env_mean_heatmap%03d.pdf"), width = 15, height = 15, onefile = FALSE )
#   print(OTU_heatmap_formatted)
#   dev.off()
#   
#   env_otu_melt_not_scaled <- melt(as.matrix(all_env_ordered_by_tree))
#   
#   environmental_plot_average_value <- function (long_data_set, 
#                                   parameter,
#                                   ind_plot_title) {
#    env_bar_plot <- ggplot(droplevels(subset(long_data_set,
#                                                 Var2 == parameter)),
#                               aes(x=Var1, y=value)) +
#     geom_bar(stat="identity")+
#     coord_flip()+
#     ggtitle(ind_plot_title)+
#     ylab(NULL)
#   
#     return(env_bar_plot)
#   }
  
#   av_temp_by_otu <- environmental_plot_average_value(env_otu_melt_not_scaled, 
#                                     "average_temp",
#                                     "Average temperature") +xlab(NULL)
#   
#   av_sal_by_otu <- environmental_plot_average_value(env_otu_melt_not_scaled, 
#                                                      "average_sal",
#                                                      "Average salinity")+ 
#    theme(axis.ticks.y = element_blank(), 
#          axis.text.y = element_blank(),
#          axis.title.y =element_blank()) +
#    ylab(NULL) 
#   
#   av_chla_by_otu <- environmental_plot_average_value(env_otu_melt_not_scaled, 
#                                                      "average_chla",
#                                                      "Average chlorophyll a")+
#    theme(axis.ticks.y = element_blank(), 
#          axis.text.y = element_blank(),
#          axis.title.y =element_blank()) +
#    ylab(NULL) 
#   
#   av_VA_by_otu <- environmental_plot_average_value(env_otu_melt_not_scaled, 
#                                                      "average_VA",
#                                                      "Average viral abundance")
#   
#   av_BA_by_otu <- environmental_plot_average_value(env_otu_melt_not_scaled, 
#                                                    "average_BA",
#                                                    "Average bacterial abundance")+
#    theme(axis.ticks.y = element_blank(),
#          axis.text.y = element_blank(),
#          axis.title.y =element_blank()) +
#    ylab(NULL)
#   
#   av_PO_by_otu <- environmental_plot_average_value(env_otu_melt_not_scaled, 
#                                                      "average_PO4",
#                                                      "Average phosphate")
#   
#   av_NO_by_otu <- environmental_plot_average_value(env_otu_melt_not_scaled, 
#                                                      "average_NO",
#                                                      "Average nitrate+nitrite")+ 
#    theme(axis.ticks.y = element_blank(),
#          axis.text.y = element_blank(),
#          axis.title.y =element_blank()) +
#    ylab(NULL) 
#   
#   av_SI_by_otu <- environmental_plot_average_value(env_otu_melt_not_scaled, 
#                                                    "average_SI",
#                                                    "Average silicate")+ theme(axis.ticks.y = element_blank(), 
#                                                                               axis.text.y = element_blank(),
#                                                                               axis.title.y =element_blank()) +
#    ylab(NULL) 
#   
#   pdf(paste0("../figures/",amplicon_name,"_tree_with_env_data_panel_chart_mean%03d.pdf"), width = 15, height = 15, onefile = FALSE )
#   grid.arrange(av_temp_by_otu,av_sal_by_otu,av_chla_by_otu,nrow = 1)
#   grid.arrange(av_VA_by_otu,av_BA_by_otu,nrow = 1)
#   grid.arrange(av_PO_by_otu,av_NO_by_otu,av_SI_by_otu,nrow = 1)
#   dev.off()
  
 ### try with range of values not just mean ####
  
  convert_otu_table_to_env_min_max <- function (OTU_table, env_param, env_col_min, env_col_max,amplicon_tree ) {
   ## add in the isolates
   OTU_table_with_isolates <- OTU_table
   sample_length <- length(rownames(OTU_table_with_isolates))
   
   #Make vector of 0s based on number of sites. 
   zeros <- rep(0, sample_length)
   
   #Insert vector in table based on number of tips in tree that are not found in the OTU table. 
   OTU_table_with_isolates  <- lapply(amplicon_tree$tip.label,function(i){
    if (is.element(i, colnames(OTU_table_with_isolates))) {
     zero_df <- c()
     #break
    }
    else {
     print("no")
     print(i)
     zero_df <- data.frame(zeros) 
     colnames(zero_df) <- paste(i)
     return(zero_df) 
    }
   })
   ## remove the empty from the lists
   OTU_table_with_isolates_non_null <- OTU_table_with_isolates[!sapply(OTU_table_with_isolates,is.null)]
   OTU_table_with_isolates1 <- do.call(cbind,OTU_table_with_isolates_non_null)
   OTU_table_with_OTUs_and_isolates <- cbind(OTU_table, OTU_table_with_isolates1)
   
   OTU_table_with_isolates_with_env <-   as.data.frame(apply(OTU_table_with_OTUs_and_isolates,
                                                             c(1,2),
                                                             function(x){
                                                              ## replace each value above 0 with the
                                                              ## value of the environmental param
                                                              replace(x,
                                                                      x > 0,
                                                                      env_param[match(rownames(OTU_table_with_OTUs_and_isolates)[x], Jericho_data$VC)])
                                                             }) )
  
  ## to get means to work on available samples
  OTU_table_with_isolates_with_env[OTU_table_with_isolates_with_env == 0] <- NA
  melt_env_otu_table <- melt(as.matrix(OTU_table_with_isolates_with_env))
  min_max_env_otu_table <- ddply(melt_env_otu_table,
                                  .(Var2),
                                  summarize,
                                  max=max(value, na.rm=TRUE),
                                  min=min(value, na.rm=TRUE))
  #### order matrix so that it is ordered like the tree ####
  row.names(min_max_env_otu_table) <- min_max_env_otu_table$Var2
  
  min_max_env_otu_table[min_max_env_otu_table == Inf] <- NA
  min_max_env_otu_table[min_max_env_otu_table == -Inf] <- NA
  
  min_max_env_otu_table_reordered <- min_max_env_otu_table[rev(amplicon_tree$tip.label),]
  names(min_max_env_otu_table_reordered) <- c("otu_name", env_col_max,env_col_min)
  return(min_max_env_otu_table_reordered)
  }
  
  
  MPL_otus_min_max_temp <- convert_otu_table_to_env_min_max(normalized_MPL_OTUs, 
                                                            Jericho_data$Temperature_YSI,
                                                            "min_temp",
                                                            "max_temp",
                                                            amplicon_tree)
  MPL_otus_min_max_sali <- convert_otu_table_to_env_min_max(normalized_MPL_OTUs,
                                                            Jericho_data$Salinity_ppt_YSI,
                                                            "min_sal",
                                                            "max_sal",
                                                            amplicon_tree)
  MPL_otus_min_max_chla <- convert_otu_table_to_env_min_max(normalized_MPL_OTUs,
                                                            Jericho_data$Average_chl_a, 
                                                            "min_chla",
                                                            "max_chla",
                                                            amplicon_tree )
  MPL_otus_min_max_VA <- convert_otu_table_to_env_min_max(normalized_MPL_OTUs,
                                                          Jericho_data$Average_viral_abundance,
                                                          "min_VA",
                                                          "max_VA",
                                                          amplicon_tree )
  MPL_otus_min_max_BA <- convert_otu_table_to_env_min_max(normalized_MPL_OTUs, Jericho_data$Average_bacterial_abundance,"min_BA", "max_BA",amplicon_tree)
  MPL_otus_min_max_PO4 <- convert_otu_table_to_env_min_max(normalized_MPL_OTUs, Jericho_data$Average_PO4, "min_PO4", "max_PO4", amplicon_tree)
  MPL_otus_min_max_NO <- convert_otu_table_to_env_min_max(normalized_MPL_OTUs, Jericho_data$Average_NO3_NO2, "min_NO", "max_NO",amplicon_tree)
  MPL_otus_min_max_SI <- convert_otu_table_to_env_min_max(normalized_MPL_OTUs, Jericho_data$Average_SiO2, "min_SI","max_SI", amplicon_tree)
  
  ### make into function and do with many others
  
  plot_crossbar_env <- function (MPL_otus_min_max_temp, min_temp, max_temp, bar_colour, ylimits) {
   MPL_otus_min_max_temp$otu_name <- factor(MPL_otus_min_max_temp$otu_name, levels = factor(amplicon_tree$tip.label))
    ggplot(MPL_otus_min_max_temp) +
     geom_crossbar(aes_string(ymin = min_temp, ymax = max_temp, x = "otu_name", y = min_temp),
                   fill = bar_colour, fatten = 5,
                   colour=bar_colour)+
     scale_y_continuous(limits = ylimits,
                        expand = c(0,0))+
     #theme_bw()+
     coord_flip()
  }
  
  pdf(paste0("../figures/",amplicon_name,"_tree_with_env_data_panel_chart_min_max%03d.pdf"), width = 15, height = 15, onefile = FALSE )
  print(plot_crossbar_env(MPL_otus_min_max_temp, "min_temp", "max_temp", "blue", c(5,25)))
  print(plot_crossbar_env(MPL_otus_min_max_sali, "min_sal", "max_sal", "green", c(6,32)))
  print(plot_crossbar_env(MPL_otus_min_max_chla, "min_chla", "max_chla", "green", c(0,10)))
  print(plot_crossbar_env(MPL_otus_min_max_VA, "min_VA", "max_VA", "yellow", c(1000000,50000000)))
  print(plot_crossbar_env(MPL_otus_min_max_BA, "min_BA", "max_BA", "purple", c(0,5000000)))
  print(plot_crossbar_env(MPL_otus_min_max_PO4, "min_PO4", "max_PO4", "lavender", c(0,3)))
  print(plot_crossbar_env(MPL_otus_min_max_NO, "min_NO", "max_NO", "pink", c(0,40)))
  print(plot_crossbar_env(MPL_otus_min_max_SI, "min_SI", "max_SI", "pink", c(0,80)))
  dev.off()
}

plot_tree_root_and_ladderize <- function (tree_file, root, plot_tree_file_name) {
 tree_amplicon <- read.tree(tree_file)
 tree_amplicon <- root(tree_amplicon, root, resolve.root=TRUE)
 tree_amplicon  <- ladderize(tree_amplicon)
 #exporting tree.
 write.tree(tree_amplicon, "../results/temp_tree.tree")
 #re-import to be able to match up the tip labels with the OTU table. 
 reordered_tree <- read.tree("../results/temp_tree.tree")
 tree_amplicon <- reordered_tree
 return(tree_amplicon)
}



### MPL ##### 

RdRp_95_ref_tree <- plot_tree_root_and_ladderize("../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.tree",
                                                 root = "Equine_rhinitis_B_virus_2_Picornaviridae_gi_15192762",
                                                 "../figures/RdRp_95_miseq_data_with_env_iso_and_ref_Fasttree.pdf")

make_env_plots_based_on_tree(RdRp_95_ref_tree,
                             normalized_MPL_OTUs,
                             "MPL_all_refs") 

###Top 100 with references ####

normalized_MPL_OTUs_top_100 <- read.delim("../data/OTU_table_Jericho_time_series_MPL_normalized_top_100.tsv",row.names=1)
## want only OTUs in normalized_mpls 100 and the refs..
## so want to get the names of OTUs not in the 100
normalized_MPL_to_remove <- normalized_MPL_OTUs[,!(colnames(normalized_MPL_OTUs) %in% colnames(normalized_MPL_OTUs_top_100))]

RdRp_95_ref_tree_top_100 <- drop.tip(RdRp_95_ref_tree, colnames(normalized_MPL_to_remove))

make_env_plots_based_on_tree(RdRp_95_ref_tree_top_100,
                             normalized_MPL_OTUs,
                             "MPL_refs_top_100") 
## Only MPL misequ ###

## so want to get the names of OTUs not in the 100
normalized_MPL_to_remove <- normalized_MPL_OTUs[,!(colnames(normalized_MPL_OTUs) %in% colnames(normalized_MPL_OTUs_top_100))]

reference_tips <- RdRp_95_ref_tree$tip.label[!(RdRp_95_ref_tree$tip.label %in% colnames(normalized_MPL_OTUs))]

RdRp_95_only_miseq <- drop.tip(RdRp_95_ref_tree, reference_tips)

# make_env_plots_based_on_tree(RdRp_95_only_miseq,
#                              normalized_MPL_OTUs,
#                              "MPL_only_miseq") 

### gp23  ####

#### With references ####

gp23_95_ref_tree <- read.tree("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.tree")
gp23_95_ref_tree <- root( gp23_95_ref_tree, "Enterobacteria_phage_T4", resolve.root=TRUE)

## having problems with memory for this tree...

gp23_95_ref_tree <- plot_tree_root_and_ladderize("../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimalled_no_colons.tree",
                                                   root = "Enterobacteria_phage_T4",
                                                   "../figures/gp23_95_miseq_data_with_env_iso_and_ref_Fasttree.pdf")
 
make_env_plots_based_on_tree(gp23_95_ref_tree,
                             normalized_gp23_OTUs,
                             "gp23_all_refs") 

###Top 100 with references ####

normalized_gp23_OTUs_top_100 <- read.delim("../data/OTU_table_Jericho_time_series_gp23_normalized_top_100.tsv",row.names=1)

## want only OTUs in normalized_gp23s 100 and the refs..
## so want to get the names of OTUs not in the 100
normalized_gp23_to_remove <- normalized_gp23_OTUs[,!(colnames(normalized_gp23_OTUs) %in% colnames(normalized_gp23_OTUs_top_100))]

gp23_95_ref_tree_top_100 <- drop.tip(gp23_95_ref_tree, colnames(normalized_gp23_to_remove))

make_env_plots_based_on_tree(gp23_95_ref_tree_top_100,
                             normalized_gp23_OTUs,
                             "gp23_refs_top_100") 

## Only gp23 misequ ###

## want only OTUs in normalized_gp23s 100 and the refs..
## so want to get the names of OTUs not in the 100
# normalized_gp23_to_remove <- normalized_gp23_OTUs[,!(colnames(normalized_gp23_OTUs) %in% colnames(normalized_gp23_OTUs_top_100))]
# 
# reference_tips <- gp23_95_ref_tree$tip.label[!(gp23_95_ref_tree$tip.label %in% colnames(normalized_gp23_OTUs))]
# 
# gp23_95_only_miseq <- drop.tip(gp23_95_ref_tree, reference_tips)

# make_env_plots_based_on_tree(gp23_95_only_miseq,
#                              normalized_gp23_OTUs,
#                              "gp23_only_miseq") 
##################################
# ### AVS #####
# 
# #### With references ####
# AVS_90_ref_tree <- plot_tree_root_and_ladderize("../results/AVS_90_miseq_data_with_env_iso_and_ref_align_clustal_edited.tree",
#                                                 root = "Ectocarpus_siliculosus_virus_1_gi_13177282",
#                                                 "../figures/AVS_90_miseq_data_with_env_iso_and_ref_Fasttree.pdf")
# 
# make_env_plots_based_on_tree(AVS_90_ref_tree,
#                              normalized_AVS_OTUs,
#                              "AVS_all_refs") 
# 
# 
# ###Top 100 with references ####
# 
# normalized_AVS_OTUs_top_100 <- read.delim("../data/OTU_table_Jericho_time_series_AVS_normalized_top_100.tsv",row.names=1)
# 
# ## want only OTUs in normalized_AVSs 100 and the refs..
# ## so want to get the names of OTUs not in the 100
# normalized_AVS_to_remove <- normalized_AVS_OTUs[,!(colnames(normalized_AVS_OTUs) %in% colnames(normalized_AVS_OTUs_top_100))]
# 
# AVS_90_ref_tree_top_100 <- drop.tip(AVS_90_ref_tree, colnames(normalized_AVS_to_remove))
# 
# make_env_plots_based_on_tree(AVS_90_ref_tree_top_100,
#                              normalized_AVS_OTUs,
#                              "AVS_refs_top_100")
# 
# ## Only AVS misequ ###
# ## want only OTUs in normalized_AVSs 100 and the refs..
# ## so want to get the names of OTUs not in the 100
# normalized_AVS_to_remove <- normalized_AVS_OTUs[,!(colnames(normalized_AVS_OTUs) %in% colnames(normalized_AVS_OTUs_top_100))]
# 
# reference_tips <- AVS_90_ref_tree$tip.label[!(AVS_90_ref_tree$tip.label %in% colnames(normalized_AVS_OTUs))]
# 
# AVS_90_only_miseq <- drop.tip(AVS_90_ref_tree, reference_tips)
# 
# make_env_plots_based_on_tree(AVS_90_only_miseq,
#                              normalized_AVS_OTUs,
#                              "AVS_only_miseq")
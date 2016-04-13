### What about Divesity over time? And OTUs shared between times??


## Like in Needham 2013, want to look at the community similarity matrices to see what it looks like over time. 

## make distance matrices over time and perform regression

library(vegan)
library(reshape2)
library(ggplot2)
library(scales)
library(gridExtra)
library(cowplot)

args <- commandArgs(TRUE)
inputFile <- args[1]

## test to see if input file is given, so I can decide whether to use this argument or the orginal one
if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
 path_colour <- "black"
 line_colour <- "black"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
 }
}

normalized_18s_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_MPL.tsv", row.names="VC_number")
#normalized_AVS_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_bloom_time_series_normalized_16S_R1.tsv", row.names="VC_number") 


## ok maybe easier to change to dates before hand
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

Env_data_for_merging <- Jericho_data [,c("Date","VC_number")]
row.names(Env_data_for_merging) <- Env_data_for_merging$VC_number

## will just merge OTUs with date. GOod for plotting then. 
get_community_similarity_plot_over_time <- function (normalized_otus) {
 new_table <- cbind(normalized_otus, Env_data_for_merging[, "Date"][match(rownames(normalized_otus), rownames(Env_data_for_merging))])
 row.names(new_table) <- new_table[,(dim(new_table)[2])]
 new_table <- new_table[,-(dim(new_table)[2])]
 
 
 diversity_new <- diversity(new_table)
 diversity_table <- as.data.frame(diversity_new)
 diversity_table$group <- "all"
  
 ggplot(diversity_table, aes(x=row.names(diversity_table),y=diversity_new)) + geom_point()
 ## plots diversity over time. 
 ggplot(diversity_table, aes(x=row.names(diversity_table),y=diversity_new, group=group)) + geom_line()
 
 
 ## what about community dissimilarity between times?
 ## want to look at the diss between following times. Time A to Time B
 
 OTU_diss <- vegdist(new_table)
 OTU_matrix <- as.matrix(OTU_diss)
 
 sorted_by_date_OTU_matrix <- OTU_matrix[sort(row.names(OTU_matrix)),sort(colnames(OTU_matrix))]
 
 community_over_time <- c("Date", "Community_diss_between_dates")
 
 ## or just a loop through all the dates using the indexing to get at them
 for (sample_date in row.names(sorted_by_date_OTU_matrix)){
  print(sample_date)
  index_time_a <- (which(row.names(sorted_by_date_OTU_matrix)==sample_date))
  ## use the next row in the sorted matrix!
  index_time_b <- index_time_a - 1
  print(index_time_b)
  ## at time 2 how does it compare to the last time
  print(sorted_by_date_OTU_matrix[index_time_a, index_time_b])
  comm_diss_between_dates <- sorted_by_date_OTU_matrix[index_time_a, index_time_b]
  row_add <- c(sample_date, comm_diss_between_dates)
  community_over_time <- rbind(community_over_time, row_add)
 }
 
 ## don't want the first few because the first is an artifact and the other is the header names
 colnames(community_over_time) <- community_over_time[1,]
 
 community_over_time <- community_over_time[-c(1:2),1:2]
 community_over_time <- as.data.frame(community_over_time)
 str(community_over_time)
 community_over_time$Date <- as.Date(community_over_time$Date)
 community_over_time$Community_diss_between_dates  <- as.numeric(as.character(community_over_time$Community_diss_between_dates))
 community_over_time$Community_sim_between_dates <- 1-community_over_time$Community_diss_between_dates
 community_over_time$group <- "all"
 str(community_over_time)

 p_all <- ggplot(community_over_time, aes(x=Date, y=Community_sim_between_dates, group=group))+ geom_line(colour=line_colour)+ geom_point(colour="grey")+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16),
        axis.line=element_line(colour = line_colour, size = 0.3),
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_line(colour="grey20"))+
  ggtitle(NULL)+ scale_x_date(limits = c(as.Date("2011-06-20"),as.Date("2011-07-06")))+
  ylim(0,1.0)
 
 return(p_all)
 
}


community_sim_over_time_18s <- get_community_similarity_plot_over_time(normalized_18s_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)
community_sim_over_time_18s
community_sim_over_time_16s <- get_community_similarity_plot_over_time(normalized_16s_OTUs)
community_sim_over_time_16s
# community_sim_over_time_AVS <- get_community_similarity_plot_over_time(normalized_AVS_OTUs)+
#  theme(axis.ticks.x = element_blank(), 
#        axis.text.x = element_blank(), 
#        axis.title.x=element_blank())+
#  xlab(NULL)
# community_sim_over_time_AVS
community_sim_over_time_MPL <- get_community_similarity_plot_over_time(normalized_MPL_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)

community_sim_over_time_MPL
community_sim_over_time_gp23 <- get_community_similarity_plot_over_time(normalized_gp23_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)
community_sim_over_time_gp23



pdf("../figures/time_series_community_similarity_all_amplicons_line_graphs%03d.pdf", width = 20, height = 15, onefile = FALSE)
plot_grid(
 #community_sim_over_time_MPL,
            # community_sim_over_time_AVS,
             community_sim_over_time_gp23,
             community_sim_over_time_16s,
             community_sim_over_time_18s, 
             ncol=1, labels = c("A", "B", "C"))
dev.off()

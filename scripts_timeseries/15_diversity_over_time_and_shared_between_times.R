
##
## make distance matrices over time and perform regression
library(vegan)
library(reshape2)
library(ggplot2)
library(scales)
library(gridExtra)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv", row.names=1)

args <- commandArgs(TRUE)
inputFile <- args[1]

## test to see if input file is given, so I can decide whether to use this argument or the orginal one. 
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

 p_all <- ggplot(community_over_time, aes(x=Date, y=Community_sim_between_dates, group=group))+ 
  season_line +
  spring_bloom_line+
  geom_line(colour=line_colour)+
  geom_point(colour="grey")+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16),
        axis.line=element_line(size = 0.3))+
  ggtitle(NULL)+ 
  date_scaling+
  ylim(0,1)
 
 return(p_all)
 
}


get_community_similarity_plot_over_time_lagged <- function (normalized_otus) {
 new_table <- cbind(normalized_otus, Env_data_for_merging[, "Date"][match(rownames(normalized_otus), rownames(Env_data_for_merging))])
 row.names(new_table) <- new_table[,(dim(new_table)[2])]
 new_table <- new_table[,-(dim(new_table)[2])]
 
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
 
 p_all <- ggplot(community_over_time, aes(x=as.Date(Date)-14, y=Community_sim_between_dates, group=group))+ 
  season_line +
  spring_bloom_line+
  geom_line(colour=line_colour)+
  geom_point(colour="grey")+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16),
        axis.line=element_line(size = 0.3))+
  ggtitle(NULL)+ 
  date_scaling+
  ylim(0,1)
 
 return(p_all)
 
}

community_sim_over_time_18s <- get_community_similarity_plot_over_time(normalized_18s_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)
community_sim_over_time_18s

## phytos
community_sim_over_time_18s_phytos <- get_community_similarity_plot_over_time(normalized_18s_OTUs_phytos)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)
community_sim_over_time_18s_phytos

community_sim_over_time_18s_hetero <- get_community_similarity_plot_over_time(normalized_18s_OTUs_hetero)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)
community_sim_over_time_18s_hetero

community_sim_over_time_16s <- get_community_similarity_plot_over_time(normalized_16s_OTUs)
community_sim_over_time_16s
community_sim_over_time_AVS <- get_community_similarity_plot_over_time(normalized_AVS_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)
community_sim_over_time_AVS

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
grid.arrange(community_sim_over_time_MPL,
             community_sim_over_time_AVS,
             community_sim_over_time_gp23,
             community_sim_over_time_18s, 
             community_sim_over_time_16s,
             ncol=1)
dev.off()

pdf("../figures/time_series_community_similarity_all_amplicons_no_AVS_line_graphs%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(community_sim_over_time_MPL,
             community_sim_over_time_18s,
             #community_sim_over_time_AVS,
             community_sim_over_time_16s,
             community_sim_over_time_gp23,
             ncol=1)
dev.off()

pdf("../figures/time_series_community_similarity_all_amplicons_no_AVS_with_phytos_and_het_line_graphs%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(community_sim_over_time_MPL,
             community_sim_over_time_18s,
             community_sim_over_time_18s_phytos,
             community_sim_over_time_18s_hetero,
             #community_sim_over_time_AVS,
             community_sim_over_time_16s,
             community_sim_over_time_gp23,
             ncol=1)
dev.off()

### plot with the viral sim moved back 2 weeks. ?

community_sim_over_time_gp23_lagged <- get_community_similarity_plot_over_time_lagged(normalized_gp23_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)
community_sim_over_time_gp23_lagged

community_sim_over_time_MPL_lagged <- get_community_similarity_plot_over_time_lagged(normalized_MPL_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)

pdf("../figures/time_series_community_similarity_lagged_virus_all_amplicons_no_AVS_line_graphs%03d.pdf", width = 20, height = 15, onefile = FALSE)
grid.arrange(community_sim_over_time_MPL_lagged,
             community_sim_over_time_18s,
             community_sim_over_time_18s_phytos,
             community_sim_over_time_18s_hetero,
             #community_sim_over_time_AVS,
             community_sim_over_time_16s,
             community_sim_over_time_gp23_lagged,
             ncol=1)
dev.off()

get_diversity_plot_over_time <- function(normalized_otus) {
 #normalized_otus <- normalized_gp23_OTUs
 ## adds in date
 new_table <- cbind(normalized_otus, Env_data_for_merging[, "Date"][match(rownames(normalized_otus), rownames(Env_data_for_merging))])
 ## rownames table by date
 row.names(new_table) <- new_table[,(dim(new_table)[2])]
 new_table <- new_table[,-(dim(new_table)[2])] ## remove extra date column
 
 diversity_new <- diversity(new_table)
 diversity_table <- as.data.frame(diversity_new)
 diversity_table$group <- "all"
 diversity_table$date <- as.Date(row.names(diversity_table))
 str(diversity_table)
 
#  diversity_points <- ggplot(diversity_table, aes(x=date,y=diversity_new))+
#   geom_point()+
#   theme_JAG_presentation()+
#   theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
#  print(diversity_points)
 ## plots diversity over time. 
 diversity_line <- ggplot(diversity_table, aes(x=date,y=diversity_new, group=group))+
  geom_line(colour=line_colour)+ 
  scale_x_date(breaks = date_breaks("month"), 
               #labels = date_format("%b"), 
               limits = c(as.Date("2010-06-01"), as.Date("2011-07-30")))+
  theme_JAG_presentation()+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
 print(diversity_line)
 return(diversity_line)
}

pdf("../figures/time_series_diversity_all_amplicons_line_graphs%03d.pdf", width = 20, height = 15, onefile = FALSE)

get_diversity_plot_over_time(normalized_MPL_OTUs)
get_diversity_plot_over_time(normalized_AVS_OTUs)
get_diversity_plot_over_time(normalized_gp23_OTUs)
get_diversity_plot_over_time(normalized_18s_OTUs)
get_diversity_plot_over_time(normalized_18s_OTUs_phytos)
get_diversity_plot_over_time(normalized_18s_OTUs_hetero)
get_diversity_plot_over_time(normalized_16s_OTUs)

dev.off()

MPL_div <- get_diversity_plot_over_time(normalized_MPL_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)

AVS_div <- get_diversity_plot_over_time(normalized_AVS_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)
gp23_div <- get_diversity_plot_over_time(normalized_gp23_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)
S18_div <- get_diversity_plot_over_time(normalized_18s_OTUs)+
 theme(axis.ticks.x = element_blank(), 
       axis.text.x = element_blank(), 
       axis.title.x=element_blank())+
 xlab(NULL)
S16_div <- get_diversity_plot_over_time(normalized_16s_OTUs)

pdf("../figures/time_series_diversity_all_amplicons_line_graphs_one_page.pdf", width = 20, height = 15, onefile = TRUE)
grid.arrange(MPL_div,
AVS_div,
gp23_div,
S18_div,
S16_div,
ncol=1)
dev.off()

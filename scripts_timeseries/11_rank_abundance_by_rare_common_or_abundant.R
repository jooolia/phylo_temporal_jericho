library(ggplot2)
library(reshape2)
library(BiodiversityR)
library(scales)
library(gridExtra)
library(plyr)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 


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


theme_overall <- theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
## will start by using this MPL

generate_frequency_by_otu_plot <- function (normalized_OTUs, name) {
  proportional_by_row <- prop.table(as.matrix(normalized_OTUs), margin=1)
  
  melted_prop_table <- melt(proportional_by_row)
  
  melted_prop_table <- droplevels(subset(melted_prop_table, value > 0))
  length_data_frame <- dim(melted_prop_table)[1]
 
  frequency_of_otu <- c()
  
  for (n in melted_prop_table$value){
   if (n >= 0.1){
    frequency_of_otu <- c(frequency_of_otu, "abundant")
   }
   else if (n >= 0.01 & n < 0.1){
    frequency_of_otu <- c(frequency_of_otu, "common")
   }
   else if (n < 0.01){
    frequency_of_otu <- c(frequency_of_otu, "rare")
   }
   else {
    print("oh no")
    print(n)
   }
  }
  
  melted_prop_table  <- cbind(frequency_of_otu, melted_prop_table)
   ## make frequency table
  tb <- table(melted_prop_table$Var2)
  
  ## rearranage the factors by the counts of the factors. Biggest first and then decreasing. 
  melted_prop_table <- within(melted_prop_table, Var2 <- factor(melted_prop_table$Var2,
                                                  levels = names(tb[order(tb, decreasing = TRUE)])))
  ## overall distribution of the types of OTUs -how often are rare, common and abundant found
  p <- ggplot(melted_prop_table, aes(x=frequency_of_otu))+ 
   geom_bar()+
   theme_JAG_presentation()
  
  pdf( paste("../figures/",name, "_frequency_of_otus_overall%03d.pdf", sep=""),width = 11, height = 17,onefile = FALSE)
  print(p)
  dev.off()
  
  p <- ggplot(melted_prop_table, aes(x=Var2, fill=frequency_of_otu))+
   geom_bar(width=0.7)+ scale_fill_manual(values = c("abundant" = "red","common" = "blue","rare" = "grey80"))+
   theme_JAG_presentation()+
   theme(panel.grid.major.y = element_blank(),
         panel.grid.major.x = element_blank(),
         axis.text.x = element_blank(), axis.title.x =element_blank())+
   xlab(NULL)

  pdf( paste("../figures/",name, "_frequency_by_otus%03d.pdf", sep=""),width = 17, height = 11,onefile = FALSE)
  print(p)
  dev.off()
  
}
# 
# ### Rank abundance by date for each amplicon
# rank_abundance_by_sample <- function (normalized_OTUs) {
#  rank_abundance_all_dates <- c()
#  for (row in row.names(normalized_OTUs)){
#   sample_for_rank_abundance <- normalized_OTUs[row,]
#   rank_abundance_sample <- rankabundance(sample_for_rank_abundance)
#   rank_abundance_sample <- as.data.frame(rank_abundance_sample)
#   rank_abundance_sample$VC_number <- row
#   rank_abundance_all_dates <- rbind(rank_abundance_all_dates, rank_abundance_sample)
#  }
#  rank_abundance_all_dates <- rank_abundance_all_dates[rank_abundance_all_dates$abundance > 0,]
#  ### now ggplot and facet!
#  rank_abundance_by_date_plot <- ggplot(rank_abundance_all_dates, aes(x=rank, y=abundance)) +geom_point(colour=line_colour)+
#   scale_y_continuous(trans=log10_trans())+
#   facet_wrap(~VC_number)+theme_JAG_presentation()
#  return(rank_abundance_by_date_plot)
# }
# 
# 
# 
# MPL <- rankabundance(normalized_MPL_OTUs)
# MPL.frame <- data.frame(MPL)
# MPL.frame <- MPL.frame[MPL.frame$abundance > 0,]
# 
# MPL_plot <- ggplot(MPL.frame, aes(x=rank, y=abundance)) + geom_point(colour=line_colour) +scale_y_continuous(trans=log10_trans())+theme_JAG_presentation()+theme(panel.grid.major = element_blank())
# plot(radfit(MPL.frame$abundance), log="y")
# 
# MPL_rank_abundance_by_date <- rank_abundance_by_sample(normalized_MPL_OTUs)
# 
# 
# AVS <- rankabundance(normalized_AVS_OTUs)
# AVS.frame <- data.frame(AVS)
# AVS.frame <- AVS.frame[AVS.frame$abundance > 0,]
# 
# AVS_plot <- ggplot(AVS.frame, aes(x=rank, y=abundance)) + geom_point(colour=line_colour) +scale_y_continuous(trans=log10_trans())+theme_JAG_presentation()+theme(panel.grid.major = element_blank())
# 
# plot(radfit(AVS.frame$abundance), log="y")
# AVS_rank_abundance_by_date <- rank_abundance_by_sample(normalized_AVS_OTUs)
# 
# gp23 <- rankabundance(normalized_gp23_OTUs)
# gp23.frame <- data.frame(gp23)
# gp23.frame <- gp23.frame[gp23.frame$abundance > 0,]
# 
# gp23_plot <- ggplot(gp23.frame, aes(x=rank, y=abundance)) + geom_point(colour=line_colour) +scale_y_continuous(trans=log10_trans())+theme_JAG_presentation()+theme(panel.grid.major = element_blank())
# plot(radfit(gp23.frame$abundance), log="y")
# gp23_rank_abundance_by_date <- rank_abundance_by_sample(normalized_gp23_OTUs)
# 
# 
# S18 <- rankabundance(normalized_18s_OTUs)
# S18.frame <- data.frame(S18)
# S18.frame <- S18.frame[S18.frame$abundance > 0,]
# 
# S18_plot <- ggplot(S18.frame, aes(x=rank, y=abundance)) + geom_point(colour=line_colour) +scale_y_continuous(trans=log10_trans())+theme_JAG_presentation()+theme(panel.grid.major = element_blank())
# plot(radfit(S18.frame$abundance), log="y")
# S18_rank_abundance_by_date <- rank_abundance_by_sample(normalized_18s_OTUs)
# 
# 
# S16 <- rankabundance(normalized_16s_OTUs)
# S16.frame <- data.frame(S16)
# S16.frame <- S16.frame[S16.frame$abundance > 0,]
# 
# S16_plot <- ggplot(S16.frame, aes(x=rank, y=abundance)) + geom_point(colour=line_colour) +scale_y_continuous(trans=log10_trans())+theme_JAG_presentation()+theme(panel.grid.major = element_blank())
# 
# plot(radfit(S16.frame$abundance), log="y")
# S16_rank_abundance_by_date <- rank_abundance_by_sample(normalized_16s_OTUs)
# 
# pdf("../figures/rank_abundance_curves.pdf",width = 20, height = 15,onefile = TRUE)
# MPL_plot
# AVS_plot
# gp23_plot
# S18_plot
# S16_plot
# dev.off()
# 
# pdf("../figures/rank_abundance_curves_one_pager.pdf",width = 20, height = 30,onefile = TRUE)
# grid.arrange(
# MPL_plot,
# AVS_plot,
# gp23_plot,
# S18_plot,
# S16_plot,
# ncol=1)
# dev.off()
# 
# 
# pdf("../figures/rank_abundance_curves_by_sample.pdf",width = 20, height = 15,onefile = TRUE)
# MPL_rank_abundance_by_date
# AVS_rank_abundance_by_date
# gp23_rank_abundance_by_date
# S18_rank_abundance_by_date
# S16_rank_abundance_by_date
# dev.off()


generate_frequency_by_otu_plot(normalized_MPL_OTUs, "MPL")
generate_frequency_by_otu_plot(normalized_gp23_OTUs, "gp23")
generate_frequency_by_otu_plot(normalized_AVS_OTUs, "AVS")
generate_frequency_by_otu_plot(normalized_18s_OTUs, "18S")
generate_frequency_by_otu_plot(normalized_16s_OTUs, "16S")

## would be nice to divide up the 16s and 18s into the major taxonomic groups
taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)

taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)

plot_bacteria_or_eularyote_otu_frequency_by_taxon <- function (normalized_table, taxonomy_table, name) {
  
 ## final result is otu number with taxonomy - could just write out this table?
  clean_taxonomy <- taxonomy_table
  ## but what if I want to remove those that are unclassified??
  proportional_by_row <- prop.table(as.matrix(normalized_table), margin=1)
  
  melted_prop_table <- melt(proportional_by_row)
  melted_prop_table <- droplevels(subset(melted_prop_table, value > 0))
  frequency_of_otu <- c()
  
  for (n in melted_prop_table$value){
   if (n >= 0.1){
    frequency_of_otu <- c(frequency_of_otu, "abundant")
   }
   else if (n >= 0.01 & n < 0.1){
    frequency_of_otu <- c(frequency_of_otu, "common")
   }
   else if (n < 0.01){
    frequency_of_otu <- c(frequency_of_otu, "rare")
   }
   else {
    print("oh no")
    print(n)
   }
  }
  
  melted_prop_table  <- cbind(frequency_of_otu, melted_prop_table)
   ## matching
  melted_prop_table$Domain <- clean_taxonomy$Domain[match(melted_prop_table$Var2, clean_taxonomy$otu_number)]
  melted_prop_table$Phylum <- clean_taxonomy$Phylum[match(melted_prop_table$Var2, clean_taxonomy$otu_number)]
  melted_prop_table$Class <- clean_taxonomy$Class[match(melted_prop_table$Var2, clean_taxonomy$otu_number)]
   
  ## make frequency table
  tb <- table(melted_prop_table$Var2)
  
  ## rearranage the factors by the counts of the factors. Biggest first and then decreasing. 
  melted_prop_table <- within(melted_prop_table, Var2 <- factor(melted_prop_table$Var2,
                                                                levels = names(tb[order(tb, decreasing = TRUE)])))

  p <- ggplot(melted_prop_table, aes(x=Var2, fill=frequency_of_otu))+
   geom_bar(width=0.7)+
   facet_wrap(~Phylum, scales = "free_x")+
   scale_fill_manual(values = c("abundant" = "red","common" = "blue","rare" = "grey80"))+
   theme_JAG_presentation()+
   theme(panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank())
  
  print(p)
  pdf( paste("../figures/",name, "_frequency_of_otus_by_Phylum%03d.pdf", sep=""),width = 17, height = 17,onefile = FALSE)
  print(p)
  #dev.off()
    
  p <- ggplot(melted_prop_table, aes(x=Var2, fill=frequency_of_otu))+
   geom_bar(width=0.7)+
   facet_wrap(~Class, scales = "free_x")+
   scale_fill_manual(values = c("abundant" = "red","common" = "blue","rare" = "grey80"))+
   theme_JAG_presentation()+
   theme(panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank())

    print(p)
 
   pdf( paste("../figures/",name, "_frequency_of_otus_by_Class%03d.pdf", sep=""),width = 17, height = 17,onefile = FALSE)
  print(p)
  dev.off()

}

test_function <-plot_bacteria_or_eularyote_otu_frequency_by_taxon(normalized_18s_OTUs, taxonomy_18s, "18S")
test_function <-plot_bacteria_or_eularyote_otu_frequency_by_taxon(normalized_16s_OTUs, taxonomy_16s, "16S")


## remove the unclassified and Chloroplasts ####

plot_bacteria_or_eularyote_otu_frequency_by_taxon_no_unclassified <- function (normalized_table, taxonomy_table, name) {

 ## final result is otu number with taxonomy - could just write out this table?
 #normalized_table <-normalized_18s_OTUs 
 clean_taxonomy <- taxonomy_table
 #clean_taxonomy <- taxonomy_18s
 unclassified_otus <- subset(clean_taxonomy, Domain == "unknown" | Phylum =="unclassified" | Class =="Chloroplast")

 ## OTUs not in unclassified
 normalized_table <- droplevels(subset(normalized_table, select=!(names(normalized_table) %in% unclassified_otus$otu_number)))
 
 proportional_by_row <- prop.table(as.matrix(normalized_table), margin=1)
 
 melted_prop_table <- melt(proportional_by_row)
 
 melted_prop_table <- droplevels(subset(melted_prop_table, value > 0))
 
  frequency_of_otu <- c()
 
 for (n in melted_prop_table$value){
  if (n >= 0.1){
   frequency_of_otu <- c(frequency_of_otu, "abundant")
  }
  else if (n >= 0.01 & n < 0.1){
   frequency_of_otu <- c(frequency_of_otu, "common")
  }
  else if (n < 0.01){
   frequency_of_otu <- c(frequency_of_otu, "rare")
  }
  else {
   print("oh no")
   print(n)
  }
 }
  melted_prop_table  <- cbind(frequency_of_otu, melted_prop_table)
  ## matching
 melted_prop_table$Domain <- clean_taxonomy$Domain[match(melted_prop_table$Var2, clean_taxonomy$otu_number)]
 melted_prop_table$Phylum <- clean_taxonomy$Phylum[match(melted_prop_table$Var2, clean_taxonomy$otu_number)]
 melted_prop_table$Class <- clean_taxonomy$Class[match(melted_prop_table$Var2, clean_taxonomy$otu_number)]
 ## make frequency table
 tb <- table(melted_prop_table$Var2)
  ## rearranage the factors by the counts of the factors. Biggest first and then decreasing. 
 melted_prop_table <- within(melted_prop_table, Var2 <- factor(melted_prop_table$Var2,
                                                               levels = names(tb[order(tb, decreasing = TRUE)])))
  
 p <- ggplot(melted_prop_table, aes(x=Var2, fill=frequency_of_otu))+
  geom_bar(width=0.7)+
  facet_wrap(~Phylum,scales = "free_x")+
  scale_fill_manual(values = c("abundant" = "red","common" = "blue","rare" = "grey80"))+
  theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank())
 
 print(p)
 pdf( paste("../figures/",name, "_frequency_of_otus_by_Phylum_no_unclassified.pdf", sep=""),width = 17, height = 17,onefile = TRUE)
 print(p)
 dev.off()
 
 p <- ggplot(melted_prop_table, aes(x=Var2, fill=frequency_of_otu))+
  geom_bar(width=0.7)+
  facet_wrap(~Class,scales = "free_x")+
  scale_fill_manual(values = c("abundant" = "red","common" = "blue","rare" = "grey80"))+
  theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank())
 print(p)
  pdf( paste("../figures/",name, "_frequency_of_otus_by_Class_no_unclassified.pdf", sep=""),width = 17, height = 17,onefile = TRUE)
 print(p)
 dev.off()
  }

 
 S18_frequency_OTUs_no_unclassified <-plot_bacteria_or_eularyote_otu_frequency_by_taxon_no_unclassified(normalized_18s_OTUs, taxonomy_18s, "18S")
S16_frequency_OTUs_no_unclassified <-plot_bacteria_or_eularyote_otu_frequency_by_taxon_no_unclassified(normalized_16s_OTUs, taxonomy_16s, "16S")
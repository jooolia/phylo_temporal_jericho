## 16_5 Comm sim vs eveness

## look at betadiversity vs. eveness
library(vegan)
library(reshape2)
library(plyr)
library(dplyr)
library(Hmisc)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv",                                  row.names=1)

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

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
## remove the winter set for now so that I can look at the lags

Env_data_for_merging <- Jericho_data [,c("Date","VC_number")]
row.names(Env_data_for_merging) <- Env_data_for_merging$VC_number

## subset data that is interesting to me:
keep <- c("Date",
          "Average_viral_abundance",
          "Average_bacterial_abundance",
          "Average_chl_a",
          "Average_PO4",
          "Average_SiO2",
          "Average_NO3_NO2",
          "Temperature_YSI",
          "Salinity_ppt_YSI",
          "Dissolved_oxygen_percent",  
          "pH")
## might be useful to do some more cleaning before
## keep only certain rows. 
params_Jericho <- droplevels(Jericho_data[,keep])



lm_eqn = function(m) {
 
 l <- list(a = format(coef(m)[1], digits = 2),
           b = format(abs(coef(m)[2]), digits = 2),
           r2 = format(summary(m)$r.squared, digits = 3));
 
 if (coef(m)[2] >= 0)  {
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
 } else {
  eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
 }
 
 as.character(as.expression(eq));                 
}


betadiversity_vs_alpha_diversity <- function (normalized_OTUs_A,normalized_OTUs_B, Env_data, title) {
# normalized_OTUs_A <- normalized_16s_OTUs
# normalized_OTUs_B <- normalized_gp23_OTUs
# Env_data <- params_Jericho
 ## add date and make new table
 
 new_table <- cbind(normalized_OTUs_A, Jericho_data[, "Date"][match(rownames(normalized_OTUs_A), Jericho_data$VC_number)])
 ## rename rownames as date
 row.names(new_table) <- new_table[,(dim(new_table)[2])]
 ## get rid of VC number column
 new_table <- new_table[,-(dim(new_table)[2])]
 
 OTU_diss <- vegdist(new_table)
 SetA_diss <- melt(as.matrix(OTU_diss))
 
 OTU_matrix <- as.matrix(OTU_diss)
 
 sorted_by_date_OTU_matrix <- OTU_matrix[sort(row.names(OTU_matrix)),sort(colnames(OTU_matrix))]
 
 community_over_time <- c("Date", "Community_diss_between_dates")
 
 ## or just a loop through all the dates using the indexing to get at them
 for (sample_date in row.names(sorted_by_date_OTU_matrix)){
  #print(sample_date)
  index_time_a <- (which(row.names(sorted_by_date_OTU_matrix)==sample_date))
  ## use the next row in the sorted matrix!
  index_time_b <- index_time_a - 1
  #print(index_time_b)
  ## at time 2 how does it compare to the last time
  #print(sorted_by_date_OTU_matrix[index_time_a, index_time_b])
  comm_diss_between_dates <- sorted_by_date_OTU_matrix[index_time_a, index_time_b]
  row_add <- c(sample_date, comm_diss_between_dates)
  community_over_time <- rbind(community_over_time, row_add)
 }
 
 colnames(community_over_time) <- community_over_time[1,]
 community_over_time <- community_over_time[-c(1:2),1:2]
 community_over_time <- as.data.frame(community_over_time)
 community_over_time$Date <- as.Date(community_over_time$Date)
 community_over_time$Community_diss_between_dates  <- as.numeric(as.character(community_over_time$Community_diss_between_dates))
 community_over_time$Community_sim_between_dates <- 1-community_over_time$Community_diss_between_dates
 ## so now compare to the eveness!!!
 
 # normalized_OTUs_B_long_otus$Date <- Jericho_data$Date[match(normalized_OTUs_B_long_otus$VC_number, Jericho_data$VC_number)]
 rownames(normalized_OTUs_B) <- Jericho_data$Date[match(rownames(normalized_OTUs_B), Jericho_data$VC_number)]
# normalized_OTUs_B$VC_number = rownames(normalized_OTUs_B )
 #normalized_OTUs_B_long_otus <- melt(normalized_OTUs_B , 
  #                                   id="VC_number", 
   #                                  variable.name="OTUid")
 ## Add in dates:
 # normalized_OTUs_B_long_otus$Date <- Jericho_data$Date[match(normalized_OTUs_B_long_otus$VC_number, Jericho_data$VC_number)]
data(BCI)
 H <- diversity(normalized_OTUs_B)
 ## Species richness (S) and Pielou's evenness (J):
 S <- specnumber(normalized_OTUs_B) ## rowSums(BCI > 0) does the same...
 J <- H/log(S)

 
  
#  
#  
#  normalized_OTUs_B_summarise_eveness_otus <- normalized_OTUs_B_long_otus %>%
#   filter(value > 0) %>%
#   ddply(~Date,summarise,eveness=length(OTUid))
 
 intersect_dates <- as.Date(intersect(as.Date(rownames(normalized_OTUs_B)), as.Date(community_over_time$Date)),  origin = "1970-01-01")
 
 community_over_time_subset <- subset(community_over_time, Date %in% intersect_dates)
 normalized_OTUs_B_eveness_subset <- subset(J, as.Date(names(J)) %in% intersect_dates)
 
 cor(community_over_time_subset$Community_diss_between_dates, normalized_OTUs_B_eveness_subset)
 
 spearman <- rcorr(community_over_time_subset$Community_diss_between_dates, normalized_OTUs_B_eveness_subset,type="spearman") 
 correlations <- spearman$r
 print(correlations)
 p <- spearman$P
 print(p)
 
 ## subset to non-blitz times
 Jericho_data_no_blitz <- Jericho_data[-c(18:22),]
 
 community_over_time_no_blitz_subset <- subset(community_over_time_subset, Date %in% Jericho_data_no_blitz$Date)
 normalized_OTUs_B_eveness_no_blitz_subset <- subset(normalized_OTUs_B_eveness_subset, as.Date(names(normalized_OTUs_B_eveness_subset)) %in% Jericho_data_no_blitz$Date)
 
 # #lag(community_over_time_no_blitz_subset$Community_sim_between_dates)
 # spearman <- rcorr(as.data.frame(lag(community_over_time_no_blitz_subset$Community_sim_between_dates, k=1)[-1]),
 #                   as.data.frame(normalized_OTUs_B_eveness_no_blitz_subset$eveness[-length(normalized_OTUs_B_eveness_no_blitz_subset$eveness)]),
 #                   type="spearman") 
 # 
#  test_beta <- merge(community_over_time_no_blitz_subset, normalized_OTUs_B_eveness_no_blitz_subset, by=c("Date", as.Date(names(normalized_OTUs_B_eveness_no_blitz_subset))))
#  
#  not_lagged_plot <- ggplot(test_beta, aes(y=eveness, x=Community_diss_between_dates))+
#   geom_point(colour=point_colour)+
#   stat_smooth(method="lm", se=TRUE)+
#   ggtitle(paste0(title, " no correlation r=", round(correlations[2,1],digits=4)," p-value ", round(p[2,1], digits=4 )))+
#   theme_JAG_presentation()+
#   geom_text(aes(y = min(eveness)+0.1, x = max(Community_sim_between_dates)-0.3), label = lm_eqn(lm(eveness~Community_sim_between_dates, test_beta_lagged)), parse = TRUE, colour=line_colour)
#  
#  pdf(file=paste0("../figures/not_lagged_comm_sim_to_eveness", title, ".pdf"))
#  print(not_lagged_plot)
#  dev.off()
#  
 
 
#  test_beta_lagged <-  mutate(test_beta, Community_diss_between_dates=lead(Community_diss_between_dates))
#  test_beta_lagged <-  mutate(test_beta_lagged, Community_sim_between_dates=lead(Community_sim_between_dates) )
#  test_beta_lagged <-  na.omit(test_beta_lagged )
#  test_beta_lagged <- as.data.frame(test_beta_lagged)
#  
#  spearman <- rcorr(lead(test_beta_lagged$Community_diss_between_dates), test_beta_lagged$eveness,type="pearson")
#  correlations <- spearman$r
#  print(correlations)
#  p <- spearman$P
#  print(p)
#  
#  lagged_plot <- ggplot(test_beta_lagged , aes(y=eveness, x=Community_diss_between_dates))+geom_point()+stat_smooth(method="lm", se=TRUE)+ggtitle(paste0(title, " with lagged correlation r=", round(correlations[2,1],digits=4)," p-value ", round(p[2,1], digits=4 )))+ theme_JAG_presentation()
#  #+ geom_text(aes(y = min(eveness)+0.1, x = max(Community_sim_between_dates)-0.3, label = lm_eqn(lm(eveness~Community_sim_between_dates, test_beta_lagged))), parse = TRUE)
#  pdf(file=paste0("../figures/lagged_comm_sim_to_eveness", title, ".pdf"))
#  print(lagged_plot)
#  dev.off()
 
 pdf(file=paste0("../figures/ccf_plot", title, ".pdf"))
 ccf(community_over_time_no_blitz_subset$Community_diss_between_dates, normalized_OTUs_B_eveness_subset)
 dev.off()
 
 #merged_data_with_lag <- merge()
}

betadiversity_vs_alpha_diversity(normalized_18s_OTUs,normalized_MPL_OTUs, params_Jericho, "S18_comm_MPL_eveness")
betadiversity_vs_alpha_diversity(normalized_18s_OTUs_hetero,normalized_MPL_OTUs, params_Jericho, "S18_hetero_comm_MPL_eveness")
betadiversity_vs_alpha_diversity(normalized_18s_OTUs_phytos,normalized_MPL_OTUs, params_Jericho,"S18_phyto_comm_MPL_eveness")

betadiversity_vs_alpha_diversity(normalized_MPL_OTUs,normalized_18s_OTUs, params_Jericho, "MPL_comm_18s_eveness")

betadiversity_vs_alpha_diversity(normalized_16s_OTUs,normalized_gp23_OTUs, params_Jericho,"S16_comm_gp23_eveness")

betadiversity_vs_alpha_diversity(normalized_18s_OTUs,normalized_gp23_OTUs, params_Jericho,"S18_comm_gp23_eveness")

betadiversity_vs_alpha_diversity(normalized_16s_OTUs,normalized_18s_OTUs, params_Jericho, "S16_comm_18s_eveness")
betadiversity_vs_alpha_diversity(normalized_18s_OTUs,normalized_16s_OTUs, params_Jericho, "S18_comm_16s_eveness")



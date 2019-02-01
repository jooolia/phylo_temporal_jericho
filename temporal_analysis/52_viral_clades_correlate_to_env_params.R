

## virus groups correlate to env params?
library(Hmisc)
library(dplyr)
library(stringr)
library(reshape2)

MPL_group_sums <- read.csv("../results/MPL_group_sums_by_site_prop.csv")
gp23_group_sums <- read.csv("../results/gp23_group_sums_by_site_prop.csv")

S18_ord_sums <- read.csv("../results/normalized_18s_summarized_by_order.csv")
S18_fam_sums <- read.csv("../results/normalized_18s_summarized_by_family.csv")

names(S18_ord_sums)[1] <- "VC"
names(S18_fam_sums)[1] <- "VC"

S16_ord_sums <- read.csv("../results/normalized_16s_summarized_by_order.csv")

names(S16_ord_sums)[1] <- "VC"

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv",
                         row.names = 1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

parameters_to_exclude <- c("season",
                           # "Average_NO3_NO2",
                           "Standard_error_NO3_NO2",
                           "Secchi_disk_disappears",
                           "Secchi_disk_reappears",
                           "VC_number",
                           "Raw_Bacteria_Rep_A",
                           "Raw_Bacteria_Rep_B",
                           "Raw_Viruses_rep_A",
                           "Raw_Viruses_rep_B",
                           "PO4_rep_A",
                           "PO4_rep_B",
                           "SiO2_rep_A",
                           "SiO2_rep_B",
                           "NO3_NO2_rep_A",
                           "NO3_NO2_rep_B",
                           "Standard_error_viral_abundance",
                           "Standard_error_bacterial_abundance",
                           "Standard_error_PO4",
                           "Standard_error_SiO2",
                           "Standard_error_chl_a")
params_Jericho <-  subset(Jericho_data,
                          select= !(colnames(Jericho_data) %in% parameters_to_exclude))

summary(params_Jericho)
#params_Jericho <- params_Jericho[,-12]
## get those in the these samples 

## RdRp

correlation_and_print_table <- function (group_sums,
                                         table_print_name) {
  
  VC_dates <- Jericho_data$Date[match(group_sums$VC,
                                      Jericho_data$VC)]
  params_Jericho_RdRp <- subset(params_Jericho, Date %in% VC_dates)
  
  spearman <- rcorr(as.matrix(group_sums[,-(1:2)]),
                    as.matrix(params_Jericho_RdRp[,-1]),
                    type="spearman") 
  correlations <- spearman$r
  p <- spearman$P
  
  melted_cor <- melt(correlations)
  melted_p <- melt(p)
  
  mystars <- ifelse(p < .001, "***",
                    ifelse(p < .01, "** ",
                           ifelse(p < .05, "* ", " ")))
  melted_mystars <- melt(mystars)
  
  melted_together <- cbind( melted_cor,
                            p_value = melted_p$value,
                            stars = melted_mystars$value)
  melted_together <- na.omit(melted_together) ## gets rid of the leftover diagonals
  names(melted_together)
  
  melted_together$cor_with_stars <- paste(round(melted_together$value, digits=3),
                                          " ",
                                          melted_together$stars)
  
  #filter for correlation about 0.6 and p value less than 0.01
  #filtered_data <- filter(melted_together, p_value <= 0.05  & abs(value) > 0.5)
  
  ## makes data frame of significance
  #filter for correlation about 0.6 and p value less than 0.01
  #filtered_data <- filter(melted_together, p_value <= 0.01  & abs(value) > 0.5)
  
  filtered_data_groups <- melted_together[str_split_fixed(melted_together$Var1,
                                                          "_",
                                                          n=2)[,1]=="sum" |str_split_fixed(melted_together$Var1,
                                                                                           "_",
                                                                                           n=2)[,1]=="sum",]
  
  #significant_cor <- dcast(filtered_data_groups, Var1 ~Var2)
  significant_cor <- dcast(melted_together,
                           Var1 ~ Var2,
                           value.var = "cor_with_stars")
  
  write.csv(significant_cor,
            file = table_print_name)
  
}

correlation_and_print_table(MPL_group_sums,
                            "../results/RdRP_groups_cor.csv")

# gp23

correlation_and_print_table(gp23_group_sums,
                            "../results/gp23_groups_cor.csv")
### 18s

correlation_and_print_table(S18_ord_sums,
                            "../results/S18_orders_cor.csv")

### compare groups to S18 and env to MPL ####
correlation_between_amplicons_and_print_table <- function (group_sums_1,
                                                           group_sums_2,
                                                           table_print_name) {
  
  group_sums_1$Date <- Jericho_data$Date[match(group_sums_1$VC,
                                               Jericho_data$VC)]
  ## want to add the S18 matrix to the viral fam one. 
  group_sums_1_with_env <- merge(as.matrix(group_sums_1),
                                 as.matrix(params_Jericho),
                                 by="Date")
  group_sums_1_with_env$Date <- as.Date(group_sums_1_with_env$Date)
  
  union_vcs <- intersect(group_sums_2$VC,
                         group_sums_1_with_env$VC)
  
  group_sums_2_jer <- droplevels(subset(group_sums_2,
                                        group_sums_2$VC %in% union_vcs))
  
  group_sums_1_with_env_jer <- droplevels(subset(group_sums_1_with_env,
                                                 group_sums_1_with_env$VC%in% union_vcs))
  
  spearman <- rcorr(as.matrix(group_sums_2_jer[,-(1:2)]),
                    as.matrix(group_sums_1_with_env_jer[,-(1:2)]),
                    type = "spearman") 
  correlations <- spearman$r
  p <- spearman$P
  
  melted_cor <- melt(correlations)
  melted_p <- melt(p)
  mystars <- ifelse(p < .001, "***",
                    ifelse(p < .01, "** ",
                           ifelse(p < .05, "* ", " ")))
  melted_mystars <- melt(mystars)
  
  melted_together <- cbind( melted_cor,
                            p_value = melted_p$value,
                            stars = melted_mystars$value)
  melted_together <- na.omit(melted_together) ## gets rid of the leftover diagonals
  names(melted_together)
  
  melted_together$cor_with_stars <- paste(round(melted_together$value,
                                                digits=3),
                                          " ",
                                          melted_together$stars)
  
  #filter for correlation about 0.6 and p value less than 0.01
  #filtered_data <- filter(melted_together, p_value <= 0.01  & abs(value) > 0.5)
  
  filtered_data_groups <- melted_together[str_split_fixed(melted_together$Var1,
                                                          "_",
                                                          n=2)[,1]=="sum" |str_split_fixed(melted_together$Var1,
                                                                                           "_",
                                                                                           n=2)[,1]=="sum",]
  
  significant_cor <- dcast(filtered_data_groups,
                           Var1 ~ Var2,
                           value.var = "cor_with_stars")
  write.csv(significant_cor,
            file = table_print_name)
}

correlation_between_amplicons_and_print_table(S18_ord_sums,
                                              MPL_group_sums,
                                              "../results/RdRp_groups_with_18S_orders_cor.csv")

### compare groups to S18 and env to gp23 ####
correlation_between_amplicons_and_print_table(S18_ord_sums,
                                              gp23_group_sums,
                                              "../results/gp23_groups_with_18S_orders_cor.csv")

### 16S and gp23
correlation_between_amplicons_and_print_table(S16_ord_sums,
                                              gp23_group_sums,
                                              "../results/gp23_groups_with_16S_orders_cor.csv")

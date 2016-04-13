## PD correlated to env?
library(psych)


RdRp_PD <- read.csv("../results/RdRp_only_miseq_phylogenetic_diversity_table.csv")
gp23_PD <- read.csv("../results/gp23_only_miseq_phylogenetic_diversity_table.csv")
AVS_PD <- read.csv("../results/AVS_only_miseq_phylogenetic_diversity_table.csv")
S18_PD <- read.csv("../results/S18_phylogenetic_diversity_table.csv")
S18_PD_phyto <- read.csv("../results/S18_phytos_phylogenetic_diversity_table.csv")
S18_PD_hetero <- read.csv("../results/S18_hetero_phylogenetic_diversity_table.csv")
S16_PD <- read.csv("../results/S16_phylogenetic_diversity_table.csv")


Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
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
params_Jericho <-  subset(Jericho_data, select= !(colnames(Jericho_data) %in% parameters_to_exclude))

summary(params_Jericho)
#params_Jericho <- params_Jericho[,-12]
## get those in the these samples 


## RdRp
corelations_of_PD_to_env <- function (env_data, amplicon_PD) {
  VC_dates <- env_data$Date[match(as.Date(amplicon_PD$Date), as.Date(params_Jericho$Date))]
  params_Jericho_amplicon <- subset(params_Jericho, Date %in% VC_dates)
  
  cor(amplicon_PD$PD, params_Jericho_amplicon[,-1])
  cor(amplicon_PD$PD, params_Jericho_amplicon[,-1], use = "pairwise")
  cor(amplicon_PD$PD, params_Jericho_amplicon[,-1], use = "pairwise", method="spearman")
  cov(amplicon_PD$PD, params_Jericho_amplicon[,-1], use = "pairwise", method="pearson")
 print( corr.test(as.data.frame(amplicon_PD[,"PD"]), params_Jericho_amplicon[,-1], method="spearman", use = "pairwise")$r)
 print( corr.test(as.data.frame(amplicon_PD[,"PD"]), params_Jericho_amplicon[,-1], method="spearman",use = "pairwise")$p)
  
 print( corr.test(as.data.frame(amplicon_PD[,"SR"]), params_Jericho_amplicon[,-1], method="spearman", use = "pairwise")$r)
  corr.test(as.data.frame(amplicon_PD[,"SR"]), params_Jericho_amplicon[,-1], method="spearman",use = "pairwise")$p
  
}

corelations_of_PD_to_env(Jericho_data, RdRp_PD)

ccf(RdRp_PD$PD, params_Jericho[,3],  na.action=na.pass)

## gp23
corelations_of_PD_to_env(Jericho_data, gp23_PD)

## 16s
corelations_of_PD_to_env(Jericho_data, S16_PD)

## 18s
corelations_of_PD_to_env(Jericho_data, S18_PD)

## 18s phytos

corelations_of_PD_to_env(Jericho_data, S18_PD_phyto)


## 18s Hetero

corelations_of_PD_to_env(Jericho_data, S18_PD_hetero)

## but would have to exclude the winter times, so that the lags were equal.
## correlations between  RdRp and 18s
correlation_between_2_amplicons <- function (env_data, Amplicon_1_PD, Amplicon_2_PD) {
  VC_dates <- env_data$Date[match(as.Date(Amplicon_1_PD$Date), as.Date(Amplicon_2_PD$Date))]
  
  union_vcs <- intersect(Amplicon_1_PD$Date, Amplicon_2_PD$Date)
  
  Jericho_Amplicon_2_PD <- droplevels(subset(Amplicon_2_PD, Date %in% union_vcs))
  
  Jericho_Amplicon_1_PD <- droplevels(subset(Amplicon_1_PD, Date %in% union_vcs))
  
  print(cor(Jericho_Amplicon_1_PD$PD, Jericho_Amplicon_2_PD$PD))
  
  print( corr.test(as.data.frame(Jericho_Amplicon_1_PD[,"PD"]), as.data.frame(Jericho_Amplicon_2_PD[,"PD"]), method="spearman", use = "pairwise")$r)
  print( corr.test(as.data.frame(Jericho_Amplicon_1_PD[,"PD"]), as.data.frame(Jericho_Amplicon_2_PD[,"PD"]), method="spearman",use = "pairwise")$p)
  
  
  ccf(Jericho_Amplicon_1_PD$PD, Jericho_Amplicon_2_PD$PD, na.action=na.pass)
}

correlation_between_2_amplicons(Jericho_data, S18_PD, RdRp_PD)
correlation_between_2_amplicons(Jericho_data, S18_PD_phyto, RdRp_PD)
correlation_between_2_amplicons(Jericho_data, S18_PD_hetero, RdRp_PD)


## correlations between  gp23 and 16s

correlation_between_2_amplicons(Jericho_data, S16_PD, gp23_PD)

## correlations between  gp23 and 18s
correlation_between_2_amplicons(Jericho_data, gp23_PD, S18_PD)

correlation_between_2_amplicons(Jericho_data, S16_PD, RdRp_PD)

amplicon_richness <- read.csv("../results/amplicon_richness_by_date.csv")
corr.test(as.data.frame(amplicon_richness[,-1]), method="spearman", use = "pairwise")$r
corr.test(as.data.frame(amplicon_richness[,-1]), method="spearman", use = "pairwise")$p

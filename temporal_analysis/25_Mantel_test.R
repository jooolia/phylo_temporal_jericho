## Perform mantel tests
library(vegan)
library(reshape2)
library(dplyr)


normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")

normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 


Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

Env_data_for_merging <- Jericho_data [,c("Date","VC_number")]
row.names(Env_data_for_merging) <- Env_data_for_merging$VC_number

Jericho_data <- read.csv("../results/Jericho_env_data_mean_imputed.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

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

## try with 18s and MPL

mantel_and_partial_mantel_tests <- function (normalized_OTUs_A,normalized_OTUs_B, Env_data, env_params_to_choose) {
#  normalized_OTUs_A <- normalized_MPL_OTUs
#  normalized_OTUs_B <- normalized_18s_OTUs
#  Env_data <- Env_data_for_merging
#  env_params_to_choose <- params_Jericho$Temperature_YSI
 
 list_of_tests <- c()
 new_table <- cbind(normalized_OTUs_A, Env_data[, "Date"][match(rownames(normalized_OTUs_A), rownames(Env_data))])
 row.names(new_table) <- new_table[,(dim(new_table)[2])]
 new_table <- new_table[,-(dim(new_table)[2])]
 
 OTU_diss <- vegdist(new_table)
 SetA_diss <- melt(as.matrix(OTU_diss))
 
 SetB_table <- cbind(normalized_OTUs_B, Env_data[, "Date"][match(rownames(normalized_OTUs_B), rownames(Env_data))])
 row.names(SetB_table) <- SetB_table[,(dim(SetB_table)[2])]
 SetB_table <- SetB_table[,-(dim(SetB_table)[2])]
 
 OTU_diss1 <- vegdist(SetB_table)
 SetB_diss <- melt(as.matrix(OTU_diss1))
 
 test_beta <- merge(SetA_diss, SetB_diss, by=c("Var1", "Var2"))
 test_beta <- filter(test_beta, value.x >0 & value.y >0)
 
 
 SetA_diss_filt <- dcast(test_beta, Var1 ~  Var2, value.var="value.x", fill=0)
 rownames(SetA_diss_filt) <- SetA_diss_filt$Var1
 SetA_diss_filt <- SetA_diss_filt[,-1]
 
 SetB_diss_filt <- dcast(test_beta, Var1 ~  Var2, value.var="value.y",fill=0)
 rownames(SetB_diss_filt) <- SetB_diss_filt$Var1
 SetB_diss_filt <- SetB_diss_filt[,-1]
 
 print(mantel_otu_tables <- mantel(as.dist(as.matrix(SetA_diss_filt)), as.dist(as.matrix(SetB_diss_filt))))
 
 params_Jericho <- subset(params_Jericho, Date %in% as.Date(rownames(SetA_diss_filt) ))
 pc <- prcomp(params_Jericho[,-1], scale = TRUE)
 pc<- scores(pc, display = "sites", choices = 1:4)
 edis <- vegdist(pc, method = "euclid")
 
 #env.dist <- dist(scale(params_Jericho[,-1]))
 
 env.dist <- dist(params_Jericho[,env_params_to_choose])
 
 print(env_otu_table_A_mantel <- mantel(env.dist,
                                        as.dist(as.matrix(SetA_diss_filt)),
                                        method="spearman"))
 print(env_otu_table_B_mantel <-mantel(env.dist,
                                       as.dist(as.matrix(SetB_diss_filt)),
                                       method="spearman"))
 
 partial_mantel <- mantel.partial(as.dist(as.matrix(SetA_diss_filt)),
                                  as.dist(as.matrix(SetB_diss_filt)),
                                  env.dist,
                                  method = "pearson",
                                  permutations = 999, 
                                  strata = NULL,
                                  na.rm = FALSE,
                                  parallel = getOption("mc.cores"))
 
 
 ## set up lagged correlations
 SetA_lead <-  mutate(test_beta,
                      value.x=lead(value.x) )
 SetA_lead <- na.omit(SetA_lead)
 SetA_diss_filt_lead <- dcast(SetA_lead,
                              Var1 ~  Var2,
                              value.var="value.x",
                              fill=0)
 rownames(SetA_diss_filt_lead) <- SetA_diss_filt_lead$Var1
 SetA_diss_filt_lead <- SetA_diss_filt_lead[,-1]
 
print( mantel_lagged_A <- mantel(as.dist(as.matrix(SetA_diss_filt_lead)),
                                 as.dist(as.matrix(SetB_diss_filt))))
 
SetB_lead <-  mutate(test_beta,
                     value.y=lead(value.y) )
SetB_lead <- na.omit(SetB_lead)
SetB_diss_filt_lead <- dcast(SetB_lead,
                             Var1 ~  Var2,
                             value.var="value.x",
                             fill=0)
rownames(SetB_diss_filt_lead) <- SetB_diss_filt_lead$Var1
SetB_diss_filt_lead <- SetB_diss_filt_lead[,-1]

print( mantel_lagged_B <- mantel(as.dist(as.matrix(SetA_diss_filt)),
                                 as.dist(as.matrix(SetB_diss_filt_lead))))

 
list_of_tests <- list(mantel_otu_tables,
                       env_otu_table_A_mantel,
                       env_otu_table_B_mantel,
                       partial_mantel,
                       mantel_lagged_A,
                       mantel_lagged_B)
 
 return(list_of_tests)
 
}

test_result <- mantel_and_partial_mantel_tests(normalized_18s_OTUs,
                                               normalized_MPL_OTUs,
                                               Env_data_for_merging)
test_result[[1]]$statistic
test_result[[1]]$signif
## want to collect all of the mantel test results and see what is happening. 
## also want to look at all of the individual ones. 

table_to_hold_all_mantle <- data.frame(env_param = colnames(subset(params_Jericho, select = !(colnames(params_Jericho)) %in% "Date")))
names_otu_tables <- c("normalized_16s_OTUs",
                      "normalized_18s_OTUs",
                      "normalized_MPL_OTUs",
                      "normalized_gp23_OTUs",
                      "normalized_16s_OTUs_lag",
                      "normalized_18s_OTUs_lag",
                      "normalized_MPL_OTUs_lag",
                      "normalized_gp23_OTUs_lag")

table_to_hold_all_mantle <- as.data.frame(c(as.character(table_to_hold_all_mantle$env_param),
                                            names_otu_tables))

colnames(table_to_hold_all_mantle)[1] <- "env_param"

rownames(table_to_hold_all_mantle) <- table_to_hold_all_mantle$env_param
table_to_hold_all_mantle[,"env_param"] <- NULL

## need to make all the possibilities?
for (x in names_otu_tables){
 print(x)
 table_to_hold_all_mantle[,x] <- cbind(rep(0, dim(table_to_hold_all_mantle)[1])) 
}

mantel_test_env_and_otu_tables <- function (normalized_A_OTUs,
                                            normalized_B_OTUs,
                                            table_to_hold_all_mantle,
                                            string_table_A,
                                            string_table_B) {
 for (x in colnames(subset(params_Jericho, select = !(colnames(params_Jericho)) %in% "Date"))){
  
  print(x)
  
  mantel_results <- mantel_and_partial_mantel_tests(normalized_A_OTUs,
                                                    normalized_B_OTUs,
                                                    Env_data_for_merging,
                                                    x )
  mystars <- ifelse(mantel_results[[1]]$signif < .001, "***",
                    ifelse(mantel_results[[1]]$signif < .01, "** ", 
                           ifelse(mantel_results[[1]]$signif < .05, "* ",
                                  " ")))
  
  cor_with_stars <- paste(round(mantel_results[[1]]$statistic, digits=2),
                          " ",
                          mystars)

  table_to_hold_all_mantle[string_table_A,string_table_B] <- cor_with_stars
  table_to_hold_all_mantle[string_table_B,string_table_A] <- cor_with_stars
  
  mystars <- ifelse(mantel_results[[2]]$signif < .001,"***", 
                    ifelse(mantel_results[[2]]$signif < .01,"** ", 
                           ifelse(mantel_results[[2]]$signif < .05, "* ",
                                  " ")))
  cor_with_stars <- paste(round(mantel_results[[2]]$statistic, digits=2), " ", mystars)  
  table_to_hold_all_mantle[x, string_table_A] <- cor_with_stars
  
  mystars <- ifelse(mantel_results[[3]]$signif < .001, "***",
                    ifelse(mantel_results[[3]]$signif < .01, "** ",
                           ifelse(mantel_results[[3]]$signif < .05, "* ",
                                  " ")))
  cor_with_stars <- paste(round(mantel_results[[3]]$statistic, digits=2), " ", mystars)    
  table_to_hold_all_mantle[x, string_table_B] <- cor_with_stars
  
  mystars <- ifelse(mantel_results[[5]]$signif < .001, "***",
                    ifelse(mantel_results[[5]]$signif < .01, "** ",
                           ifelse(mantel_results[[5]]$signif < .05, "* ",
                                  " ")))
  cor_with_stars <- paste(round(mantel_results[[5]]$statistic, digits=2), " ", mystars)  
  
  table_to_hold_all_mantle[paste0(string_table_B, "_lag"), string_table_A] <- cor_with_stars
  
  
  mystars <- ifelse(mantel_results[[6]]$signif < .001, "***",
                    ifelse(mantel_results[[6]]$signif < .01, "** ",
                           ifelse(mantel_results[[6]]$signif < .05, "* ",
                                  " ")))
  cor_with_stars <- paste(round(mantel_results[[6]]$statistic, digits=2), " ", mystars)  
  
  table_to_hold_all_mantle[paste0(string_table_A, "_lag"), string_table_B] <- cor_with_stars
 
  
  ## want to add these to a data frame....
  ## for each combination of tables and env parameters want to make a table
 }
 return(table_to_hold_all_mantle)
}


table_to_hold_all_mantle <- mantel_test_env_and_otu_tables(normalized_18s_OTUs,
                                                           normalized_MPL_OTUs,
                                                           table_to_hold_all_mantle,
                                                           "normalized_18s_OTUs",
                                                           "normalized_MPL_OTUs")

table_to_hold_all_mantle <- mantel_test_env_and_otu_tables(normalized_16s_OTUs,
                                                           normalized_gp23_OTUs,
                                                           table_to_hold_all_mantle,
                                                           "normalized_16s_OTUs",
                                                           "normalized_gp23_OTUs")

table_to_hold_all_mantle <- mantel_test_env_and_otu_tables(normalized_18s_OTUs,
                                                           normalized_16s_OTUs,
                                                           table_to_hold_all_mantle,
                                                           "normalized_18s_OTUs",
                                                           "normalized_16s_OTUs")

table_to_hold_all_mantle <- mantel_test_env_and_otu_tables(normalized_18s_OTUs,
                                                           normalized_gp23_OTUs,
                                                           table_to_hold_all_mantle,
                                                           "normalized_18s_OTUs",
                                                           "normalized_gp23_OTUs")

table_to_hold_all_mantle <- mantel_test_env_and_otu_tables(normalized_MPL_OTUs,
                                                           normalized_gp23_OTUs,
                                                           table_to_hold_all_mantle,
                                                           "normalized_MPL_OTUs",
                                                           "normalized_gp23_OTUs")

table_to_hold_all_mantle <- mantel_test_env_and_otu_tables(normalized_MPL_OTUs,
                                                           normalized_16s_OTUs,
                                                           table_to_hold_all_mantle,
                                                           "normalized_MPL_OTUs",
                                                           "normalized_16s_OTUs")


table_test <- table_to_hold_all_mantle

write.csv(table_to_hold_all_mantle, file="../results/mantel_tests_otus_and_env.csv")

write.table(table_test, file="../results/mantel_tests_otus_and_env.txt")

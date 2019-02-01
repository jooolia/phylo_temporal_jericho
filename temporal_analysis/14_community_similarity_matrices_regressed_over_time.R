## Like in Needham 2013, want to look at the community similarity matrices to see what it looks like over time. 
## make distance matrices over time and perform regression

library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv",
                                  row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",
                                   row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv",
                                  row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv",
                                  row.names="VC_number") 

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
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv",
                         row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)
Env_data_for_merging <- Jericho_data [,c("Date","VC_number")]
row.names(Env_data_for_merging) <- Env_data_for_merging$VC_number

lm_eqn = function(m) {
eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                 list(a = format(coef(m)[1], digits = 2), 
                      b = format(coef(m)[2], digits = 2), 
                      r2 = format(summary(m)$r.squared, digits = 3)))
as.character(as.expression(eq))  
}

plot_normalized_OTUs_by_date_and_community_dissimilarity <- function (normalized_otus, title) {
 new_table <- cbind(normalized_otus,
                    Env_data_for_merging[, "Date"][match(rownames(normalized_otus), 
                                                         rownames(Env_data_for_merging))])
 row.names(new_table) <- new_table[,(dim(new_table)[2])]
 new_table <- new_table[,-(dim(new_table)[2])]
 
 OTU_diss <- vegdist(new_table)
 ## get all the distances by VC
 
 OTU_matrix <- as.matrix(OTU_diss)
 
 OTU_matrix[upper.tri(OTU_matrix )] <- NA
 ## need to replace all the VC numbers with dates
 melted_matrix <- melt(OTU_matrix, 
                       value.name = "Bray_curtis_dissimilarity")
 
 melted_matrix <- na.omit(melted_matrix) ## gets rid of the leftover diagonals
 ## get rid of the ones where the date is the same:
 melted_matrix <- filter(melted_matrix, 
                         Var1 !=Var2)
 
 melted_matrix$Var1 <- as.Date(melted_matrix$Var1)
 melted_matrix$Var2 <- as.Date(melted_matrix$Var2)
 
 melted_matrix$date_diff <- as.numeric(abs(difftime(melted_matrix$Var1,
                                                    melted_matrix$Var2,
                                                    units = "days")))
 
 p <-  ggplot(melted_matrix,
              aes(x=date_diff,
                  y=Bray_curtis_dissimilarity))+
  geom_point()+
  ggtitle(title)+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold",
                                    size=20),
        axis.text.x  = element_text(size=16),
        axis.line=element_line(size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour="grey20"))
 return(p)
}

pdf("../figures/Community_dissimilarity_matrix_regressed_over_time%03d.pdf",
    width = 15,
    height = 11,
    onefile = FALSE)
plot_normalized_OTUs_by_date_and_community_dissimilarity(normalized_18s_OTUs, 
                                                         "18S by date")
plot_normalized_OTUs_by_date_and_community_dissimilarity(normalized_18s_OTUs,
                                                         "18S by date")+
 geom_smooth(method = "lm",
             se=FALSE,
             formula = y ~ x)+
 geom_text(aes(x = 200,
               y = 0.25,
               label = lm_eqn(lm(Bray_curtis_dissimilarity ~ date_diff))),
           parse = TRUE)

plot_normalized_OTUs_by_date_and_community_dissimilarity(normalized_gp23_OTUs,
                                                         "gp23 by date")

plot_normalized_OTUs_by_date_and_community_dissimilarity(normalized_gp23_OTUs,
                                                         "gp23 by date")+
 geom_smooth(method = "lm",
             se=FALSE, 
             formula = y ~ x)+
 geom_text(aes(x = 200,
               y = 0.25,
               label = lm_eqn(lm(Bray_curtis_dissimilarity ~ date_diff))),
           parse = TRUE)

plot_normalized_OTUs_by_date_and_community_dissimilarity(normalized_MPL_OTUs,
                                                         "MPL by date")

plot_normalized_OTUs_by_date_and_community_dissimilarity(normalized_MPL_OTUs,
                                                         "MPL by date")+
 geom_smooth(method = "lm",
             se=FALSE,
             formula = y ~ x)+
 geom_text(aes(x = 200,
               y = 0.25,
               label = lm_eqn(lm(Bray_curtis_dissimilarity ~ date_diff))),
           parse = TRUE)

plot_normalized_OTUs_by_date_and_community_dissimilarity(normalized_16s_OTUs,
                                                         "16s by date")

plot_normalized_OTUs_by_date_and_community_dissimilarity(normalized_16s_OTUs,
                                                         "16s by date")+ 
 geom_smooth(method = "lm",
             se=FALSE,
             formula = y ~ x)+
 geom_text(aes(x = 200,
               y = 0.25,
               label = lm_eqn(lm(Bray_curtis_dissimilarity ~ date_diff))),
           parse = TRUE)
dev.off()


### Do it with comunity similarity

plot_normalized_OTUs_by_date_and_community_similarity <- function (normalized_otus,
                                                                   title) {
 
 new_table <- cbind(normalized_otus,
                    Env_data_for_merging[, "Date"][match(rownames(normalized_otus), rownames(Env_data_for_merging))])
 row.names(new_table) <- new_table[,(dim(new_table)[2])]
 new_table <- new_table[,-(dim(new_table)[2])]
 
 OTU_diss <- vegdist(new_table)
 ## get all the distances by VC
 
 OTU_matrix <- as.matrix(OTU_diss)
 
 OTU_matrix[upper.tri(OTU_matrix )] <- NA
 
 ## need to replace all the VC numbers with dates
 melted_matrix <- melt(OTU_matrix,
                       value.name = "Bray_curtis_dissimilarity")
 melted_matrix <- na.omit(melted_matrix) ## gets rid of the leftover diagonals
 ## get rid of the ones where the date is the same:
 melted_matrix <- filter(melted_matrix,
                         Var1 != Var2)
 melted_matrix$Var1 <- as.Date(melted_matrix$Var1)
 melted_matrix$Var2 <- as.Date(melted_matrix$Var2)
 melted_matrix$Bray_curtis_similarity <- 1-melted_matrix$Bray_curtis_dissimilarity
 
 melted_matrix$date_diff <- as.numeric(abs(difftime(melted_matrix$Var1,
                                                    melted_matrix$Var2,
                                                    units = "days")))
 
 p <-  ggplot(melted_matrix,
              aes(x=date_diff,
                  y=Bray_curtis_similarity))+
  geom_point(colour=point_colour)+
  ggtitle(title)+
  theme_JAG_presentation()+
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=16),
        axis.line=element_line(size = 0.3))+
  ylim(0,1)+
  xlim(0,400)
 return(p)
}

pdf("../figures/Community_similarity_matrix_regressed_over_time%03d.pdf",
    width = 15,
    height = 11,
    onefile = FALSE)
plot_normalized_OTUs_by_date_and_community_similarity(normalized_18s_OTUs,
                                                      "18S by date") 

plot_normalized_OTUs_by_date_and_community_similarity(normalized_gp23_OTUs,
                                                      "gp23 by date")
plot_normalized_OTUs_by_date_and_community_similarity(normalized_MPL_OTUs,
                                                      "MPL by date")
plot_normalized_OTUs_by_date_and_community_similarity(normalized_16s_OTUs,
                                                      "16s by date")
dev.off()

png("../figures/Community_similarity_matrix_regressed_over_time_one_pager.png",
    width = 1500,
    height = 1500)

grid.arrange( 
 plot_normalized_OTUs_by_date_and_community_similarity(normalized_MPL_OTUs,
                                                       "MPL by date")+
               theme(axis.ticks.x = element_blank(), 
                     axis.text.x = element_blank(), 
                     axis.title.x=element_blank())+
               xlab(NULL),
 plot_normalized_OTUs_by_date_and_community_similarity(normalized_gp23_OTUs,
                                                       "gp23 by date")+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x=element_blank())+
  xlab(NULL),
 plot_normalized_OTUs_by_date_and_community_similarity(normalized_18s_OTUs,
                                                       "18S by date")+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x=element_blank())+
  xlab(NULL) ,
 plot_normalized_OTUs_by_date_and_community_similarity(normalized_16s_OTUs,
                                                       "16s by date"),
  ncol=1)
dev.off()


## Beta diversity compared: ####

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


compare_beta_diversity <- function (normalized_A_OTUs,
                                    normalized_B_OTUs) {
#  normalized_A_OTUs <- normalized_18s_OTUs
 # normalized_B_OTUs <- normalized_MPL_OTUs
  new_table <- cbind(normalized_A_OTUs,
                     Env_data_for_merging[, "Date"][match(rownames(normalized_A_OTUs),
                                                          rownames(Env_data_for_merging))])
  row.names(new_table) <- new_table[,(dim(new_table)[2])]
  new_table <- new_table[,-(dim(new_table)[2])]
  
  OTU_diss <- vegdist(new_table)
  A_diss <- melt(as.matrix(OTU_diss))
  
  B_table <- cbind(normalized_B_OTUs,
                     Env_data_for_merging[, "Date"][match(rownames(normalized_B_OTUs),
                                                          rownames(Env_data_for_merging))])
  row.names(B_table) <- B_table[,(dim(B_table)[2])]
  B_table <- B_table[,-(dim(B_table)[2])]
  
  OTU_diss1 <- vegdist(B_table)
  B_diss <- melt(as.matrix(OTU_diss1))
  
  test_beta <- merge(A_diss,
                     B_diss,
                     by=c("Var1", "Var2"))
 test_beta <- dplyr::filter(test_beta,
                            value.x >0 & value.y >0)
  str(test_beta)

 p1 <-  ggplot(test_beta,
               aes(x=value.x,
                   y=value.y))+
   geom_point(colour=point_colour)
  print(p1)
  
  print(lm_comparison <- lm(value.y ~ value.x, test_beta))
  
 p2 <-  ggplot(test_beta,
               aes(x=value.x,
                   y=value.y))+
   geom_point(colour=point_colour)+
   stat_smooth(method="lm", se=FALSE)+ 
  ## problem with this in the function. If it is run  beforehand it is evaluated, but it always takes the old version. 
  geom_text(aes(x = min(value.x)+0.1,
                y = max(value.y)-0.1), 
            label = lm_eqn(lm_comparison),
            parse = TRUE)+
  theme_bw()
  print(p2)
}

compare_beta_diversity(normalized_18s_OTUs,
                       normalized_MPL_OTUs)

compare_beta_diversity(normalized_16s_OTUs,
                       normalized_gp23_OTUs)

compare_beta_diversity(normalized_16s_OTUs,
                       normalized_MPL_OTUs)

compare_beta_diversity(normalized_18s_OTUs,
                       normalized_MPL_OTUs)

compare_beta_diversity(normalized_18s_OTUs,
                       normalized_16s_OTUs)

compare_beta_diversity(normalized_MPL_OTUs,
                       normalized_gp23_OTUs)

## Calculate PCA and RDA for all the OTU tables

library(vegan)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(packfor)
library(geosphere)
library("ggvegan")
library("grid")
library(gridExtra)
library(venneuler)


args <- commandArgs(TRUE)
inputFile <- args[1]

## test to see if input file is given, so I can decide whether to use this argument or the orginal one. 
if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
}


### Load in data ###
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv", row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
#normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_concat.tsv", row.names="VC_number")
normalized_AVS <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R1.tsv", row.names="VC_number")
#normalized_AVS_R2_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS_R2.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 



original_jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
original_jericho_data$Date <- as.Date(original_jericho_data$Date)

## from 001 script
## don't think should be using transformed
# Jericho_data <- read.csv("../results/Jericho_env_data_mean_imputed_transformed.csv", row.names=1)
# Jericho_data$Date <- as.Date(Jericho_data$Date)

Jericho_data <- read.csv("../results/Jericho_env_data_mean_imputed.csv", row.names=1)
 Jericho_data$Date <- as.Date(Jericho_data$Date)

## amplicon richness calculated in earlier script.
amplicon_richness <- read.csv("../results/amplicon_richness_by_date.csv")
amplicon_richness$Date <- as.Date(amplicon_richness$Date)

## exclude season and NO2
parameters_to_exclude <- c("season",
                           "Average_NO3_NO2",
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
params_Jericho <- params_Jericho[,-12]

VC_dates <- Jericho_data$Date[match(rownames(normalized_18s_OTUs), original_jericho_data$VC_number)]
params_Jericho_18s <- subset(params_Jericho, Date %in% VC_dates)
## leave out phosphate


amplicon_richness_for_study <- subset(amplicon_richness, Date %in% VC_dates )
## temporary fix
amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0

#Physiogrpahy
envbio <- params_Jericho_18s[,c("Average_viral_abundance","Average_bacterial_abundance", "Average_chl_a")]
names(envbio)

envchem <- params_Jericho_18s[,c("Average_PO4", "Average_SiO2",
                                # "Average_NO3_NO2",
                                 "Temperature_YSI", "Salinity_ppt_YSI", "Dissolved_oxygen_percent" , "pH", "Tide_height")]
names(envchem)
envtemp <- params_Jericho_18s[,c("day_length", "month_number")]
envrich <- amplicon_richness_for_study[,c("richness.AVS","richness.MPL","richness.16S","richness.gp23")]


########################

## Now analyse the community matrix data

#Now I think that really I should try by removing all the metazoa, fungi and plants. Need to work on the OTU table and then retry to the analysis
taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)
## so want to exclude those with Order Metazoa, Fungi, 
#taxonomy_18s_no_met_fun <- subset(taxonomy_18s, !(Order %in% c("Metazoa","Fungi")))
taxonomy_18s_no_met_fun <- subset(taxonomy_18s,!(Phylum %in% "Opisthokonta"))

colnames(normalized_18s_OTUs) <- gsub(".size.*.", "", colnames(normalized_18s_OTUs))
## so now pull out OTUs from that data frame

test_hel <- decostand(normalized_18s_OTUs, "hellinger")
spe.rda <- rda(test_hel ~ ., params_Jericho_18s[,-1])

normalized_18s_OTUs_no_met_fun <- subset(normalized_18s_OTUs, select=(colnames(normalized_18s_OTUs) %in% taxonomy_18s_no_met_fun$otu_number ))

RDA_of_OTU_table <- function (OTU_table,OTU_table_name, transformation, env_data, env_data_name) {
#   OTU_table <- normalized_18s_OTUs_no_met_fun
#   transformation <- "hellinger"
#   env_data <- params_Jericho_18s[,-1]
 
  #scaled_env <- as.data.frame(scale(env_data))
  spe.hel <- decostand(OTU_table, transformation)
  spe.rda <- rda(spe.hel ~ ., env_data)
  head(summary(spe.rda))
  coef(spe.rda) ## gives canonical coefficients
  ## Unadjusted R^2 retrieved form the rda result
  R2 <- RsquareAdj(spe.rda)$r.squared
  
  # Adjusted R^2 retrieved from the rda object
  R2adj <- RsquareAdj(spe.rda)$adj.r.squared
  print(anova_rda <- anova.cca(spe.rda, step=1000))
  print(anova_rda <- anova.cca(spe.rda, by="axis", step=1000))
  #It can be also negative - it means that real variables explain less variation than would explain the same number of randomly generated variables.
  pdf(paste("../figures/Triplot_RDA_",OTU_table_name,"_with_", env_data_name,"_plots.pdf",sep=""), width = 15, height = 15, onefile = TRUE)
  ## Triplots of the rda results
  ## Scaling 1: distance triplot
 
  plot(spe.rda, scaling=1, main=paste("Triplot RDA \n", OTU_table_name," with ", transformation," ~ ", env_data_name, "\n -scaling 1 - wa scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value ", anova_rda$"Pr(>F)"[1],collapse="\n"))
  ## add in arrows for speices
  ## note you have to specify the scaling
  spe.sc <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
  arrows(0,0, spe.sc[,1], spe.sc[,2], length=0,lty=1, col="red")
 
    ## scaling 2
  plot(spe.rda, scaling=2, main=paste("Triplot RDA \n", OTU_table_name," with ", transformation," ~ ", env_data_name, "\n -scaling 2 - wa scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value ", anova_rda$"Pr(>F)"[1],collapse="\n"))
  ## add in arrows for speices
  spe2.sc <- scores(spe.rda, choices=1:2, scaling=2, display="sp")
  arrows(0,0, spe2.sc[,1], spe2.sc[,2], length=0,lty=1, col="red")
  
  ## Scaling 1
  plot(spe.rda, scaling=1, display=c("sp", "lc", "cn"),main=paste("Triplot RDA \n", OTU_table_name," with ", transformation," ~ ", env_data_name, "\n -scaling 1 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
  arrows(0,0, spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
  
  #Scaling 2
  plot(spe.rda,  display=c("sp", "lc", "cn"),main=paste("Triplot RDA \n", OTU_table_name," with ", transformation," ~ ", env_data_name, "\n -scaling 2 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
  arrows(0,0, spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col="red")
  
  ## testing out side by side plot
  layout(matrix(1:4, ncol = 2))
  plot(spe.rda, display = "sp",type = "n", scaling = 1)
  orditorp(spe.rda, display = "sites", scaling = 1)
 # ordiarrows(spe.rda, "species")
  arrows(0,0, spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
  
  ordiplot (spe.rda, display = 'sp', type = 'n', scaling = 1)
  orditorp (spe.rda, display = 'sp',  scaling = 1)
  
  
  plot(spe.rda, display = "sites",type = "n", scaling = 1)
  orditorp(spe.rda, display = "sites", scaling = 1)
  # ordiarrows(spe.rda, "species")
  arrows(0,0, spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
  
  ordiplot (spe.rda, display = 'sites', type = 'n', scaling = 1)
  orditorp (spe.rda, display = 'sp',  scaling = 1)
  
  layout(1)

  ## use fortify method to extract scores in ggplot-friendly format
  scrs <- fortify(spe.rda, scaling = 1)
  #print(scrs)
  ## take only site scores for this
  sites <- with(scrs, scrs[Score == "sites", ])
  species <- with(scrs, scrs[Score == "species", ])
  biplot <- with(scrs, scrs[Score == "biplot",])
  ## add in something to plot for the arrows

  ## add in arrows
  plt <- ggplot(species, aes(x = Dim1, y = Dim2, label=Label)) + 
   geom_segment(aes(x=0, y=0, xend=Dim1, yend=Dim2), colour="purple", arrow = arrow(length = unit(0.5, "cm")))+
   coord_fixed(xlim = c(-1,1), ylim = c(-1,1))+
   #geom_text(hjust=1.2)+
   geom_point(data=sites,aes(x = Dim1, y = Dim2, label=Label), colour="blue")+
   theme_JAG_presentation()
#print(plt)
  
  plt2 <- ggplot(species, aes(x = Dim1, y = Dim2, label=Label)) + 
   geom_point(data=sites,aes(x = Dim1, y = Dim2, label=Label), colour="blue")+
   coord_fixed(xlim = c(-1,1), ylim = c(-1,1))+
   #geom_text(hjust=1.2)+
   geom_segment(data=biplot, aes(x=0, y=0, xend=Dim1, yend=Dim2), colour="green",arrow = arrow(length = unit(0.5, "cm")))+
   geom_text(data=biplot, aes(x=Dim1, y=Dim2), colour="black", size=3)+
   theme_JAG_presentation()
  #print(plt2)
  
  plt3 <- ggplot(species, aes(x = Dim1, y = Dim2, label=Label)) + 
  # geom_point() +
   coord_fixed(xlim = c(-1,1), ylim = c(-1,1))+
   geom_text(hjust=1.2, colour="purple",size=3)+
   geom_segment(data=biplot, aes(x=0, y=0, xend=Dim1, yend=Dim2), colour="green",arrow = arrow(length = unit(0.5, "cm")))+
   geom_text(data=biplot, aes(x=Dim1, y=Dim2), colour="black", size=3 )+
   theme_JAG_presentation()
  #print(plt2)
  
  plt4 <- ggplot(species, aes(x = Dim1, y = Dim2, label=Label)) + 
 #  geom_point() +
   coord_fixed(xlim = c(-1,1), ylim = c(-1,1))+
   geom_text(hjust=1.2, colour="purple", size=3)+
   geom_segment(data=biplot, aes(x=0, y=0, xend=Dim1, yend=Dim2), colour="green",arrow = arrow(length = unit(0.5, "cm")))+
   geom_text(data=biplot, aes(x=Dim1, y=Dim2), colour="black", size=3)+
   geom_point(data=sites,aes(x = Dim1, y = Dim2, label=Label), colour="blue")+
   theme_JAG_presentation()
  #print(plt2)
 
  grid.arrange(plt, plt2,plt3,plt4, ncol=2)
  
  ## use fortify method to extract scores in ggplot-friendly format
  scrs <- fortify(spe.rda, scaling = 2)
 # print(scrs)
  ## take only site scores for this
  sites <- with(scrs, scrs[Score == "sites", ])
  species <- with(scrs, scrs[Score == "species", ])
  biplot <- with(scrs, scrs[Score == "biplot",])
  ## add in something to plot for the arrows
  
  ## add in arrows
  plt <- ggplot(species, aes(x = Dim1, y = Dim2, label=Label)) + 
   geom_segment(aes(x=0, y=0, xend=Dim1, yend=Dim2), colour="purple", arrow = arrow(length = unit(0.2, "cm")))+
   coord_fixed(xlim = c(-1,1), ylim = c(-1,1))+
   #geom_text(hjust=1.2)+
   geom_point(data=sites,aes(x = Dim1, y = Dim2, label=Label), colour="blue")+
   theme_JAG_presentation()
  #print(plt)
  
  plt2 <- ggplot(species, aes(x = Dim1, y = Dim2, label=Label)) + 
   geom_point(data=sites,aes(x = Dim1, y = Dim2, label=Label), colour="blue")+
   coord_fixed(xlim = c(-1,1), ylim = c(-1,1))+
   #geom_text(hjust=1.2)+
   geom_segment(data=biplot, aes(x=0, y=0, xend=Dim1, yend=Dim2), colour="green",arrow = arrow(length = unit(0.5, "cm")))+
   geom_text(data=biplot, aes(x=Dim1, y=Dim2), colour="black", size=3)+
   theme_JAG_presentation()

  
  plt3 <- ggplot(species, aes(x = Dim1, y = Dim2, label=Label)) + 
   # geom_point() +
   coord_fixed(xlim = c(-1,1), ylim = c(-1,1))+
   geom_text(hjust=1.2, colour="purple", size=3)+
   geom_segment(data=biplot, aes(x=0, y=0, xend=Dim1, yend=Dim2), colour="green",arrow = arrow(length = unit(0.5, "cm")))+
   geom_text(data=biplot, aes(x=Dim1, y=Dim2), colour="black", size=5)+
   theme_JAG_presentation()

  
  plt4 <- ggplot(species, aes(x = Dim1, y = Dim2, label=Label)) + 
   #  geom_point() +
   coord_fixed(xlim = c(-1,1), ylim = c(-1,1))+
   geom_text(hjust=1.2, colour="purple", size=3)+
   geom_segment(data=biplot, aes(x=0, y=0, xend=Dim1, yend=Dim2), colour="green",arrow = arrow(length = unit(0.5, "cm")))+
   geom_text(data=biplot, aes(x=Dim1, y=Dim2), colour="black", size=3)+
   geom_point(data=sites,aes(x = Dim1, y = Dim2, label=Label), colour="blue")+
   theme_JAG_presentation()

  grid.arrange(plt, plt2,plt3,plt4, ncol=2)
  
  
    dev.off()
}


## RDA not significant so therefore commented out. 
# RDA_18s_no_metazoan_no_fungi <- RDA_of_OTU_table(normalized_18s_OTUs_no_met_fun,"normalized-18s-OTUs-no-met-fun","hellinger", params_Jericho_18s[,-1], "params-Jericho-no-date")

## test out top 100 or 50
find_top50_OTUs <- function (long_otus) {
 summarised_by_OTU <- ddply(long_otus,~Var2,summarise,relative_abundance=sum(value))
 OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
 top_50 <- droplevels(head(OTU_sorted_by_rel_abun, n=50))
 return(top_50)
} 

melted_otus <- melt(as.matrix(normalized_18s_OTUs_no_met_fun))

top50_18s <- find_top50_OTUs(melted_otus)

top50_normalized_18s_no_met_fun <- subset(normalized_18s_OTUs_no_met_fun, select= colnames(normalized_18s_OTUs_no_met_fun) %in% top50_18s$Var2)

## not significant overall. 
#RDA_18s_top_50 <- RDA_of_OTU_table(top50_normalized_18s_no_met_fun,"top50_normalized-18s-OTUs-fam-sum-otus","hellinger", params_Jericho_18s[,-1], "params-Jericho-no-date")

find_top20_OTUs <- function (long_otus) {
 summarised_by_OTU <- ddply(long_otus,~Var2,summarise,relative_abundance=sum(value))
 OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
 top_20 <- droplevels(head(OTU_sorted_by_rel_abun, n=20))
 return(top_20)
} 

top20_18s <- find_top20_OTUs(melted_otus)

top20_normalized_18s_no_met_fun <- subset(normalized_18s_OTUs_no_met_fun, select= colnames(normalized_18s_OTUs_no_met_fun) %in% top20_18s$Var2)

## not significant overall
# RDA_18s_top_20 <- RDA_of_OTU_table(top20_normalized_18s_no_met_fun,"top20_normalized-18s-OTUs-fam-sum-otus","hellinger", params_Jericho_18s[,-1], "params-Jericho-no-date")

### Summarise 18s OTUs to the Family level ####
melted_otus$Family <- taxonomy_18s_no_met_fun$Family[match(melted_otus$Var2, taxonomy_18s_no_met_fun$otu_number)]
## remove unclassified
melted_otus <- subset(melted_otus, Family != "unclassified"  )

fam_melted <- melted_otus %>%
 group_by(Var1,Family) %>%
 summarise(total=sum(value))

## Ok need to put it back together again.
## Families on top
fam_sum_otus <- dcast(fam_melted, Var1 ~  Family)
## need to make them rownames again
row.names(fam_sum_otus) <- fam_sum_otus$Var1
fam_sum_otus <- fam_sum_otus[,-1]

### Summarise 18s OTUs to the order level ####
melted_otus <- melt(as.matrix(normalized_18s_OTUs_no_met_fun))
melted_otus$Order <- taxonomy_18s_no_met_fun$Order[match(melted_otus$Var2, taxonomy_18s_no_met_fun$otu_number)]
## remove unclassified
melted_otus <- subset(melted_otus, Order != "unclassified"  )

ord_melted <- melted_otus %>%
 group_by(Var1,Order) %>%
 summarise(total=sum(value))

## Ok need to put it back together again.
## Families on top
ord_sum_otus <- dcast(ord_melted, Var1 ~  Order)
## need to make them rownames again
row.names(ord_sum_otus) <- ord_sum_otus$Var1
ord_sum_otus <- ord_sum_otus[,-1]


## overall model significant
## axies are not significant except for the fourth one. 
#RDA_18s_summarised_by_fam <- RDA_of_OTU_table(fam_sum_otus,"normalized-18s-OTUs-fam-sum-otus","hellinger", params_Jericho_18s[,-1], "params-Jericho-no-date")

test_hel <- decostand(fam_sum_otus, "hellinger")
test_cca <- cca(test_hel ~ .,data=params_Jericho_18s[,-1]) 
test_cca
summary(test_cca)
print(anova_cca <- anova.cca(test_cca, step=1000)) # significant
print(anova_rda <- anova.cca(test_cca, by="axis", step=1000)) #5 axes show some significance
plot(test_cca)

## overall model significatn
## axes are not significat except fro 9th one. 
#RDA_18s_summarised_by_ord <- RDA_of_OTU_table(ord_sum_otus,"normalized-18s-OTUs-ord-sum-otus","hellinger", params_Jericho_18s[,-1], "params-Jericho-no-date")

test_hel <- decostand(ord_sum_otus, "hellinger")
test_cca <- cca(test_hel ~ .,data=params_Jericho_18s[,-1]) 
test_cca
summary(test_cca)
print(anova_cca <- anova.cca(test_cca, step=1000))
print(anova_rda <- anova.cca(test_cca, by="axis", step=1000)) ## first 2 significant below 0.005, 5 significant axes
plot(test_cca)

## Partial RDA: effect of water chemistry, holding physiography constanct
#3 X and W may be seprarte tables of quant variables
partial_RDAs <- function (OTU_table,OTU_table_name, transformation,env_data_name,envchem,envbio, envtemp, envrich) {
 spe.hel <- decostand(OTU_table, transformation)
 spechem.bio <- rda(spe.hel, envchem, envbio)
 print(spechem.bio)
 
 R2 <- RsquareAdj(spechem.bio)$r.squared
 
 # Adjusted R^2 retrieved from the rda object
 R2adj <- RsquareAdj(spechem.bio)$adj.r.squared
 print(anova_rda <- anova.cca(spechem.bio, step=1000))
 #print(anova_rda_axes <- anova.cca(spechem.bio, by="axis", step=1000))
 
 pdf(paste("../figures/Partial_RDA_",OTU_table_name,"_with_", env_data_name,"_plots.pdf",sep=""), width = 15, height = 15, onefile = TRUE)
 #Scaling 1
 plot(spechem.bio, scaling=1, display=c("sp", "lc", "cn"), 
      main=paste("Partial RDA spe.hel ~chem | bio Top =scaling1 -1c scores\n", OTU_table_name," with ", transformation," ~ ","\n -scaling 1 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe3.sc <- scores(spechem.bio, choices=1:2, scaling=1, display="sp")
 arrows(0,0, spe3.sc[,1], spe3.sc[,2], length=0, lty=1, col="red")
 
 ## Scaling 2
 plot(spechem.bio, scaling=2, display=c("sp", "lc", "cn"), 
      main=paste("Partial RDA spe.hel ~chem | bio Top =scaling1 -1c scores\n", OTU_table_name," with ", transformation," ~ ","\n -scaling 2 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe4.sc <- scores(spechem.bio, choices=1:2, display="sp")
 arrows(0,0, spe4.sc[,1], spe4.sc[,2], length=0, lty=1, col="red")
 #vif.cca(spe.rda)
 vif.cca(spechem.bio)
 
 spebio.chem <- rda(spe.hel, envbio, envchem)
 R2 <- RsquareAdj(spebio.chem)$r.squared
 
 # Adjusted R^2 retrieved from the rda object
 R2adj <- RsquareAdj(spebio.chem)$adj.r.squared
 print(anova_rda <- anova.cca(spebio.chem, step=1000))
 #print(anova_rda_axes <- anova.cca(spebio.chem, by="axis", step=1000))
 
 #Scaling 1
 plot(spebio.chem, scaling=1, display=c("sp", "lc", "cn"), 
      main=paste("Partial RDA spe.hel ~bio | chem Top =scaling1 -1c scores\n", OTU_table_name," with ", transformation," ~ ","\n -scaling 1 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe3.sc <- scores(spebio.chem, choices=1:2, scaling=1, display="sp")
 arrows(0,0, spe3.sc[,1], spe3.sc[,2], length=0, lty=1, col="red")
 
 ## Scaling 2
 plot(spebio.chem, scaling=2, display=c("sp", "lc", "cn"), 
      main=paste("Partial RDA spe.hel ~bio | chem Top =scaling1 -1c scores\n", OTU_table_name," with ", transformation," ~ ","\n -scaling 2 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe4.sc <- scores(spebio.chem, choices=1:2, display="sp")
 arrows(0,0, spe4.sc[,1], spe4.sc[,2], length=0, lty=1, col="red") 
 
 vif.cca(spebio.chem) 
 
 spebio.temp <- rda(spe.hel, envbio, envtemp)
 R2 <- RsquareAdj(spebio.temp)$r.squared
 
 # Adjusted R^2 retrieved from the rda object
 R2adj <- RsquareAdj(spebio.temp)$adj.r.squared
 print(anova_rda <- anova.cca(spebio.temp, step=1000))
 #print(anova_rda_axes <- anova.cca(spebio.temp, by="axis", step=1000))
 
  #Scaling 1
 plot(spebio.temp, scaling=1, display=c("sp", "lc", "cn"), 
      main=paste("Partial RDA spe.hel ~~bio | temp  Top =scaling1 -1c scores\n", OTU_table_name," with ", transformation," ~ ","\n -scaling 1 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe3.sc <- scores(spebio.temp, choices=1:2, scaling=1, display="sp")
 arrows(0,0, spe3.sc[,1], spe3.sc[,2], length=0, lty=1, col="red")
 
 ## Scaling 2
 plot(spebio.temp, scaling=2, display=c("sp", "lc", "cn"), 
      main=paste("Partial RDA spe.hel ~~bio | temp  Top =scaling1 -1c scores\n", OTU_table_name," with ", transformation," ~ ","\n -scaling 2 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe4.sc <- scores(spebio.temp, choices=1:2, display="sp")
 arrows(0,0, spe4.sc[,1], spe4.sc[,2], length=0, lty=1, col="red") 
 
 vif.cca(spebio.temp)
 
 spebio.rich <- rda(spe.hel, envbio, envrich)
 R2 <- RsquareAdj(spebio.rich)$r.squared
 
 # Adjusted R^2 retrieved from the rda object
 R2adj <- RsquareAdj(spebio.rich)$adj.r.squared
 print(anova_rda <- anova.cca(spebio.rich, step=1000))
 #print(anova_rda_axes <- anova.cca(spebio.rich, by="axis", step=1000))
 
  plot(spebio.rich, scaling=1, display=c("sp", "lc", "cn"), 
      main=paste("Partial RDA spe.hel ~bio | richness   Top =scaling1 -1c scores\n", OTU_table_name," with ", transformation," ~ ","\n -scaling 1 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe3.sc <- scores(spebio.rich, choices=1:2, scaling=1, display="sp")
 arrows(0,0, spe3.sc[,1], spe3.sc[,2], length=0, lty=1, col="red")
 
 ## Scaling 2
 plot(spebio.rich, scaling=2, display=c("sp", "lc", "cn"), 
      main=paste("Partial RDA spe.hel ~bio | richness   Top =scaling1 -1c scores\n", OTU_table_name," with ", transformation," ~ ","\n -scaling 2 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe4.sc <- scores(spebio.rich, choices=1:2, display="sp")
 arrows(0,0, spe4.sc[,1], spe4.sc[,2], length=0, lty=1, col="red") 
 
 vif.cca(spebio.rich)
 
 sperich.bio <- rda(spe.hel, envrich, envbio)
 R2 <- RsquareAdj(sperich.bio)$r.squared
 
 # Adjusted R^2 retrieved from the rda object
 R2adj <- RsquareAdj(sperich.bio)$adj.r.squared
 print(anova_rda <- anova.cca(sperich.bio, step=1000))
 #print(anova_rda_axes <- anova.cca(sperich.bio, by="axis", step=1000))
 
 plot(sperich.bio, scaling=1, display=c("sp", "lc", "cn"), 
      main=paste("Partial RDA spe.hel ~richness  | bio  Top =scaling1 -1c scores\n", OTU_table_name," with ", transformation," ~ ","\n -scaling 1 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe3.sc <- scores(sperich.bio, choices=1:2, scaling=1, display="sp")
 arrows(0,0, spe3.sc[,1], spe3.sc[,2], length=0, lty=1, col="red")
 
 ## Scaling 2
 plot(sperich.bio, scaling=2, display=c("sp", "lc", "cn"), 
      main=paste("Partial RDA spe.hel ~richness  | bio  Top =scaling1 -1c scores\n", OTU_table_name," with ", transformation," ~ ","\n -scaling 2 - 1c scores", " adjusted R2 ",R2adj,"\n Anova with 1000 steps p-value", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe4.sc <- scores(sperich.bio, choices=1:2, display="sp")
 arrows(0,0, spe4.sc[,1], spe4.sc[,2], length=0, lty=1, col="red") 
 
 vif.cca(sperich.bio)
 dev.off()
}

RDA_with_forward_selection <- function (OTU_table,OTU_table_name, transformation, env_data, env_data_name) {
 spe.hel <- decostand(OTU_table, transformation)
 ## does not allow factor variables
 spe.rda.all <- rda(spe.hel ~ .,data=env_data)
 #global adjusted R2
 (R2a.all <- RsquareAdj(spe.rda.all)$adj.r.squared)
 # Forward selection using packfor's forward.sel()
 ##how to get at these things:
 forward_select_var <- forward.sel(spe.hel, env_data, adjR2thresh=R2a.all)
 forward_select_var$variables # gives the forward selected variables

 spe.rda.forward <- rda(as.formula(paste("spe.hel  ~ ",paste(forward_select_var$variables, collapse="+"),sep = "")),data=env_data)
 print(spe.rda.forward)
 (R2a.forward <- RsquareAdj(spe.rda.forward)$adj.r.squared)
 print(anova_rda <- anova.cca(spe.rda.forward, step=1000))  
 print(anova_rda_axes <- anova.cca(spe.rda.forward, by="axis", step=1000))
 vif.cca(spe.rda.forward)
 
 pdf(paste("../figures/Triplot_RDA_",OTU_table_name,"_with_forward_selection", env_data_name,"_plots.pdf",sep=""), width = 15, height = 15, onefile = TRUE) 
 #Scaling 1
 plot(spe.rda.forward, scaling=1, display=c("sp", "lc", "cn"), main=paste("Triplot RDA \n", OTU_table_name," with ", transformation," ~ forward selection \n -scaling 1 - 1c scores", " adjusted R2 ",R2a.all,"\n Anova with 1000 steps p-value ", anova_rda$"Pr(>F)"[1],collapse="\n"))
 
 spe3.sc <- scores(spe.rda.forward, choices=1:2, scaling=1, display="sp")
 arrows(0,0, spe3.sc[,1], spe3.sc[,2], length=0, lty=1, col="red")
 
 ## Scaling 2
 plot(spe.rda.forward, scaling=2, display=c("sp", "lc", "cn"), main=paste("Triplot RDA \n", OTU_table_name," with ", transformation," ~ forward selection \n -scaling 2 - 1c scores", " adjusted R2 ",R2a.all,"\n Anova with 1000 steps p-value ", anova_rda$"Pr(>F)"[1],collapse="\n"))
 spe4.sc <- scores(spe.rda.forward, choices=1:2, display="sp")
 arrows(0,0, spe4.sc[,1], spe4.sc[,2], length=0, lty=1, col="red")
 dev.off()
 return(spe.rda.forward)
}

## overall significant
## individual axes not significant
#RDA_18s_with_forward_selection <- RDA_with_forward_selection(fam_sum_otus,"normalized-18s-OTUs-fam-sum-otus","hellinger", params_Jericho_18s[,-1], "params-Jericho-no-date")

## axes are not significant.... overall is significant.
Partial_RDA_18s <- partial_RDAs(fam_sum_otus,"normalized-18s-OTUs-fam-sum-otus","hellinger", "params-Jericho-no-date", envchem, envbio, envtemp, envrich)

#### Now try to variation partitioning ####


## variation partitionign with all explanatory variables
OTU_table <- fam_sum_otus
transformation <- "hellinger"
spe.hel <- decostand(OTU_table, transformation)
spe.part.bio_chem <- varpart(spe.hel, envchem, envbio)
spe.part.bio_chem
plot(spe.part.bio_chem, digits=2)

v <- venneuler(c(envchem=spe.part.bio_chem$part$indfract$Adj.R.squared[1], envbio=spe.part.bio_chem$part$indfract$Adj.R.squared[3], "envchem&envbio"=spe.part.bio_chem$part$indfract$Adj.R.squared[2]))

# vd <- venneuler(c(A=0.3, B=0.3, C=1.1, "A&B"=0.1, "A&C"=0.2, "B&C"=0.1 ,"A&B&C"=0.1))

v$labels<- c(
 paste("envchem\n",round(spe.part.bio_chem$part$fract$Adj.R.squared[1], digits=2)),
 paste("envbio\n",round(spe.part.bio_chem$part$fract$Adj.R.squared[2], digits=2))
)

pdf("../figures/Variation_partitioning_18s_normalized_fam_sum.pdf")
plot(v, main=paste0("Variation Partitioning of fam_sum_otus \n residuals:",
                    round(spe.part.bio_chem$part$indfract$Adj.R.squared[4], digits=2),
                    "   only chem: ",
                    round(spe.part.bio_chem$part$indfract$Adj.R.squared[1], digits=2),
                    "\n only bio: ",
                    round(spe.part.bio_chem$part$indfract$Adj.R.squared[3], digits=2),
                    "   overlap: ",
                    round(spe.part.bio_chem$part$indfract$Adj.R.squared[2], digits=2)))
dev.off()


## so it seems that the variables are intercorrelated so a good reason to do parsimony and combine variation partiionaing iwth forward selection.

# spe.part.all <- varpart(spe.hel, envchem, envbio,envtemp)
# spe.part.all
# plot(spe.part.all, digits=2)
# 
# spe.part.all <- varpart(spe.hel, envchem, envbio,envrich)
# spe.part.all
# plot(spe.part.all, digits=2)
# 
# spe.part.all <- varpart(spe.hel, envchem, envbio,envrich, envtemp)
# spe.part.all
# plot(spe.part.all, digits=2)

anova.cca(rda(spe.hel,envchem),step=1000) ## significant
anova.cca(rda(spe.hel, envbio),step=1000) ## significant
anova.cca(rda(spe.hel, envtemp),step=1000) ## 0.006 so < 0.005
anova.cca(rda(spe.hel, envrich),step=1000) ## 0.099 not significant



## Separate forward selection in each subset of environmental variables
spe.chem <- rda(spe.hel, envchem)
R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
forward.sel(spe.hel, envchem, adjR2thresh = R2a.all.chem, nperm=9999)

spe.bio <- rda(spe.hel, envbio)
R2a.all.bio <- RsquareAdj(spe.bio)$adj.r.squared
forward.sel(spe.hel, envbio, adjR2thresh=R2a.all.bio, nperm=9999)

spe.temp <- rda(spe.hel, envtemp)
R2a.all.temp <- RsquareAdj(spe.temp)$adj.r.squared
forward.sel(spe.hel, envtemp, adjR2thresh=R2a.all.temp, nperm=9999)

spe.rich <- rda(spe.hel, envrich)
R2a.all.rich <- RsquareAdj(spe.rich)$adj.r.squared
forward.sel(spe.hel, envrich, adjR2thresh=R2a.all.rich, nperm=9999)

#parsimonioius subsets of explanatory variables(based on forward selections)
envchem.pars <- envchem[,c(1,3,4)]
envbio.pars <- envbio[,c(1,3)] ## but this doesn't quite make sense because it was only 1. I will redo with only 1
envtemp.pars <- envtemp[,c(1)]
envrich.pars <- envrich[,c(2)]
## Variation partitioning
(spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envrich.pars))
plot(spe.part, digits=2)

(spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envrich.pars, envtemp.pars))
plot(spe.part, digits=2)


anova.cca(rda(spe.hel, envchem.pars),step=1000) ## significant
anova.cca(rda(spe.hel, envbio.pars),step=1000) ## significant
anova.cca(rda(spe.hel, envtemp.pars),step=1000) ## significant
anova.cca(rda(spe.hel, envrich.pars),step=1000) ## not significant

## using those that are significant
(spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envtemp.pars))
plot(spe.part, digits=2)


#### Ok try with removing VC 1252 ####

## remove VC 1252 from env param and from OTUs
fam_sum_otus_no_1252 <-fam_sum_otus[-25,]
params_Jericho_without_date_no_1252 <- params_Jericho_18s[-25,-1]

## overall significant
## axes not significant
#RDA_18s_summarised_by_fam_no_1252 <- RDA_of_OTU_table(fam_sum_otus_no_1252,"normalized-18s-OTUs-no-1252-fam-sum-otus","hellinger", params_Jericho_without_date_no_1252 , "params-Jericho-no-date")

## remove VC 1252 from env param and from OTUs
ord_sum_otus_no_1252 <-ord_sum_otus[-25,]

## overall significant
## axes not significant
#RDA_18s_summarised_by_fam_no_1252 <- RDA_of_OTU_table(ord_sum_otus_no_1252,"normalized-18s-OTUs-no-1252-ord-sum-otus","hellinger", params_Jericho_without_date_no_1252 , "params-Jericho-no-date")



## Interpretation of the triplot. This plot has an R2 adjusted of 0.3 and when testing the overall plot with anova we have p-value of 0.001.
## Looking at scaling 1
## The constrained variance is 0.609 (unadjusted. ), so to find the amount explained by first axis it is the R^2 adjusted multiplied by that axis. RDA1=0.3737 * 0.301 and then RDA2 0.2315 * 0.301 So we see 11% explained by the first axis and 7% explained by the second.
0.3727*0.301
0.2315*0.301
## Is this really low?
## It seems that PO4, Salinity, DO, day length, Average viral abundance, Chla, bacterial abundance, pH explain the first axis. While month number and silicate explain the 2nd axis and temperature seems to influence both.
## Looking at scaling 2
## Noctilucales seems correlated with high silicate. Dinophycae seems to be correlated with high temperature. 
##Diatomea might be dominating the 1235,1234 and then 1203 and 1225 sites. Are these bloom times? Would be useful to have dates on here. 
## Some of the species with shorter arrows might be seen all the time or are less strongly correlated. 
## 

## nothing very interesting. 
#Partial_RDA_18s <- partial_RDAs(fam_sum_otus_no_1252,"normalized-18s-OTUs-fam-sum-otus-no-1252","hellinger", "params-Jericho-no-date", envchem[-25,], envbio[-25,], envtemp[-25,], envrich[-25,])

#RDA_18s_with_forward_selection_no_1252 <- RDA_with_forward_selection(fam_sum_otus_no_1252,"normalized-18s-OTUs-no-1252-fam-sum-otus","hellinger", params_Jericho_without_date_no_1252 , "params-Jericho-no-date")


#### Now try to variation partitioning ####
## variation partitionign with all explanatory variables
spe.hel <- decostand(fam_sum_otus_no_1252, "hellinger")
spe.part.all <- varpart(spe.hel, envchem[-25,], envbio[-25,])
spe.part.all
plot(spe.part.all, digits=2)
str(spe.part.all)
str(spe.part.all$part)
spe.part.all$part$SS.Y
spe.part.all$part$fract$Adj.R.squared
spe.part.all$part$fract
spe.part.all$part$indfract$Adj.R.squared

# fractions [a+b]:
rda.chem <- rda (spe.hel ~ ., envchem[-25,])
# fractions [b+c]:
rda.bio <- rda (spe.hel ~ ., envbio[-25,])
anova(rda.chem)
anova(rda.bio)

v <- venneuler(c(envchem=spe.part.all$part$indfract$Adj.R.squared[1], envbio=spe.part.all$part$indfract$Adj.R.squared[3], "envchem&envbio"=spe.part.all$part$indfract$Adj.R.squared[2]))

# vd <- venneuler(c(A=0.3, B=0.3, C=1.1, "A&B"=0.1, "A&C"=0.2, "B&C"=0.1 ,"A&B&C"=0.1))
 
v$labels<- c(
 paste("envchem\n",round(spe.part.all$part$fract$Adj.R.squared[1], digits=2)),
 paste("envbio\n",round(spe.part.all$part$fract$Adj.R.squared[2], digits=2))
 )

# plot(vd)

pdf("../figures/Variation_partitioning_18s_normalized_fam_sum_no_1252.pdf")
plot(v, main=paste0("Variation Partitioning of fam_sum_otus_no_1252 \n residuals:",
                    round(spe.part.all$part$indfract$Adj.R.squared[4], digits=2),
                    "   only chem: ",
                    round(spe.part.all$part$indfract$Adj.R.squared[1], digits=2),
                    "\n only bio: ",
                    round(spe.part.all$part$indfract$Adj.R.squared[3], digits=2),
                    "   overlap: ",
                    round(spe.part.all$part$indfract$Adj.R.squared[2], digits=2)))
dev.off()
## so it seems that the variables are intercorrelated so a good reason to do parsimony and combine variation partiionaing iwth forward selection.


##### Now try with 16s ####
VC_dates <- Jericho_data$Date[match(rownames(normalized_16s_OTUs), original_jericho_data$VC_number)]
params_Jericho_16s <- subset(params_Jericho, Date %in% VC_dates )
amplicon_richness_for_study <- subset(amplicon_richness, Date %in% VC_dates )
## temporary fix
amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0
#Physiogrpahy
envbio <- params_Jericho_16s[,c("Average_viral_abundance","Average_bacterial_abundance", "Average_chl_a")]
names(envbio)

envchem <- params_Jericho_16s[,c("Average_PO4", "Average_SiO2",
                                # "Average_NO3_NO2", 
                                 "Temperature_YSI", "Salinity_ppt_YSI", "Dissolved_oxygen_percent" , "pH", "Tide_height")]
names(envchem)
envtemp <- params_Jericho_16s[,c("day_length", "month_number")]
envrich <- amplicon_richness_for_study[,c("richness.AVS","richness.MPL","richness.18S","richness.gp23")]

## overall not siginificant
## axes not significant
#RDA_16s <- RDA_of_OTU_table(normalized_16s_OTUs,"normalized-16s-OTUs","hellinger", params_Jericho_16s[,-1] , "params-Jericho-no-date")

## nothing very interesting
#Partial_RDA_16s <- partial_RDAs(normalized_16s_OTUs,"normalized-16s-OTUs","hellinger", "params-Jericho-no-date", envchem, envbio, envtemp, envrich)

## overall significant
# axes not significant
#RDA_16s_with_forward_selection <- RDA_with_forward_selection(normalized_16s_OTUs,"normalized-16s-OTUs","hellinger", params_Jericho_16s[,-1] , "params-Jericho-no-date")


### Try same approach by family with bacteria too ###
colnames(normalized_16s_OTUs) <- gsub(".size.*.", "", colnames(normalized_16s_OTUs))
melted_otus <- melt(as.matrix(normalized_16s_OTUs))
taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)
melted_otus$Family <- taxonomy_16s$Family[match(melted_otus$Var2, taxonomy_16s$otu_number)]
## remove unclassified
melted_otus <- subset(melted_otus, Family != "unclassified"  )

fam_melted <- melted_otus %>%
 group_by(Var1,Family) %>%
 summarise(total=sum(value))

## Ok need to put it back together again.
## Families on top
fam_sum_otus <- dcast(fam_melted, Var1 ~  Family)
## need to make them rownames again
row.names(fam_sum_otus) <- fam_sum_otus$Var1
fam_sum_otus <- fam_sum_otus[,-1]


## not significant overall
#RDA_16s_fam <- RDA_of_OTU_table(fam_sum_otus,"normalized-16s-OTUs-summarised-by-fam","hellinger", params_Jericho_16s[,-1] , "params-Jericho-no-date")

melted_otus <- melt(as.matrix(normalized_16s_OTUs))
melted_otus$Order <- taxonomy_16s$Order[match(melted_otus$Var2, taxonomy_16s$otu_number)]
## remove unclassified
melted_otus <- subset(melted_otus, Order != "unclassified"  )

org_melted <- melted_otus %>%
 group_by(Var1,Order) %>%
 summarise(total=sum(value))

## Ok need to put it back together again.
## Families on top
ord_sum_otus <- dcast(ord_melted, Var1 ~  Order)
## need to make them rownames again
row.names(ord_sum_otus) <- ord_sum_otus$Var1
ord_sum_otus <- ord_sum_otus[,-1]


#RDA_16s_ord <- RDA_of_OTU_table(ord_sum_otus,"normalized-16s-OTUs-summarised-by-ord","hellinger", params_Jericho_16s[,-1] , "params-Jericho-no-date")

## nothing interesting of significant
#Partial_RDA_16s <- partial_RDAs(fam_sum_otus,"normalized-16s-OTUs-summarised-by-fam","hellinger", "params-Jericho-no-date", envchem, envbio, envtemp, envrich)

## no variables could be selected with such a shitty model. 
#RDA_16s_with_forward_selection <- RDA_with_forward_selection(fam_sum_otus,"normalized-16s-OTUs-summarised-by-fam","hellinger", params_Jericho_16s[,-1] , "params-Jericho")



##### Now try with MPL ####
VC_dates <- Jericho_data$Date[match(rownames(normalized_MPL_OTUs), original_jericho_data$VC_number)]
params_Jericho_MPL <- subset(params_Jericho, Date %in% VC_dates )
amplicon_richness_for_study <- subset(amplicon_richness, Date %in% VC_dates )
## temporary fix
amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0
#Physiogrpahy
envbio <- params_Jericho_MPL[,c("Average_viral_abundance","Average_bacterial_abundance", "Average_chl_a")]
names(envbio)

envchem <- params_Jericho_MPL[,c("Average_PO4", "Average_SiO2",
                                # "Average_NO3_NO2", 
                                 "Temperature_YSI", "Salinity_ppt_YSI", "Dissolved_oxygen_percent" , "pH", "Tide_height")]
names(envchem)
envtemp <- params_Jericho_MPL[,c("day_length", "month_number")]
envrich <- amplicon_richness_for_study[,c("richness.AVS","richness.16S","richness.18S","richness.gp23")]

## overall not significant
#RDA_MPL <- RDA_of_OTU_table(normalized_MPL_OTUs,"normalized-MPL-OTUs","hellinger", params_Jericho_MPL[,-1] , "params-Jericho-no-date")

## nothing significant or interesting
#Partial_RDA_MPL <- partial_RDAs(normalized_MPL_OTUs,"normalized-MPL-OTUs","hellinger", "params-Jericho-no-date", envchem, envbio, envtemp, envrich)

##nothing selected
#RDA_MPL_with_forward_selection <- try(RDA_with_forward_selection(normalized_MPL_OTUs,"normalized-MPL-OTUs","hellinger", params_Jericho_MPL[,-1] , "params-Jericho-no-date"))

# ## try with cca
# test_hel <- decostand(normalized_MPL_OTUs, "hellinger")
# test_cca <- cca(test_hel ~ .,data=params_Jericho_MPL[,-1]) 
# test_cca
# summary(test_cca)
# print(anova_cca <- anova.cca(test_cca, step=1000)) # not significant
# print(anova_rda <- anova.cca(test_cca, by="axis", step=1000)) ## first 2 significant below 0.005, 5 significant axes
# 
# melted_mpl <- melt(as.matrix(normalized_MPL_OTUs))
# top20_MPL <- find_top20_OTUs(melted_mpl)
# top20_normalized_MPL <- subset(normalized_MPL_OTUs, select= colnames(normalized_MPL_OTUs) %in% top20_MPL$Var2)
# 
# test_hel <- decostand(top20_normalized_MPL, "hellinger")
# test_cca <- cca(test_hel ~ .,data=params_Jericho_MPL[,-1]) 
# test_cca
# summary(test_cca)
# print(anova_cca <- anova.cca(test_cca, step=1000)) # not significant
# #print(anova_rda <- anova.cca(test_cca, by="axis", step=1000)) 


##### Now try with gp23 ####
VC_dates <- Jericho_data$Date[match(rownames(normalized_gp23_OTUs), original_jericho_data$VC_number)]
params_Jericho_gp23 <- subset(params_Jericho, Date %in% VC_dates )
amplicon_richness_for_study <- subset(amplicon_richness, Date %in% VC_dates )
## temporary fix
amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0
#Physiogrpahy
envbio <- params_Jericho_gp23[,c("Average_viral_abundance","Average_bacterial_abundance", "Average_chl_a")]
names(envbio)

envchem <- params_Jericho_gp23[,c("Average_PO4", "Average_SiO2",
                                 # "Average_NO3_NO2",
                                  "Temperature_YSI", "Salinity_ppt_YSI", "Dissolved_oxygen_percent" , "pH", "Tide_height")]
names(envchem)
envtemp <- params_Jericho_gp23[,c("day_length", "month_number")]
envrich <- amplicon_richness_for_study[,c("richness.AVS","richness.16S","richness.18S","richness.MPL")]

## not significant overall
#RDA_gp23 <- RDA_of_OTU_table(normalized_gp23_OTUs,"normalized-gp23-OTUs","hellinger", params_Jericho_gp23[,-1] , "params-Jericho-no-date")

## neither significant nor interesting.
#Partial_RDA_gp23 <- partial_RDAs(normalized_gp23_OTUs,"normalized-gp23-OTUs","hellinger", "params-Jericho-no-date", envchem, envbio, envtemp, envrich)

## only picked silicate.
# RDA_gp23_with_forward_selection <- RDA_with_forward_selection(normalized_gp23_OTUs,"normalized-gp23-OTUs","hellinger", params_Jericho_gp23[,-1] , "params-Jericho-no-date")

## try with cca
# test_hel <- decostand(normalized_gp23_OTUs, "hellinger")
# test_cca <- cca(test_hel ~ .,data=params_Jericho_gp23[,-1]) 
# test_cca
# summary(test_cca)
# print(anova_cca <- anova.cca(test_cca, step=1000)) # not significant
# print(anova_rda <- anova.cca(test_cca, by="axis", step=1000)) ## first 2 significant below 0.005, 5 significant axes
# 
# melted_gp23 <- melt(as.matrix(normalized_gp23_OTUs))
# top20_gp23 <- find_top20_OTUs(melted_gp23)
# top20_normalized_gp23 <- subset(normalized_gp23_OTUs, select= colnames(normalized_gp23_OTUs) %in% top20_gp23$Var2)
# 
# test_hel <- decostand(top20_normalized_gp23, "hellinger")
# test_cca <- cca(test_hel ~ .,data=params_Jericho_gp23[,-1]) 
# test_cca
# summary(test_cca)
# print(anova_cca <- anova.cca(test_cca, step=1000)) # not significant
# #print(anova_rda <- anova.cca(test_cca, by="axis", step=1000)) 


##### Now try with AVS ####
VC_dates <- Jericho_data$Date[match(rownames(normalized_AVS_OTUs), original_jericho_data$VC_number)]
params_Jericho_AVS <- subset(params_Jericho, Date %in% VC_dates )
amplicon_richness_for_study <- subset(amplicon_richness, Date %in% VC_dates )
## temporary fix
amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0
#Physiogrpahy
envbio <- params_Jericho_AVS[,c("Average_viral_abundance","Average_bacterial_abundance", "Average_chl_a")]
names(envbio)

envchem <- params_Jericho_AVS[,c("Average_PO4", "Average_SiO2",
                                 #"Average_NO3_NO2", 
                                 "Temperature_YSI", "Salinity_ppt_YSI", "Dissolved_oxygen_percent" , "pH", "Tide_height")]
names(envchem)
envtemp <- params_Jericho_AVS[,c("day_length", "month_number")]
envrich <- amplicon_richness_for_study[,c("richness.AVS","richness.16S","richness.18S","richness.MPL")]

## not significant overall
#RDA_AVS <- RDA_of_OTU_table(normalized_AVS_OTUs,"normalized-AVS-OTUs","hellinger", params_Jericho_AVS[,-1] , "params-Jericho-no-date")

## not significant overall
#Partial_RDA_AVS <- partial_RDAs(normalized_AVS_OTUs,"normalized-AVS-OTUs","hellinger", "params-Jericho-no-date", envchem, envbio, envtemp, envrich)

## nothing selected
#RDA_AVS_with_forward_selection <- try(RDA_with_forward_selection(normalized_AVS_OTUs,"normalized-AVS-OTUs","hellinger", params_Jericho_AVS[,-1] , "params-Jericho-no-date"))

# test_hel <- decostand(normalized_AVS_OTUs, "hellinger")
# test_cca <- cca(test_hel ~ .,data=params_Jericho_AVS[,-1]) 
# test_cca
# summary(test_cca)
# print(anova_cca <- anova.cca(test_cca, step=1000)) # not significant
# print(anova_rda <- anova.cca(test_cca, by="axis", step=1000)) ## first 2 significant below 0.005, 5 significant axes
# 
# melted_AVS <- melt(as.matrix(normalized_AVS_OTUs))
# top20_AVS <- find_top20_OTUs(melted_AVS)
# top20_normalized_AVS <- subset(normalized_AVS_OTUs, select= colnames(normalized_AVS_OTUs) %in% top20_AVS$Var2)
# 
# test_hel <- decostand(top20_normalized_AVS, "hellinger")
# test_cca <- cca(test_hel ~ .,data=params_Jericho_AVS[,-1]) 
# test_cca
# summary(test_cca)
# print(anova_cca <- anova.cca(test_cca, step=1000)) # not significant
# #print(anova_rda <- anova.cca(test_cca, by="axis", step=1000))

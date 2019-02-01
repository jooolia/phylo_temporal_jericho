
library(vegan)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
#library(packfor)
library(geosphere)
library(ggvegan)
library(grid)
library(gridExtra)
library(cowplot)
library(venneuler)
library(VennDiagram)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
## tidy up some OTU sizes in the the names
colnames(normalized_18s_OTUs) <- gsub(".size.*.", "", colnames(normalized_18s_OTUs))

#normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
#normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv",                                  row.names=1)
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv", row.names="VC_number")

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")

normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 
## tidy up some OTU sizes in the the names
colnames(normalized_16s_OTUs) <- gsub(".size.*.", "", colnames(normalized_16s_OTUs))

original_jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
original_jericho_data$Date <- as.Date(original_jericho_data$Date)

## missing data filled in by mean imputation
Jericho_data <- read.csv("../results/Jericho_env_data_mean_imputed.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

## amplicon richness 
#amplicon_richness <- read.csv("../results/amplicon_richness_by_date.csv")
#amplicon_richness$Date <- as.Date(amplicon_richness$Date)

## exclude season and NO2
parameters_to_exclude <- c("season",
                           #"Average_NO3_NO2",
                           "Standard_error_NO3_NO2",
                           "Secchi_disk_disappears",
                           "Secchi_disk_reappears",
                           #"Dissolved_oxygen_percent",
                           "VC_number",
                           "Tide_height",
                           #"month_number",
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
                           #"Average_viral_abundance",
                           "Standard_error_viral_abundance",
                           "Standard_error_bacterial_abundance",
                           "Standard_error_PO4",
                           "Standard_error_SiO2",
                           "Standard_error_chl_a")

params_Jericho <-  subset(Jericho_data, select= !(colnames(Jericho_data) %in% parameters_to_exclude))
head(params_Jericho)


library(psych)
corr.test(params_Jericho[,-1])$r



### First look at 18s
VC_dates <- Jericho_data$Date[match(rownames(normalized_18s_OTUs), original_jericho_data$VC_number)]

params_Jericho_18s <- subset(params_Jericho, Date %in% VC_dates)

#amplicon_richness_for_study <- subset(amplicon_richness, Date %in% VC_dates )

#amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0

## split up the params for later use: 

envbio <- params_Jericho_18s[,c("Average_viral_abundance",
                                "Average_bacterial_abundance", 
                                "Average_chl_a")]
envchem <- params_Jericho_18s[,c("Average_PO4",
                                 "Average_SiO2",
                                 "Average_NO3_NO2",
                                 "Temperature_YSI", 
                                 "Salinity_ppt_YSI",
                                 "pH")]
envtemp <- params_Jericho_18s[,c("month_number",
                                 "day_length")]

taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)

hell_transform_and_rda <- function (otu_table, env_table) {
  
  hel_trans <- decostand(otu_table, "hellinger")

  print(spe.rda <- rda(hel_trans ~ ., env_table))
  
  print(R2 <- RsquareAdj(spe.rda)$r.squared)
  
  # Adjusted R^2 retrieved from the rda object
  print(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)

  return(spe.rda)
}

check_sig_of_rda <- function (population_rda) {
 try(anova_rda <- anova.cca(population_rda, step=1000))
 print(anova_rda)
}

#Plot of the RDA scaling 1
plot_scaling_1 <- function (population_rda) {
 
 # Total variation of the species dataset (sum of all eigenvalues)
 population_rda$tot.chi
 anova_rda <- anova.cca(population_rda, step=1000)
 
 # Variation explained by each canonical axis (%) sing adjusted R2
 adjust_var_axe <- 100 * population_rda$CCA$eig / RsquareAdj(population_rda)$adj.r.squared
 
 R2 <- RsquareAdj(population_rda)$r.squared
 # Adjusted R^2 retrieved from the rda object
 R2adj <- RsquareAdj(population_rda)$adj.r.squared
 
 scrs <- fortify(population_rda, scaling = 1)
 
 ## take only site scores for this
 sites <- with(scrs, scrs[Score == "sites", ])
 species <- with(scrs, scrs[Score == "species", ])
 biplot <- with(scrs, scrs[Score == "biplot",])
 
 
 sites$Label <- Jericho_data$Date[match(sites$Label, as.character(Jericho_data$VC_number))]

 overall_plot <-   ggplot(species, 
                          aes(x = Dim1,
                              y = Dim2,
                              label=Label))+
  coord_fixed(xlim = c(-1.5,1.5),
                ylim = c(-1.5,1.5))+
  xlab(paste0( "RDA1 (percent explained ",
               round(adjust_var_axe, digits=2), ")"))+
  ylab(paste0( "RDA2 (percent explained ",
               round(adjust_var_axe[2], digits=2), ")"))
 
 ## add in arrows
 plt <- overall_plot+
 geom_segment(aes(x=0, y=0, xend=Dim1, yend=Dim2),
              colour="purple", 
              arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(data=sites,
            aes(x = Dim1, y = Dim2, label=Label),
            hjust=1.2,
            size=3)+
  geom_point(data=sites,
             aes(x = Dim1, y = Dim2, label=Label),
             colour="blue")+
  annotate("text", x=-1, y=-1,
           label=paste0("Adjusted r2 =",
                        round(R2adj,digits=2), 
                        "\nOverall significance=",
                        anova_rda$"Pr(>F)"[1]),
           hjust = 0 ,
           size=4)

 plt2 <-  overall_plot+ 
  geom_point(data=sites,
             aes(x = Dim1, y = Dim2, label=Label), 
             colour="blue")+
  geom_segment(data=biplot,
               aes(x=0, y=0, xend=Dim1, yend=Dim2), 
               colour="green",
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(data=biplot, 
            aes(x=Dim1, y=Dim2), 
            colour="black", size=3)

 plt3 <-  overall_plot + 
  geom_text(hjust=1.2, 
            colour="purple",
            size=3)+
  geom_segment(data=biplot,
               aes(x=0, y=0, xend=Dim1, yend=Dim2),
               colour="green",
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(data=biplot, 
            aes(x=Dim1, y=Dim2), 
            colour="black", 
            size=3 )

 
 plt4 <-  plt3 +
  geom_point(data=sites,
             aes(x = Dim1, y = Dim2, label=Label), 
             colour="blue")
 

 sites$Season <- Jericho_data$season[match(sites$Label, Jericho_data$Date)]

 plt4 <-  plt3 +
  geom_point(data=sites,
             aes(x = Dim1, y = Dim2, label=Label, size=2.1))+
  geom_point(data=sites,
             aes(x = Dim1, y = Dim2, label=Label, colour=Season))
 
#  p <- ggplot(data=data.scores,aes(x=Axis.1,y=Axis.2))+ 
#   geom_point(aes_string(colour=environmental_parameter),
#              size=10) + 

 plot_grid(plt, plt2, plt3, plt4, ncol=2, labels=c("A","B", "C", "D"))
}

#Scaling 2 with the same data
plot_scaling_2 <- function (population_rda) {
 population_rda$tot.chi
 anova_rda <- anova.cca(population_rda, step=1000)
 
 # Variation explained by each canonical axis (%) sing adjusted R2
 adjust_var_axe <- 100 * population_rda$CCA$eig / RsquareAdj(population_rda)$adj.r.squared
 
 R2 <- RsquareAdj(population_rda)$r.squared
 # Adjusted R^2 retrieved from the rda object
 R2adj <- RsquareAdj(population_rda)$adj.r.squared
 
scrs <- fortify(population_rda, scaling = 2)

sites <- with(scrs, scrs[Score == "sites", ])
species <- with(scrs, scrs[Score == "species", ])
biplot <- with(scrs, scrs[Score == "biplot",])
## add in something to plot for the arrows

sites$Label <- Jericho_data$Date[match(sites$Label, as.character(Jericho_data$VC_number))]
sites$Season <- Jericho_data$season[match(sites$Label, Jericho_data$Date)]

overall_plot <-   ggplot(species, 
                         aes(x = Dim1,
                             y = Dim2,
                             label=Label))+
 coord_fixed(xlim = c(-1.5,1.5),
             ylim = c(-1.5,1.5))+
 xlab(paste0( "RDA1 (percent explained ",
              round(adjust_var_axe, digits=2), ")"))+
 ylab(paste0( "RDA2 (percent explained ",
              round(adjust_var_axe[2], digits=2), ")"))
## add in arrows
 plt <- overall_plot+
 geom_segment(aes(x=0, y=0, xend=Dim1, yend=Dim2),
              colour="purple", 
              arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(data=sites,
            aes(x = Dim1, y = Dim2, label=Label),
            hjust=1.2,
            size=1)+
  geom_point(data=sites,
             aes(x = Dim1, y = Dim2, label=Label),
             colour="blue")+
  annotate("text", x=-1, y=-1,
           label=paste0("Adjusted r2 =",
                        round(R2adj,digits=2), 
                        "\nOverall significance=",
                        anova_rda$"Pr(>F)"[1]),
           hjust = 0 ,
           size=4)

 plt2 <-  overall_plot+ 
  geom_point(data=sites,
             aes(x = Dim1, y = Dim2, label=Label), 
             colour="blue")+
  geom_segment(data=biplot,
               aes(x=0, y=0, xend=Dim1, yend=Dim2), 
               colour="green",
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(data=biplot, 
            aes(x=Dim1, y=Dim2), 
            colour="black", size=3)

 plt3 <-  overall_plot + 
  geom_text(hjust=1.2, 
            colour="purple",
            size=3)+
  geom_segment(data=biplot,
               aes(x=0, y=0, xend=Dim1, yend=Dim2),
               colour="green",
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(data=biplot, 
            aes(x=Dim1, y=Dim2), 
            colour="black", 
            size=3 )

 plt4 <-  plt3 +
  geom_point(data=sites,
             aes(x = Dim1, y = Dim2, label=Label, size=2.1))+
  geom_point(data=sites,
             aes(x = Dim1, y = Dim2, label=Label, colour=Season))

 plot_grid(plt, plt2, plt3, plt4, ncol=2, labels=c("A","B", "C", "D"))
 return(list(plt, plt2, plt3, plt4))
}


taxonomy_18s_no_met_fun <- subset(taxonomy_18s,!(Phylum %in% "Opisthokonta"))

normalized_18s_OTUs_no_met_fun <- subset(normalized_18s_OTUs, select=(colnames(normalized_18s_OTUs) %in% taxonomy_18s_no_met_fun$otu_number ))

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

write.csv(ord_sum_otus, file="../results/normalized_18s_summarized_by_order.csv")

melted_otus <- melt(as.matrix(normalized_18s_OTUs_no_met_fun))
melted_otus$Family <- taxonomy_18s$Family[match(melted_otus$Var2, taxonomy_18s$otu_number)]
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

write.csv(fam_sum_otus, file="../results/normalized_18s_summarized_by_family.csv")



# 
# normalized_18s_rda_ord_sum <- hell_transform_and_rda(ord_sum_otus  ,params_Jericho_18s[,-1] ) 
# vif.cca(normalized_18s_rda_ord_sum)
# ## try to decrease that. Take temp or salinity
# env_selected <- params_Jericho_18s[,c("Average_PO4",
#                                  "Average_SiO2",
#                                  "Temperature_YSI", 
#                                  "Average_bacterial_abundance",
#                                  "Average_chl_a",
#                                  "day_length")]
# 
# normalized_18s_rda_ord_sum <- hell_transform_and_rda(ord_sum_otus  ,env_selected ) 
# vif.cca(normalized_18s_rda_ord_sum)
# 
# 
# check_sig_of_rda(normalized_18s_rda_ord_sum) 
# ## check the significance of the axes.
# anova.cca(rda(decostand(ord_sum_otus, "hellinger") ~.,params_Jericho_18s[,-1]),by="axis", step=1000)
# 
# pdf("../figures/RDA_18s_ord_sum%03d.pdf", onefile=FALSE, width=20, height = 20)
# S18_scaling1 <- plot_scaling_1(normalized_18s_rda_ord_sum)
# S18_scaling1 
# S18_scaling2 <- plot_scaling_2(normalized_18s_rda_ord_sum)
# S18_scaling2
# dev.off()
# 
# ## so see that the variables are too correlated, so will do forward selection of the global variables. 
# RDA_with_forward_selection <- function (OTU_table,
#                                         env_data){
#  spe.hel <- decostand(OTU_table, "hellinger")
#  spe.rda.all <- rda(spe.hel ~ .,data=env_data)
#  #global adjusted R2
#  (R2a.all <- RsquareAdj(spe.rda.all)$adj.r.squared)
#  # Forward selection using packfor's forward.sel()
#  try(forward_select_var <- forward.sel(spe.hel, env_data, adjR2thresh=R2a.all))
#  
#  try(spe.rda.forward <- rda(as.formula(paste("spe.hel  ~ ",paste(forward_select_var$variables, collapse="+"),sep = "")),
#                             data=env_data))
#  print(spe.rda.forward)
#  try(R2a.forward <- RsquareAdj(spe.rda.forward)$adj.r.squared)
#  try(return(spe.rda.forward))
# }
# 
# 
# rda_forward_sel <- RDA_with_forward_selection(ord_sum_otus, params_Jericho_18s[,-1])
# check_sig_of_rda(rda_forward_sel) 
# 
# vif.cca(rda_forward_sel)
# save(rda_forward_sel, file ="../results/normalized_18s_rda_ord_sum_for_sel.rda")
# 
# 
# ## ------------------------------------------------------------------------
# pdf("../figures/RDA_18s_ord_sum_with_forward_selection%03d.pdf", onefile=FALSE, width=20, height = 20)
# S18_for_sel_scaling1 <- plot_scaling_1(rda_forward_sel)
# S18_for_sel_scaling1
# S18_for_sel_scaling2 <- plot_scaling_2(rda_forward_sel)
# S18_for_sel_scaling2
# dev.off()
# 
# partial_RDA <- function(OTU_table, table1, table2){
#  spe.hel <- decostand(OTU_table, "hellinger")
#  rda_partial <- rda(spe.hel, table1, table2)
#  print(rda_partial)
#  R2 <- RsquareAdj(rda_partial)$r.squared
#  # Adjusted R^2 retrieved from the rda object
#  R2adj <- RsquareAdj(rda_partial)$adj.r.squared
#  return(rda_partial)
# } 
# 
# ## if hold biology constact
# partial_ord_sum_18s <- partial_RDA(ord_sum_otus, envchem, envbio)
# plot_scaling_1(partial_ord_sum_18s)
# 
# ## if hold chemistry 
# partial_ord_sum_18s <- partial_RDA(ord_sum_otus, envbio,envchem)
# plot_scaling_1(partial_ord_sum_18s)
# ## not significant...
# 
# OTU_table <- ord_sum_otus
# transformation <- "hellinger"
# spe.hel <- decostand(OTU_table, transformation)
# spe.part.bio_chem <- varpart(spe.hel, envchem, envbio)
# spe.part.bio_chem
# plot(spe.part.bio_chem, digits=2)
# 
# 
# ## ------------------------------------------------------------------------
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
# 
# anova.cca(rda(spe.hel,envchem),step=1000) ## significant
# anova.cca(rda(spe.hel, envbio),step=1000) ## significant
# anova.cca(rda(spe.hel, envtemp),step=1000) ## 0.006 so > 0.005
# anova.cca(rda(spe.hel, envrich),step=1000) ## 0.099 not significant
# 
# 
# pdf("../figures/Variation_partitioning_18s_ord_sum_chem_bio_time_richness.pdf")
# plot(spe.part.all, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()
# 
# ## ------------------------------------------------------------------------
# spe.chem <- rda(spe.hel, envchem)
# R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
# forward.sel(spe.hel, envchem, adjR2thresh = R2a.all.chem, nperm=9999)
# 
# spe.bio <- rda(spe.hel, envbio)
# R2a.all.bio <- RsquareAdj(spe.bio)$adj.r.squared
# forward.sel(spe.hel, envbio, adjR2thresh=R2a.all.bio, nperm=9999)
# 
# spe.temp <- rda(spe.hel, envtemp)
# R2a.all.temp <- RsquareAdj(spe.temp)$adj.r.squared
# forward.sel(spe.hel, envtemp, adjR2thresh=R2a.all.temp, nperm=9999)
# 
# spe.rich <- rda(spe.hel, envrich)
# R2a.all.rich <- RsquareAdj(spe.rich)$adj.r.squared
# forward.sel(spe.hel, envrich, adjR2thresh=R2a.all.rich, nperm=9999)
# 
# #parsimonioius subsets of explanatory variables(based on forward selections)
# envchem.pars <- envchem[,c("Average_PO4","Temperature_YSI","Salinity_ppt_YSI")]
# envbio.pars <- envbio[,c("Average_chl_a","Average_bacterial_abundance")] ## but this doesn't quite make sense because it was only 1. I will redo with only 1
# envtemp.pars <- envtemp[,c("day_length")]
# envrich.pars <- envrich[,c("richness.gp23","richness.MPL")]
# ## Variation partitioning
# 
# anova.cca(rda(spe.hel, envchem.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envbio.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envtemp.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envrich.pars),step=1000) ## not significant
# 
# ## using those that are significant
# (spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envtemp.pars))
# plot(spe.part, digits=2)

#(spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envtemp.pars, envrich.pars))

# pdf("../figures/Variation_partitioning_18s_ord_sum_chem_bio_time_richness_forward_selection.pdf")
# plot(spe.part, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()


# spe.part.bio_chem <- varpart(spe.hel, envchem, envbio)
# spe.part.bio_chem
# plot(spe.part.bio_chem, digits=2)
# 
# res.venn <- spe.part$part$indfract[-8,3]
# tot.venn <- spe.part$part$fract[-8,3]
#names(res.venn) <- c("envchem","envbio","envtemp","envchem&envbio","envbio&envtemp","envchem&envtemp","envchem&envbio&envtemp")

# v <- venneuler(res.venn)
# plot(v)
# 
# v$labels<- c(
#  paste("envchem\n",round(spe.part$part$fract$Adj.R.square[1], digits=2)),
#  paste("envbio\n",round(spe.part$part$fract$Adj.R.square[2], digits=2)),
#  paste("envtemp\n", round(spe.part$part$indfract$Adj.R.square[3], digits=2))
# )
# 
# plot(v)


# draw.triple.venn(65, 75, 85,
#                  35, 15, 25, 5, c("First", "Second", "Third"))
# 
# draw.triple.venn(20, 40, 60, 2, 2, 2, 1,
#                  c("First", "Second", "Third"), sep.dist = 0.1, rotation.degree = 30)
# 
# vennDiagram_18s <- draw.triple.venn(area1 = round(tot.venn[1]*100),
#                  area2 = round(tot.venn[2]*100),
#                  area3 = round(tot.venn[3]*100),
#                  n12 = round((res.venn[4]+res.venn[7]) * 100),
#                  n23 = round((res.venn[5]+res.venn[7]) *100),
#                  n13 = round((res.venn[6]+res.venn[7]) * 100), 
#                  n123 = round((res.venn[7]) *100),
#                  category = c("envchem","envbio", "envtemp" ),
#                  fill = c("skyblue", "pink1", "mediumorchid"), euler.d = TRUE,
#                  scaled = TRUE)
# 
# # plot(v, main=paste0("Variation Partitioning of fam_sum_otus \n residuals:",
# #                     round(spe.part$part$indfract$Adj.R.square[4], digits=2),
# #                     "   only chem: ",
# #                     round(spe.part$part$indfract$Adj.R.square[1], digits=2),
# #                     "\n only bio: ",
# #                     round(spe.part$part$indfract$Adj.R.square[3], digits=2),
# #                     "   overlap: ",
# #                     round(spe.part$part$indfract$Adj.R.square[2], digits=2)))
# 
# 
# 
# ## Do all this on 18s without the VC 1252 ####
# ## ------------------------------------------------------------------------
# 
# ## remove VC 1252 from env param and from OTUs
# ord_sum_otus_no_1252 <-ord_sum_otus[!(rownames(ord_sum_otus) %in% "1252"),]
# params_Jericho_without_date_no_1252 <- params_Jericho_18s[!(rownames(ord_sum_otus) %in% "1252"),-1]
# 
# envbio <- params_Jericho_without_date_no_1252[,c("Average_viral_abundance",
#                                 "Average_bacterial_abundance", 
#                                 "Average_chl_a")]
# envchem <- params_Jericho_without_date_no_1252[,c("Average_PO4",
#                                  "Average_SiO2",
#                                  "Average_NO3_NO2",
#                                  "Temperature_YSI", 
#                                  "Salinity_ppt_YSI",
#                                  "pH")]
# envtemp <- params_Jericho_without_date_no_1252[,c("month_number",
#                                  "day_length")]
# # envrich <- amplicon_richness_for_study[!(amplicon_richness_for_study$Date %in% as.Date("2011-07-05")),c("richness.MPL",
# #                                           "richness.16S",
# #                                           "richness.gp23")]
# 
# normalized_18s_rda_ord_sum_no_1252 <- hell_transform_and_rda(ord_sum_otus_no_1252 ,params_Jericho_without_date_no_1252 ) 
# 
# 
# check_sig_of_rda(normalized_18s_rda_ord_sum_no_1252) 
# anova.cca(rda(decostand(ord_sum_otus_no_1252, "hellinger") ~.,params_Jericho_without_date_no_1252),by="axis", step=1000)
# 
# plot_scaling_1(normalized_18s_rda_ord_sum_no_1252)
# plot_scaling_2(normalized_18s_rda_ord_sum_no_1252)
# 
# vif.cca(normalized_18s_rda_ord_sum_no_1252) ## vifs are too high!!
# check_sig_of_rda(normalized_18s_rda_ord_sum_no_1252) 
# ## check the significance of the axes.
# anova.cca(rda(decostand(ord_sum_otus, "hellinger") ~.,params_Jericho_18s[,-1]),by="axis", step=1000)
# 
# pdf("../figures/RDA_18s_ord_sum_no_1252%03d.pdf", onefile=FALSE, width=20, height = 20)
# plot_scaling_1(normalized_18s_rda_ord_sum_no_1252)
# plot_scaling_2(normalized_18s_rda_ord_sum_no_1252)
# dev.off()
# 
# rda_forward_sel <- RDA_with_forward_selection(ord_sum_otus_no_1252, params_Jericho_without_date_no_1252)
# check_sig_of_rda(rda_forward_sel) 
# 
# vif.cca(rda_forward_sel)
# 
# pdf("../figures/RDA_18s_ord_sum_no_1252_with_forward_selection%03d.pdf", onefile=FALSE, width=20, height = 20)
# plot_scaling_1(rda_forward_sel)
# plot_scaling_2(rda_forward_sel)
# dev.off()
# 
# ## if hold biology constact
# partial_ord_sum_18s <- partial_RDA(ord_sum_otus_no_1252, envchem, envbio)
# plot_scaling_1(partial_ord_sum_18s)
# 
# ## if hold chemistry 
# partial_ord_sum_18s <- partial_RDA(ord_sum_otus_no_1252, envbio,envchem)
# plot_scaling_1(partial_ord_sum_18s)
# ## not significant...
# 
# OTU_table <- ord_sum_otus_no_1252
# transformation <- "hellinger"
# spe.hel <- decostand(OTU_table, transformation)
# spe.part.bio_chem <- varpart(spe.hel, envchem, envbio)
# spe.part.bio_chem
# plot(spe.part.bio_chem, digits=2)
# 
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
# 
# anova.cca(rda(spe.hel,envchem),step=1000) ## significant
# anova.cca(rda(spe.hel, envbio),step=1000) ## not significant
# anova.cca(rda(spe.hel, envtemp),step=1000) ## 0.006 so > 0.005
# anova.cca(rda(spe.hel, envrich),step=1000) ## 0.099 not significant
# 
# pdf("../figures/Variation_partitioning_18s_no_1252_ord_sum_chem_bio_time_richness.pdf")
# spe.part.all <- varpart(spe.hel, envchem,envtemp,envrich)
# plot(spe.part.all, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()
# 
# spe.chem <- rda(spe.hel, envchem)
# R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
# forward.sel(spe.hel, envchem, adjR2thresh = R2a.all.chem, nperm=9999)
# 
# spe.bio <- rda(spe.hel, envbio)
# R2a.all.bio <- RsquareAdj(spe.bio)$adj.r.squared
# forward.sel(spe.hel, envbio, adjR2thresh=R2a.all.bio, nperm=9999)
# 
# spe.temp <- rda(spe.hel, envtemp)
# R2a.all.temp <- RsquareAdj(spe.temp)$adj.r.squared
# forward.sel(spe.hel, envtemp, adjR2thresh=R2a.all.temp, nperm=9999)
# 
# spe.rich <- rda(spe.hel, envrich)
# R2a.all.rich <- RsquareAdj(spe.rich)$adj.r.squared
# forward.sel(spe.hel, envrich, adjR2thresh=R2a.all.rich, nperm=9999)
# 
# #parsimonioius subsets of explanatory variables(based on forward selections)
# envchem.pars <- envchem[,c("Average_PO4","Temperature_YSI")]
# envbio.pars <- envbio[,c("Average_chl_a")] ## but this doesn't quite make sense because it was only 1. I will redo with only 1
# envtemp.pars <- envtemp[,c("day_length")]
# envrich.pars <- envrich[,c("richness.gp23","richness.MPL")]
# ## Variation partitioning
# 
# anova.cca(rda(spe.hel, envchem.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envbio.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envtemp.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envrich.pars),step=1000) ## not significant
# 
# ## using those that are significant
# (spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envtemp.pars))
# plot(spe.part, digits=2)
# 
# (spe.part <- varpart(spe.hel, envchem.pars,envtemp.pars, envrich.pars))
# 
# pdf("../figures/Variation_partitioning_18s_no_1252_ord_sum_chem_bio_time_richness_forward_selection.pdf")
# plot(spe.part, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()
# 
# 
# 
# ## ------------------------------------------------------------------------
# ## RDA with phytoplankton
# 
# phytoplankton_fam <- c("Diatomea",
#                        "Dinophyceae",
#                        "Mamiellophyceae",
#                        "Chrysophyceae",
#                        "Chlorophyceae",
#                        "Euglenozoa",
#                        "Dictyochophyceae",
#                        "Chlorarachniophyta",
#                        "Raphidophyceae")
# phytoplankton_class <- c("Rhodophyceae",
#                          "Chloroplastida")
# phytoplantkon_phylum <- c("Haptophyta",
#                           "Cryptophyceae")
# phytoplankton_order <- c("Dinoflagellata",
#                          "Ochrophyta") ## not sure about this, are all dinos phytoplankton?
# 
# ## algae but maybe not phytoplankton,
# algae_class <- c("Chloroplastida")
# 
# ## masts we don't know...
# 
# phytoplankton_tax <- taxonomy_18s %>% 
#  filter(Family %in% phytoplankton_fam | Class %in% phytoplankton_class | Phylum %in% phytoplantkon_phylum| Order %in% phytoplankton_order)
# phytoplankton_otus <- droplevels(subset(ord_sum_otus, select= colnames(ord_sum_otus) %in% phytoplankton_tax$Order))
# 
# 
# normalized_18s_rda_ord_sum_phyto <- hell_transform_and_rda(phytoplankton_otus ,params_Jericho_18s[,-1] ) 
# 
# check_sig_of_rda(normalized_18s_rda_ord_sum_phyto) 
# ## overall significant.
# ## check the significance of the axes.
# anova.cca(rda(decostand(phytoplankton_otus, "hellinger") ~.,params_Jericho_18s[,-1]),by="axis", step=1000)
# vif.cca(normalized_18s_rda_ord_sum_phyto)
# 
# pdf("../figures/RDA_18s_ord_sum_phytos%03d.pdf", onefile=FALSE, width=20, height = 20)
# plot_scaling_1(normalized_18s_rda_ord_sum_phyto)
# plot_scaling_2(normalized_18s_rda_ord_sum_phyto)
# dev.off()
# 
# 
# OTU_table <- phytoplankton_otus
# transformation <- "hellinger"
# spe.hel <- decostand(OTU_table, transformation)
# 
# envbio <- params_Jericho_18s[,c("Average_viral_abundance",
#                                 "Average_bacterial_abundance", 
#                                 "Average_chl_a")]
# envchem <- params_Jericho_18s[,c("Average_PO4",
#                                  "Average_SiO2",
#                                  "Average_NO3_NO2",
#                                  "Temperature_YSI", 
#                                  "Salinity_ppt_YSI",
#                                  "pH")]
# envtemp <- params_Jericho_18s[,c("month_number",
#                                  "day_length")]
# 
# spe.part.bio_chem <- varpart(spe.hel, envchem, envbio)
# spe.part.bio_chem
# plot(spe.part.bio_chem, digits=2)
# 
# 
# rda_forward_sel <- RDA_with_forward_selection(phytoplankton_otus, params_Jericho_18s[,-1])
# check_sig_of_rda(rda_forward_sel) 
# 
# vif.cca(rda_forward_sel)
# 
# ## ------------------------------------------------------------------------
# pdf("../figures/RDA_18s_ord_sum_phytos_with_forward_selection%03d.pdf", onefile=FALSE, width=20, height = 20)
# plot_scaling_1(rda_forward_sel)
# plot_scaling_2(rda_forward_sel)
# dev.off()
# 
# ## if hold biology constact
# partial_ord_sum_18s <- partial_RDA(phytoplankton_otus, envchem, envbio)
# plot_scaling_1(partial_ord_sum_18s)
# 
# ## if hold chemistry 
# partial_ord_sum_18s <- partial_RDA(phytoplankton_otus, envbio,envchem)
# plot_scaling_1(partial_ord_sum_18s)
# ## not significant...
# 
# OTU_table <- phytoplankton_otus
# transformation <- "hellinger"
# spe.hel <- decostand(OTU_table, transformation)
# spe.part.bio_chem <- varpart(spe.hel, envchem, envbio)
# spe.part.bio_chem
# plot(spe.part.bio_chem, digits=2)
# 
# 
# ## ------------------------------------------------------------------------
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
# 
# anova.cca(rda(spe.hel,envchem),step=1000) ## significant
# anova.cca(rda(spe.hel, envbio),step=1000) ## not significant
# anova.cca(rda(spe.hel, envtemp),step=1000) ## 0.006 so > 0.005
# anova.cca(rda(spe.hel, envrich),step=1000) ## 0.099 not significant
# 
# pdf("../figures/Variation_partitioning_18s_phytoplanktok_ord_sum_chem_bio_time.pdf")
# spe.part.all <- varpart(spe.hel,envchem,envbio, envtemp)
# plot(spe.part.all, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()
# 
# ## ------------------------------------------------------------------------
# spe.chem <- rda(spe.hel, envchem)
# R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
# forward.sel(spe.hel, envchem, adjR2thresh = R2a.all.chem, nperm=9999)
# 
# spe.bio <- rda(spe.hel, envbio)
# R2a.all.bio <- RsquareAdj(spe.bio)$adj.r.squared
# forward.sel(spe.hel, envbio, adjR2thresh=R2a.all.bio, nperm=9999)
# 
# spe.temp <- rda(spe.hel, envtemp)
# R2a.all.temp <- RsquareAdj(spe.temp)$adj.r.squared
# forward.sel(spe.hel, envtemp, adjR2thresh=R2a.all.temp, nperm=9999)
# 
# spe.rich <- rda(spe.hel, envrich)
# R2a.all.rich <- RsquareAdj(spe.rich)$adj.r.squared
# forward.sel(spe.hel, envrich, adjR2thresh=R2a.all.rich, nperm=9999)
# 
# #parsimonioius subsets of explanatory variables(based on forward selections)
# envchem.pars <- envchem[,c("Average_PO4","Temperature_YSI")]
# envbio.pars <- envbio[,c("Average_chl_a","Average_viral_abundance")] ## but this doesn't quite make sense because it was only 1. I will redo with only 1
# envtemp.pars <- envtemp[,c("day_length")]
# envrich.pars <- envrich[,c("richness.gp23")]
# ## Variation partitioning
# 
# anova.cca(rda(spe.hel, envchem.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envbio.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envtemp.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envrich.pars),step=1000) ## not significant
# 
# ## using those that are significant
# (spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envtemp.pars, envrich.pars))
# plot(spe.part, digits=2)
# 
# (spe.part <- varpart(spe.hel, envchem.pars,envbio.pars, envtemp.pars, envrich.pars))
# 
# pdf("../figures/Variation_partitioning_18s_phytos_ord_sum_chem_bio_time_richness_forward_selection.pdf")
# plot(spe.part, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()
# 
# 

## RDA on heterotrophs 

## or coudl try approach that heterotrophs are those that are not phytoplankton
# heterotroph_tax <- taxonomy_18s %>% 
#  filter(!(otu_number %in% phytoplankton_tax$otu_number))
# heterotroph_otus <- droplevels(subset(ord_sum_otus, select = colnames(ord_sum_otus) %in% heterotroph_tax$Order)) 
# 
# normalized_18s_rda_ord_sum_het <- hell_transform_and_rda(heterotroph_otus ,params_Jericho_18s[,-1] ) 
# check_sig_of_rda(normalized_18s_rda_ord_sum_het) 
# ## overall significant.
# ## check the significance of the axes.
# anova.cca(rda(decostand(heterotroph_otus, "hellinger") ~.,params_Jericho_18s[,-1]),by="axis", step=1000)
# vif.cca(normalized_18s_rda_ord_sum_het)
# 
# pdf("../figures/RDA_18s_ord_sum_hets%03d.pdf", onefile=FALSE, width=20, height = 20)
# plot_scaling_1(normalized_18s_rda_ord_sum_het)
# plot_scaling_2(normalized_18s_rda_ord_sum_het)
# dev.off()
# 
# envbio <- params_Jericho_18s[,c("Average_viral_abundance",
#                                 "Average_bacterial_abundance", 
#                                 "Average_chl_a")]
# envchem <- params_Jericho_18s[,c("Average_PO4",
#                                  "Average_SiO2",
#                                  "Average_NO3_NO2",
#                                  "Temperature_YSI", 
#                                  "Salinity_ppt_YSI",
#                                  "pH")]
# envtemp <- params_Jericho_18s[,c("month_number",
#                                  "day_length")]
# 
# rda_forward_sel <- RDA_with_forward_selection(heterotroph_otus, params_Jericho_18s[,-1])
# check_sig_of_rda(rda_forward_sel) 
# 
# vif.cca(rda_forward_sel)
# 
# pdf("../figures/RDA_18s_ord_sum_hets_with_forward_selection%03d.pdf", onefile=FALSE, width=20, height = 20)
# plot_scaling_1(rda_forward_sel)
# plot_scaling_2(rda_forward_sel)
# dev.off()
# 
# ## if hold biology constact
# partial_ord_sum_18s <- partial_RDA(heterotroph_otus, envchem, envbio)
# plot_scaling_1(partial_ord_sum_18s)
# 
# ## if hold chemistry 
# partial_ord_sum_18s <- partial_RDA(heterotroph_otus, envbio,envchem)
# plot_scaling_1(partial_ord_sum_18s)
# ## not significant...
# 
# OTU_table <- heterotroph_otus
# transformation <- "hellinger"
# spe.hel <- decostand(OTU_table, transformation)
# spe.part.bio_chem <- varpart(spe.hel, envchem, envbio)
# spe.part.bio_chem
# plot(spe.part.bio_chem, digits=2)
# 
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
# 
# anova.cca(rda(spe.hel,envchem),step=1000) 
# anova.cca(rda(spe.hel, envbio),step=1000) ## significant
# anova.cca(rda(spe.hel, envtemp),step=1000)
# anova.cca(rda(spe.hel, envrich),step=1000) 
# 
# pdf("../figures/Variation_partitioning_18s_heterotrophs_ord_chem_bio.pdf")
# spe.part.all <- varpart(spe.hel,envchem,envbio)
# plot(spe.part.all, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()
# 
# spe.chem <- rda(spe.hel, envchem)
# R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
# forward.sel(spe.hel, envchem, adjR2thresh = R2a.all.chem, nperm=9999)
# 
# spe.bio <- rda(spe.hel, envbio)
# R2a.all.bio <- RsquareAdj(spe.bio)$adj.r.squared
# forward.sel(spe.hel, envbio, adjR2thresh=R2a.all.bio, nperm=9999)
# 
# spe.temp <- rda(spe.hel, envtemp)
# R2a.all.temp <- RsquareAdj(spe.temp)$adj.r.squared
# forward.sel(spe.hel, envtemp, adjR2thresh=R2a.all.temp, nperm=9999)
# 
# spe.rich <- rda(spe.hel, envrich)
# R2a.all.rich <- RsquareAdj(spe.rich)$adj.r.squared
# try(forward.sel(spe.hel, envrich, adjR2thresh=R2a.all.rich, nperm=9999))
# #try(forward_select_var <- forward.sel(spe.hel, env_data, adjR2thresh=R2a.all))
# 
# 
# #parsimonioius subsets of explanatory variables(based on forward selections)
# envchem.pars <- envchem[,c("Average_PO4","Temperature_YSI")]
# envbio.pars <- envbio[,c("Average_chl_a","Average_bacterial_abundance")] ## but this doesn't quite make sense because it was only 1. I will redo with only 1
# envtemp.pars <- envtemp[,c("day_length")]
# envrich.pars <- envrich[,c("richness.gp23")]
# ## Variation partitioning
# 
# anova.cca(rda(spe.hel, envchem.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envbio.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envtemp.pars),step=1000) ## significant
# #anova.cca(rda(spe.hel, envrich.pars),step=1000) ## not significant
# 
# ## using those that are significant
# (spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envtemp.pars))
# plot(spe.part, digits=2)
# 
# (spe.part <- varpart(spe.hel, envchem.pars,envbio.pars, envtemp.pars))
# 
# pdf("../figures/Variation_partitioning_18s_hets_ord_sum_chem_bio_time_richness_forward_selection.pdf")
# plot(spe.part, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()
# 



## ------------------------------------------------------------------------
## Now do with 16s

VC_dates <- Jericho_data$Date[match(rownames(normalized_16s_OTUs), original_jericho_data$VC_number)]
params_Jericho_16s <- subset(params_Jericho, Date %in% VC_dates )
# amplicon_richness_for_study <- subset(amplicon_richness, Date %in% VC_dates )
# ## temporary fix
# amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0


## ------------------------------------------------------------------------
#Physiogrpahy
# envbio <- params_Jericho_16s[,c( "Average_viral_abundance",
#                                  "Average_bacterial_abundance",
#                                 "Average_chl_a")]
# envchem <- params_Jericho_16s[,c("Average_PO4",
#                                  "Average_SiO2",
#                                  "Average_NO3_NO2",
#                                  "Temperature_YSI", 
#                                  "Salinity_ppt_YSI",
#                                  "pH")]
# envtemp <- params_Jericho_16s[,c("month_number", "day_length")]
# envrich <- amplicon_richness_for_study[,c("richness.MPL",
#                                           "richness.18S",
#                                           "richness.gp23")]

# ## ------------------------------------------------------------------------
# melted_otus <- melt(as.matrix(normalized_16s_OTUs))
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

write.csv(fam_sum_otus, file="../results/normalized_16s_summarized_by_family.csv")

## ------------------------------------------------------------------------
melted_otus <- melt(as.matrix(normalized_16s_OTUs))
melted_otus$Var2 <- gsub(".size.*.", "", melted_otus$Var2)
melted_otus$Order <- taxonomy_16s$Order[match(melted_otus$Var2, taxonomy_16s$otu_number)]## remove unclassified
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

write.csv(ord_sum_otus, file="../results/normalized_16s_summarized_by_order.csv")


# ## ------------------------------------------------------------------------
melted_otus <- melt(as.matrix(normalized_16s_OTUs))
melted_otus$Var2 <- gsub(".size.*.", "", melted_otus$Var2)
melted_otus$Class <- taxonomy_16s$Class[match(melted_otus$Var2, taxonomy_16s$otu_number)]## remove unclassified
melted_otus <- subset(melted_otus, Class != "unclassified"  )

class_melted <- melted_otus %>%
 group_by(Var1,Class) %>%
 summarise(total=sum(value))
# 
# ## Ok need to put it back together again.
# ## Families on top
class_sum_otus <- dcast(class_melted, Var1 ~  Class)
## need to make them rownames again
row.names(class_sum_otus) <- class_sum_otus$Var1
class_sum_otus <- class_sum_otus[,-1]

write.csv(class_sum_otus, file="../results/normalized_16s_summarized_by_class.csv")


# ## ------------------------------------------------------------------------
# normalized_16s_rda_ord_sum <- hell_transform_and_rda(ord_sum_otus ,params_Jericho_16s[,-1] ) 
# 
# ## ------------------------------------------------------------------------
# check_sig_of_rda(normalized_16s_rda_ord_sum) 
# anova.cca(rda(decostand(ord_sum_otus, "hellinger") ~.,params_Jericho_16s[,-1] ),by="axis", step=1000)
# vif.cca(normalized_16s_rda_ord_sum)
# save(normalized_16s_rda_ord_sum, file ="../results/normalized_16s_rda_ord_sum.rda")
# 
# 
# 
# env_selected <- params_Jericho_16s[,c("Average_PO4",
#                                       "Average_SiO2",
#                                       "Temperature_YSI", 
#                                       "Average_bacterial_abundance",
#                                       "Average_chl_a",
#                                       "day_length")]
# 
# normalized_16s_rda_ord_sum <- hell_transform_and_rda(ord_sum_otus ,env_selected) 
# 
# ## ------------------------------------------------------------------------
# check_sig_of_rda(normalized_16s_rda_ord_sum) 
# anova.cca(rda(decostand(ord_sum_otus, "hellinger") ~.,params_Jericho_16s[,-1] ),by="axis", step=1000)
# vif.cca(normalized_16s_rda_ord_sum)
# 
# ## ------------------------------------------------------------------------
# pdf("../figures/RDA_normalized_16s_rda_ord_sum%03d.pdf", onefile=FALSE)
# S16_scaling1 <- plot_scaling_1(normalized_16s_rda_ord_sum)
# S16_scaling1
# S16_scaling2 <- plot_scaling_2(normalized_16s_rda_ord_sum)
# S16_scaling2
# dev.off()
# 
# ## ------------------------------------------------------------------------
#  rda_forward_sel <- RDA_with_forward_selection(ord_sum_otus, params_Jericho_16s[,-1])
#  check_sig_of_rda(rda_forward_sel) 
#  vif.cca(rda_forward_sel)
#  save(rda_forward_sel, file ="../results/normalized_16s_rda_ord_sum_for_sel.rda")
#  
#  pdf("../figures/RDA_16s_ord_sum_with_forward_selection%03d.pdf", onefile=FALSE, width=20, height = 20)
#  S16_for_sel_scaling1 <- plot_scaling_1(rda_forward_sel)
#  S16_for_sel_scaling1
#  
#  S16_for_sel_scaling2 <- plot_scaling_2(rda_forward_sel)
#  S16_for_sel_scaling2
#  dev.off()
#  
#  
#  ## if hold biology constact
#  partial_ord_sum_16s <- partial_RDA(ord_sum_otus, envchem, envbio)
#  plot_scaling_1(partial_ord_sum_16s)
#  
#  ## if hold chemistry 
#  partial_ord_sum_16s <- partial_RDA(ord_sum_otus, envbio,envchem)
#  plot_scaling_1(partial_ord_sum_16s)
#  ##  this one is significant.....if hold chemistry constant....interesting diffs form 18s. 
#  
#  OTU_table <- ord_sum_otus
#  transformation <- "hellinger"
#  spe.hel <- decostand(OTU_table, transformation)
#  spe.part.bio_chem <- varpart(spe.hel, envchem, envbio)
#  spe.part.bio_chem
#  plot(spe.part.bio_chem, digits=2)
#  
#  spe.part.all <- varpart(spe.hel, envchem, envbio,envtemp)
#  spe.part.all
#  plot(spe.part.all, digits=2)
#  
#  spe.part.all <- varpart(spe.hel, envchem, envbio,envrich)
#  spe.part.all
#  plot(spe.part.all, digits=2)
#  
#  spe.part.all <- varpart(spe.hel, envchem, envbio,envrich, envtemp)
#  spe.part.all
#  plot(spe.part.all, digits=2)
#  
#  anova.cca(rda(spe.hel,envchem),step=1000) 
#  anova.cca(rda(spe.hel, envbio),step=1000) ## significant
#  anova.cca(rda(spe.hel, envtemp),step=1000)
#  anova.cca(rda(spe.hel, envrich),step=1000) 
#  
#  pdf("../figures/Variation_partitioning_16s_ord_chem_bio.pdf")
#  spe.part.all <- varpart(spe.hel,envchem,envbio)
#  plot(spe.part.all, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
#  dev.off()
#  
#  spe.chem <- rda(spe.hel, envchem)
#  R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
#  forward.sel(spe.hel, envchem, adjR2thresh = R2a.all.chem, nperm=9999)
#  
#  spe.bio <- rda(spe.hel, envbio)
#  R2a.all.bio <- RsquareAdj(spe.bio)$adj.r.squared
#  forward.sel(spe.hel, envbio, adjR2thresh=R2a.all.bio, nperm=9999)
#  
#  spe.temp <- rda(spe.hel, envtemp)
#  R2a.all.temp <- RsquareAdj(spe.temp)$adj.r.squared
#  forward.sel(spe.hel, envtemp, adjR2thresh=R2a.all.temp, nperm=9999)
#  
#  spe.rich <- rda(spe.hel, envrich)
#  R2a.all.rich <- RsquareAdj(spe.rich)$adj.r.squared
#  forward.sel(spe.hel, envrich, adjR2thresh=R2a.all.rich, nperm=9999)
#  
#  #parsimonioius subsets of explanatory variables(based on forward selections)
#  envchem.pars <- envchem[,c("Salinity_ppt_YSI","Temperature_YSI")]
#  envbio.pars <- envbio[,c("Average_chl_a","Average_bacterial_abundance")] ## but this doesn't quite make sense because it was only 1. I will redo with only 1
#  envtemp.pars <- envtemp[,c("month_number")]
#  envrich.pars <- envrich[,c("richness.18S")]
#  ## Variation partitioning
#  
#  anova.cca(rda(spe.hel, envchem.pars),step=1000) ## significant
#  anova.cca(rda(spe.hel, envbio.pars),step=1000) ## significant
#  anova.cca(rda(spe.hel, envtemp.pars),step=1000) ## significant
#  anova.cca(rda(spe.hel, envrich.pars),step=1000) ## not significant
#  
#  ## using those that are significant
#  (spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envtemp.pars))
# plot(spe.part, digits=2,Xnames=c("Chem", "Bio", "Temporal"))
# 
# 
# 
#  (spe.part <- varpart(spe.hel, envchem.pars,envbio.pars, envtemp.pars))
#  
#  pdf("../figures/Variation_partitioning_16s_ord_sum_chem_bio_time_richness_forward_selection.pdf")
#  plot(spe.part, digits=2, bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
#  dev.off()
#  
#  res.venn <- spe.part$part$indfract[-8,3]
#  tot.venn <- spe.part$part$fract[-8,3]
#  vennDiagram_16s <- draw.triple.venn(area1 = round(tot.venn[1]*100),
#                                      area2 = round(tot.venn[2]*100),
#                                      area3 = round(tot.venn[3]*100),
#                                      n12 = 0.1,
#                                      n23 = 0.1,
#                                      n13 = 0.1 , 
#                                      n123 = round((res.venn[7]) *100),
#                                      category = c("envchem","envbio", "envtemp" ),
#                                      fill = c("skyblue", "pink1", "mediumorchid"), euler.d = TRUE,
#                                      scaled = TRUE)
 
 

 ## test out top 50 OTUs
 find_top50_OTUs <- function (long_otus) {
  summarised_by_OTU <- ddply(long_otus, 
                             ~Var2,
                             summarise,
                             relative_abundance=sum(value))
  OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
  top_50 <- droplevels(head(OTU_sorted_by_rel_abun, n=50))
  return(top_50)
 } 
 
 find_top20_OTUs <- function (long_otus) {
  summarised_by_OTU <- ddply(long_otus, 
                             ~Var2,
                             summarise,
                             relative_abundance=sum(value))
  OTU_sorted_by_rel_abun <- summarised_by_OTU[order(-summarised_by_OTU$relative_abundance),]
  top_20 <- droplevels(head(OTU_sorted_by_rel_abun, n=20))
  return(top_20)
 }
 

## ------------------------------------------------------------------------
## MPL now.. 
VC_dates <- Jericho_data$Date[match(rownames(normalized_MPL_OTUs), original_jericho_data$VC_number)]
params_Jericho_MPL <- subset(params_Jericho, Date %in% VC_dates )
# amplicon_richness_for_study <- subset(amplicon_richness, Date %in% VC_dates )
# ## temporary fix
# amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0
#Physiogrpahy
# envbio <- params_Jericho_MPL[,c( "Average_viral_abundance",
#                                  "Average_bacterial_abundance",
#                                  "Average_chl_a")]
# envchem <- params_Jericho_MPL[,c("Average_PO4",
#                                  "Average_SiO2",
#                                  "Average_NO3_NO2",
#                                  "Temperature_YSI", 
#                                  "Salinity_ppt_YSI",
#                                  "pH")]
# envtemp <- params_Jericho_MPL[,c("month_number", "day_length")]
# envrich <- amplicon_richness_for_study[,c("richness.16S",
#                                           "richness.18S",
#                                           "richness.gp23")]

## ------------------------------------------------------------------------
normalized_MPL_rda <- hell_transform_and_rda(normalized_MPL_OTUs,
                                             params_Jericho_MPL[,-1] ) 
check_sig_of_rda(normalized_MPL_rda) 
## nope


## Try top 20
melted_otus <- melt(as.matrix(normalized_MPL_OTUs))
top20_otus_mpl <- find_top20_OTUs(melted_otus)

top20_normalized_MPL <- subset(normalized_MPL_OTUs,
                                          select=colnames(normalized_MPL_OTUs) %in% top20_otus_mpl$Var2)

normalized_MPL_top_20_rda <- hell_transform_and_rda(top20_normalized_MPL,
                                             params_Jericho_MPL[,-1] ) 
check_sig_of_rda(normalized_MPL_top_20_rda) 
## noope

## ------------------------------------------------------------------------
MPL_group_sums <- read.csv("../results/MPL_group_sums_by_site_prop.csv", row.names = 1)
rownames(MPL_group_sums) <- MPL_group_sums[,"VC"]
MPL_group_sums <- MPL_group_sums[,!(colnames(MPL_group_sums) %in% "VC")] ## remove VC column because now rownames

group_sums_MPL_rda <- hell_transform_and_rda(MPL_group_sums,
                                             params_Jericho_MPL[,-1] ) 
check_sig_of_rda(group_sums_MPL_rda) 
anova.cca(rda(decostand(MPL_group_sums, "hellinger") ~.,params_Jericho_MPL[,-1] ),by="axis", step=1000)
vif.cca(group_sums_MPL_rda)



env_selected <- params_Jericho_MPL[,c("Average_PO4",
                                      "Average_SiO2",
                                      "Temperature_YSI", 
                                      "Average_bacterial_abundance",
                                      "Average_chl_a",
                                      "day_length")]
group_sums_MPL_rda <- hell_transform_and_rda(MPL_group_sums,
                                             env_selected ) 
check_sig_of_rda(group_sums_MPL_rda) 
anova.cca(rda(decostand(MPL_group_sums, "hellinger") ~.,params_Jericho_MPL[,-1] ),by="axis", step=1000)
vif.cca(group_sums_MPL_rda)


## ------------------------------------------------------------------------
pdf("../figures/RDA_normalized_MPL_group_sums_rda_%03d.pdf", onefile=FALSE)
MPL_scaling1 <- plot_scaling_1(group_sums_MPL_rda)
MPL_scaling1
MPL_scaling2 <- plot_scaling_2(group_sums_MPL_rda)
MPL_scaling2
dev.off()

# ## ------------------------------------------------------------------------
# rda_forward_sel <- RDA_with_forward_selection(MPL_group_sums, params_Jericho_MPL[,-1])
# check_sig_of_rda(rda_forward_sel) 
# vif.cca(rda_forward_sel)
# save(rda_forward_sel, file ="../results/normalized_MPL_rda_for_sel.rda")
# 
# pdf("../figures/RDA_MPL_group_sum_with_forward_selection%03d.pdf", onefile=FALSE, width=20, height = 20)
# MPL_for_sel_scaling1 <- plot_scaling_1(rda_forward_sel)
# MPL_for_sel_scaling1
# MPL_for_sel_scaling2 <- plot_scaling_2(rda_forward_sel)
# MPL_for_sel_scaling2
# dev.off()
# 
# 
# ## if hold biology constact
# partial_group_sum_MPL <- partial_RDA(MPL_group_sums, envchem, envbio)
# plot_scaling_1(partial_group_sum_MPL)
# 
# ## if hold chemistry 
# partial_group_sum_MPL <- partial_RDA(MPL_group_sums, envbio,envchem)
# plot_scaling_1(partial_group_sum_MPL)
# 
# 
# OTU_table <- MPL_group_sums
# transformation <- "hellinger"
# spe.hel <- decostand(OTU_table, transformation)
# spe.part.bio_chem <- varpart(spe.hel, envchem, envbio)
# spe.part.bio_chem
# plot(spe.part.bio_chem, digits=2,bg=c("skyblue", "pink1", "mediumorchid"))
# 
# 
# ## ------------------------------------------------------------------------
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
# 
# anova.cca(rda(spe.hel,envchem),step=1000) 
# anova.cca(rda(spe.hel, envbio),step=1000) ## significant
# anova.cca(rda(spe.hel, envtemp),step=1000)
# anova.cca(rda(spe.hel, envrich),step=1000) 

# pdf("../figures/Variation_partitioning_MPL_group_sum_chem_bio.pdf")
# spe.part.all <- varpart(spe.hel,envchem,envbio)
# plot(spe.part.all, digits=2,bg=c("skyblue", "pink1", "mediumorchid"), Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()

## ------------------------------------------------------------------------
# spe.chem <- rda(spe.hel, envchem)
# R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
# forward.sel(spe.hel, envchem, adjR2thresh = R2a.all.chem, nperm=9999)
# 
# spe.bio <- rda(spe.hel, envbio)
# R2a.all.bio <- RsquareAdj(spe.bio)$adj.r.squared
# forward.sel(spe.hel, envbio, adjR2thresh=R2a.all.bio, nperm=9999)
# 
# spe.temp <- rda(spe.hel, envtemp)
# R2a.all.temp <- RsquareAdj(spe.temp)$adj.r.squared
# forward.sel(spe.hel, envtemp, adjR2thresh=R2a.all.temp, nperm=9999)
# 
# spe.rich <- rda(spe.hel, envrich)
# R2a.all.rich <- RsquareAdj(spe.rich)$adj.r.squared
# forward.sel(spe.hel, envrich, adjR2thresh=R2a.all.rich, nperm=9999)

#parsimonioius subsets of explanatory variables(based on forward selections)
# envchem.pars <- envchem[,c("Salinity_ppt_YSI","Temperature_YSI","Average_NO3_NO2", "Average_PO4")]
# envbio.pars <- envbio[,c("Average_chl_a","Average_bacterial_abundance")] ## but this doesn't quite make sense because it was only 1. I will redo with only 1
# envtemp.pars <- envtemp[,c("day_length")]
# envrich.pars <- envrich[,c("richness.18S")]
# ## Variation partitioning
# 
# anova.cca(rda(spe.hel, envchem.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envbio.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envtemp.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envrich.pars),step=1000) ## not significant
# 
# ## using those that are significant
# (spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envtemp.pars))
# plot(spe.part, digits=2)


# pdf("../figures/Variation_partitioning_MPL_group_sum_chem_bio_time_forward_selection.pdf")
# plot(spe.part, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()

# res.venn <- spe.part$part$indfract[-8,3]
# tot.venn <- spe.part$part$fract[-8,3]
# vennDiagram_MPL <- draw.triple.venn(area1 = round(tot.venn[1]*100),
#                                     area2 = round(tot.venn[2]*100),
#                                     area3 = round(tot.venn[3]*100),
#                                     n12 = round((res.venn[4]+res.venn[7]) * 100),
#                                     n23 = round((res.venn[5]+res.venn[7]) *100),
#                                     n13 = round((res.venn[6]+res.venn[7]) * 100), 
#                                     n123 = round((res.venn[7]) *100),
#                                     category = c("envchem","envbio", "envtemp" ),
#                                     fill = c("skyblue", "pink1", "mediumorchid"), euler.d = TRUE,
#                                     scaled = TRUE)


## ------------------------------------------------------------------------
## gp23 time ####


VC_dates <- Jericho_data$Date[match(rownames(normalized_gp23_OTUs), original_jericho_data$VC_number)]
params_Jericho_gp23 <- subset(params_Jericho, Date %in% VC_dates )
# amplicon_richness_for_study <- subset(amplicon_richness, Date %in% VC_dates )
# ## temporary fix
# amplicon_richness_for_study[is.na(amplicon_richness_for_study)] <- 0

## Physiogrpahy
# envbio <- params_Jericho_gp23[,c( "Average_viral_abundance",
#                                  "Average_bacterial_abundance",
#                                  "Average_chl_a")]
# envchem <- params_Jericho_gp23[,c("Average_PO4",
#                                  "Average_SiO2",
#                                  "Average_NO3_NO2",
#                                  "Temperature_YSI", 
#                                  "Salinity_ppt_YSI",
#                                  "pH")]
# envtemp <- params_Jericho_gp23[,c("month_number", "day_length")]
# envrich <- amplicon_richness_for_study[,c("richness.MPL",
#                                           "richness.18S",
#                                           "richness.16S")]


normalized_gp23_rda <- hell_transform_and_rda(normalized_gp23_OTUs,
                                              params_Jericho_gp23[,-1] ) 
check_sig_of_rda(normalized_gp23_rda) 
## nope


## Try top 20
melted_otus <- melt(as.matrix(normalized_gp23_OTUs))
top20_otus_gp23 <- find_top20_OTUs(melted_otus)

top20_normalized_gp23 <- subset(normalized_gp23_OTUs,
                               select=colnames(normalized_gp23_OTUs) %in% top20_otus_gp23$Var2)

normalized_gp23_top_20_rda <- hell_transform_and_rda(top20_normalized_gp23,
                                                    params_Jericho_gp23[,-1] ) 
check_sig_of_rda(normalized_gp23_top_20_rda) 
## noope


gp23_group_sums <- read.csv("../results/gp23_group_sums_by_site_prop.csv",row.names = 1)
rownames(gp23_group_sums) <- gp23_group_sums[,"VC"]
gp23_group_sums <- gp23_group_sums[,!(colnames(gp23_group_sums) %in% "VC")] ## remove VC column because now rownames

group_sums_gp23_rda <- hell_transform_and_rda(gp23_group_sums,
                                              params_Jericho_gp23[,-1] ) 

check_sig_of_rda(group_sums_gp23_rda) 
anova.cca(rda(decostand(gp23_group_sums, "hellinger") ~.,params_Jericho_gp23[,-1] ),by="axis", step=1000)
vif.cca(group_sums_gp23_rda)


env_selected <- params_Jericho_gp23[,c("Average_PO4",
                                      "Average_SiO2",
                                      "Temperature_YSI", 
                                      "Average_bacterial_abundance",
                                      "Average_chl_a",
                                      "day_length")]

group_sums_gp23_rda <- hell_transform_and_rda(gp23_group_sums,
                                              env_selected  ) 

check_sig_of_rda(group_sums_gp23_rda) 
anova.cca(rda(decostand(gp23_group_sums, "hellinger") ~.,params_Jericho_gp23[,-1] ),by="axis", step=1000)
vif.cca(group_sums_gp23_rda)


# pdf("../figures/RDA_normalized_gp23_group_sums_rda_%03d.pdf", onefile=FALSE)
# gp23_scaling1 <- plot_scaling_1(group_sums_gp23_rda)
# gp23_scaling1
# gp23_scaling2 <- plot_scaling_2(group_sums_gp23_rda)
# gp23_scaling2 
# dev.off()

## With forward selection ####
# rda_forward_sel <- RDA_with_forward_selection(gp23_group_sums, params_Jericho_gp23[,-1])
# check_sig_of_rda(rda_forward_sel) 
# vif.cca(rda_forward_sel)
# save(rda_forward_sel, file ="../results/normalized_gp23_rda_for_sel.rda")

# pdf("../figures/RDA_gp23_group_sum_with_forward_selection%03d.pdf", onefile=FALSE, width=20, height = 20)
# gp23_for_sel_scaling1 <- plot_scaling_1(rda_forward_sel)
# gp23_for_sel_scaling1
# 
# gp23_for_sel_scaling2 <- plot_scaling_2(rda_forward_sel)
# gp23_for_sel_scaling2
# dev.off()

## Partial RDA ####
## if hold biology constact
# partial_group_sum_gp23 <- partial_RDA(gp23_group_sums, envchem, envbio)
# plot_scaling_1(partial_group_sum_gp23)
# 
# ## if hold chemistry 
# partial_group_sum_gp23 <- partial_RDA(gp23_group_sums, envbio,envchem)
# plot_scaling_1(partial_group_sum_gp23)
# # nope
# 
# 
# ## Variation partitioning
# OTU_table <- gp23_group_sums
# transformation <- "hellinger"
# spe.hel <- decostand(OTU_table, transformation)
# spe.part.bio_chem <- varpart(spe.hel, envchem, envbio)
# spe.part.bio_chem
# plot(spe.part.bio_chem, digits=2)
# 
# 
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
# 
# anova.cca(rda(spe.hel,envchem),step=1000) 
# anova.cca(rda(spe.hel, envbio),step=1000) ## significant
# anova.cca(rda(spe.hel, envtemp),step=1000)
# anova.cca(rda(spe.hel, envrich),step=1000) 
# 
# pdf("../figures/Variation_partitioning_gp23_group_sum_chem_bio_time.pdf")
# spe.part.all <- varpart(spe.hel,envchem,envbio, envtemp)
# plot(spe.part.all, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()
# 
# ## ------------------------------------------------------------------------
# ## Variation partitioning with forward selection to get parsimonious variables
# spe.chem <- rda(spe.hel, envchem)
# R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
# try(forward.sel(spe.hel, envchem, adjR2thresh = R2a.all.chem, nperm=9999))
# 
# spe.bio <- rda(spe.hel, envbio)
# R2a.all.bio <- RsquareAdj(spe.bio)$adj.r.squared
# try(forward.sel(spe.hel, envbio, adjR2thresh=R2a.all.bio, nperm=9999))
# 
# spe.temp <- rda(spe.hel, envtemp)
# R2a.all.temp <- RsquareAdj(spe.temp)$adj.r.squared
# try(forward.sel(spe.hel, envtemp, adjR2thresh=R2a.all.temp, nperm=9999))
# 
# spe.rich <- rda(spe.hel, envrich)
# R2a.all.rich <- RsquareAdj(spe.rich)$adj.r.squared
# try(forward.sel(spe.hel, envrich, adjR2thresh=R2a.all.rich, nperm=9999))
# 
# #parsimonioius subsets of explanatory variables(based on forward selections)
# envchem.pars <- envchem[,c("Temperature_YSI","Average_SiO2")]
# envbio.pars <- envbio[,c("Average_viral_abundance")] ## but this doesn't quite make sense because it was only 1. I will redo with only 1
# envtemp.pars <- envtemp[,c("day_length", "month_number")]
# #envrich.pars <- envrich[,c("richness.18S")]
# ## Variation partitioning
# 
# anova.cca(rda(spe.hel, envchem.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envbio.pars),step=1000) ## significant
# anova.cca(rda(spe.hel, envtemp.pars),step=1000) ## significant
# #anova.cca(rda(spe.hel, envrich.pars),step=1000) ## not significant

## using those that are significant
# (spe.part <- varpart(spe.hel, envchem.pars, envbio.pars, envtemp.pars))
# plot(spe.part, digits=2)
# 
# 
# pdf("../figures/Variation_partitioning_gp23_group_sum_chem_bio_time_forward_selection.pdf")
# plot(spe.part, digits=2,bg=c("skyblue", "pink1", "mediumorchid"),Xnames=c("Chem", "Bio", "Temporal"))
# dev.off()

# res.venn <- spe.part$part$indfract[-8,3]
# tot.venn <- spe.part$part$fract[-8,3]
# vennDiagram_gp23 <- draw.triple.venn(area1 = round(tot.venn[1]*100),
#                                     area2 = round(tot.venn[2]*100),
#                                     area3 = round(tot.venn[3]*100),
#                                     n12 = round((res.venn[4]+res.venn[7]) * 100),
#                                     n23 = round((res.venn[5]+res.venn[7]) *100),
#                                     n13 = round((res.venn[6]+res.venn[7]) * 100), 
#                                     n123 = round((res.venn[7]) *100),
#                                     category = c("envchem","envbio", "envtemp" ),
#                                     fill = c("skyblue", "pink1", "mediumorchid"), euler.d = TRUE,
#                                     scaled = TRUE)


# ### Print all together
# pdf("../figures/variation_partitioning_%03d.pdf", onefile=FALSE, width=20, height = 20)
# plot_grid(vennDiagram_18s,
#           vennDiagram_16s,
#           vennDiagram_MPL,
#           vennDiagram_gp23,
#           ncol=2,
#           labels = c("A", "B", "C", "D"),
#           align ="hv"
# )
# dev.off()


### Print all together
# pdf("../figures/RDA_scaling_2_all_communities%03d.pdf", onefile=FALSE, width=20, height = 20)
# plot_grid(S18_scaling2[[4]],
#           S16_scaling2[[4]],
#           MPL_scaling2[[4]],
#           gp23_scaling2[[4]],
#           ncol=2,
#           labels = c("A", "B", "C", "D"),
#           align ="hv"
# )
# dev.off()


# ### Print all together
# pdf("../figures/RDA_for_selected_scaling_2_all_communities%03d.pdf", onefile=FALSE, width=20, height = 20)
# plot_grid(S18_for_sel_scaling2[[4]],
#           S16_for_sel_scaling2[[4]],
#           MPL_for_sel_scaling2[[4]],
#           gp23_for_sel_scaling2[[4]],
#           ncol=2,
#           labels = c("A", "B", "C", "D"),
#           align ="hv"
#           )
# dev.off()

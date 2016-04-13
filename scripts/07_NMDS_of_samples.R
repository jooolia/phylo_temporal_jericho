## NMDS plots of the communities of the OTUs. 

library(ggplot2)
library(gridExtra)
library(vegan)
library(reshape2)
library(phyloseq)
library(dplyr)
library(ape)
library(cowplot)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_AVS_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_AVS1.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 

normalized_18s_OTUs_phytos <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Phytoplankton.tsv", row.names=1)
normalized_18s_OTUs_hetero <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_Heterotrophs.tsv",                                  row.names=1)


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

### get relevant VC numbers
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

## don't need the standard error:
enviro_keep <- c("VC_number",
                 "Date",
                 "season",
                 "pH",
                 "Temperature_YSI",
                 "Salinity_ppt_YSI",
                 "Dissolved_oxygen_percent",
                 "Average_viral_abundance",
                 "Average_bacterial_abundance",
                 "Average_chl_a",
                 "Average_PO4",
                 "Average_SiO2",
                 "Average_NO3_NO2"
) 

Jericho_data <- Jericho_data[ , names(Jericho_data) %in% enviro_keep]

## want to add the date and place of sequencing to the environmental data as a way to colour the nmds. replace the pool numbers with the VC numbers
Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)

new_Jericho <- Jericho_data
row.names(new_Jericho) <- new_Jericho$VC_number

## want to add the date and place of sequencing to the environmental data 
## as a way to colour the nmds. replace the pool numbers with the VC numbers
Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)

## add in the seq facility by amplicon
new_Jericho$seq_facility_gp23 <- na.omit(Library_metadata$Sequencing_facility[match(Library_metadata$gp23_VC_number, new_Jericho$VC_number)])
new_Jericho$seq_facility_AVS <- na.omit(Library_metadata$Sequencing_facility[match(Library_metadata$AVS_VC_number, new_Jericho$VC_number)])
new_Jericho$seq_facility_MPL <- na.omit(Library_metadata$Sequencing_facility[match(Library_metadata$MPL_VC_number, new_Jericho$VC_number)])
new_Jericho$seq_facility_16s <- na.omit(Library_metadata$Sequencing_facility[match(Library_metadata$X16s_VC_number, new_Jericho$VC_number)])
new_Jericho$seq_facility_18s <- na.omit(Library_metadata$Sequencing_facility[match(Library_metadata$X18s_VC_number, new_Jericho$VC_number)])

## row names are the VC number so this is redundant
new_Jericho <- subset(new_Jericho, select=-(VC_number))
env_data <- sample_data(new_Jericho)


#### not printing right...fix. 

plot_pcoa_OTUs <- function (normalized_OTUs, env_data, name) {
  otu_table <- otu_table(normalized_OTUs, taxa_are_rows=FALSE)
  physeq <- phyloseq(otu_table, env_data)
  chao_dist <- ordinate(physeq, method="PCoA", distance = "bray")
  p <- plot_ordination(physeq,chao_dist, color="season", label="Date")+
   geom_point(size = 3)
  # Add title to each plot
  p2 <- p + ggtitle(paste("Ordination", name," using distance method Chao with PCoA")) + theme_bw()
 # return(p)
  pdf( paste("../figures/OTU_table_PCoA_", name,"%03d.pdf", sep=""),width = 17, height = 17,onefile = FALSE) 
  print(p2)
  dev.off()
}

plot_pcoa_OTUs(normalized_18s_OTUs, env_data, "18s")
plot_pcoa_OTUs(normalized_18s_OTUs_phytos, env_data, "18s_phytos")
plot_pcoa_OTUs(normalized_18s_OTUs_hetero, env_data, "18s_hetero")
plot_pcoa_OTUs(normalized_16s_OTUs, env_data, "16s")
plot_pcoa_OTUs(normalized_gp23_OTUs, env_data, "gp23")
plot_pcoa_OTUs(normalized_AVS_OTUs, env_data, "AVS")
plot_pcoa_OTUs(normalized_MPL_OTUs, env_data, "MPL")

## need to plot points as different colours for environmental factors. 
## how to trace line through ordination plot?

NMDS_of_amplicons <- function (normalized_OTUs, name) {
 #testNMDS <- metaMDS(normalized_OTUs,distance="bray", autotransform=FALSE, trymax = 200)
 #normalized_OTUs <- normalized_AVS_OTUs
 normalized_otus_bray <- vegdist(normalized_OTUs, method="bray")
 PcoaBray <- pcoa(normalized_otus_bray )
 #plot ordination
 cords <- PcoaBray$vectors[,c(1:2)]
 df <-data.frame(as.character(rownames(cords)),cords)
 str(df)
 colnames(df)[1] = "Var2"
 str(df)
 print(df$Var2)
 
 p1 <-ggplot(data = df, aes(x = Axis.1, y = Axis.2)) + 
  geom_text(aes(label = Var2),size = 3,hjust = 0,vjust = 0, colour=line_colour) + 
  theme_bw() + 
  theme(legend.position = "none")
 
 pdf(paste("../figures/", name, "_OTU_table_NMDS%03d.pdf", sep=""), width = 11, height = 11, onefile = FALSE)
 print(p1)
dev.off()
 return(PcoaBray)

}


MPL_nmds <-NMDS_of_amplicons(normalized_MPL_OTUs, "MPL")
AVS_nmds <-NMDS_of_amplicons(normalized_AVS_OTUs,"AVS")
gp23_nmds <-NMDS_of_amplicons(normalized_gp23_OTUs, "gp23")
S18_nmds <-NMDS_of_amplicons(normalized_18s_OTUs, "18S")
S18_nmds_phyto <-NMDS_of_amplicons(normalized_18s_OTUs_phytos, "18S_phyto")
S18_nmds_hetero <-NMDS_of_amplicons(normalized_18s_OTUs_hetero, "18S_hetero")
S16_nmds <-NMDS_of_amplicons(normalized_16s_OTUs, "16S")

## Colour the NMDS by different environmental parameters ======

convert_vegan_NMDS_to_ggplot_coloured_by_enviro <- function (normalized_OTUs,
                                                             new_Jericho,
                                                             environmental_parameter) {
 
#  normalized_OTUs <- normalized_MPL_OTUs
#  environmental_parameter <- "pH"
 normalized_otus_bray <- vegdist(normalized_OTUs, method="bray")
 PcoaBray <- pcoa(normalized_otus_bray )
 #plot ordination
 cords <- PcoaBray$vectors[,c(1:2)]

#testNMDS <- metaMDS(normalized_OTUs,distance="bray")
Jericho_OTU <- droplevels(subset(new_Jericho, row.names(new_Jericho) %in% row.names(normalized_OTUs)))

data.scores  <-data.frame(as.character(rownames(cords)),cords)
#convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data

subsetted_jericho <- subset(Jericho_OTU, select=environmental_parameter)
subsetted_jericho$site <- rownames(subsetted_jericho)
data.scores <- merge(data.scores,subsetted_jericho, by=c("site" ))

### plotting ####
head(data.scores)
names(data.scores)
# print(data.scores$environmental_parameter)

p <- ggplot(data=data.scores,aes_string(x="Axis.1",y="Axis.2",group=environmental_parameter))+
 geom_point(aes_string(colour=environmental_parameter),
            size=3) +  # add the site labels
 scale_colour_gradient(low = "blue", high="red")+ # add the point markers
 geom_text(aes(label=site),size=6,vjust=0, colour="white") +
coord_equal() +
 theme_JAG_presentation()
print(p)
}

## take out the env parameters that aren't continuous and seq facility stuff for that function. Plot them individually
new_jericho_no_data <- subset(new_Jericho, select=-c(Date, pH,season,seq_facility_gp23,seq_facility_AVS,seq_facility_MPL,seq_facility_16s,seq_facility_18s ))

# ## Loop through all the environmental params ====
# for (x in names(new_jericho_no_data)){
#  print(x)
#  ## MPL =====
#  NMDS_param <- convert_vegan_NMDS_to_ggplot_coloured_by_enviro(normalized_MPL_OTUs,
#                                                                new_jericho_no_data, 
#                                                                x)
#  pdf(paste("../figures/MPL_OTU_table_NMDS_coloured_by_", x,".pdf", sep=""),width = 11, height = 11, onefile = TRUE)
#  print(NMDS_param)
#  dev.off()
#  
#  ## gp23 =====
#  NMDS_param <- convert_vegan_NMDS_to_ggplot_coloured_by_enviro(normalized_gp23_OTUs,
#                                                                new_jericho_no_data, 
#                                                                x)
#  pdf(paste("../figures/gp23_OTU_table_NMDS_coloured_by_", x,".pdf", sep=""),width = 11, height = 11, onefile = TRUE)
#  print(NMDS_param)
#  dev.off()
#  
#  ## AVS =====
#  NMDS_param <- convert_vegan_NMDS_to_ggplot_coloured_by_enviro(normalized_AVS_OTUs,
#                                                                new_jericho_no_data, 
#                                                                x)
#  pdf(paste("../figures/AVS_OTU_table_NMDS_coloured_by_", x,".pdf", sep=""),width = 11, height = 11, onefile = TRUE)
#  print(NMDS_param)
#  dev.off()
#  
#  ## 18s =====
#  NMDS_param <- convert_vegan_NMDS_to_ggplot_coloured_by_enviro(normalized_18s_OTUs,
#                                                                new_jericho_no_data, 
#                                                                x)
#  pdf(paste("../figures/18s_OTU_table_NMDS_coloured_by_", x,".pdf", sep=""),width = 11, height = 11, onefile = TRUE)
#  print(NMDS_param)
#  dev.off()
#  
#  ## 16s =====
#  NMDS_param <- convert_vegan_NMDS_to_ggplot_coloured_by_enviro(normalized_16s_OTUs,
#                                                                new_jericho_no_data, 
#                                                                x)
#  pdf(paste("../figures/16s_OTU_table_NMDS_coloured_by_", x,".pdf", sep=""),width = 11, height = 11, onefile = TRUE)
#  print(NMDS_param)
#  dev.off()
# }
# 

## Colour the NMDS by sequencing facility ======
convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility <-
 function (normalized_OTUs,
           new_Jericho,
           environmental_parameter) {
  
  normalized_otus_bray <- vegdist(normalized_OTUs, method="bray")
  PcoaBray <- pcoa(normalized_otus_bray )
  #plot ordination
  cords <- PcoaBray$vectors[,c(1:2)]

  Jericho_OTU <- droplevels(subset(new_Jericho, row.names(new_Jericho) %in% row.names(normalized_OTUs)))
  
  data.scores  <-data.frame(as.character(rownames(cords)),cords)
  
 # testNMDS <- metaMDS(normalized_OTUs,distance = "bray")
  #Jericho_OTU <-droplevels(subset(new_Jericho, row.names(new_Jericho) %in% row.names(normalized_OTUs)))
  #data.scores <- as.data.frame(scores(testNMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
  data.scores$site <-   rownames(data.scores)  # create a column of site names, from the rownames of data.scores
  #data.scores$grp <- grp  #  add the grp variable created earlier
  
  subsetted_jericho <- subset(Jericho_OTU, select = environmental_parameter)
  subsetted_jericho$site <- rownames(subsetted_jericho)
  data.scores <- merge(data.scores,subsetted_jericho, by = c("site"))
  
  ### plotting ####
  p <-  ggplot(data = data.scores,aes_string(x = "Axis.1",y = "Axis.2",group = environmental_parameter)) +
   geom_point(aes_string(colour = environmental_parameter),size = 3) +  # add the site labels
   #  scale_colour_gradient(low = "blue", high="red")+ # add the point markers
   geom_text(aes(label = site),size = 6,vjust = 0, colour=line_colour) +
   coord_equal() +
 #  annotate("text", x = min(data.scores$Axis.1), y = max(data.scores$Axis.2),label = paste0("Stress: ", round(testNMDS$stress, digits =2)), hjust=1)+
   theme_JAG_presentation()
  print(p)
 }

# pdf("../figures/16s_OTU_table_NMDS_coloured_by_sequencing_facility.pdf",width = 11, height = 11, onefile = TRUE)
# convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility(normalized_16s_OTUs,
#                                                                    new_Jericho,
#                                                       "seq_facility_16s")
# dev.off()
# pdf("../figures/18s_OTU_table_NMDS_coloured_by_sequencing_facility.pdf",width = 11, height = 11, onefile = TRUE)
# convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility(normalized_18s_OTUs,
#                                                       new_Jericho,
#                                                       "seq_facility_18s")
# dev.off()
# pdf("../figures/AVS_OTU_table_NMDS_coloured_by_sequencing_facility.pdf",width = 11, height = 11, onefile = TRUE)
# convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility(normalized_AVS_OTUs,
#                                                       new_Jericho,
#                                                       "seq_facility_AVS")
# dev.off()
# pdf("../figures/gp23_OTU_table_NMDS_coloured_by_sequencing_facility.pdf",width = 11, height = 11, onefile = TRUE)
# convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility(normalized_gp23_OTUs,
#                                                       new_Jericho,
#                                                       "seq_facility_gp23")
# dev.off()
# pdf("../figures/MPL_OTU_table_NMDS_coloured_by_sequencing_facility.pdf",width = 11, height = 11, onefile = TRUE)
# convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility(normalized_MPL_OTUs,
#                                                       new_Jericho,
#                                                       "seq_facility_MPL")
# dev.off()
# ## by season


# pdf("../figures/MPL_OTU_table_NMDS_coloured_by_season.pdf",width = 11, height = 11, onefile = TRUE)
# convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility(normalized_MPL_OTUs,
#                                                       new_Jericho,
#                                                       "season")
# dev.off()

# pdf("../figures/gp23_OTU_table_NMDS_coloured_by_season.pdf",width = 11, height = 11, onefile = TRUE)
# convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility(normalized_gp23_OTUs,
#                                                       new_Jericho,
#                                                       "season")
# dev.off()
# pdf("../figures/16S_OTU_table_NMDS_coloured_by_season.pdf",width = 11, height = 11, onefile = TRUE)
# convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility(normalized_16s_OTUs,
#                                                       new_Jericho,
#                                                       "season")
# dev.off()
# pdf("../figures/18S_OTU_table_NMDS_coloured_by_season.pdf",width = 11, height = 11, onefile = TRUE)
# convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility(normalized_18s_OTUs,
#                                                       new_Jericho,
#                                                       "season")
# dev.off()
# pdf("../figures/AVS_OTU_table_NMDS_coloured_by_season.pdf",width = 11, height = 11, onefile = TRUE)
# convert_vegan_NMDS_to_ggplot_coloured_by_seq_facility(normalized_AVS_OTUs,
#                                                       new_Jericho,
#                                                       "season")
# dev.off()


## try season with a path 

convert_vegan_NMDS_to_ggplot_coloured_by_season_with_lines <- function (normalized_OTUs, new_Jericho,environmental_parameter) {
 
#  normalized_OTUs <- normalized_MPL_OTUs
#  environmental_parameter <- "season"
 
 normalized_otus_bray <- vegdist(normalized_OTUs, method="bray")
 PcoaBray <- pcoa(normalized_otus_bray )
 #plot ordination
 cords <- PcoaBray$vectors[,c(1:2)]
 
 #testNMDS <- metaMDS(normalized_OTUs,distance="bray")
 Jericho_OTU <- droplevels(subset(new_Jericho, row.names(new_Jericho) %in% row.names(normalized_OTUs)))
 
 data.scores  <-data.frame(as.character(rownames(cords)),cords)
 
#  testNMDS <- metaMDS(normalized_OTUs,distance="bray", trymax=100)
#  Jericho_OTU <- droplevels(subset(new_Jericho, row.names(new_Jericho) %in% row.names(normalized_OTUs)))
#  data.scores <- as.data.frame(scores(testNMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
 data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores

 subsetted_jericho <- subset(Jericho_OTU, select=environmental_parameter)
 subsetted_jericho$site <- rownames(subsetted_jericho)
 data.scores <- merge(data.scores,subsetted_jericho, by=c("site" ))
 
 ## match the dates...
 data.scores$date <- Jericho_OTU$Date[match(data.scores$site, rownames(Jericho_OTU))]
 
 ## need to annotate the start and end
 start_date <- filter(data.scores, date==min(data.scores$date))
 end_date <- filter(data.scores, date==max(data.scores$date))
 
 ## in between these points..

 spring_bloom <- as.Date("2011-04-08")
 #spring_bloom_date <- filter(data.scores, )
test_date_1 <-  which(abs(data.scores$date-spring_bloom) == min(abs(data.scores$date-spring_bloom)))
test_date_next <-  which(abs(data.scores$date-spring_bloom) == min(abs(data.scores$date-spring_bloom)))+1

spring_bloom_point_1 <- filter(data.scores, date==data.scores$date[test_date_1])
spring_bloom_point_next <- filter(data.scores, date==data.scores$date[test_date_next])
spring_bloom_point <- data.frame(Axis.1=as.numeric(spring_bloom_point_1$Axis.1-((spring_bloom_point_1$Axis.1 - spring_bloom_point_next$Axis.1)/2)),
                                 Axis.2=as.numeric(spring_bloom_point_1$Axis.2-((spring_bloom_point_1$Axis.2 - spring_bloom_point_next$Axis.2)/2)))


 p <- ggplot(data=data.scores,aes(x=Axis.1,y=Axis.2))+ 
  geom_point(data=start_date, aes(x=Axis.1,y=Axis.2),
             size=12,
             colour=line_colour,
             fill=NA)+
  geom_point(data=end_date,
             aes(x=Axis.1,y=Axis.2),
             size=12,
             colour=line_colour,
             fill=NA)+
  geom_point(aes_string(colour=environmental_parameter),
             size=10) +  # add the site labels
  geom_path(aes(x=Axis.1,y=Axis.2),
            size=0.2,
            linetype=1,
            colour=line_colour)+
  annotate("text",
           label = "start",
           x=start_date$Axis.1,
           y=start_date$Axis.2,
           size = 6,
           hjust=1,
           vjust=-1.5,
           colour=line_colour)+
 annotate("text",
          label = "end",
          x=end_date$Axis.1,
          y=end_date$Axis.2,
          size = 6,
          hjust=1,
          vjust=-1.5,
          colour=line_colour)+
  annotate("text",
           label = "spring bloom",
           x=spring_bloom_point$Axis.1,
           y=spring_bloom_point$Axis.2,
           size = 3,
           colour=line_colour)+
  coord_equal() +
 #annotate("text", x = (min(data.scores$Axis.1))+0.5, y = (max(data.scores$Axis.2)),label = paste0("Stress: ", round(testNMDS$stress, digits =2))) +
 theme_JAG_presentation()
 print(p)
 #stressplot(testNMDS)
 return(p)
}

pdf("../figures/MPL_OTU_table_NMDS_coloured_by_season_with_lines%03d.pdf",width = 11, height = 11, onefile = FALSE)
MPL_NMDS <- convert_vegan_NMDS_to_ggplot_coloured_by_season_with_lines(normalized_MPL_OTUs,new_Jericho, "season")
MPL_NMDS
dev.off()

# pdf("../figures/AVS_OTU_table_NMDS_coloured_by_season_with_lines.pdf",width = 11, height = 11, onefile = TRUE)
# AVS_NMDS <-convert_vegan_NMDS_to_ggplot_coloured_by_season_with_lines(normalized_AVS_OTUs,new_Jericho,"season")
# AVS_NMDS
# dev.off()

pdf("../figures/18s_OTU_table_NMDS_coloured_by_season_with_lines%03d.pdf",width = 11, height = 11, onefile = FALSE)
S18_NMDS <- convert_vegan_NMDS_to_ggplot_coloured_by_season_with_lines(normalized_18s_OTUs,new_Jericho, "season")
S18_NMDS
dev.off()

pdf("../figures/18s_OTU_table_phytos_NMDS_coloured_by_season_with_lines%03d.pdf",width = 11, height = 11, onefile = FALSE)
S18_NMDS <- convert_vegan_NMDS_to_ggplot_coloured_by_season_with_lines(normalized_18s_OTUs_phytos, new_Jericho,"season")
S18_NMDS
dev.off()

pdf("../figures/18s_OTU_table_hetero_NMDS_coloured_by_season_with_lines%03d.pdf",width = 11, height = 11, onefile = FALSE)
S18_NMDS <- convert_vegan_NMDS_to_ggplot_coloured_by_season_with_lines(normalized_18s_OTUs_hetero, new_Jericho,"season")
S18_NMDS
dev.off()


pdf("../figures/gp23_OTU_table_NMDS_coloured_by_season_with_lines%03d.pdf",width = 11, height = 11, onefile = FALSE)
gp23_NMDS <-convert_vegan_NMDS_to_ggplot_coloured_by_season_with_lines(normalized_gp23_OTUs, new_Jericho, "season")
gp23_NMDS
dev.off()

pdf("../figures/16s_OTU_table_NMDS_coloured_by_season_with_lines%03d.pdf",width = 11, height = 11, onefile = FALSE)
S16_NMDS <- convert_vegan_NMDS_to_ggplot_coloured_by_season_with_lines(normalized_16s_OTUs,new_Jericho,"season")
S16_NMDS
dev.off()


pdf("../figures/all_NMDS_coloured_by_season_with_lines%03d.pdf",width = 20, height = 20, onefile = FALSE)
plot_grid(MPL_NMDS,gp23_NMDS,S18_NMDS,S16_NMDS, 
          ncol=2,
          labels = c("A", "B", "C", "D"),
          label_size = 20,
          align = c("v", "h"))

dev.off()



procrustes_of_pcoa <- function (normalized_A_OTUs, normalized_B_OTUs) {
  ## need same dates and then compare
  union_vcs_A_and_B <- intersect(row.names(normalized_A_OTUs), row.names(normalized_B_OTUs))
  
  Jericho_A <- droplevels(subset(normalized_A_OTUs, row.names(normalized_A_OTUs) %in% union_vcs_A_and_B))
  
  Jericho_B <- droplevels(subset(normalized_B_OTUs, row.names(normalized_B_OTUs) %in% union_vcs_A_and_B))
  
  testNMDS_A <- metaMDS(Jericho_A ,distance="bray", trymax=100)
  plot(testNMDS_A)
  
  testNMDS_B <- metaMDS(Jericho_B ,distance="bray", trymax=100)
  plot(testNMDS_B)
  
  normalized_otus_bray <- vegdist(Jericho_B, method="bray")
  PcoaBray_B <- cmdscale(normalized_otus_bray )
  
  normalized_otus_bray <- vegdist(Jericho_A, method="bray")
  PcoaBray_A <- cmdscale(normalized_otus_bray )
  
  proc_A_B <- procrustes(PcoaBray_B,PcoaBray_A)
  print(summary(proc_A_B))
  plot(proc_A_B)
}
## S18 and MPL
proc_MPL_and_18s <- procrustes_of_pcoa(normalized_MPL_OTUs, normalized_18s_OTUs)
## try with AVS and 18s
proc_AVS_and_18s <- procrustes_of_pcoa(normalized_AVS_OTUs, normalized_18s_OTUs)
## try with gp23 and 16s
proc_gp23_and_16s <- procrustes_of_pcoa(normalized_gp23_OTUs, normalized_16s_OTUs)

proc_gp23_and_18s <- procrustes_of_pcoa(normalized_gp23_OTUs, normalized_18s_OTUs)

proc_MPL_and_16s <- procrustes_of_pcoa(normalized_MPL_OTUs, normalized_16s_OTUs)

proc_MPL_and_18s_phyto <- procrustes_of_pcoa(normalized_MPL_OTUs, normalized_18s_OTUs_phytos)

proc_MPL_and_18s_hetero <- procrustes_of_pcoa(normalized_MPL_OTUs, normalized_18s_OTUs_hetero)

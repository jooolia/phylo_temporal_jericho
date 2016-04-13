## find number of edges between seasons
## find number of edges between sites X days apart


library(igraph)
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(cowplot)


args <- commandArgs(TRUE)
inputFile <- args[1]
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

graph_file <- read.graph("../results/LSA_tables/network_LSA_overall.graphml",format="graphml")


edges <- as.data.frame(ends(graph_file, E(graph_file)))

weight <- E(graph_file)$weight_LS

date_2010.06.22 <- E(graph_file)$t_2010.06.22
date_2010.07.06 <- E(graph_file)$t_2010.07.06
date_2010.07.20 <- E(graph_file)$t_2010.07.20
date_2010.08.05 <- E(graph_file)$t_2010.08.05
date_2010.08.17 <- E(graph_file)$t_2010.08.17
date_2010.08.31 <- E(graph_file)$t_2010.08.31
date_2010.09.14 <- E(graph_file)$t_2010.09.14
date_2010.09.28 <- E(graph_file)$t_2010.09.28
date_2010.10.12 <- E(graph_file)$t_2010.10.12
date_2010.10.26 <- E(graph_file)$t_2010.10.26
date_2010.11.09 <- E(graph_file)$t_2010.11.09
date_2010.11.23 <- E(graph_file)$t_2010.11.23
date_2010.12.06 <- E(graph_file)$t_2010.12.06
date_2010.12.20 <- E(graph_file)$t_2010.12.20
date_2011.01.10 <- E(graph_file)$t_2011.01.10
date_2011.01.18 <- E(graph_file)$t_2011.01.18
date_2011.01.29 <- E(graph_file)$t_2011.01.29
date_2011.01.31 <- E(graph_file)$t_2011.01.31
date_2011.02.02 <- E(graph_file)$t_2011.02.02
date_2011.02.04 <- E(graph_file)$t_2011.02.04
date_2011.02.06 <- E(graph_file)$t_2011.02.06
date_2011.02.08 <- E(graph_file)$t_2011.02.08
date_2011.02.10 <- E(graph_file)$t_2011.02.10
date_2011.02.24 <- E(graph_file)$t_2011.02.24
date_2011.03.10 <- E(graph_file)$t_2011.03.10
date_2011.03.24 <- E(graph_file)$t_2011.03.24
date_2011.04.07 <- E(graph_file)$t_2011.04.07
date_2011.04.21 <- E(graph_file)$t_2011.04.21
date_2011.05.05 <- E(graph_file)$t_2011.05.05
date_2011.05.18 <- E(graph_file)$t_2011.05.18
date_2011.06.07 <- E(graph_file)$t_2011.06.07
date_2011.06.21 <- E(graph_file)$t_2011.06.21
date_2011.07.05 <- E(graph_file)$t_2011.07.05
date_2011.07.19 <- E(graph_file)$t_2011.07.19


edges <- cbind(edges, weight,date_2010.06.22,
               date_2010.07.06,
               date_2010.07.20,
               date_2010.08.05,
               date_2010.08.17,
               date_2010.08.31,
               date_2010.09.14,
               date_2010.09.28 ,
               date_2010.10.12 ,
               date_2010.10.26 ,
               date_2010.11.09 ,
               date_2010.11.23 ,
               date_2010.12.06 ,
               date_2010.12.20 ,
               date_2011.01.10 ,
               date_2011.01.18 ,
               date_2011.01.29 ,
               date_2011.01.31 ,
               date_2011.02.02 ,
               date_2011.02.04 ,
               date_2011.02.06 ,
               date_2011.02.08 ,
               date_2011.02.10 ,
               date_2011.02.24 ,
               date_2011.03.10 ,
               date_2011.03.24 ,
               date_2011.04.07 ,
               date_2011.04.21 ,
               date_2011.05.05 ,
               date_2011.05.18 ,
               date_2011.06.07 ,
               date_2011.06.21 ,
               date_2011.07.05 ,
               date_2011.07.19)

## ok so know which edges appear when.

## so if edge occurs in summer we count that occurence, if in fall, count that. 
spring_dates <- c(seq(as.Date("2010-03-22"), 
                      as.Date("2010-06-22"), "days"),
                  seq(as.Date("2011-03-22"),
                      as.Date("2011-06-22"), "days"))
summer_dates <- c(seq(as.Date("2010-06-22"), 
                    as.Date("2010-09-22"), "days"),
                  seq(as.Date("2011-06-22"), 
                      as.Date("2011-09-22"), "days"))
fall_dates <- seq(as.Date("2010-09-22"), 
                  as.Date("2010-12-22"), "days")
winter_dates <- seq(as.Date("2010-12-22"), 
                    as.Date("2011-03-22"), "days")
 


##fix naming
colnames(edges ) <- gsub("date_([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})","\\1-\\2-\\3" , colnames(edges ))
names(edges)

spring_presence <- edges[colnames(edges) %in% c("V1", "V2", "weight",as.character(spring_dates))]
spring_presence[colnames(spring_presence) %in% c(weight, as.character(spring_dates))] <- sapply(spring_presence[colnames(spring_presence) %in% c(weight, as.character(spring_dates))], as.character)
spring_presence[colnames(spring_presence) %in% c(weight, as.character(spring_dates))] <- sapply(spring_presence[colnames(spring_presence) %in% c(weight, as.character(spring_dates))], as.numeric)



summer_presence <- edges[colnames(edges) %in% c("V1", "V2", "weight",as.character(summer_dates))]
## convert appropriate columns to character and then numeric?
summer_presence[colnames(summer_presence) %in% c(weight, as.character(summer_dates))] <- sapply(summer_presence[colnames(summer_presence) %in% c(weight, as.character(summer_dates))], as.character)
summer_presence[colnames(summer_presence) %in% c(weight, as.character(summer_dates))] <- sapply(summer_presence[colnames(summer_presence) %in% c(weight, as.character(summer_dates))], as.numeric)

fall_presence <- edges[colnames(edges) %in% c("V1", "V2", "weight",as.character(fall_dates))]
fall_presence[colnames(fall_presence) %in% c(weight, as.character(fall_dates))] <- sapply(fall_presence[colnames(fall_presence) %in% c(weight, as.character(fall_dates))], as.character)
fall_presence[colnames(fall_presence) %in% c(weight, as.character(fall_dates))] <- sapply(fall_presence[colnames(fall_presence) %in% c(weight, as.character(fall_dates))], as.numeric)


winter_presence <- edges[colnames(edges) %in% c("V1", "V2", "weight",as.character(winter_dates))]
winter_presence[colnames(winter_presence) %in% c(weight, as.character(winter_dates))] <- sapply(winter_presence[colnames(winter_presence) %in% c(weight, as.character(winter_dates))], as.character)
winter_presence[colnames(winter_presence) %in% c(weight, as.character(winter_dates))] <- sapply(winter_presence[colnames(winter_presence) %in% c(weight, as.character(winter_dates))], as.numeric)


## rows present -so any date above 0
## so sum columns in date columns and see if above 0
spring_presence_above_0 <- filter(spring_presence, rowSums(spring_presence[colnames(spring_presence) %in% as.character(spring_dates)], na.rm=TRUE) > 0)

summer_presence_above_0 <- filter(summer_presence, rowSums(summer_presence[colnames(summer_presence) %in% as.character(summer_dates)], na.rm=TRUE) > 0)

fall_presence_above_0 <- filter(fall_presence, rowSums(fall_presence[colnames(fall_presence) %in% as.character(fall_dates)], na.rm=TRUE) > 0)

winter_presence_above_0 <- filter(winter_presence, rowSums(winter_presence[colnames(winter_presence) %in% as.character(winter_dates)], na.rm=TRUE) > 0)

## how to combine to find edges present in all?
## want to see if edge in both ones. intersection

spring_summer_edges <- merge(spring_presence_above_0, summer_presence_above_0, by=c("V1", "V2", "weight"))

fall_winter_edges <- merge(fall_presence_above_0, winter_presence_above_0, by=c("V1", "V2","weight"))

spring_fall_edges <- merge(spring_presence_above_0, fall_presence_above_0, by=c("V1", "V2","weight"))

summer_winter_edges <- merge(summer_presence_above_0, winter_presence_above_0, by=c("V1", "V2","weight"))

spring_winter_edges <- merge(spring_presence_above_0, winter_presence_above_0, by=c("V1", "V2","weight"))

summer_fall_edges <- merge(summer_presence_above_0, fall_presence_above_0, by=c("V1", "V2","weight"))


present_all_season_edges <- merge(spring_summer_edges, fall_winter_edges, by=c("V1", "V2","weight"))


number_nodes_in_season <- function (season_presence_above_0) {
  summer_otusV1 <- as.character(unique(season_presence_above_0$V1))
  summer_otusV2 <- as.character(unique(season_presence_above_0$V2))
  summer_otus_V1_V2 <- c(summer_otusV1, summer_otusV2)
  summer_otus <- unique(summer_otus_V1_V2)
  return(summer_otus)
}

spring_otus <- number_nodes_in_season(spring_presence_above_0)
summer_otus <- number_nodes_in_season(summer_presence_above_0)
fall_otus <- number_nodes_in_season(fall_presence_above_0)
winter_otus <- number_nodes_in_season(winter_presence_above_0)

summer_spring_otus <- number_nodes_in_season(spring_summer_edges)
fall_winter_otus <- number_nodes_in_season(fall_winter_edges)
spring_fall_otus <- number_nodes_in_season(spring_fall_edges)
summer_winter_otus <- number_nodes_in_season(summer_winter_edges)
spring_winter_otus <- number_nodes_in_season(spring_winter_edges)
summer_fall_otus <- number_nodes_in_season(summer_fall_edges)
present_all_season_nodes <- number_nodes_in_season(present_all_season_edges)


## if divided by number of dates??
spring_date_names <- colnames(edges)[colnames(edges) %in% c(as.character(spring_dates))]
summer_date_names <- colnames(edges)[colnames(edges) %in% c(as.character(summer_dates))]
fall_date_names <- colnames(edges)[colnames(edges) %in% c(as.character(fall_dates))]
winter_date_names <- colnames(edges)[colnames(edges) %in% c(as.character(winter_dates))]


count_edges_summer <- c(pos=(dim(filter(summer_presence_above_0, weight > 0))[1])/length(summer_date_names),
                        neg=(dim(filter(summer_presence_above_0, weight < 0))[1])/length(summer_date_names),
                        nodes = length(summer_otus))
count_edges_spring <- c(pos=(dim(filter(spring_presence_above_0, weight > 0))[1])/length(spring_date_names),
                        neg=(dim(filter(spring_presence_above_0, weight < 0))[1])/length(spring_date_names),
                        nodes = length(spring_otus))

count_edges_fall <- c(pos=(dim(filter(fall_presence_above_0, weight > 0))[1])/length(fall_date_names),
                      neg=(dim(filter(fall_presence_above_0, weight < 0))[1])/length(fall_date_names),
                      nodes = length(fall_otus))

count_edges_winter <- c(pos=(dim(filter(winter_presence_above_0, weight > 0))[1])/length(winter_date_names),
                        neg=(dim(filter(winter_presence_above_0, weight < 0))[1])/length(winter_date_names),
                        nodes = length(winter_otus))

count_edges_spring_summer <- c(pos=(dim(filter(spring_summer_edges, weight > 0))[1])/length(c(summer_date_names, spring_date_names)),
                        neg=(dim(filter(spring_summer_edges, weight < 0))[1])/length(c(summer_date_names, spring_date_names)),
                        nodes = length(summer_spring_otus))

count_edges_fall_winter <- c(pos=(dim(filter(fall_winter_edges, weight > 0))[1])/length(c(fall_date_names,winter_date_names)),
                               neg=(dim(filter(fall_winter_edges, weight < 0))[1])/length(c(fall_date_names,winter_date_names)),
                               nodes = length(fall_winter_otus))

count_edges_spring_fall <- c(pos=(dim(filter(spring_fall_edges, weight > 0))[1])/length(c(fall_date_names, spring_date_names)),
                               neg=(dim(filter(spring_fall_edges, weight < 0))[1])/length(c(fall_date_names, spring_date_names)),
                               nodes = length(spring_fall_otus))

count_edges_summer_winter <- c(pos=(dim(filter(summer_winter_edges, weight > 0))[1])/length(c(summer_date_names,winter_date_names)),
                             neg=(dim(filter(summer_winter_edges, weight < 0))[1])/length(c(summer_date_names,winter_date_names)),
                             nodes = length(summer_winter_otus))


count_edges_spring_winter <- c(pos=(dim(filter(spring_winter_edges, weight > 0))[1])/length(c(spring_date_names,winter_date_names)),
                               neg=(dim(filter(spring_winter_edges, weight < 0))[1])/length(c(spring_date_names,winter_date_names)),
                               nodes = length(spring_winter_otus))

count_edges_summer_fall <- c(pos=(dim(filter(summer_fall_edges, weight > 0))[1])/length(c(summer_date_names,fall_date_names)),
                               neg=(dim(filter(summer_fall_edges, weight < 0))[1])/length(c(summer_date_names,fall_date_names)),
                               nodes = length(summer_fall_otus))

count_edges_always_present <- c(pos=(dim(filter(present_all_season_edges, weight > 0))[1])/length(c(spring_date_names,summer_date_names,fall_date_names,winter_date_names)),
                               neg=(dim(filter(present_all_season_edges, weight < 0))[1])/length(c(spring_date_names,summer_date_names,fall_date_names,winter_date_names)),
                               nodes = length(present_all_season_nodes))



counts_for_barplot <- rbind(count_edges_spring,
                            count_edges_summer,
                            count_edges_fall,
                            count_edges_winter,
                            count_edges_spring_summer,
                            count_edges_fall_winter,
                            count_edges_spring_fall,
                            count_edges_spring_summer,
                            count_edges_summer_winter,
                            count_edges_spring_winter,
                            count_edges_summer_fall,
                            count_edges_always_present
) 

melted_counts <- melt(counts_for_barplot, value.name="Number_of_Edges_Div_by_Date")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"
melted_counts$Between_which_communities <- gsub("count_edges_","", melted_counts$Between_which_communities)

melted_counts$Between_which_communities <- ordered(melted_counts$Between_which_communities, levels = c("summer","fall","winter", "spring","summer_fall", "fall_winter" , "spring_winter", "spring_summer"  ,"summer_winter", "spring_fall","always_present" ))


## need to standardize this effort...by potential sample

pdf("../figures/overall_network_barplot_of_edges_by_season.pdf", width = 15, height = 11)




edges_by_season <- ggplot(subset(melted_counts, Sign_of_relationship %in% c("neg", "pos")), aes(x=Between_which_communities, y=Number_of_Edges_Div_by_Date, fill=Sign_of_relationship))+
 geom_bar(stat="identity")+
 
 # geom_text(aes(label=paste0("nodes:\n",
 #                            subset(melted_counts, Sign_of_relationship %in% "nodes")$Number_of_Edges_Div_by_Date)), vjust=1.5,
 #           colour="white")+
 annotate("text",
          x=subset(melted_counts, Sign_of_relationship %in% "nodes")$Between_which_communities,
          y=75,
          # label="test")
          label=paste0("nodes:\n",
                       subset(melted_counts, Sign_of_relationship %in% "nodes")$Number_of_Edges_Div_by_Date))+
 ggtitle("")+
 theme_JAG_presentation() +
 theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))+
 theme(axis.title.y = element_text(face="bold", size=14),
       legend.text = element_text(size = 20, face = "bold"))+
 ylab("Edges divided by \n# of samples in group")+
 xlab("")+
 scale_fill_discrete(name="",
                     labels=c("Positive edges", "Negative edges"))
edges_by_season

dev.off()


### so to combine

## Make plot of Number of edges vs. the tupe of relationship

# only_different <- delete.edges(graph_file, E(graph_file) [str_split(ends(graph_file, 1), "_")[[1]][1]== str_split(ends(graph_file, 1), "_", n=2)[[2]][1]])
## so want to find edges -euk to euk S18 to S18, Bac to euk, euk ot vir,

edges <- as.data.frame(ends(graph_file, E(graph_file)))

weight <- E(graph_file)$weight_LS

edges <- cbind(edges, weight)

## between same types
edges_between_like <-
 edges[str_split_fixed(edges[,1], "_", n = 2)[,1] == str_split_fixed(edges[,2], "_", n =
                                                                      2)[,1],]

euk_to_euk <-
 edges_between_like [str_split_fixed(edges_between_like [,1], "_", n = 2)[,1] ==
                      "S18",]
bac_to_bac <-
 edges_between_like [str_split_fixed(edges_between_like [,1], "_", n = 2)[,1] ==
                      "S16",]
gp23_to_gp23 <-
 edges_between_like [str_split_fixed(edges_between_like [,1], "_", n = 2)[,1] ==
                      "gp23OTU",]
MPL_to_MPL <-
 edges_between_like [str_split_fixed(edges_between_like [,1], "_", n = 2)[,1] ==
                      "MPLOTU",]


edges_between_unlike <-
 edges[str_split_fixed(edges[,1], "_", n = 2)[,1] != str_split_fixed(edges[,2], "_", n =
                                                                      2)[,1],]

## find the edges either way: either in column 1 or in column 2
euk_to_bac <-
 edges_between_unlike[str_split_fixed(edges_between_unlike [,1], "_", n =
                                       2)[,1] == "S18" |
                       str_split_fixed(edges_between_unlike [,1], "_", n = 2)[,1] == "S16" ,]
euk_to_bac <-
 euk_to_bac[str_split_fixed(euk_to_bac [,2], "_", n = 2)[,1] == "S18" |
             str_split_fixed(euk_to_bac [,2], "_", n = 2)[,1] == "S16" ,]

# euk_to_AVS <- edges_between_unlike[str_split_fixed(edges_between_unlike [,1], "_", n=2)[,1]=="S18" |str_split_fixed(edges_between_unlike [,1], "_", n=2)[,1]=="AVSOTU" ,]
# euk_to_AVS <- euk_to_AVS[str_split_fixed(euk_to_AVS[,2], "_", n=2)[,1]=="S18" |str_split_fixed(euk_to_AVS[,2], "_", n=2)[,1]=="AVSOTU" ,]

euk_to_MPL <-
 edges_between_unlike[str_split_fixed(edges_between_unlike [,1], "_", n =
                                       2)[,1] == "S18" |
                       str_split_fixed(edges_between_unlike [,1], "_", n = 2)[,1] == "MPLOTU" ,]
euk_to_MPL <-
 euk_to_MPL[str_split_fixed(euk_to_MPL[,2], "_", n = 2)[,1] == "S18" |
             str_split_fixed(euk_to_MPL[,2], "_", n = 2)[,1] == "MPLOTU" ,]

bac_to_gp23 <-
 edges_between_unlike[str_split_fixed(edges_between_unlike [,1], "_", n =
                                       2)[,1] == "S16" |
                       str_split_fixed(edges_between_unlike [,1], "_", n = 2)[,1] == "gp23OTU" ,]
bac_to_gp23 <-
 bac_to_gp23[str_split_fixed(bac_to_gp23[,2], "_", n = 2)[,1] == "S16" |
              str_split_fixed(bac_to_gp23[,2], "_", n = 2)[,1] == "gp23OTU" ,]

taxonomy_18s <-
 read.csv("../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names =
           1)
taxonomy_16s <-
 read.csv("../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names =
           1)


## want to see the main taxonomic groups in the bac to gp23

unique(taxonomy_16s$Phylum[match(gsub("S16_","OTU_", bac_to_gp23$V1), taxonomy_16s$otu_number)])
unique(taxonomy_16s$Class[match(gsub("S16_","OTU_", bac_to_gp23$V1), taxonomy_16s$otu_number)])
unique(taxonomy_16s$Order[match(gsub("S16_","OTU_", bac_to_gp23$V1), taxonomy_16s$otu_number)])


## what about for 18s

unique(taxonomy_18s$Phylum[match(gsub("S18_","OTU_", euk_to_MPL$V1), taxonomy_18s$otu_number)])
unique(taxonomy_18s$Class[match(gsub("S18_","OTU_", euk_to_MPL$V1), taxonomy_18s$otu_number)])
unique(taxonomy_18s$Order[match(gsub("S18_","OTU_", euk_to_MPL$V1), taxonomy_18s$otu_number)])
unique(taxonomy_18s$Family[match(gsub("S18_","OTU_", euk_to_MPL$V1), taxonomy_18s$otu_number)])




number_nodes_in_season <- function (season_presence_above_0) {
 summer_otusV1 <- as.character(unique(season_presence_above_0$V1))
 summer_otusV2 <- as.character(unique(season_presence_above_0$V2))
 summer_otus_V1_V2 <- c(summer_otusV1, summer_otusV2)
 summer_otus <- unique(summer_otus_V1_V2)
 return(summer_otus)
}

nodes_between_like <- number_nodes_in_season(edges_between_like)
nodes_between_unlike <- number_nodes_in_season(edges_between_unlike)
nodes_between_euk_to_euk <- number_nodes_in_season(euk_to_euk)
nodes_between_bac_to_bac <- number_nodes_in_season(bac_to_bac)
nodes_between_euk_to_bac <- number_nodes_in_season(euk_to_bac)
nodes_between_euk_to_MPL <- number_nodes_in_season(euk_to_MPL)
nodes_between_bac_to_gp23 <- number_nodes_in_season(bac_to_gp23)
nodes_between_MPL_to_MPL <- number_nodes_in_season(MPL_to_MPL)
nodes_between_gp23_to_gp23 <- number_nodes_in_season(gp23_to_gp23)

## count negative vs. positive and do a histogram
count_edges_between_like <-
 c(
  pos = dim(filter(edges_between_like, weight > 0))[1],
  neg = dim(filter(edges_between_like, weight < 0))[1],
  nodes = length(nodes_between_like)
 )

count_edges_between_unlike <-
 c(
  pos = dim(filter(edges_between_unlike, weight > 0))[1],
  neg = dim(filter(edges_between_unlike, weight < 0))[1],
  nodes = length(nodes_between_unlike)
 )

count_euk_to_euk <- c(
 pos = dim(filter(euk_to_euk, weight > 0))[1],
 neg = dim(filter(euk_to_euk, weight < 0))[1],
 nodes = length(nodes_between_euk_to_euk)
)

count_bac_to_bac <- c(
 pos = dim(filter(bac_to_bac, weight > 0))[1],
 neg = dim(filter(bac_to_bac, weight < 0))[1],
 nodes = length(nodes_between_bac_to_bac)
)

count_euk_to_bac <- c(
 pos = dim(filter(euk_to_bac, weight > 0))[1],
 neg = dim(filter(euk_to_bac, weight < 0))[1],
 nodes = length(nodes_between_euk_to_bac)
)
#
# count_euk_to_AVS <- c(pos=dim(filter(euk_to_AVS, weight > 0))[1],
#                       neg=dim(filter(euk_to_AVS, weight < 0))[1])

count_euk_to_MPL <- c(
 pos = dim(filter(euk_to_MPL, weight > 0))[1],
 neg = dim(filter(euk_to_MPL, weight < 0))[1],
 nodes = length(nodes_between_euk_to_MPL)
)

count_bac_to_gp23 <- c(
 pos = dim(filter(bac_to_gp23, weight > 0))[1],
 neg = dim(filter(bac_to_gp23, weight < 0))[1],
 nodes = length(nodes_between_bac_to_gp23)
)

count_MPL_to_MPL <- c(
 pos = dim(filter(MPL_to_MPL, weight > 0))[1],
 neg = dim(filter(MPL_to_MPL, weight < 0))[1],
 nodes = length(nodes_between_MPL_to_MPL)
)

count_gp23_to_gp23 <-
 c(
  pos = dim(filter(gp23_to_gp23, weight > 0))[1],
  neg = dim(filter(gp23_to_gp23, weight < 0))[1],
  nodes = length(nodes_between_gp23_to_gp23)
 )

counts_for_barplot <- rbind(
 count_edges_between_like,
 count_edges_between_unlike,
 count_euk_to_euk,
 count_bac_to_bac,
 count_euk_to_bac,
 #count_euk_to_AVS,
 count_euk_to_MPL,
 count_bac_to_gp23,
 count_MPL_to_MPL,
 count_gp23_to_gp23
)
## test some negs
#counts_for_barplot[,2] <- c(30,50,5,10,20,30,50,60)

melted_counts <- melt(counts_for_barplot, value.name = "Edges")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"
melted_counts$Between_which_communities <-
 gsub("count_edges_", "", melted_counts$Between_which_communities)
melted_counts$Between_which_communities <-
 gsub("count_", "", melted_counts$Between_which_communities)
melted_counts$Between_which_communities <-
 gsub("_", " ", melted_counts$Between_which_communities)
melted_counts$Between_which_communities <-
 ordered(
  melted_counts$Between_which_communities, levels = c(
   "between like",
   "between unlike",
   "euk to euk",
   "bac to bac",
   "MPL to MPL",
   "gp23 to gp23",
   "bac to gp23",
   "euk to bac" ,
   "euk to MPL"
  )
 )



pdf(
 "../figures/overall_network_barplot_of_edges.pdf", width = 11, height = 11
)
overall_edges_by_type <- ggplot(
 subset(melted_counts, Sign_of_relationship %in% c("neg", "pos")), aes(x =
                                                                        Between_which_communities, y = Edges, fill = Sign_of_relationship)
) +
 geom_bar(stat = "identity") +
 annotate(
  "text",
  x = subset(melted_counts, Sign_of_relationship %in% "nodes")$Between_which_communities,
  y = 500,
  label = paste0(
   "nodes:\n",
   subset(melted_counts, Sign_of_relationship %in% "nodes")$Edges
  )
 ) +
 theme_JAG_presentation() +
 theme(axis.title.y = element_text(face="bold",
                                   size=14),
       axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
       legend.text = element_text(size = 20, face = "bold")
       )+
 ylab("Edges divided by \n# of samples in group") +
 xlab("") +
 scale_fill_discrete(name = "",
                     labels = c("Positive edges", "Negative edges"))+
scale_x_discrete(labels=c("Between like",
                             "Between unlike",
                             "Between eukaryotes",
                             "Between bacteria",
                             "Between marine \npicorna-like viruses",
                             "Between T4-like myoviruses",
                             "Bacteria to T4-like myovirus",
                             "Eukaryotes to bacteria" ,
                             "Eukaryotes to marine\n picorna-like viruses"))
overall_edges_by_type
dev.off()


pdf("../figures/overall_network_barplots%1d.pdf", width = 15, height = 11, onefile=FALSE)
plot_grid(edges_by_season,
          overall_edges_by_type,
          ncol=1,
          labels = c("A", "B"),
          align = "v")
dev.off()

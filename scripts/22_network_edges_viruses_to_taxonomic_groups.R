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
 if (inputFile == "../../JAG_black_presentation.R") {
  path_colour <- "white"
  line_colour <- "white"
  point_colour <- "white"
 }
}



graph_file <-
 read.graph("../results/LSA_tables/network_LSA_overall.graphml",format =
             "graphml")


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
 "../figures/overall_network_barplot_of_edges.pdf", width = 15, height = 11
)
ggplot(
 subset(melted_counts, Sign_of_relationship %in% c("neg", "pos")), aes(x =
                                                                        Between_which_communities, y = Edges, fill = Sign_of_relationship)
) +
 geom_bar(stat = "identity") +
 annotate(
  "text",
  x = subset(melted_counts, Sign_of_relationship %in% "nodes")$Between_which_communities,
  y = 300,
  label = paste0(
   "nodes:\n",
   subset(melted_counts, Sign_of_relationship %in% "nodes")$Edges
  )
 ) +
 theme_JAG_presentation() +
 ylab("Edges divided by # of samples in group") +
 xlab("") +
 scale_fill_discrete(name = "",
                     labels = c("Positive edges", "Negative edges"))

dev.off()


## look at it between the environmental data
vector_organisms <- c("S16", "S18", "MPLOTU", "gp23OTU")

### Temperature to orgs ####


## find the edges either way: either in column 1 or in column 2
get_edges_between_unlike_all_orgs <-
 function (matrix_edges_between_unlike, string_env_variable) {
  temp_to_euk <-
   matrix_edges_between_unlike[str_split_fixed(matrix_edges_between_unlike [,1], "_", n =
                                                2)[,1] == "S18" |
                                str_split_fixed(matrix_edges_between_unlike [,1], "_", n = 2)[,1] == string_env_variable ,]
  temp_to_euk <-
   temp_to_euk[str_split_fixed(temp_to_euk [,2], "_", n = 2)[,1] == string_env_variable |
                str_split_fixed(temp_to_euk [,2], "_", n = 2)[,1] == "S18" ,]
  
  temp_to_bac <-
   matrix_edges_between_unlike[str_split_fixed(matrix_edges_between_unlike [,1], "_", n =
                                                2)[,1] == "S16" |
                                str_split_fixed(matrix_edges_between_unlike [,1], "_", n = 2)[,1] == string_env_variable ,]
  temp_to_bac <-
   temp_to_bac[str_split_fixed(temp_to_bac [,2], "_", n = 2)[,1] == string_env_variable |
                str_split_fixed(temp_to_bac [,2], "_", n = 2)[,1] == "S16" ,]
  
  
  temp_to_MPL <-
   matrix_edges_between_unlike[str_split_fixed(matrix_edges_between_unlike [,1], "_", n =
                                                2)[,1] == "MPLOTU" |
                                str_split_fixed(matrix_edges_between_unlike [,1], "_", n = 2)[,1] == string_env_variable ,]
  temp_to_MPL <-
   temp_to_MPL[str_split_fixed(temp_to_MPL [,2], "_", n = 2)[,1] == string_env_variable |
                str_split_fixed(temp_to_MPL [,2], "_", n = 2)[,1] == "MPLOTU" ,]
  
  
  temp_to_gp23 <-
   matrix_edges_between_unlike[str_split_fixed(matrix_edges_between_unlike [,1], "_", n =
                                                2)[,1] == "gp23OTU" |
                                str_split_fixed(matrix_edges_between_unlike [,1], "_", n = 2)[,1] == string_env_variable ,]
  temp_to_gp23 <-
   temp_to_gp23[str_split_fixed(temp_to_gp23 [,2], "_", n = 2)[,1] == string_env_variable |
                 str_split_fixed(temp_to_gp23 [,2], "_", n = 2)[,1] == "gp23OTU" ,]
  
  return(list(temp_to_euk, temp_to_bac, temp_to_MPL, temp_to_gp23))
 }

list_temp_to_orgs <-
 get_edges_between_unlike_all_orgs(edges_between_unlike, "Temperature")

temp_to_orgs <-
 rbind(list_temp_to_orgs [[1]], list_temp_to_orgs [[2]], list_temp_to_orgs [[3]], list_temp_to_orgs [[4]])

##### Salinity ####

list_sal_to_orgs <-
 get_edges_between_unlike_all_orgs(edges_between_unlike, "Salinity")

sal_to_orgs <-
 rbind(list_sal_to_orgs [[1]], list_sal_to_orgs [[2]], list_sal_to_orgs [[3]], list_sal_to_orgs [[4]])


###### Dissolved oxygen ####

list_dO_to_orgs <-
 get_edges_between_unlike_all_orgs(edges_between_unlike, "Dissolved")

dO_to_orgs <-
 rbind(list_dO_to_orgs [[1]], list_dO_to_orgs [[2]], list_dO_to_orgs [[3]], list_dO_to_orgs [[4]])


###### pH ####

list_pH_to_orgs <-
 get_edges_between_unlike_all_orgs(edges_between_unlike, "pH")

pH_to_orgs <-
 rbind(list_pH_to_orgs [[1]], list_pH_to_orgs [[2]], list_pH_to_orgs [[3]], list_pH_to_orgs [[4]])

#### Average bacterial abundance ####



get_edges_between_unlike_all_orgs_two_part_var <-
 function (matrix_edges_between_unlike, string_env_variable) {
  BA_to_euk <-
   matrix_edges_between_unlike[str_split_fixed(matrix_edges_between_unlike [,1], "_", n =
                                                2)[,1] == "S18" |
                                matrix_edges_between_unlike[,1] == string_env_variable ,]
  BA_to_euk <-
   BA_to_euk[BA_to_euk[,2] == string_env_variable |
              str_split_fixed(BA_to_euk [,2], "_", n = 2)[,1] == "S18" ,]
  
  
  BA_to_bac <-
   matrix_edges_between_unlike[str_split_fixed(matrix_edges_between_unlike [,1], "_", n =
                                                2)[,1] == "S16" |
                                matrix_edges_between_unlike[,1] == string_env_variable ,]
  BA_to_bac <-
   BA_to_bac[BA_to_bac[,2] == string_env_variable  |
              str_split_fixed(BA_to_bac [,2], "_", n = 2)[,1] == "S16" ,]
  
  
  BA_to_MPL <-
   matrix_edges_between_unlike[str_split_fixed(matrix_edges_between_unlike [,1], "_", n =
                                                2)[,1] == "MPLOTU" |
                                matrix_edges_between_unlike[,1] == string_env_variable ,]
  BA_to_MPL <-
   BA_to_MPL[BA_to_MPL[,2] == string_env_variable |
              str_split_fixed(BA_to_MPL [,2], "_", n = 2)[,1] == "MPLOTU" ,]
  
  
  BA_to_gp23 <-
   matrix_edges_between_unlike[str_split_fixed(matrix_edges_between_unlike [,1], "_", n =
                                                2)[,1] == "gp23OTU" |
                                matrix_edges_between_unlike[,1] == string_env_variable ,]
  BA_to_gp23 <-
   BA_to_gp23[BA_to_gp23[,2] == string_env_variable |
               str_split_fixed(BA_to_gp23 [,2], "_", n = 2)[,1] == "gp23OTU" ,]
  
  return(list(BA_to_euk, BA_to_bac, BA_to_MPL, BA_to_gp23))
 }

list_BA_to_orgs <-
 get_edges_between_unlike_all_orgs_two_part_var(edges_between_unlike, "Average_bacterial_abundance")

BA_to_orgs <-
 rbind(list_BA_to_orgs [[1]], list_BA_to_orgs [[2]], list_BA_to_orgs [[3]], list_BA_to_orgs [[4]])


#### Average PO4 ####


list_PO4_to_orgs <-
 get_edges_between_unlike_all_orgs_two_part_var(edges_between_unlike, "Average_PO4")

PO4_to_orgs <-
 rbind(list_PO4_to_orgs [[1]], list_PO4_to_orgs [[2]], list_PO4_to_orgs [[3]], list_PO4_to_orgs [[4]])

##### Silicate ####

list_sil_to_orgs <-
 get_edges_between_unlike_all_orgs_two_part_var(edges_between_unlike, "Average_SiO2")

sil_to_orgs <-
 rbind(list_sil_to_orgs [[1]], list_sil_to_orgs [[2]], list_sil_to_orgs [[3]], list_sil_to_orgs [[4]])



#### Average NO3_NO2 ####

list_NO3_NO2_to_orgs <-
 get_edges_between_unlike_all_orgs_two_part_var(edges_between_unlike, "Average_NO3_NO2")

NO3_NO2_to_orgs <-
 rbind(
  list_NO3_NO2_to_orgs [[1]], list_NO3_NO2_to_orgs [[2]], list_NO3_NO2_to_orgs [[3]], list_NO3_NO2_to_orgs [[4]]
 )


## count negative vs. positive and do a histogram

temp_to_orgs <- c(pos = dim(filter(temp_to_orgs, weight > 0))[1],
                  neg = dim(filter(temp_to_orgs, weight < 0))[1])

sal_to_orgs <- c(pos = dim(filter(sal_to_orgs, weight > 0))[1],
                 neg = dim(filter(sal_to_orgs, weight < 0))[1])

dO_to_orgs <- c(pos = dim(filter(dO_to_orgs, weight > 0))[1],
                neg = dim(filter(dO_to_orgs, weight < 0))[1])

pH_to_orgs <- c(pos = dim(filter(pH_to_orgs, weight > 0))[1],
                neg = dim(filter(pH_to_orgs, weight < 0))[1])

BA_to_orgs <- c(pos = dim(filter(BA_to_orgs, weight > 0))[1],
                neg = dim(filter(BA_to_orgs, weight < 0))[1])

PO4_to_orgs <- c(pos = dim(filter(PO4_to_orgs, weight > 0))[1],
                 neg = dim(filter(PO4_to_orgs, weight < 0))[1])

NO3_NO2_to_orgs <-
 c(pos = dim(filter(NO3_NO2_to_orgs, weight > 0))[1],
   neg = dim(filter(NO3_NO2_to_orgs, weight < 0))[1])

SI_to_orgs <- c(pos = dim(filter(sil_to_orgs, weight > 0))[1],
                neg = dim(filter(sil_to_orgs, weight < 0))[1])


counts_for_barplot <- rbind(
 temp_to_orgs,
 sal_to_orgs,
 dO_to_orgs,
 pH_to_orgs,
 BA_to_orgs,
 PO4_to_orgs,
 NO3_NO2_to_orgs,
 SI_to_orgs
)

## test some negs
#counts_for_barplot[,2] <- c(30,50,5,10,20,30,50,60)

melted_counts <- melt(counts_for_barplot, value.name = "Edges")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"

pdf(
 "../figures/overall_network_barplot_of_environmental_to_orgs.pdf", width = 15, height = 11
)
ggplot(
 melted_counts, aes(x = Between_which_communities, y = Edges, fill = Sign_of_relationship)
) +
 geom_bar(stat = "identity") +
 theme_JAG_presentation() +
 theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
 scale_fill_discrete(name = "",
                     labels = c("Positive edges", "Negative edges"))

dev.off()

### then do for individual commmunities!!!

## Eukaryotes

temp_to_euk <-
 c(pos = dim(filter(list_temp_to_orgs[[1]], weight > 0))[1],
   neg = dim(filter(list_temp_to_orgs[[1]], weight < 0))[1])

sal_to_euk <-
 c(pos = dim(filter(list_sal_to_orgs[[1]], weight > 0))[1],
   neg = dim(filter(list_sal_to_orgs[[1]], weight < 0))[1])

dO_to_euk <- c(pos = dim(filter(list_dO_to_orgs[[1]], weight > 0))[1],
               neg = dim(filter(list_dO_to_orgs[[1]], weight < 0))[1])

pH_to_euk <- c(pos = dim(filter(list_pH_to_orgs[[1]], weight > 0))[1],
               neg = dim(filter(list_pH_to_orgs[[1]], weight < 0))[1])

BA_to_euk <- c(pos = dim(filter(list_BA_to_orgs[[1]], weight > 0))[1],
               neg = dim(filter(list_BA_to_orgs[[1]], weight < 0))[1])

PO4_to_euk <-
 c(pos = dim(filter(list_PO4_to_orgs[[1]], weight > 0))[1],
   neg = dim(filter(list_PO4_to_orgs[[1]], weight < 0))[1])

NO3_NO2_to_euk <-
 c(pos = dim(filter(list_NO3_NO2_to_orgs[[1]], weight > 0))[1],
   neg = dim(filter(list_NO3_NO2_to_orgs[[1]], weight < 0))[1])

SI_to_euk <-
 c(pos = dim(filter(list_sil_to_orgs[[1]], weight > 0))[1],
   neg = dim(filter(list_sil_to_orgs[[1]], weight < 0))[1])


counts_for_barplot <- rbind(
 temp_to_euk,
 sal_to_euk,
 dO_to_euk,
 pH_to_euk,
 BA_to_euk,
 PO4_to_euk,
 NO3_NO2_to_euk,
 SI_to_euk
)



melted_counts <- melt(counts_for_barplot, value.name = "Edges")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"

pdf(
 "../figures/overall_network_barplot_of_environmental_to_euks.pdf", width = 15, height = 11
)
env_to_euk <-
 ggplot(
  melted_counts, aes(x = Between_which_communities, y = Edges, fill = Sign_of_relationship)
 ) +
 geom_bar(stat = "identity") +
 theme_JAG_presentation() +
 theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
 xlab(NULL) +
 scale_fill_discrete(name = "",
                     labels = c("Positive edges", "Negative edges"))
env_to_euk

dev.off()


## Bacteria
temp_to_bac <-
 c(pos = dim(filter(list_temp_to_orgs[[2]], weight > 0))[1],
   neg = dim(filter(list_temp_to_orgs[[2]], weight < 0))[1])

sal_to_bac <-
 c(pos = dim(filter(list_sal_to_orgs[[2]], weight > 0))[1],
   neg = dim(filter(list_sal_to_orgs[[2]], weight < 0))[1])

dO_to_bac <- c(pos = dim(filter(list_dO_to_orgs[[2]], weight > 0))[1],
               neg = dim(filter(list_dO_to_orgs[[2]], weight < 0))[1])

pH_to_bac <- c(pos = dim(filter(list_pH_to_orgs[[2]], weight > 0))[1],
               neg = dim(filter(list_pH_to_orgs[[2]], weight < 0))[1])

BA_to_bac <- c(pos = dim(filter(list_BA_to_orgs[[2]], weight > 0))[1],
               neg = dim(filter(list_BA_to_orgs[[2]], weight < 0))[1])

PO4_to_bac <-
 c(pos = dim(filter(list_PO4_to_orgs[[2]], weight > 0))[1],
   neg = dim(filter(list_PO4_to_orgs[[2]], weight < 0))[1])

NO3_NO2_to_bac <-
 c(pos = dim(filter(list_NO3_NO2_to_orgs[[2]], weight > 0))[1],
   neg = dim(filter(list_NO3_NO2_to_orgs[[2]], weight < 0))[1])

SI_to_bac <-
 c(pos = dim(filter(list_sil_to_orgs[[2]], weight > 0))[1],
   neg = dim(filter(list_sil_to_orgs[[2]], weight < 0))[1])


counts_for_barplot <- rbind(
 temp_to_bac,
 sal_to_bac,
 dO_to_bac,
 pH_to_bac,
 BA_to_bac,
 PO4_to_bac,
 NO3_NO2_to_bac,
 SI_to_bac
)



melted_counts <- melt(counts_for_barplot, value.name = "Edges")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"

pdf(
 "../figures/overall_network_barplot_of_environmental_to_bacs.pdf", width = 15, height = 11
)
env_to_bac <-
 ggplot(
  melted_counts, aes(x = Between_which_communities, y = Edges, fill = Sign_of_relationship)
 ) +
 geom_bar(stat = "identity") +
 theme_JAG_presentation() +
 theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
 xlab(NULL) +
 scale_fill_discrete(name = "",
                     labels = c("Positive edges", "Negative edges"))
env_to_bac
dev.off()


## MPL

temp_to_MPL <-
 c(pos = dim(filter(list_temp_to_orgs[[3]], weight > 0))[1],
   neg = dim(filter(list_temp_to_orgs[[3]], weight < 0))[1])

sal_to_MPL <-
 c(pos = dim(filter(list_sal_to_orgs[[3]], weight > 0))[1],
   neg = dim(filter(list_sal_to_orgs[[3]], weight < 0))[1])

dO_to_MPL <- c(pos = dim(filter(list_dO_to_orgs[[3]], weight > 0))[1],
               neg = dim(filter(list_dO_to_orgs[[3]], weight < 0))[1])

pH_to_MPL <- c(pos = dim(filter(list_pH_to_orgs[[3]], weight > 0))[1],
               neg = dim(filter(list_pH_to_orgs[[3]], weight < 0))[1])

BA_to_MPL <- c(pos = dim(filter(list_BA_to_orgs[[3]], weight > 0))[1],
               neg = dim(filter(list_BA_to_orgs[[3]], weight < 0))[1])

PO4_to_MPL <-
 c(pos = dim(filter(list_PO4_to_orgs[[3]], weight > 0))[1],
   neg = dim(filter(list_PO4_to_orgs[[3]], weight < 0))[1])

NO3_NO2_to_MPL <-
 c(pos = dim(filter(list_NO3_NO2_to_orgs[[3]], weight > 0))[1],
   neg = dim(filter(list_NO3_NO2_to_orgs[[3]], weight < 0))[1])

SI_to_MPL <-
 c(pos = dim(filter(list_sil_to_orgs[[3]], weight > 0))[1],
   neg = dim(filter(list_sil_to_orgs[[3]], weight < 0))[1])



counts_for_barplot <- rbind(
 temp_to_MPL,
 sal_to_MPL,
 dO_to_MPL,
 pH_to_MPL,
 BA_to_MPL,
 PO4_to_MPL,
 NO3_NO2_to_MPL,
 SI_to_MPL
)



melted_counts <- melt(counts_for_barplot, value.name = "Edges")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"

pdf(
 "../figures/overall_network_barplot_of_environmental_to_MPLs.pdf", width = 15, height = 11
)
env_to_MPL <-
 ggplot(
  melted_counts,
  aes(x = Between_which_communities,
      y = Edges,
      fill = Sign_of_relationship)
 ) +
 geom_bar(stat = "identity") +
 theme_JAG_presentation() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
 xlab(NULL) +
 scale_fill_discrete(name = "",
                     labels = c("Positive edges",
                                "Negative edges"))
env_to_MPL

dev.off()


## gp23

temp_to_gp23 <-
 c(pos = dim(filter(list_temp_to_orgs[[4]], weight > 0))[1],
   neg = dim(filter(list_temp_to_orgs[[4]], weight < 0))[1])

sal_to_gp23 <-
 c(pos = dim(filter(list_sal_to_orgs[[4]], weight > 0))[1],
   neg = dim(filter(list_sal_to_orgs[[4]], weight < 0))[1])

dO_to_gp23 <-
 c(pos = dim(filter(list_dO_to_orgs[[4]], weight > 0))[1],
   neg = dim(filter(list_dO_to_orgs[[4]], weight < 0))[1])

pH_to_gp23 <-
 c(pos = dim(filter(list_pH_to_orgs[[4]], weight > 0))[1],
   neg = dim(filter(list_pH_to_orgs[[4]], weight < 0))[1])

BA_to_gp23 <-
 c(pos = dim(filter(list_BA_to_orgs[[4]], weight > 0))[1],
   neg = dim(filter(list_BA_to_orgs[[4]], weight < 0))[1])

PO4_to_gp23 <-
 c(pos = dim(filter(list_PO4_to_orgs[[4]], weight > 0))[1],
   neg = dim(filter(list_PO4_to_orgs[[4]], weight < 0))[1])

NO3_NO2_to_gp23 <-
 c(pos = dim(filter(list_NO3_NO2_to_orgs[[4]], weight > 0))[1],
   neg = dim(filter(list_NO3_NO2_to_orgs[[4]], weight < 0))[1])

SI_to_gp23 <-
 c(pos = dim(filter(list_sil_to_orgs[[4]], weight > 0))[1],
   neg = dim(filter(list_sil_to_orgs[[4]], weight < 0))[1])


counts_for_barplot <- rbind(
 temp_to_gp23,
 sal_to_gp23,
 dO_to_gp23,
 pH_to_gp23,
 BA_to_gp23,
 PO4_to_gp23,
 NO3_NO2_to_gp23,
 SI_to_gp23
)



melted_counts <- melt(counts_for_barplot, value.name = "Edges")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"

pdf(
 "../figures/overall_network_barplot_of_environmental_to_gp23s.pdf", width = 15, height = 11
)
env_to_gp23 <-
 ggplot(
  melted_counts, aes(x = Between_which_communities, y = Edges, fill = Sign_of_relationship)
 ) +
 geom_bar(stat = "identity") +
 theme_JAG_presentation() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
 xlab(NULL) +
 scale_fill_discrete(name = "",
                     labels = c("Positive edges", "Negative edges"))
env_to_gp23
dev.off()



# arrange the three plots in a single row
prow <- plot_grid(
 env_to_gp23 + theme(legend.position = "none"),
 env_to_MPL + theme(legend.position = "none"),
 env_to_euk + theme(legend.position = "none"),
 env_to_bac + theme(legend.position = "none"), labels = c("A", "B", "C", "D")
)
prow
# extract the legend from one of the plots
# (clearly the whole thing only makes sense if all plots
# have the same legend, so we can arbitrarily pick one.)
grobs <- ggplotGrob(env_to_gp23)$grobs
legend <-
 grobs[[which(sapply(grobs, function(x)
  x$name) == "guide-box")]]

# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
p <- plot_grid(prow, legend, rel_widths = c(3, .3))
p


pdf(
 "../figures/overall_network_barplot_of_environmental_to_orgs_labelled.pdf", width = 15, height = 11
)
p
dev.off()
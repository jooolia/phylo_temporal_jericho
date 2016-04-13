## histogram of degree distribution

library(igraph)
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)


graph_file <- read.graph("../results/LSA_tables/network_LSA_overall.graphml",format="graphml")


Jericho_data <- read.csv("../results/Jericho_env_data_mean_imputed.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

# ## degree distribution
# 
# name_degree_date <-aV(graph_file)$degree)
# names(degree_df) <- "degree"
# 
# ggplot(degree_df, aes(x=degree))+geom_histogram(binwidth=2)

## so we have date, just need node and degreee

generate_network_stats <- function (graph_file, graph_name, type) {
 
 degree <-V(graph_file)$degree
 name <-V(graph_file)$name
 amplicon_type <-V(graph_file)$type
 betweenness <- V(graph_file)$betweenness
 table_of_stats <- cbind(type,
                         graph_name,
                         degree,
                         name, 
                         amplicon_type,
                         betweenness
)
 
 return(table_of_stats)

}


all_network_stats <- c()
graphml_for_cytoscape <- list.files(path = "../results/LSA_tables", pattern = "time_LSA_overallt.*.graphml", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


overall_table_dates <- c()
## cytoscape needs to be open
for (i in 1:length(graphml_for_cytoscape)){
 graph_file <- read.graph(paste0("../results/LSA_tables/",graphml_for_cytoscape[i]),format="graphml")
 print(graphml_for_cytoscape[i])
 date <- gsub("time_LSA_overallt_([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})\\.graphml","\\1-\\2-\\3", graphml_for_cytoscape[i])
 stats_overall <- generate_network_stats(graph_file, paste0(date), "overall")
 overall_table_dates <- rbind(overall_table_dates, stats_overall)
 
}



overall_network_stats <- as.data.frame(overall_table_dates)
overall_network_stats$Temperature <- Jericho_data$Temperature_YSI[match(as.Date(overall_network_stats$graph_name), Jericho_data$Date)]
overall_network_stats$pH <- Jericho_data$pH[match(as.Date(overall_network_stats$graph_name), Jericho_data$Date)]
overall_network_stats$VA <- Jericho_data$Average_viral_abundance[match(as.Date(overall_network_stats$graph_name), Jericho_data$Date)] 

## Temperature
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=Temperature))+geom_point()+geom_smooth()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=Temperature))+geom_point(alpha = 1/10)+geom_smooth()+facet_wrap(~amplicon_type)

ggplot(overall_network_stats, aes(y=as.numeric(degree), x=Temperature))+geom_density2d()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=Temperature))+geom_bin2d()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=Temperature))+geom_hex()+facet_wrap(~amplicon_type)

## with betweeness
ggplot(overall_network_stats, aes(y=as.numeric(betweenness), x=Temperature))+geom_point()+geom_smooth()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(betweenness), x=Temperature))+geom_point(alpha = 1/10)+geom_smooth()+facet_wrap(~amplicon_type)

ggplot(overall_network_stats, aes(y=as.numeric(betweenness), x=Temperature))+geom_density2d()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(betweenness), x=Temperature))+geom_bin2d()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(betweenness), x=Temperature))+geom_hex()+facet_wrap(~amplicon_type)




## pH
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=pH))+geom_point()+geom_smooth()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=pH))+geom_point(alpha = 1/10)+geom_smooth()+facet_wrap(~amplicon_type)

ggplot(overall_network_stats, aes(y=as.numeric(degree), x=pH))+geom_density2d()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=pH))+geom_bin2d()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=pH))+geom_hex()+facet_wrap(~amplicon_type)


## VA
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=VA))+geom_point()+geom_smooth()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=VA))+geom_point(alVAa = 1/10)+geom_smooth()+facet_wrap(~amplicon_type)

ggplot(overall_network_stats, aes(y=as.numeric(degree), x=VA))+geom_density2d()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=VA))+geom_bin2d()+facet_wrap(~amplicon_type)
ggplot(overall_network_stats, aes(y=as.numeric(degree), x=VA))+geom_hex()+facet_wrap(~amplicon_type)


library(igraph)
library(ggplot2)
library(reshape2)
library(scales)
library(cowplot)
library(RColorBrewer)

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

# adding in the seasons and the spring bloom for the ggplots
season_line <- geom_vline(xintercept = as.numeric(c(as.Date("2010-03-22"), 
                                                    as.Date("2010-06-22"), 
                                                    as.Date("2010-09-22"),
                                                    as.Date("2010-12-22"),
                                                    as.Date("2011-03-22"),
                                                    as.Date("2011-06-22"))),
                          colour="grey",
                          size=1.5)

mid_summer_x <- as.Date("2010-09-22")-as.numeric(as.Date("2010-09-22")-as.Date("2010-06-22"))/2
mid_fall_x <- as.Date("2010-12-22")-as.numeric(as.Date("2010-12-22")-as.Date("2010-09-22"))/2
mid_winter_x <- as.Date("2011-03-22")-as.numeric(as.Date("2011-03-22")-as.Date("2010-12-22"))/2
mid_spring_x <- as.Date("2011-06-22")-as.numeric(as.Date("2011-06-22")-as.Date("2011-03-22"))/2

season_text <- annotate("text",x=c(mid_summer_x, mid_fall_x, mid_winter_x, mid_spring_x), y=1.05, label=c("Summer", "Fall", "Winter", "Spring"))

spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=1)
date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")))

### so want to read in the graphs make a dataframe of all of the different traits

generate_network_stats <- function (graph_file, graph_name, type) {
 density <- graph_file$density
 modularity <- graph_file$modularity
 edges <- ecount(graph_file)
 nodes <- vcount(graph_file)
 diameter <- diameter(graph_file) ## largest distance between any two nodes.
 farthest_nodes <- farthest.nodes(graph_file)
 largest_cliques <-largest.cliques(graph_file)
 
 average_nearest_neighbour_degree <- mean(graph.knn(graph_file)$knn) #average nearest neighbour degree -average of this
 average_degree <- mean(V(graph_file)$degree)
 median_degree <- median(V(graph_file)$degree)
 min_degree <- min(V(graph_file)$degree)
 max_degree <- max(V(graph_file)$degree)
 average_closeness <- mean(V(graph_file)$closeness)
 average_betweeness <- mean(V(graph_file)$betweenness)
 average_page_rank <- mean(V(graph_file)$page_rank) ## google page rank algorithm
 
 average_path_length <- average.path.length(graph_file) # average distance between two nodes in a graph. 
 
 ## positive edges
 positive_edges <- length(E(graph_file)$weight_LS[E(graph_file)$weight_LS > 0])
 percent_positive_edges <- (positive_edges/ecount(graph_file))*100
 
 ## negative
 negative_edges <- length(E(graph_file)$weight_LS[E(graph_file)$weight_LS < 0])
 percent_negative_edges <- (negative_edges/ecount(graph_file))*100
 
 ## delay in relationship
 delay_0 <- length(E(graph_file)$delay[E(graph_file)$delay == 0])
 percent_delay_0 <- (delay_0/ecount(graph_file))*100
 delay_pos <- length(E(graph_file)$delay[E(graph_file)$delay > 0])
 percent_delay_pos <- (delay_pos/ecount(graph_file))*100
 delay_neg <- length(E(graph_file)$delay[E(graph_file)$delay < 0])
 percent_delay_neg <- (delay_neg/ecount(graph_file))*100
 
 ## do I have to generate random graphs here too: 
g1 <- erdos.renyi.game(nodes, edges,type="gnm") # random nodes and then edges
#plot(g1) 
## checking these random graphs
random_edges <- ecount(g1)
random_nodes <- vcount(g1)

random_diameter <- diameter(g1)
wc <- cluster_walktrap(g1)
V(g1)$membership <- wc$membership
random_modularity <- modularity(wc)
random_density = graph.density(g1)
random_average_path_length <- average.path.length(g1) 
random_average_degree <- mean(degree(g1))
random_median_degree <- median(degree(g1))

g2 <- barabasi.game(nodes, m=edges/nodes) # scale free
#plot(g2)
scale_free_edges <- ecount(g2)
scale_free_nodes <- vcount(g2)
scale_free_diameter <- diameter(g2)

scale_free_diameter <- diameter(g2)
wc <- cluster_walktrap(g2)
V(g2)$membership <- wc$membership
scale_free_modularity <- modularity(wc)
scale_free_density = graph.density(g2)
scale_free_average_path_length <- average.path.length(g2) 
scale_free_average_degree <- mean(degree(g2))
scale_free_median_degree <- median(degree(g2))
 
 table_of_stats <- cbind(type,
                         graph_name,
                         #deparse(substitute(graph_name)),
                         nodes,
                         edges,
                         diameter,
                         density,
                         modularity,
                         random_density,
                         random_modularity,
                         random_average_path_length,
                         random_average_degree,
                         random_median_degree,
                         # random_edges, 
                         # random_nodes,
                         random_diameter,
                         # scale_free_edges,
                         # scale_free_nodes,
                         scale_free_diameter,
                         scale_free_modularity,
                         scale_free_density,
                         scale_free_average_path_length,
                         scale_free_average_degree,
                         scale_free_median_degree,
                         # largest_cliques,
                         average_nearest_neighbour_degree,
                         average_degree,
                         median_degree,
                         # min_degree,
                         max_degree,
                         average_path_length,
                         average_closeness,
                         average_betweeness,
                         # positive_edges,
                         percent_positive_edges,
                         # negative_edges,
                         percent_negative_edges,
                         # delay_0,
                         percent_delay_0,
                         # delay_pos,
                         percent_delay_pos,
                         # delay_neg,
                         percent_delay_neg)
 
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

melted_network_stats <- reshape2::melt(overall_network_stats, id=c("graph_name", "type"))
all_network_stats <-rbind(all_network_stats , melted_network_stats)

ggplot(overall_network_stats, aes(x=as.Date(graph_name),
                                  y=as.numeric(average_degree)))+
 geom_line()

ggplot(melted_network_stats, aes(x=as.Date(graph_name), y=as.numeric(value)))+
 geom_line()+
 facet_wrap(~ variable, scales="free_y")



plot_stats_with_generated_networks <- function (melted_network_stats) {
  ## just diameter
 x_label <- "Date"
 
 nodes_edges <- c("nodes","edges" )
 nodes_edges_plot <- ggplot(subset(melted_network_stats, variable %in% nodes_edges),
                         aes(x=as.Date(graph_name),
                             y=as.numeric(value), 
                             group = variable))+
  season_line +
  #season_text+
  spring_bloom_line+
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank())+
  # geom_line(aes(colour=variable))+
  geom_line(aes(linetype=variable)) +
  scale_linetype_manual(values=c(nodes="solid",
                                 edges="longdash"),
                        name="",
                        breaks=nodes_edges,
                        labels=c("Nodes", "Edges"))+
  xlab(NULL)+
  ylab("Count")
 

  diameters <- c("diameter","random_diameter", "scale_free_diameter" )
  diameter_plot <- ggplot(subset(melted_network_stats, variable %in% diameters),
         aes(x=as.Date(graph_name),
             y=as.numeric(value), 
             group = variable))+
   season_line +
   #season_text+
   spring_bloom_line+
   date_scaling+
   theme_JAG_presentation()+
   theme(panel.grid.major = element_blank())+
   theme(axis.text.x = element_blank())+
   # geom_line(aes(colour=variable))+
   geom_line(aes(linetype=variable)) +
   scale_linetype_manual(values=c(diameter="solid",
                                  random_diameter="dotdash",
                                  scale_free_diameter="dotted"),
                         name="",
                         breaks=diameters,
                         labels=c("Network", "RN", "SF"))+
  xlab(NULL)+
   ylab("Diameter")
 
  modularitys <- c("modularity","random_modularity", "scale_free_modularity" )
  modularity_plot <- ggplot(subset(melted_network_stats, variable %in% modularitys),
                          aes(x=as.Date(graph_name),
                              y=as.numeric(value), 
                              group = variable))+
   season_line +
   #season_text+
   spring_bloom_line+
   date_scaling+
   theme_JAG_presentation()+
   theme(panel.grid.major = element_blank())+
   theme(axis.text.x = element_blank())+
   # geom_line(aes(colour=variable))+
   geom_line(aes(linetype=variable)) +
   scale_linetype_manual(values=c(modularity="solid",
                                  random_modularity="dotdash",
                                  scale_free_modularity="dotted"),
                         name="",
                         breaks=modularitys,
                         labels=c("Network", "RN", "SF"))+
   xlab(NULL)+
   ylab("Modularity")
   
  density <- c("density","random_density", "scale_free_density" )
  density_plot <- ggplot(subset(melted_network_stats, variable %in% density),
         aes(x=as.Date(graph_name),
             y=as.numeric(value),
             group = variable))+
   season_line +
   # season_text+
   spring_bloom_line+
   date_scaling+
   theme_JAG_presentation()+
   theme(panel.grid.major = element_blank())+
   theme(axis.text.x = element_blank())+
   # geom_line(aes(colour=variable))+
   geom_line(aes(linetype=variable)) +
   scale_linetype_manual(values=c(density="solid",random_density="dotdash", scale_free_density="dotted"),
                         name="",
                         breaks=density,
                         labels=c("Network", "RN", "SF"))+
   xlab(NULL)+
   ylab("Density")
  
  
  path_length <- c("average_path_length",
                   "random_average_path_length",
                   "scale_free_average_path_length")
  path_length_plot <- ggplot(subset(melted_network_stats, variable %in% path_length),
         aes(x=as.Date(graph_name),
             y=as.numeric(value),
             group = variable))+
   season_line +
   # season_text+
   spring_bloom_line+
   date_scaling+
   theme_JAG_presentation()+
   theme(panel.grid.major = element_blank())+
   theme(axis.text.x = element_blank())+
   #geom_line(aes(colour=variable))+
   geom_line(aes(linetype=variable)) +
   scale_linetype_manual(values=c(average_path_length="solid",random_average_path_length="dotdash", scale_free_average_path_length="dotted"),
                         name="",
                         breaks=path_length,
                         labels=c("Network", "RN", "SF"))+
   xlab(NULL)+
   ylab("Average \npath length")
  
  ## average degree from network and random average degree are the same.Maybe median in better. 
  median_degree <- c("median_degree", 
                      "random_median_degree",
                      "scale_free_median_degree")
  median_degree_plot <- ggplot(subset(melted_network_stats,
                                       variable %in% median_degree),
         aes(x=as.Date(graph_name), 
             y=as.numeric(value), 
             group = variable))+
   season_line +
   # season_text+
   spring_bloom_line+
   date_scaling+
   theme_JAG_presentation()+
   theme(panel.grid.major = element_blank())+
   geom_line(aes(linetype=variable)) +
   scale_linetype_manual(values=c(median_degree="solid",random_median_degree="dotdash", scale_free_median_degree="dotted"),
                         name="",
                         breaks=median_degree,
                         labels=c("Network", "RN", "SF"))+
   xlab(x_label)+
   ylab("Median \ndegree")

  
  
  
#   prow <- plot_grid( nodes_edges_plot,
#                      diameter_plot+ theme(legend.position="none"),
#                      density_plot+ theme(legend.position="none"),
#                      modularity_plot+ theme(legend.position="none"),
#                      path_length_plot+ theme(legend.position="none"),
#                      median_degree_plot+ theme(legend.position="none"),ncol=1, align="v"
#   )
#   
#   prow
#   # extract the legend from one of the plots
#   # (clearly the whole thing only makes sense if all plots
#   # have the same legend, so we can arbitrarily pick one.)
#   grobs <- ggplotGrob(nodes_edges_plot)$grobs
#   legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
#   
#   # add the legend to the row we made earlier. Give it one-third of the width
#   # of one plot (via rel_widths).
#   p <- plot_grid( prow, legend, rel_widths = c(3, .3))
#   p
  
  
  plot_grid(nodes_edges_plot,
            diameter_plot,
            density_plot,
            modularity_plot,
            path_length_plot,
            median_degree_plot,ncol=1, align="v")
  
}


pdf(file="../figures/plots_network_stats_overall%01d.pdf", height=15, width=11, 
    onefile=FALSE)
plot_stats_with_generated_networks(melted_network_stats)
dev.off()


## bac to gp23 

## need to think about this.
graphml_for_cytoscape <- list.files(path = "../results/LSA_tables", pattern = "time_LSA_bac_to_gp23t.*.graphml", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


overall_table_dates <- c()
## cytoscape needs to be open
for (i in 1:length(graphml_for_cytoscape)){
 graph_file <- read.graph(paste0("../results/LSA_tables/",graphml_for_cytoscape[i]),format="graphml")
 print(graphml_for_cytoscape[i])
 date <- gsub("time_LSA_bac_to_gp23t_([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})\\.graphml","\\1-\\2-\\3", graphml_for_cytoscape[i])
 stats_overall <- generate_network_stats(graph_file, paste0(date), "bac_to_gp23")
 overall_table_dates <- rbind(overall_table_dates, stats_overall)
}


overall_network_stats <- as.data.frame(overall_table_dates)
melted_network_stats <- reshape2::melt(overall_network_stats,id=c("graph_name", "type"))
all_network_stats <-rbind(all_network_stats , melted_network_stats)

ggplot(melted_network_stats, aes(x=as.Date(graph_name), y=as.numeric(value)))+
 geom_line()+
 facet_wrap(~ variable, scales="free_y")

plot_stats_with_generated_networks(melted_network_stats)
pdf(file="../figures/plots_network_stats_bac_to_gp23%01d.pdf", height=11, width=11,  
    onefile=FALSE)
plot_stats_with_generated_networks(melted_network_stats)
dev.off()


## only bac

## need to think about this.
graphml_for_cytoscape <- list.files(path = "../results/LSA_tables", pattern = "time_LSA_bac_to_bact.*.graphml", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

overall_table_dates <- c()
## cytoscape needs to be open
for (i in 1:length(graphml_for_cytoscape)){
 graph_file <- read.graph(paste0("../results/LSA_tables/",graphml_for_cytoscape[i]),format="graphml")
 print(graphml_for_cytoscape[i])
 date <- gsub("time_LSA_bac_to_bact_([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})\\.graphml","\\1-\\2-\\3", graphml_for_cytoscape[i])
 stats_overall <- generate_network_stats(graph_file, paste0(date), "bac_to_bac")
 overall_table_dates <- rbind(overall_table_dates, stats_overall)
}

overall_network_stats <- as.data.frame(overall_table_dates)
melted_network_stats <- reshape2::melt(overall_network_stats,id=c("graph_name", "type"))
all_network_stats <-rbind(all_network_stats , melted_network_stats)

ggplot(melted_network_stats, aes(x=as.Date(graph_name), y=as.numeric(value)))+
 geom_line()+
 facet_wrap(~ variable, scales="free_y")

plot_stats_with_generated_networks(melted_network_stats)
pdf(file="../figures/plots_network_stats_bac_to_bac%01d.pdf", height=11, width=11,
    onefile = FALSE)
plot_stats_with_generated_networks(melted_network_stats)
dev.off()

## euk to euk
graphml_for_cytoscape <- list.files(path = "../results/LSA_tables", pattern = "time_LSA_euk_to_eukt.*.graphml", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

overall_table_dates <- c()
## cytoscape needs to be open
for (i in 1:length(graphml_for_cytoscape)){
 graph_file <- read.graph(paste0("../results/LSA_tables/",graphml_for_cytoscape[i]),format="graphml")
 print(graphml_for_cytoscape[i])
 date <- gsub("time_LSA_euk_to_eukt_([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})\\.graphml","\\1-\\2-\\3", graphml_for_cytoscape[i])
 stats_overall <- generate_network_stats(graph_file, paste0(date), "euk_to_euk")
 overall_table_dates <- rbind(overall_table_dates, stats_overall)
}

overall_network_stats <- as.data.frame(overall_table_dates)
melted_network_stats <- reshape2::melt(overall_network_stats, id=c("graph_name", "type"))
all_network_stats <-rbind(all_network_stats , melted_network_stats)

ggplot(melted_network_stats, aes(x=as.Date(graph_name), y=as.numeric(value)))+
 geom_line()+
 facet_wrap(~ variable, scales="free_y")

plot_stats_with_generated_networks(melted_network_stats)
pdf(file="../figures/plots_network_stats_euk_to_euk%01d.pdf", height=11, width=11, 
    onefile=FALSE)
plot_stats_with_generated_networks(melted_network_stats)
dev.off()

## MPl MPL
graphml_for_cytoscape <- list.files(path = "../results/LSA_tables", pattern = "time_LSA_MPL_to_MPLt.*.graphml", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

overall_table_dates <- c()
## cytoscape needs to be open
for (i in 1:length(graphml_for_cytoscape)){
 graph_file <- read.graph(paste0("../results/LSA_tables/",graphml_for_cytoscape[i]),format="graphml")
 print(graphml_for_cytoscape[i])
 date <- gsub("time_LSA_MPL_to_MPLt_([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})\\.graphml","\\1-\\2-\\3", graphml_for_cytoscape[i])
 stats_overall <- generate_network_stats(graph_file, paste0(date), "MPL_to_MPL")
 overall_table_dates <- rbind(overall_table_dates, stats_overall)
}

overall_network_stats <- as.data.frame(overall_table_dates)
melted_network_stats <- reshape2::melt(overall_network_stats, id=c("graph_name", "type"))
all_network_stats <-rbind(all_network_stats , melted_network_stats)

ggplot(melted_network_stats, aes(x=as.Date(graph_name), y=as.numeric(value)))+
 geom_line()+
 facet_wrap(~ variable, scales="free_y")

plot_stats_with_generated_networks(melted_network_stats)
pdf(file="../figures/plots_network_stats_MPL_to_MPL%01d.pdf",height=11, width=11, 
    onefile=FALSE)
plot_stats_with_generated_networks(melted_network_stats)
dev.off()

## gp23 to gp23
graphml_for_cytoscape <- list.files(path = "../results/LSA_tables", pattern = "time_LSA_gp23_to_gp23t.*.graphml", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

overall_table_dates <- c()
## cytoscape needs to be open
for (i in 1:length(graphml_for_cytoscape)){
 graph_file <- read.graph(paste0("../results/LSA_tables/",graphml_for_cytoscape[i]),format="graphml")
 print(graphml_for_cytoscape[i])
 date <- gsub("time_LSA_gp23_to_gp23t_([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})\\.graphml","\\1-\\2-\\3", graphml_for_cytoscape[i])
 stats_overall <- generate_network_stats(graph_file, paste0(date), "gp23_to_gp23")
 overall_table_dates <- rbind(overall_table_dates, stats_overall)
}

overall_network_stats <- as.data.frame(overall_table_dates)
melted_network_stats <- reshape2::melt(overall_network_stats, id=c("graph_name", "type"))
all_network_stats <-rbind(all_network_stats , melted_network_stats)

ggplot(melted_network_stats, aes(x=as.Date(graph_name), y=as.numeric(value)))+
 geom_line()+
 facet_wrap(~ variable, scales="free_y")

plot_stats_with_generated_networks(melted_network_stats)
pdf(file="../figures/plots_network_stats_gp23_to_gp23%01d.pdf", height=11, width=11, 
    onefile=FALSE)
plot_stats_with_generated_networks(melted_network_stats)
dev.off()


## bac to euk
graphml_for_cytoscape <- list.files(path = "../results/LSA_tables", pattern = "time_LSA_euk_to_bact.*.graphml", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

overall_table_dates <- c()
## cytoscape needs to be open
for (i in 1:length(graphml_for_cytoscape)){
 graph_file <- read.graph(paste0("../results/LSA_tables/",graphml_for_cytoscape[i]),format="graphml")
 print(graphml_for_cytoscape[i])
 date <- gsub("time_LSA_euk_to_bact_([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})\\.graphml","\\1-\\2-\\3", graphml_for_cytoscape[i])
 stats_overall <- generate_network_stats(graph_file, paste0(date), "euk_to_bac")
 overall_table_dates <- rbind(overall_table_dates, stats_overall)
}

overall_network_stats <- as.data.frame(overall_table_dates)
melted_network_stats <- reshape2::melt(overall_network_stats, id=c("graph_name", "type"))
all_network_stats <-rbind(all_network_stats , melted_network_stats)

ggplot(melted_network_stats, aes(x=as.Date(graph_name), y=as.numeric(value)))+
 geom_line()+
 facet_wrap(~ variable, scales="free_y")

plot_stats_with_generated_networks(melted_network_stats)
pdf(file="../figures/plots_network_stats_euk_to_bac%01d.pdf", height=11, width=11, 
    onefile=FALSE)
plot_stats_with_generated_networks(melted_network_stats)
dev.off()


## mpl to euk
graphml_for_cytoscape <- list.files(path = "../results/LSA_tables", pattern = "time_LSA_euk_to_MPLt.*.graphml", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

overall_table_dates <- c()
## cytoscape needs to be open
for (i in 1:length(graphml_for_cytoscape)){
 graph_file <- read.graph(paste0("../results/LSA_tables/",graphml_for_cytoscape[i]),format="graphml")
 print(graphml_for_cytoscape[i])
 date <- gsub("time_LSA_euk_to_MPLt_([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})\\.graphml","\\1-\\2-\\3", graphml_for_cytoscape[i])
 stats_overall <- generate_network_stats(graph_file, paste0(date), "euk_to_MPL")
 overall_table_dates <- rbind(overall_table_dates, stats_overall)
}

overall_network_stats <- as.data.frame(overall_table_dates)
melted_network_stats <- reshape2::melt(overall_network_stats, id=c("graph_name", "type"))
all_network_stats <-rbind(all_network_stats , melted_network_stats)

ggplot(melted_network_stats, aes(x=as.Date(graph_name), y=as.numeric(value)))+
 geom_line()+
 facet_wrap(~ variable, scales="free_y")


pdf(file="../figures/plots_network_stats_euk_to_MPL%01d.pdf",
    height=11,
    width=11,
    onefile=FALSE)
plot_stats_with_generated_networks(melted_network_stats)
dev.off()


## want to pull all together into one table and then do a few overall plots.
## make colours for different 
colour_names <- c("overall","bac_to_gp23","euk_to_bac","euk_to_MPL", "overall", "bac_to_bac","euk_to_euk","MPL_to_MPL","gp23_to_gp23")
length(colour_names)
## want overall to be black
palette_plots_colours <- c("black", brewer.pal(length(colour_names)-1, "Set1"))
names(palette_plots_colours) <- colour_names

plot_stats_together <- function (melted_network_stats) {
 ## just diameter
 x_label <- "Date"
 
 nodes_edges_plot <- ggplot(subset(melted_network_stats, variable %in% "nodes"),
                            aes(x=as.Date(graph_name),
                                y=as.numeric(value), 
                                group = type))+
  season_line +
  #season_text+
  spring_bloom_line+
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())+
  # geom_line(aes(colour=variable))+
  geom_line(aes(colour=type), size =1.5, alpha= 0.8) +
  scale_colour_manual(values = palette_plots_colours)+
  theme(axis.text.x = element_blank(),
        legend.text = element_text(size = 20))+
#   scale_linetype_manual(values=c(nodes="solid",
#                                  edges="longdash"),
#                         name="",
#                         breaks=nodes_edges,
#                         labels=c("Nodes", "Edges"))+

  xlab(NULL)+
  ylab("Node \nCount\n")+ theme(legend.title = element_blank())
 

edges_plot <- ggplot(subset(melted_network_stats, variable %in% "edges"),
                            aes(x=as.Date(graph_name),
                                y=as.numeric(value), 
                                group = type))+
  season_line +
  #season_text+
  spring_bloom_line+
  date_scaling+
  theme_JAG_presentation()+
 theme(panel.grid.major = element_blank())+
 theme(axis.text.x = element_blank())+
  # geom_line(aes(colour=variable))+
  geom_line(aes(colour=type), size =1.5, alpha= 0.8) +
 scale_colour_manual(values = palette_plots_colours)+
  #   scale_linetype_manual(values=c(nodes="solid",
  #                                  edges="longdash"),
  #                         name="",
  #                         breaks=nodes_edges,
  #                         labels=c("Nodes", "Edges"))+
  xlab(NULL)+
  ylab("Edge \nCount\n")+ theme(legend.title = element_blank())
 
 diameter_plot <- ggplot(subset(melted_network_stats, variable == "diameter"),
                         aes(x=as.Date(graph_name),
                             y=as.numeric(value), 
                             group = type))+
  season_line +
  #season_text+
  spring_bloom_line+
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank())+
  # geom_line(aes(colour=variable))+
  geom_line(aes(colour=type), size =1.5, alpha= 0.8) +
  scale_colour_manual(values = palette_plots_colours)+
  xlab(NULL)+
  ylab("Diameter\n\n")+ theme(legend.title = element_blank())
 
 modularity_plot <- ggplot(subset(melted_network_stats, variable == "modularity"),
                           aes(x=as.Date(graph_name),
                               y=as.numeric(value), 
                               group = type))+
  season_line +
  #season_text+
  spring_bloom_line+
  date_scaling+
  theme_JAG_presentation()+
  theme(axis.text.x = element_blank())+
  # geom_line(aes(colour=variable))+
  geom_line(aes(colour=type), size =1.5, alpha= 0.8) +
  scale_colour_manual(values = palette_plots_colours)+
  xlab(NULL)+
  ylab("Modularity\n\n")+ theme(legend.title = element_blank())
 

 density_plot <- ggplot(subset(melted_network_stats, variable == "density"),
                        aes(x=as.Date(graph_name),
                            y=as.numeric(value),
                            group = type))+
  season_line +
  # season_text+
  spring_bloom_line+
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank())+
  # geom_line(aes(colour=variable))+
  geom_line(aes(colour=type), size =1.5, alpha= 0.8) +
  scale_colour_manual(values = palette_plots_colours)+
  xlab(NULL)+
  ylab("Density\n\n")+ theme(legend.title = element_blank())
 
 path_length_plot <- ggplot(subset(melted_network_stats, variable == "average_path_length"),
                            aes(x=as.Date(graph_name),
                                y=as.numeric(value),
                                group = type))+
  season_line +
  # season_text+
  spring_bloom_line+
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank())+
  #geom_line(aes(colour=variable))+
  geom_line(aes(colour=type), size =1.5, alpha= 0.8) +
  scale_colour_manual(values = palette_plots_colours)+
  xlab(NULL)+
  ylab("Average \npath \nlength")+ theme(legend.title = element_blank())
 
 ## average degree from network and random average degree are the same.Maybe median in better. 
 median_degree_plot <- ggplot(subset(melted_network_stats,
                                     variable == "median_degree"),
                              aes(x=as.Date(graph_name), 
                                  y=as.numeric(value), 
                                  group = type))+
  season_line +
  # season_text+
  spring_bloom_line+
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())+
  geom_line(aes(colour=type), size =1.5, alpha= 0.8) +
  scale_colour_manual(values = palette_plots_colours)+
  xlab(x_label)+
  ylab("Median \ndegree\n")+ theme(legend.title = element_blank(),
                                   axis.text.x = element_text(size = 13),
                                   axis.title.x = element_text(vjust = -2))
 
 plot_grid(nodes_edges_plot,
           edges_plot,
           diameter_plot,
           density_plot,
           modularity_plot,
           path_length_plot,
           median_degree_plot,ncol=1, align="v")
 
 
 prow <- plot_grid(nodes_edges_plot+ theme(legend.position="none"),
                   edges_plot+ theme(legend.position="none"),
                   diameter_plot+ theme(legend.position="none"),
                   density_plot+ theme(legend.position="none"),
                   modularity_plot+ theme(legend.position="none"),
                   path_length_plot+ theme(legend.position="none"),
                   median_degree_plot+ theme(legend.position="none"),
                   ncol=1,
                   align="v")
 
 prow
 
 grobs <- ggplotGrob(nodes_edges_plot)$grobs
 legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
 
 # add the legend to the row we made earlier. Give it one-third of the width
 # of one plot (via rel_widths).
 p <- plot_grid( prow, legend, rel_widths = c(2.50, 1))
 p
 
}


pdf(file="../figures/plots_network_stats_all_together%01d.pdf", height=15, width=11,  onefile=FALSE)
plot_stats_together(all_network_stats)
dev.off()

## divide into only single and only together. still compare to everyting


overall_and_two_types_stats <- c("overall","bac_to_gp23","euk_to_bac","euk_to_MPL") 
overall_and_within_types_stats <- c("overall", "bac_to_bac","euk_to_euk","MPL_to_MPL","gp23_to_gp23")

pdf(file="../figures/plots_network_stats_two_types%01d.pdf", height=14, width=9,  onefile=FALSE)
plot_stats_together(subset(all_network_stats, type %in% overall_and_two_types_stats))
dev.off()


pdf(file="../figures/plots_network_stats_within_type%01d.pdf", height=14, width=9,  onefile=FALSE)
plot_stats_together(subset(all_network_stats, type %in% overall_and_within_types_stats))
dev.off()




# ## testing correlations between network stats and env parameters, richenss....
# library(Hmisc)
# Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
# Jericho_data$Date <- as.Date(Jericho_data$Date)
# 
# Dates_in_common <- intersect(as.Date(overall_network_stats$graph_name), as.Date(Jericho_data$Date))
# 
# subset_Jericho <- subset(Jericho_data, Date %in% Dates_in_common)
# 
# overall_network_stats$nodes <- as.numeric(as.character(overall_network_stats$nodes))
# overall_network_stats$edges<- as.numeric(as.character(overall_network_stats$edges))
# overall_network_stats$diameter <- as.numeric(as.character(overall_network_stats$diameter))
# overall_network_stats$density <- as.numeric(as.character(overall_network_stats$density))
# overall_network_stats$modularity <- as.numeric(as.character(overall_network_stats$modularity))
# 
# cor.test(overall_network_stats[3:7] , subset_Jericho[2:4])
# 
# rcorr(as.matrix(overall_network_stats[3:7]) , as.matrix(subset_Jericho[2:4]))
# 
# richness <- read.csv("../results/amplicon_richness_by_date.csv")
# subset_richness <- richness[!(is.na(richness$richness.MPL)),]
# 
# Dates_in_common <- intersect(as.Date(subset_richness$Date),as.Date(overall_network_stats$graph_name))
# 
# subset_richness_MPL<- subset(subset_richness, as.Date(Date) %in% Dates_in_common)
# 
# rcorr(as.matrix(overall_network_stats[3:7]) ,as.matrix(subset_richness_MPL[c(-1,-2)]))

library(igraph)
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(scales)
library(grid)

args <- commandArgs(TRUE)
inputFile <- args[1]
if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
 path_colour <- "black"
 line_colour <- "black"
 point_colour <- "black"
 figures_dir <- "../figures/"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
  point_colour <- "white"
  figures_dir <- "../figures_pres/"
 }
}




## make graphs for all of these
counts_by_triplets <- list.files(path = "../results/LSA_tables", pattern = "counts_triplets_associated_.*with.*csv", all.files = FALSE,full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

## exclude the triplets that have less than 6 edges

counts_by_triplets <- counts_by_triplets[!(counts_by_triplets == "counts_triplets_associated_withpH.csv")]
counts_by_triplets <- counts_by_triplets[!(counts_by_triplets == "counts_triplets_associated_withAverage_chl_a.csv")]
counts_by_triplets <- counts_by_triplets[!(counts_by_triplets == "counts_triplets_associated_withmonth_number.csv")]


list_of_plots <- list()

mypalette<-brewer.pal(4,"Set2")

for (i in 1:length(counts_by_triplets)){

 count_triplet_df <- read.csv(paste0("../results/LSA_tables/",counts_by_triplets[i]))

 title.csv <- gsub("counts_triplets_associated_with", "", counts_by_triplets[i])
 title <- gsub(".csv", "", title.csv)
 title <- gsub("Average", "", title)
 title <- gsub("YSI", "", title)
 title <- gsub("percent", "", title)
 title <- gsub("ppt", "", title)
 title <- gsub("_", " ", title)
 
 melted_counts_triplets <- melt(count_triplet_df)
 ## remove the counts and positive and negative
 melted_counts_triplets <- droplevels(subset(melted_counts_triplets, !(variable %in% c("counts", "positive", "negative", "pos_neg_env_pos_otus", "pos_neg_env_neg_otus"))))


  melted_counts_triplets$X <- gsub("_", " to ", melted_counts_triplets$X)
 
  p <- ggplot(melted_counts_triplets, aes(x= X, y=value, group=variable))+
   geom_bar(aes(fill=variable),stat="identity")+
   scale_fill_manual(values = mypalette,
    name="",
    breaks=c("pos_env_pos_otus","neg_env_pos_otus","pos_env_neg_otus",  "neg_env_neg_otus"),
    labels=  c("+ to env,\n+ between OTUs", 
               "- to env,\n+ between OTUs",   "+ to env,\n- between OTUs", 
               "- to env,\n- between OTUs"))+
   ggtitle(title)+
   scale_x_discrete(labels=c("Eukaryotes\nto \nbacteria" ,
                             "Bacteria \nto \nT4-like\n myovirus",
                              "Eukaryotes \nto \nmarine\n picorna-like \nviruses"))+
   xlab(NULL)+
   ylab("Count")+
   theme_JAG_presentation()+
   theme(plot.margin = unit(c(6,0,6,0), "pt"),
         legend.text = element_text(size = 20),
         legend.key.size = unit(2, "cm"))
  
  list_of_plots[[i]] <-p 
 


}

# arrange the three plots in a single row
prow <- plot_grid( list_of_plots[[1]]+ theme(legend.position="none",axis.text.x=element_blank()),
                   list_of_plots[[2]]+
                    ylab(NULL)+
                    theme(legend.position="none",axis.text.x=element_blank()),
                   list_of_plots[[3]]+
                    ylab(NULL)+ 
                    theme(legend.position="none",axis.text.x=element_blank()),
                   list_of_plots[[4]]+ theme(legend.position="none",axis.text.x=element_blank()),
                   list_of_plots[[5]]+ 
                    ylab(NULL)+
                    theme(legend.position="none",axis.text.x=element_blank()),
                   list_of_plots[[6]]+ 
                    ylab(NULL)+
                    theme(legend.position="none",axis.text.x=element_blank()),
                   list_of_plots[[7]]+ theme(legend.position="none"),
                   list_of_plots[[8]]+ 
                    ylab(NULL)+
                    theme(legend.position="none"),
                   list_of_plots[[9]]+
                    ylab(NULL)+ 
                    theme(legend.position="none"),
                   align = 'vh',
                  # labels = c("A", "B", "C"),
                  # hjust = -1,
                   ncol = 3
)
prow

grobs <- ggplotGrob(list_of_plots[[1]])$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).


pdf(paste0(figures_dir,"overall_counts_of_triplets_by_env%03d.pdf"), width = 19, height = 11, onefile=FALSE)

p <- plot_grid( prow, legend, rel_widths = c(3, .75))
p
dev.off()
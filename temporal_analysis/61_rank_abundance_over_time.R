#'
#' @source "Ordering categories within ggplot2 Facets" by Tyler Rinker:
#' \url{https://trinkerrstuff.wordpress.com/2016/12/23/ordering-categories-within-ggplot2-facets/}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}


library(reshape2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",row.names="VC_number")
normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv", row.names="VC_number") 


RdRp_OTU_annotations <- read.csv("../results/RdRp_groups_with_OTUs.csv")
gp23_OTU_annotations <- read.csv("../results/gp23_groups_with_OTUs.csv")



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


theme_overall <- theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
## will start by using this MPL

### Get relevant VC numbers ====
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv",
                         row.names = 1)
Jericho_data$Date <- as.Date(Jericho_data$Date)
Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv",
                             nrows=61)

## first rename the column names in the OTU table to something better
add_date_to_OTU_table_rtrn_long <- function (otu_table, 
                                             Jericho_data) {
  
  proportional_by_row <- as.data.frame(prop.table(as.matrix(otu_table), margin=1))
  proportional_by_row$VC_number = rownames(proportional_by_row)
  long_otus <- melt(proportional_by_row,
                    id = "VC_number",
                    variable.name="OTUid")
  ## Add in dates:
  long_otus$Date <- Jericho_data$Date[match(long_otus$VC_number,
                                            Jericho_data$VC_number)]
  return(long_otus)
}


generate_frequency_by_otu_plot <- function (normalized_OTUs, name, group_file, colour_palette) {
  long_otus <- add_date_to_OTU_table_rtrn_long(normalized_OTUs, Jericho_data)
  
#  long_otus <- dplyr::filter(long_otus, value > 0)

  length_data_frame <- dim(long_otus)[1]

  tb <- table(long_otus$OTUid)
  
  ## rearranage the factors by the counts of the factors. Biggest first and then decreasing. 
  long_otus <- within(long_otus, OTUid <- factor(long_otus$OTUid,
                                                                levels = names(tb[order(tb, decreasing = TRUE)])))
  
  long_otus$Group <- group_file$group[match(long_otus$OTUid,group_file$values)]
  
  long_otus_top_50 <-  long_otus %>% 
    group_by(Date) %>% 
    top_n(n = 50, value)
  
  ## subset by first 50 OTUs...
  #long_otus_top_50 <- droplevels(long_otus[long_otus$OTUid %in% levels(long_otus$OTUid)[1:50],])
  
print(levels(long_otus_top_50$OTUid))



## shouldn't they be reordered each time?


## also want to add in the group info for each OTU. 

# pdf( paste("../figures/",name, "_frequency_by_otus%03d.pdf", sep=""),width = 17, height = 11,onefile = FALSE)


# ggplot(long_otus_top_50 , aes(reorder_within(OTUid, -value, Date),y = value*100, fill = Group))+
ggplot(long_otus_top_50 , aes(reorder_within(OTUid, -value, Date),y = value*100, colour = Group))+
  #geom_bar(stat = "identity")+
  geom_point()+
  scale_x_reordered()+
  #scale_y_log10()+
  facet_wrap(.~Date, scales = "free")+
  scale_y_continuous(trans="log10", breaks = c(0, 0.0001, 0.001, 0.01, 0.1,1,10))+
  scale_colour_manual(values = colour_palette)+
    #theme_JAG_presentation()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(), axis.title.x =element_blank())+
  xlab(NULL)

  #dev.off()
  
}


MPL_colours <-c("grey",brewer.pal(7,
                      "Dark2"))

generate_frequency_by_otu_plot(normalized_MPL_OTUs, "MPL", RdRp_OTU_annotations,MPL_colours )

gp23_colours <- c("grey",brewer.pal(9,
                      "Set3"))

generate_frequency_by_otu_plot(normalized_gp23_OTUs, "gp23", gp23_OTU_annotations, gp23_colours)

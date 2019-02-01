

library(ggplot2)
library(plyr)
library(reshape2)
library(dplyr)
library(scales)
library(cowplot)

args <- commandArgs(TRUE)
inputFile <- args[1]

## test to see if input file is given, so I can decide whether to use this argument or the orginal one.
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
  if (inputFile == "../../JAG_black_presentation.R") {
    path_colour <- "white"
    line_colour <- "white"
    point_colour <- "white"
    figures_dir <- "../figures_pres/"
  }
}

# adding in the seasons and the spring bloom for the ggplots
season_line <-
  geom_vline(xintercept = as.numeric(c(
    as.Date("2010-03-22"),
    as.Date("2010-06-22"),
    as.Date("2010-09-22"),
    as.Date("2010-12-22"),
    as.Date("2011-03-22"),
    as.Date("2011-06-22")
  )),
  colour = "grey",
  size = 1.5)

date_scaling <-   scale_x_date(
  breaks = date_breaks("month"),
  labels = date_format("%b"),
  limits = c(as.Date("2010-06-15"),
             as.Date("2011-07-25"))
)

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv",
                                  row.names = "VC_number")
normalized_gp23_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv",
                                   row.names = "VC_number")
normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv",
                                  row.names = "VC_number")
normalized_16s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv",
                                  row.names = "VC_number")


## first plot the different stuff
high_res_vcs <- c(1198,
                  1199,
                  1200,
                  1201,
                  1202)

## remove the high res time points
normalized_MPL_OTUs <- subset(normalized_MPL_OTUs,
                              !(rownames(normalized_MPL_OTUs) %in% high_res_vcs))
normalized_gp23_OTUs <- subset(normalized_gp23_OTUs,
                               !(rownames(normalized_gp23_OTUs) %in% high_res_vcs))
normalized_18s_OTUs<- subset(normalized_18s_OTUs,
                             !(rownames(normalized_18s_OTUs) %in% high_res_vcs))
normalized_16s_OTUs <- subset(normalized_16s_OTUs,
                              !(rownames(normalized_16s_OTUs) %in% high_res_vcs))



proportional_18s <- as.data.frame(prop.table(as.matrix(normalized_18s_OTUs),
                                             margin = 1)) * 100
proportional_16s <- as.data.frame(prop.table(as.matrix(normalized_16s_OTUs),
                                             margin = 1)) * 100

taxonomy_18s <-read.csv("../results/cleaned_up_18s_taxonomy_Jericho.csv",
                        row.names = 1)
taxonomy_16s <-read.csv("../results/cleaned_up_16s_taxonomy_Jericho.csv",
                        row.names = 1)






##########################
## Try to separate by het and phototrophic 18s
## like fig 2 in de Vargas et al 2015

raphidophyte_fam <- c("Raphidophyceae")

diatom_fam <- c("Diatomea") ## not sure about this, are all dinos phytoplankton?

raphidophyte_tax <- taxonomy_18s %>%
  filter(Family %in% raphidophyte_fam)

diatom_tax <- taxonomy_18s %>%
  filter(Family %in% diatom_fam)

colnames(proportional_18s) <- gsub(".size.*.",
                                   "",
                                   colnames(proportional_18s))

raphidophyte_otu_table <- droplevels(subset(
  proportional_18s,
  select = colnames(proportional_18s) %in% raphidophyte_tax$otu_number
))

## remove high res
raphidophyte_otu_table<- subset(raphidophyte_otu_table,
                              !(rownames(raphidophyte_otu_table) %in% high_res_vcs))

diatom_otu_table <- droplevels(subset(
  proportional_18s,
  select = colnames(proportional_18s) %in% diatom_tax$otu_number
))

diatom_otu_table <- subset(diatom_otu_table ,
                                !(rownames(diatom_otu_table ) %in% high_res_vcs))

Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv",
                         row.names = 1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

Jericho_data <- subset(Jericho_data,
                              !(Jericho_data$VC_number %in% high_res_vcs))

## RdRp

proportional_MPL_groupA <- read.csv("../results/proportional_MPL_groupA.csv",
                                    row.names = 1)*100
proportional_MPL_groupB <- read.csv("../results/proportional_MPL_groupB.csv",
                                    row.names = 1)*100
#proportional_MPL_groupC <- read.csv("../results/proportional_MPL_groupC.csv", row.names = 1)
proportional_MPL_groupD <- read.csv("../results/proportional_MPL_groupD.csv",
                                    row.names = 1)*100
proportional_MPL_groupE <- read.csv("../results/proportional_MPL_groupE.csv",
                                    row.names = 1)*100
#proportional_MPL_groupF <- read.csv("../results/proportional_MPL_groupF.csv", row.names = 1)

proportional_MPL_groupG <- read.csv("../results/proportional_MPL_groupG.csv",
                                    row.names = 1)*100
proportional_MPL_groupH <- read.csv("../results/proportional_MPL_groupH.csv",
                                    row.names = 1)*100





proportional_MPL_groupA <- subset(proportional_MPL_groupA,
                                  !(rownames(proportional_MPL_groupA) %in% high_res_vcs))
proportional_MPL_groupB <- subset(proportional_MPL_groupB,
                                  !(rownames(proportional_MPL_groupB) %in% high_res_vcs))
proportional_MPL_groupD <- subset(proportional_MPL_groupD,
                                  !(rownames(proportional_MPL_groupD) %in% high_res_vcs))
proportional_MPL_groupE <- subset(proportional_MPL_groupE,
                                  !(rownames(proportional_MPL_groupE) %in% high_res_vcs))
proportional_MPL_groupG <- subset(proportional_MPL_groupG,
                                  !(rownames(proportional_MPL_groupG) %in% high_res_vcs))
proportional_MPL_groupH <- subset(proportional_MPL_groupH,
                                  !(rownames(proportional_MPL_groupH) %in% high_res_vcs))

## want Jericho_dates that are not in MPL table and high res

Jericho_no_high_res_vcs <- Jericho_data$VC_number
## want to annotate missing data days.
## want to annotate those days using star
missing_MPL_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_MPL_OTUs))]
missing_MPL_dates <- Jericho_data$Date[match(missing_MPL_samples,
                                             Jericho_data$VC_number)]

missing_gp23_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_gp23_OTUs))]
missing_gp23_dates <- Jericho_data$Date[match(missing_gp23_samples,
                                              Jericho_data$VC_number)]

missing_18s_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_18s_OTUs))]
missing_18s_dates <- Jericho_data$Date[match(missing_18s_samples,
                                             Jericho_data$VC_number)]

missing_16s_samples <- Jericho_no_high_res_vcs[!(Jericho_no_high_res_vcs %in% rownames(normalized_16s_OTUs))]
missing_16s_dates <- Jericho_data$Date[match(missing_16s_samples,
                                             Jericho_data$VC_number)]

## first rename the column names in the OTU table to something better
add_date_to_OTU_table_rtrn_long <- function(otu_table,
                                            Jericho_data) {
  otu_table$VC_number = rownames(otu_table)
  long_otus <- melt(otu_table,
                    id = "VC_number",
                    variable.name = "OTUid")
  ## Add in dates:
  long_otus$Date <- Jericho_data$Date[match(long_otus$VC_number,
                                            Jericho_data$VC_number)]
  return(long_otus)
}

## Read in OTU table and get long format with Date
long_rdrp_groupA_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupA,
                                                         Jericho_data)
long_rdrp_groupB_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupB,
                                                         Jericho_data)
long_rdrp_groupE_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupE,
                                                         Jericho_data)
long_rdrp_groupG_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupG,
                                                         Jericho_data)
long_rdrp_groupH_otus <- add_date_to_OTU_table_rtrn_long(proportional_MPL_groupH,
                                                         Jericho_data)
long_proportional_MPL <- add_date_to_OTU_table_rtrn_long(normalized_MPL_OTUs,
                                                         Jericho_data)

long_raphidophyte_otus <- add_date_to_OTU_table_rtrn_long(raphidophyte_otu_table,
                                                          Jericho_data)
long_diatom_otus <- add_date_to_OTU_table_rtrn_long(diatom_otu_table,
                                                    Jericho_data)
long_proportional_18s <- add_date_to_OTU_table_rtrn_long(normalized_18s_OTUs,
                                                         Jericho_data)

## Plot single OTUs over time
plot_richness_over_time <- function(long_otus,
                                    title,
                                    type_OTU,
                                    missing_dates,
                                    all_long_otus) {
  ## need to group long otus together by data
  # long_otus_sum <- ddply(all_long_otus,
  #                        "Date",
  #                        summarise,
  #                        rel_abun = sum(value))
  
  
  all_long_otus_summarised <- all_long_otus %>%
    group_by(Date) %>%
    filter(value > 0) %>%
    dplyr::summarise(richness = n(),
                     sum_rel_abun = sum(value, na.rm = TRUE))
  
  p_line <- ggplot(droplevels( all_long_otus_summarised ),
                   aes(x = Date,
                       y = (richness))) +
    season_line +
    #geom_line(colour = line_colour) +
    geom_bar(stat = "identity")+
    date_scaling +
    theme_JAG_presentation(base_size = 7) +
    ylab("Richness\nTotal\n") +
    xlab("\nDate\n") +
     theme(
       axis.title.x = element_blank(),
       axis.text.x = element_blank())
  number_OTUs <- length(unique(long_otus$OTUid))
  
  long_otus_summarised <- long_otus %>%
    group_by(Date) %>%
    filter(value > 0) %>%
    dplyr::summarise(richness = n(),
                     sum_rel_abun = sum(value, na.rm = TRUE))
  
  ## for each missing dates need to have for each group
  ## inelegant solution, btu it wokrs for now.
  missing_dates_df <- c()
  
  for (date in missing_dates) {
    # cat(date)
    df <- data.frame(
      Date = as.Date(date, origin = "1970-01-01"),
      richness = NA,
      sum_rel_abun = NA
    )
    missing_dates_df <- rbind(missing_dates_df, df)
  }
  
  long_otus_summarised <- bind_rows(long_otus_summarised, missing_dates_df)
  
  ## need to add in dates where it was zero
  
  dates_to_add <- Jericho_data$Date[!Jericho_data$Date %in% long_otus_summarised$Date]
  
  zero_date_rows <- do.call(rbind, (lapply(dates_to_add, function(x) {
    date_df <- data.frame(
      Date = as.Date(x, origin = "1970-01-01"),
      richness = 0,
      sum_rel_abun = 0
    )
  })))
  
  long_otus_summarised <- bind_rows(long_otus_summarised, zero_date_rows)
  
  p <- ggplot(long_otus_summarised ,
              aes(x = Date,
                  y = richness))+
    season_line+
  ## I feel this shows the missing data too much
      geom_line(colour = "grey")+
    geom_point(aes(colour = sum_rel_abun), size = 3)+
  # geom_bar(aes(colour = sum_rel_abun), stat = "identity")+
     annotate(
      "text",
      x = missing_dates,
      y = -0.1,
      label = "x",
      size = 3,
      colour = "grey"
    )+
    # facet_grid(group~., scales = "free_y")+
    date_scaling+
    scale_colour_viridis_c(name = "Relative\nAbundance\n(%)", direction = -1)+
   # scale_color_distiller(palette = "YlGnBu", name = "Relative\nAbundance(%)", direction = -1)+
    theme_JAG_presentation(base_size = 7)+
    theme(strip.background = element_rect(fill = "white"))+
    theme(legend.position = c(0.95, 0.5), legend.background = element_rect(colour = NA, fill = NA))+
    ylab("Number of OTUs observed")
  
  new_vect <- list(p_line, p)
  return(new_vect)
}

filtered_long_raphidophyte <- long_raphidophyte_otus %>%
  filter(Date < as.Date("2011-06-01"))

filtered_long_diatom <- long_diatom_otus %>%
  filter(Date < as.Date("2011-06-01"))

group_A_rdrp <-
  plot_richness_over_time(long_rdrp_groupA_otus,
                          "MPL group A",
                          "virus",
                          missing_MPL_dates,
                          long_proportional_MPL)

group_B_rdrp <-
  plot_richness_over_time(long_rdrp_groupB_otus,
                          "MPL group B",
                          "virus",
                          missing_MPL_dates,
                          long_proportional_MPL)
group_E_rdrp <-
  plot_richness_over_time(long_rdrp_groupE_otus,
                          "MPL group E",
                          "virus",
                          missing_MPL_dates,
                          long_proportional_MPL)
group_G_rdrp <-
  plot_richness_over_time(long_rdrp_groupG_otus,
                          "MPL group G",
                          "virus",
                          missing_MPL_dates,
                          long_proportional_MPL)
group_H_rdrp <-
  plot_richness_over_time(long_rdrp_groupH_otus,
                          "MPL group H",
                          "virus",
                          missing_MPL_dates,
                          long_proportional_MPL)

raphidophyte_over_time_matching_rdrp<- plot_richness_over_time(filtered_long_raphidophyte,
                                                               "Raphidophyte family",
                                                               "host",
                                                               missing_18s_dates,
                                                               long_proportional_18s)
######################

diatom_over_time_matching_rdrp <-
  plot_richness_over_time(filtered_long_diatom,
                          "Diatom family",
                          "host",
                          missing_18s_dates,
                          long_proportional_18s)


pdf(
  paste0(figures_dir, "group_B_with_Raphidophytes.pdf"),
  width = 6,
  ##height = 14,
  onefile = FALSE
)
print(
  plot_grid(
    raphidophyte_over_time_matching_rdrp[[1]],
    raphidophyte_over_time_matching_rdrp[[2]],
    group_B_rdrp[[1]],
    group_B_rdrp[[2]],
    align = "v",
    rel_heights = c(1, 4, 1, 4),
    ncol = 1,
    labels = c("A", "", "B", ""),
    label_size = 20
  )
)
dev.off()


pdf(
  paste0(figures_dir, "group_E_with_Diatoms.pdf"),
  width = 6,
 ## height = 14,
  onefile = FALSE
)
print(
  plot_grid(
    diatom_over_time_matching_rdrp[[1]],
    diatom_over_time_matching_rdrp[[2]],
    group_E_rdrp[[1]],
    group_E_rdrp[[2]],
    align = "v",
    rel_heights = c(1, 4, 1, 4),
    ncol = 1,
    labels = c("A", "", "B", ""),
    label_size = 20
  )
)
dev.off()


pdf(
  paste0(figures_dir, "group_G_with_Diatoms.pdf"),
  width = 6,
 ## height = 14,
  onefile = FALSE
)
print(
  plot_grid(
    diatom_over_time_matching_rdrp[[1]],
    diatom_over_time_matching_rdrp[[2]],
    group_G_rdrp[[1]],
    group_G_rdrp[[2]],
    align = "v",
    rel_heights = c(1, 4, 1, 4),
    ncol = 1,
    labels = c("A", "", "B", ""),
    label_size = 20
  )
)
dev.off()


pdf(
  paste0(figures_dir, "group_H_with_Diatoms.pdf"),
  width = 6,
 ## height = 14,
  onefile = FALSE
)
print(
  plot_grid(
    diatom_over_time_matching_rdrp[[1]],
    diatom_over_time_matching_rdrp[[2]],
    group_H_rdrp[[1]],
    group_H_rdrp[[2]],
    align = "v",
    rel_heights = c(1, 4, 1, 4),
    ncol = 1,
    labels = c("A", "", "B", ""),
    label_size = 20
  )
)
dev.off()

pdf(
  paste0(figures_dir, "group_A_with_Raphidophytes.pdf"),
  width = 6,
  ##height = 10
)
print(
  plot_grid(
    raphidophyte_over_time_matching_rdrp[[1]],
    raphidophyte_over_time_matching_rdrp[[2]] + xlab("") + ylim(c(-0.5, 20)),
    group_A_rdrp[[1]],
    group_A_rdrp[[2]] ,
    ncol = 1,
    rel_heights = c(1, 6, 1, 6),
    align = "v",
    labels = c("A", "", "B", ""),
    label_size = 20
  )
)
dev.off()

#################################################
### what about with gp23 group I and cyanobacteria

## gp23

proportional_gp23_groupA <- read.csv("../results/proportional_gp23_groupA.csv",
                                     row.names = 1)*100
proportional_gp23_groupB <- read.csv("../results/proportional_gp23_groupB.csv",
                                     row.names = 1)*100
proportional_gp23_groupC <- read.csv("../results/proportional_gp23_groupC.csv",
                                     row.names = 1)*100
proportional_gp23_groupD <- read.csv("../results/proportional_gp23_groupD.csv",
                                     row.names = 1)*100
proportional_gp23_groupE <- read.csv("../results/proportional_gp23_groupE.csv",
                                     row.names = 1)*100
proportional_gp23_groupF <- read.csv("../results/proportional_gp23_groupF.csv",
                                     row.names = 1)*100
proportional_gp23_groupG <- read.csv("../results/proportional_gp23_groupG.csv",
                                     row.names = 1)*100
proportional_gp23_groupH <- read.csv("../results/proportional_gp23_groupH.csv",
                                     row.names = 1)*100
proportional_gp23_groupI <- read.csv("../results/proportional_gp23_groupI.csv",
                                     row.names = 1)*100
#proportional_gp23_groupJ <- read.csv("../results/proportional_gp23_groupJ.csv", row.names = 1)


proportional_gp23_groupA <- subset(proportional_gp23_groupA,
                                  !(rownames(proportional_gp23_groupA) %in% high_res_vcs))
proportional_gp23_groupB <- subset(proportional_gp23_groupB,
                                   !(rownames(proportional_gp23_groupB) %in% high_res_vcs))
proportional_gp23_groupC <- subset(proportional_gp23_groupC,
                                   !(rownames(proportional_gp23_groupC) %in% high_res_vcs))
proportional_gp23_groupD <- subset(proportional_gp23_groupD,
                                   !(rownames(proportional_gp23_groupD) %in% high_res_vcs))
proportional_gp23_groupE <- subset(proportional_gp23_groupE,
                                   !(rownames(proportional_gp23_groupE) %in% high_res_vcs))
proportional_gp23_groupF <- subset(proportional_gp23_groupF,
                                   !(rownames(proportional_gp23_groupF) %in% high_res_vcs))
proportional_gp23_groupG <- subset(proportional_gp23_groupG,
                                   !(rownames(proportional_gp23_groupG) %in% high_res_vcs))
proportional_gp23_groupH <- subset(proportional_gp23_groupH,
                                   !(rownames(proportional_gp23_groupH) %in% high_res_vcs))
proportional_gp23_groupI <- subset(proportional_gp23_groupI,
                                   !(rownames(proportional_gp23_groupI) %in% high_res_vcs))

cyanobacteria_phylum <- c("Cyanobacteria")
cyanobacteria_class <- c("Cyanobacteria")

cyano_tax <- taxonomy_16s %>%
  filter(Class %in% cyanobacteria_class)

colnames(proportional_16s) <- gsub(".size.*.",
                                   "",
                                   colnames(proportional_16s))

cyanobacteria_otu_table <- droplevels(subset(proportional_16s,
                                             select = colnames(proportional_16s) %in% cyano_tax$otu_number
))


## try phylum

cyano_tax_phyla <- taxonomy_16s %>%
  filter(Phylum %in% cyanobacteria_phylum)

cyanobacteria_otu_table_phyla <- droplevels(subset(proportional_16s,
                                                   select = colnames(proportional_16s) %in% cyano_tax_phyla$otu_number))
write.csv(cyanobacteria_otu_table_phyla,
          "../results/normalized_16s_cyano_otus.csv")

## Read in OTU table and get long format with Date
long_gp23_groupI_otus <-add_date_to_OTU_table_rtrn_long(proportional_gp23_groupI,
                                                        Jericho_data)

## Read in OTU table and get long format with Date
long_gp23_all_otus <-add_date_to_OTU_table_rtrn_long(normalized_gp23_OTUs,
                                                        Jericho_data)

long_cyano_otus <-add_date_to_OTU_table_rtrn_long(cyanobacteria_otu_table,
                                                  Jericho_data)

long_cyano_otus_phyla <-add_date_to_OTU_table_rtrn_long(cyanobacteria_otu_table_phyla,
                                                        Jericho_data)

long_all_16s <-add_date_to_OTU_table_rtrn_long(proportional_16s,
                                                        Jericho_data)

group_I_gp23 <-plot_richness_over_time(long_gp23_groupI_otus,
                                       "T4-like myoviral group I",
                                       "virus",
                                       missing_gp23_dates,
                                       long_gp23_all_otus)


cyano_over_time_matching_gp23 <-
  plot_richness_over_time(long_cyano_otus,
                          "Class cyanobacteria",
                          "host",
                          missing_16s_dates,
                          long_all_16s)


cyano_over_time_matching_gp23_phyla <-
  plot_richness_over_time(long_cyano_otus_phyla,
                          "phylum Cyanobacteria",
                          "host",
                          missing_16s_dates,
                          long_all_16s)

#dev.off()

pdf(
  paste0(figures_dir, "group_I_with_Cyanos.pdf"),
  width = 6,
  ##height = 14,
  onefile = FALSE
)
print(
  plot_grid(
    cyano_over_time_matching_gp23[[1]],
    cyano_over_time_matching_gp23[[2]],
    group_I_gp23[[1]],
    group_I_gp23[[2]],
    align = "v",
    rel_heights = c(1, 6, 1, 6),
    ncol = 1,
    labels = c("A", "", "B", ""),
    label_size = 20
  )
)
dev.off()


pdf(
  paste0(figures_dir, "group_I_with_Cyanos_phyla.pdf"),
  width = 6,
  ##height = 10,
  onefile = FALSE)
print(
  plot_grid(
    cyano_over_time_matching_gp23_phyla[[1]],
    cyano_over_time_matching_gp23_phyla[[2]],
    group_I_gp23[[1]],
    group_I_gp23[[2]],
    align = "v",
    rel_heights = c(1, 4, 1, 4),
    ncol = 1,
    labels = c("A", "", "B", ""),
    label_size = 20
  )
)
dev.off()

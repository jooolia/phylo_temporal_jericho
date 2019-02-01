## Purpose: Create a panel chart from Jericho Pier data that I collected over 1 year.
## Author: Julia Gustavsen
## Date Created: 2 September 2013
## Date modified: 9 May 2016
library(ggplot2)
library(scales)
library(reshape2)
library(gridExtra)
library(dplyr)

args <- commandArgs(TRUE)
inputFile <- args[1]

## test to see if input file is given, so I can decide whether to use this argument or the orginal one
if (!file_test("-f", inputFile)) {
  print("input theme not defined, using orginal one for manuscript.")
  source("../../JAG_manuscript_figure.R")
  path_colour <- "black"
  line_colour <- "black"
  figures_dir <- "../figures/"
  
} else {
  print("Cool you passed a nice theme file to this script")
  source(inputFile)
  if (inputFile == "../../JAG_black_presentation.R"){
    path_colour <- "white"
    line_colour <- "white"
    figures_dir <- "../figures_pres/"
  }
}


Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv",
                         row.names = 1)
Jericho_data$Date <- as.Date(Jericho_data$Date)

irradiation_data <- read.csv("../data/ibis_shortwave_irradiance_2010_2011_date_fixed.csv",
                             row.names = 1)

irradiation_data$Date <- as.Date(irradiation_data$Date)

Jericho_data <- merge(Jericho_data,
                      irradiation_data[, c("Mean","Date")],
                      by = c("Date", "Date"))

colnames(Jericho_data)[dim(Jericho_data)[2]]<- "Mean_irradiation"

#Reshape the data into a long format so that I can use the facet wrapping on it.
Jericho_enviro_and_bio_long <- melt(Jericho_data,
                                    id = ("Date"),
                                    na.rm = TRUE)

## option so that viral and bacterial plots have better scales
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e+", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  ## plus needs to be escaped with \\
  l <- gsub("\\+", "", l)
  # return this as an expression
  parse(text = l)
}

### Plots of all variables for full data set. ###
# adding in the seasons and the spring bloom for the ggplots
season_line <- geom_vline(xintercept = as.numeric(c(as.Date("2010-03-22"), 
                                                    as.Date("2010-06-22"), 
                                                    as.Date("2010-09-22"),
                                                    as.Date("2010-12-22"),
                                                    as.Date("2011-03-22"),
                                                    as.Date("2011-06-22")
)
),
colour="grey",
size=1.5)
spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=1)

date_scaling <-   scale_x_date(breaks = date_breaks("month"), 
                               labels = date_format("%b'%y"),
                               limits = c(as.Date("2010-06-15"),
                                          as.Date("2011-07-25")
                               )
)


environmental_plot <- function (long_data_set, 
                                parameter,
                                ylimits,
                                ylabel) {
  
  env_plot <- ggplot(droplevels(subset(long_data_set,
                                       variable == parameter)),
                     aes(Date,as.numeric(value))) +
    season_line + 
  #  spring_bloom_line +
    geom_line(colour=line_colour, size=1.5) +
    geom_point(colour = "grey50") +
    date_scaling +
    ylab(ylabel) +
    scale_y_continuous(limits = ylimits,
                       expand = c(0,0)) +
    ggtitle(NULL) +
    theme_JAG_presentation()+
    theme(panel.grid.major = element_blank())
}

temperature_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long,
                                               "Temperature_YSI",
                                               c(0,25), 
                                               "Temperature \n (C)\n\n") + 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x =element_blank()) +
  xlab(NULL)

salinity_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long,
                                            "Salinity_ppt_YSI",
                                            c(0,30),
                                            "Salinity \n(psu)\n\n") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x =element_blank()) +
  xlab(NULL)

viral_abundance_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long,
                                                   "Average_viral_abundance",
                                                   c(1,60000000),
                                                   "Viral abundance \n(viruses/mL)\n") +
  scale_y_continuous(labels = fancy_scientific,
                     limits = c(1000000,50000000),
                     expand = c(0,0)) + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x =element_blank()) +
  xlab(NULL)

bacterial_abundance_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long,
                                                       "Average_bacterial_abundance",
                                                       c(1,60000000), 
                                                       "Bacterial abundance\n (cells/mL)\n") + 
  scale_y_continuous(labels = fancy_scientific,
                     limits = c(0,5000000),
                     expand = c(0,0)) + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x =element_blank()) +
  xlab(NULL)

chlorophyll_a_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long,
                                                 "Average_chl_a",
                                                 c(0,50),
                                                 "Chlorophyll a \n(ug/L)\n\n") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x =element_blank()) + 
  xlab(NULL)

chlorophyll_a_full_dataset_under_10 <- environmental_plot(Jericho_enviro_and_bio_long,
                                                          "Average_chl_a",
                                                          c(0,6), 
                                                          "Chlorophyll a \n(ug/L)\n\n") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank()) + 
  xlab(NULL)

phosphate_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long,
                                             "Average_PO4",
                                             c(0,3),
                                             "Phosphate\n (um/L)\n\n") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.title.x =element_blank()) +
  xlab(NULL)

silicate_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, 
                                            "Average_SiO2", 
                                            c(0,80),
                                            "Silicate \n(um/L)\n\n") + 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) + 
  xlab(NULL)

nitrate_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, 
                                           "Average_NO3_NO2",
                                           c(0,80), 
                                           "Nitrate+Nitrite \n(um/L)\n\n") +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,30))

irradiation_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long,
                                               "Mean_irradiation", 
                                               c(0,350),
                                               "Mean_irradiation units?") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,350))

# Get the widths
gA <- ggplot_gtable(ggplot_build(temperature_full_dataset))
gB <- ggplot_gtable(ggplot_build(salinity_full_dataset))
gC <- ggplot_gtable(ggplot_build(viral_abundance_full_dataset))
gD <- ggplot_gtable(ggplot_build(bacterial_abundance_full_dataset))
gE <- ggplot_gtable(ggplot_build(chlorophyll_a_full_dataset))
gF <- ggplot_gtable(ggplot_build(silicate_full_dataset))
gG <- ggplot_gtable(ggplot_build(phosphate_full_dataset))
gH <- ggplot_gtable(ggplot_build(nitrate_full_dataset))
gI <- ggplot_gtable(ggplot_build(irradiation_full_dataset))

maxWidth = grid::unit.pmax(
  gA$widths[2:3], 
  gB$widths[2:3],
  gC$widths[2:3], 
  gD$widths[2:3],
  gE$widths[2:3], 
  gF$widths[2:3],
  gG$widths[2:3], 
  gH$widths[2:3], 
  gI$widths[2:3]
)

# Set the widths
gA$widths[2:3] <- maxWidth
gB$widths[2:3] <- maxWidth
gC$widths[2:3] <- maxWidth
gD$widths[2:3] <- maxWidth
gE$widths[2:3] <- maxWidth
gF$widths[2:3] <- maxWidth
gG$widths[2:3] <- maxWidth
gH$widths[2:3] <- maxWidth
gI$widths[2:3] <- maxWidth


### Print graph to file ###

pdf(paste0(figures_dir,
           "Jericho_full_panel_chart.pdf"),
    width = 15,
    height = 15, 
    onefile = FALSE
)

grid.arrange(gA, 
             gB,
             gC,
             gD,
             gE,
             gF,
             gG,
             gH,
             #gI,
             ncol = 1)
dev.off()


### Plots With standard error of the mean ####

environmental_plot_with_se <- function(data_frame, 
                                       data_frame_mean,
                                       data_frame_se,
                                       ylimits,
                                       ylabel) {
  
  param_with_se <- data.frame(Sample_name = data_frame[,1],
                              Means = data_frame_mean,
                              param_se = data_frame_se)
  
  param_plot_line <- ggplot(param_with_se,
                            aes(Sample_name,
                                Means)) +
    season_line+
  #  spring_bloom_line+
    geom_line(colour=line_colour, size=1.5) +
    geom_point(colour = "grey50") +
    date_scaling +
    ylab(ylabel) +
    scale_y_continuous(limits = ylimits,
                       expand = c(0,0)) +
    ggtitle(NULL) +
    theme_JAG_presentation()+
    # theme(axis.title.y=element_text(hjust=0.2))+
    theme(panel.grid.major = element_blank())
  
  param_plot_line
  param_for_graph_with_se <-
    aes(ymax = Means + param_se,
        ymin = Means - param_se)
  param_plot_line_error <- param_plot_line +
    geom_errorbar(param_for_graph_with_se,
                  position = "dodge",
                  width = 0.25) 
  return(param_plot_line_error)
}


viral_abundance_full_dataset_with_se <-
  environmental_plot_with_se(
    Jericho_data,
    Jericho_data[,"Average_viral_abundance"],
    Jericho_data[,"Standard_error_viral_abundance"],
    c(1,60000000),
    "Viral \nabundance \n(viruses/mL)\n"
  ) + scale_y_continuous(labels = fancy_scientific,
                         limits = c(1000000,50000000),
                         expand = c(0,0)) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
  ) + xlab(NULL)

bacterial_abundance_full_dataset_with_se <- environmental_plot_with_se(Jericho_data,
                                                                       Jericho_data[,"Average_bacterial_abundance"],
                                                                       Jericho_data[,"Standard_error_bacterial_abundance"],
                                                                       c(1,60000000),
                                                                       "Bacterial \nabundance\n (cells/mL)\n"
)+ 
  scale_y_continuous(labels = fancy_scientific,
                     limits = c(0,5000000),
                     expand = c(0,0))+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank())+
  xlab(NULL)


chlorophyll_a_full_dataset_with_se <- environmental_plot_with_se(Jericho_data,
                                                                 Jericho_data[,"Average_chl_a"],
                                                                 Jericho_data[,"Standard_error_chl_a"],
                                                                 c(0,50), 
                                                                 "Chlorophyll a \n(ug/L)\n") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x =element_blank()) +
  xlab(NULL)

chlorophyll_a_full_dataset_under_10_with_se <- environmental_plot_with_se(Jericho_data,
                                                                          Jericho_data[,"Average_chl_a"],
                                                                          Jericho_data[,"Standard_error_chl_a"],
                                                                          c(0,6),
                                                                          "Chlorophyll a \n(ug/L)\n") +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank()) +
  xlab(NULL)

phosphate_full_dataset_with_se <- environmental_plot_with_se(Jericho_data,
                                                             Jericho_data[,"Average_PO4"],
                                                             Jericho_data[,"Standard_error_PO4"],
                                                             c(0,3), 
                                                             "Phosphate\n (um/L)\n") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.title.x =element_blank()) + 
  xlab(NULL)

silicate_full_dataset_with_se <- environmental_plot_with_se(Jericho_data,
                                                            Jericho_data[,"Average_SiO2"],
                                                            Jericho_data[,"Standard_error_SiO2"],
                                                            c(0,80),
                                                            "Silicate \n(um/L)\n") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x =element_blank()) +
  xlab(NULL)

nitrate_full_dataset_with_se <- environmental_plot_with_se(Jericho_data,
                                                           Jericho_data[,"Average_NO3_NO2"],
                                                           Jericho_data[,"Standard_error_NO3_NO2"],
                                                           c(0,80), 
                                                           "Nitrate+\nNitrite \n(um/L)") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,30))+
  xlab("\nDate")+
  theme(axis.title.x =element_text(vjust=0.5, size=25),
        axis.text.x = element_text(size=15))


gA <- ggplot_gtable(ggplot_build(temperature_full_dataset))
gB <- ggplot_gtable(ggplot_build(salinity_full_dataset))
gC <- ggplot_gtable(ggplot_build(viral_abundance_full_dataset_with_se))
gD <- ggplot_gtable(ggplot_build(bacterial_abundance_full_dataset_with_se))
gE <- ggplot_gtable(ggplot_build(chlorophyll_a_full_dataset_with_se))
gF <- ggplot_gtable(ggplot_build(silicate_full_dataset_with_se))
gG <- ggplot_gtable(ggplot_build(phosphate_full_dataset_with_se))
gH <- ggplot_gtable(ggplot_build(nitrate_full_dataset_with_se))
gI <- ggplot_gtable(ggplot_build(irradiation_full_dataset))

maxWidth = grid::unit.pmax(
  gA$widths[2:3], 
  gB$widths[2:3],
  gC$widths[2:3],
  gD$widths[2:3],
  gE$widths[2:3],
  gF$widths[2:3],
  gG$widths[2:3],
  gH$widths[2:3],
  gI$widths[2:3]
)

# Set the widths
gA$widths[2:3] <- maxWidth
gB$widths[2:3] <- maxWidth
gC$widths[2:3] <- maxWidth
gD$widths[2:3] <- maxWidth
gE$widths[2:3] <- maxWidth
gF$widths[2:3] <- maxWidth
gG$widths[2:3] <- maxWidth
gH$widths[2:3] <- maxWidth
gI$widths[2:3] <- maxWidth

### Print graph to file ###

pdf(paste0(figures_dir,
           "Jericho_full_panel_chart_with_se.pdf"),
    width = 15,
    height = 15,
    onefile = FALSE
)
grid.arrange(gA, 
             gB,
             gC,
             gD,
             gE,
             gF,
             gG,
             gH,
             #gI,
             ncol = 1)
dev.off()


pdf(paste0(figures_dir,"Jericho_env_with_se%03d.pdf"),
    width = 20,
    height = 15,
    onefile = FALSE
)
print(temperature_full_dataset)
print(salinity_full_dataset)
print(viral_abundance_full_dataset_with_se)
print(bacterial_abundance_full_dataset_with_se)
print(chlorophyll_a_full_dataset_with_se)
print(silicate_full_dataset_with_se)
print(phosphate_full_dataset_with_se)
print(nitrate_full_dataset_with_se)
print(irradiation_full_dataset)
dev.off()

### Want to compare the environmental data when I sampled to the whole datasets. ====

irradiation_data$group <- "all"
p <- ggplot(irradiation_data, 
            aes(x = Date, 
                y = Mean,
                group = group)) +
  season_line+
 # spring_bloom_line+
  geom_line(colour=line_colour,
            size=1.5)
irradiation_with_sampling <- p +
  geom_point(data = Jericho_data,
             aes(x = Date,
                 y = Mean_irradiation,
                 group = "test"),
             colour ="red",
             size = 3) +
  date_scaling+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())

pdf( paste0(figures_dir,"Jericho_irradiance_with_samples.pdf"), 
     width = 15,
     height = 11,
     onefile = FALSE)

irradiation_with_sampling
dev.off()

## Plot nutrients on one graph together ####
target <- c("Average_PO4",
            "Average_SiO2",
            "Average_NO3_NO2")

target_se <- c("Standard_error_PO4", 
               "Standard_error_SiO2", 
               "Standard_error_NO3_NO2")

nuts_data <- filter(Jericho_enviro_and_bio_long,
                    variable %in% target)


nutrients_plot_line <- ggplot(nuts_data,
                              aes(Date, 
                                  as.numeric(value),
                                  group = variable,
                                  colour = variable)) +
  season_line+
 # spring_bloom_line+
  geom_line(size=1.5) +
  geom_point(colour = "grey50") +
  date_scaling+
  scale_y_continuous(limits = c(0,80),
                     expand = c(0,0))+
  theme(legend.position="right")+
  ggtitle(NULL)+
  ylab("Nutrients (um/L)")+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())


pdf(paste0(figures_dir,"Jericho_nutrients.pdf"),
    width = 15,
    height = 11,
    onefile = FALSE
)
nutrients_plot_line
dev.off()


## Plot chlorophyll as a small subsection ####
chla_data <- filter(Jericho_enviro_and_bio_long,
                    variable == "Average_chl_a")

chla_plot_full_line <- ggplot(chla_data,
                              aes(Date,
                                  as.numeric(value),
                                  group = variable)) +
  season_line+
  # spring_bloom_line+
  geom_line(size=1.5) +
  geom_point(colour = "grey50") +
  date_scaling+
  scale_y_continuous(limits = c(0,50),
                     expand = c(0,0))+ 
  ggtitle(NULL)+
  xlab(NULL)+
  ylab(expression(mu * "g/L"))+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8,
                                   colour = "black"),
        plot.background = element_rect(size=2,
                                       linetype="solid",
                                       color="grey50"),
        axis.text.y = element_text(size = 8,
                                   colour = "black"))


chla_plot_mini_line <- ggplot(chla_data,
                              aes(Date,
                                  as.numeric(value),
                                  group = variable)) +
  season_line+
  # spring_bloom_line+
  geom_line(size=1.5) +
  geom_point(colour = "grey50",
             size=2) +
  date_scaling+
  #scale_y_continuous(limits = c(0,6),
  #                   expand = c(0,0))+
  coord_cartesian(ylim=c(0,6))+
  ggtitle(NULL)+
  ylab(expression("Chlorophyll a (" * mu * "g/L)"))+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())+
  xlab("\nDate")+
  theme(axis.title.x =element_text(vjust=0.5,
                                   size=25))


vp <- grid::viewport(width = 0.5, 
                     height = 0.4,
                     x = 0.63,
                     y = 0.575, 
                     just = c("right",
                              "bottom"))

print(chla_plot_mini_line) 
print(chla_plot_full_line,
      vp = vp)




pdf(paste0(figures_dir,"Jericho_chl_a%03d.pdf"),
    width = 15,
    height = 11,
    onefile = FALSE
)
#grid.arrange(chla_plot_full_line,
#             chla_plot_mini_line,
#             ncol = 1)
print(chla_plot_mini_line) 
print(chla_plot_full_line,
      vp = vp)

dev.off()

pdf(paste0(figures_dir,"Jericho_secchi.pdf"),
    width = 15,
    height = 11,
    onefile = FALSE
)

sechi_plot_full_line <- ggplot(droplevels(subset(Jericho_enviro_and_bio_long,
                                                 (variable == "Secchi_disk_disappears"))),
                               aes(x = Date,
                                   y = as.numeric(value),
                                   group = variable))+
  season_line+
  # spring_bloom_line+
  geom_point(position = position_jitter(width = .1),
             size = 6) +
  date_scaling+
  scale_y_reverse(limits = c(5,0)) +
  ggtitle(NULL) +
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())
sechi_plot_full_line

dev.off()


## Dissolved oxygen
pdf(paste0(figures_dir,"Jericho_dO.pdf"),
    width = 15, 
    height = 11,
    onefile = FALSE
)

DO_plot_full_line <- ggplot(droplevels(subset(Jericho_enviro_and_bio_long,
                                              (variable == "Dissolved_oxygen_percent"))),
                            aes(x = Date,
                                y = as.numeric(value),
                                group = variable))+
  season_line+
  # spring_bloom_line+
  geom_line(colour = path_colour,
            size=1.5) +
  geom_point(colour = line_colour) +
  date_scaling+
  #scale_y_reverse(limits = c(5,0)) +
  ggtitle(NULL) +
  ylab("Dissolved Oxygen (%)")+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank())
DO_plot_full_line 

dev.off()


## Generate proper environmental table for Jericho
##
## Export only useful environmental data, so that I don't have to keep using the same to import and reimport in different scripts.
## Author: Julia Gustavsen (j.gustavsen@gmail.com)
library(plotrix)
library(lubridate)
library(geosphere) # daylength package
library(reshape2)
library(e1071)
library(MASS)

Jericho_Pier_all_data <- read.csv("../data/JerichoDataCollected.csv", header=TRUE)

# Convert the date to something more sensible and so that it might help with the panel chart. 
Jericho_Pier_all_data$Date <- as.Date(Jericho_Pier_all_data$Date,
                                      format='%d-%b-%y')

## Jericho spring bloom occured in April 6-8 2010. 

#adding season to table
Jericho_Pier_all_data$quarter <- quarter(Jericho_Pier_all_data$Date)
Jericho_Pier_all_data$season <- NA
Jericho_Pier_all_data <- within(Jericho_Pier_all_data, 
                                {season[Date >= "2010-03-22" & Date <= "2010-06-21"] <- "spr"
                                season[Date >= "2010-06-22" & Date <= "2010-09-21"] <- "sum"
                                season[Date >= "2010-09-22" & Date <= "2010-12-21"] <- "fall"
                                season[Date >= "2010-12-22" & Date <= "2011-03-21"] <- "win"
                                season[Date >= "2011-03-22" & Date <= "2011-06-21"] <- "spr"
                                season[Date >= "2011-06-22" & Date <= "2011-09-21"] <- "sum"
})


## Remove the bloom data for the analysis

Jericho_Pier_all_data <- Jericho_Pier_all_data[-c(33:38),]

################# Averages and standard errors ###############
## First want to do the averages of the viral abundance data, bacterial, chla, nutrients, etc and calculate the standard error of the mean. 

Jericho_Pier_all_data$Average_viral_abundance <- cbind(rowMeans(Jericho_Pier_all_data[,c("Raw_Viruses_rep_A","Raw_Viruses_rep_B")],na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_viral_abundance <- apply(Jericho_Pier_all_data[,c("Raw_Viruses_rep_A", "Raw_Viruses_rep_B")], 1, std.error, na.rm = TRUE) 
Jericho_Pier_all_data$Average_bacterial_abundance <- cbind(rowMeans(Jericho_Pier_all_data[,c("Raw_Bacteria_Rep_A", "Raw_Bacteria_Rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_bacterial_abundance <- apply(Jericho_Pier_all_data[,c("Raw_Bacteria_Rep_A", "Raw_Bacteria_Rep_B")], 1, std.error, na.rm = TRUE) 

Jericho_Pier_all_data$Average_chl_a <- cbind(rowMeans(Jericho_Pier_all_data[,c("Chl_a_rep_A", "Chl_a_rep_B", "Chl_a_rep_C")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_chl_a <- apply(Jericho_Pier_all_data[,c("Chl_a_rep_A", "Chl_a_rep_B", "Chl_a_rep_C")], 1, std.error, na.rm = TRUE) 

Jericho_Pier_all_data$Average_PO4 <- cbind(rowMeans(Jericho_Pier_all_data[,c("PO4_rep_A", "PO4_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_PO4 <- apply(Jericho_Pier_all_data[,c("PO4_rep_A", "PO4_rep_B")], 1, std.error, na.rm = TRUE) 
Jericho_Pier_all_data$Average_SiO2 <-cbind(rowMeans(Jericho_Pier_all_data[,c("SiO2_rep_A", "SiO2_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_SiO2 <- apply(Jericho_Pier_all_data[,c("SiO2_rep_A", "SiO2_rep_B")], 1, std.error, na.rm = TRUE)
Jericho_Pier_all_data$Average_NO3_NO2 <- cbind(rowMeans(Jericho_Pier_all_data[,c("NO3_NO2_rep_A", "NO3_NO2_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_NO3_NO2 <- apply(Jericho_Pier_all_data[,c("NO3_NO2_rep_A", "NO3_NO2_rep_B")], 1, std.error, na.rm = TRUE)

############# Missing Data ###############
## Deal with missing data by using the duplicated temperature data from the YSI compared to the thermometer and the refractometer.

## Fill in temp for YSI where it was missing from the thermometer data. 

find_na <- is.na(Jericho_Pier_all_data$Temperature_YSI)

Jericho_Pier_all_data$Temperature_YSI[is.na(Jericho_Pier_all_data$Temperature_YSI)] <- Jericho_Pier_all_data$Temperature_thermometer[find_na]

### Do with salinity
## use relationship to fill in missing values
test_equation <- line(Jericho_Pier_all_data$Salinity_ppt_refractometer,Jericho_Pier_all_data$Salinity_ppt_YSI)

plot(Jericho_Pier_all_data$Salinity_ppt_refractometer,Jericho_Pier_all_data$Salinity_ppt_YSI)
abline(coef(test_equation))


transformed_refractometer_data <- Jericho_Pier_all_data$Salinity_ppt_refractometer*coef(test_equation)[2] + coef(test_equation)[2]

find_na <- is.na(Jericho_Pier_all_data$Salinity_ppt_YSI)

Jericho_Pier_all_data$Salinity_ppt_YSI[find_na] <- transformed_refractometer_data[find_na]

############## Write out new table #############

columns_useful_for_analyses <- c("Date",
# "Samples_collected_by.",
# "Time_start",
# "Time_end",
# "Weather",
# "Air_Temperature",
# "Water_Depth",
# "Temperature_thermometer",
"Temperature_YSI",
# "Salinity_ppt_refractometer",
"Salinity_ppt_YSI",
# "Conductivity_mS",
"Dissolved_oxygen_percent",
# "Dissolved_oxgen_mg_per_L",
 "Secchi_disk_disappears",
 "Secchi_disk_reappears",
# "Chl_a_rep_A",
# "Chl_a_rep_B",                    
# "Chl_a_rep_C",
"pH",
"Tide_height",
# "wind_speed",
"VC_number",
# "Prefiltration_date",
# "Prefiltration_Time_began",
# "Prefiltration_Time_ended",
# "Glass_fiber_type",
# "Glass_fiber_filter_number",
# "Glass_fiber_filter_observations",
# "PVDF_type",
# "PVDF_number",
# "PVDF_observations",
# "Ultrafiltration_date",
# "VC_initial_vol",
# "VC_final_vol",
# "Cartridge_name",
# "VC_time_began",
# "VC_time_ended",
# "Pump_type",
# "pump_speed",
# "Cartridge_back_pressure",
# "VC_comments",
# "UF_virus",
 "Raw_Bacteria_Rep_A",
 "Raw_Bacteria_Rep_B",
# "High_Fluo_Raw_Bacteria_Rep_A",
# "High_Fluo_Raw_Bacteria_rep_B",
# "Low_Fluo_Raw_Bacteria_rep_A",
# "Low_Fluo_Raw_Bacteria_rep_B",
 "Raw_Viruses_rep_A",
 "Raw_Viruses_rep_B",
# "High_Fluo_Raw_Viruses_Rep_A",
# "High_Fluo_Raw_Viruses_Rep_B",
# "Low_Fluo_Raw_Viruses_rep_A",     
# "Low_Fluo_Raw_Viruses_rep_B",
# "Filtrate_Viruses_rep_A",         
# "Filtrate_viruses_rep_B",
# "Filtrate_Bacteria_rep_A",
# "Filtrate_Bacteria_rep_B",
# "VC_viruses_rep_A",               
# "VC_viruses_rep_B",
 "PO4_rep_A",
 "PO4_rep_B",
 "SiO2_rep_A",                     
 "SiO2_rep_B",
 "NO3_NO2_rep_A",                  
 "NO3_NO2_rep_B",
"season",
# "quarter",
"Average_viral_abundance",
"Standard_error_viral_abundance", 
"Average_bacterial_abundance",
"Standard_error_bacterial_abundance",
"Average_chl_a",
"Standard_error_chl_a",
"Average_PO4",
"Standard_error_PO4",
"Average_SiO2",
"Standard_error_SiO2",
"Average_NO3_NO2",
"Standard_error_NO3_NO2")   

filtered_Jericho_Pier <- Jericho_Pier_all_data[,columns_useful_for_analyses]



############ Exclude values that are problematic #################

## Viral abundance and bacterial abundance 
## If standard error is greater 50% of the average then value is recoded as NA.

for (x in seq(nrow(filtered_Jericho_Pier))){
 if (is.na(filtered_Jericho_Pier$Standard_error_viral_abundance[x]) ){
  next
 }
 else if (filtered_Jericho_Pier$Standard_error_viral_abundance[x]/filtered_Jericho_Pier$Average_viral_abundance[x] > 0.5 ){
  filtered_Jericho_Pier$Average_viral_abundance[x] <- NA
  filtered_Jericho_Pier$Standard_error_viral_abundance[x] <- NA
 }
}

for (x in seq(nrow(filtered_Jericho_Pier))){
 if (is.na(filtered_Jericho_Pier$Standard_error_bacterial_abundance[x]) ){
  next
 }
 else if (filtered_Jericho_Pier$Standard_error_bacterial_abundance[x]/filtered_Jericho_Pier$Average_bacterial_abundance[x] > 0.5 ){
  filtered_Jericho_Pier$Average_bacterial_abundance[x] <- NA
  filtered_Jericho_Pier$Standard_error_bacterial_abundance[x] <- NA
 }
}

## fix pH
## imported as factor because of text in column
filtered_Jericho_Pier$pH <- as.character(filtered_Jericho_Pier$pH)
filtered_Jericho_Pier$pH[13]  <- "6.5"
filtered_Jericho_Pier$pH <- as.numeric(filtered_Jericho_Pier$pH)

## add in Daylength uisng Vancouver's latidue and the year month date
filtered_Jericho_Pier$day_length <- daylength(49.2827, filtered_Jericho_Pier$Date)
filtered_Jericho_Pier$month_number <- as.numeric(format(filtered_Jericho_Pier$Date, "%m"))

## fix secchi disk
## imported as factor because of text in column
filtered_Jericho_Pier$Secchi_disk_disappears<- as.character(filtered_Jericho_Pier$Secchi_disk_disappears)
filtered_Jericho_Pier$Secchi_disk_disappears[7]  <- "4.9"
filtered_Jericho_Pier$Secchi_disk_disappears <- as.numeric(filtered_Jericho_Pier$Secchi_disk_disappears)
## was still visible at the bottom
filtered_Jericho_Pier$Secchi_disk_reappears[7] <- 4.9

summary(filtered_Jericho_Pier)
write.csv(filtered_Jericho_Pier, "../results/Jericho_data_for_env_analysis.csv")

########## Make mean imputed data for analysis #####
### Do mean imputation ###

mean_imputed_Jericho_Pier <- filtered_Jericho_Pier

summary(mean_imputed_Jericho_Pier)

mean_imputed_Jericho_Pier$Average_viral_abundance[is.na(mean_imputed_Jericho_Pier$Average_viral_abundance)] = mean(mean_imputed_Jericho_Pier$Average_viral_abundance, na.rm=TRUE)

mean_imputed_Jericho_Pier$Average_bacterial_abundance[is.na(mean_imputed_Jericho_Pier$Average_bacterial_abundance)] = mean(mean_imputed_Jericho_Pier$Average_bacterial_abundance, na.rm=TRUE)
mean_imputed_Jericho_Pier$Average_chl_a[is.na(mean_imputed_Jericho_Pier$Average_chl_a)] = mean(mean_imputed_Jericho_Pier$Average_chl_a, na.rm=TRUE)
mean_imputed_Jericho_Pier$Average_PO4[is.na(mean_imputed_Jericho_Pier$Average_PO4)] = mean(mean_imputed_Jericho_Pier$Average_PO4, na.rm=TRUE)
mean_imputed_Jericho_Pier$Average_SiO2[is.na(mean_imputed_Jericho_Pier$Average_SiO2)] = mean(mean_imputed_Jericho_Pier$Average_SiO2, na.rm=TRUE)
mean_imputed_Jericho_Pier$Average_NO3_NO2[is.na(mean_imputed_Jericho_Pier$Average_NO3_NO2)] = mean(mean_imputed_Jericho_Pier$Average_NO3_NO2, na.rm=TRUE)
mean_imputed_Jericho_Pier$Temperature_YSI[is.na(mean_imputed_Jericho_Pier$Temperature_YSI)] = mean(mean_imputed_Jericho_Pier$Temperature_YSI, na.rm=TRUE)
mean_imputed_Jericho_Pier$Salinity_ppt_YSI[is.na(mean_imputed_Jericho_Pier$Salinity_ppt_YSI)] = mean(mean_imputed_Jericho_Pier$Salinity_ppt_YSI, na.rm=TRUE)
mean_imputed_Jericho_Pier$Dissolved_oxygen_percent[is.na(mean_imputed_Jericho_Pier$Dissolved_oxygen_percent)] = mean(mean_imputed_Jericho_Pier$Dissolved_oxygen_percent, na.rm=TRUE)
mean_imputed_Jericho_Pier$pH[is.na(mean_imputed_Jericho_Pier$pH)] = mean(mean_imputed_Jericho_Pier$pH, na.rm=TRUE)
summary(mean_imputed_Jericho_Pier)

write.csv(mean_imputed_Jericho_Pier, "../results/Jericho_env_data_mean_imputed.csv")


######### Transform data for analysis ######

## keep certain data to standardise
to_keep <- c("Average_viral_abundance",
             "Average_bacterial_abundance",
             "Average_chl_a", 
             "Average_PO4",
             "Average_SiO2", 
             "Average_NO3_NO2",
             "Temperature_YSI",
             "Salinity_ppt_YSI", 
             "Dissolved_oxygen_percent",
             "pH",
             "Date",
             "season",
             "day_length",
             "month_number")

subset_to_transform <- subset(mean_imputed_Jericho_Pier,
                              select = colnames(mean_imputed_Jericho_Pier) %in% to_keep)

## check for transformations
long_Jericho<- melt(subset_to_transform,id.vars="Date")

## skewness is 0

## Salinity - log seems to work well
hist(subset_to_transform$Salinity_ppt_YSI)
skewness(subset_to_transform$Salinity_ppt_YSI)
skewness(log(subset_to_transform$Salinity_ppt_YSI))
hist(log(subset_to_transform$Salinity_ppt_YSI))
skewness(sqrt(subset_to_transform$Salinity_ppt_YSI))

subset_to_transform$Salinity_ppt_YSI <- log(subset_to_transform$Salinity_ppt_YSI)

## Temperature, no transformations seem to improve the distribution
hist(subset_to_transform$Temperature_YSI)
skewness(subset_to_transform$Temperature_YSI)
skewness(log(subset_to_transform$Temperature_YSI))
hist(log(subset_to_transform$Temperature_YSI))
skewness(sqrt(subset_to_transform$Temperature_YSI))
skewness(1/(subset_to_transform$Temperature_YSI))
skewness(1/sqrt(subset_to_transform$Temperature_YSI))
boxcox(subset_to_transform$Temperature_YSI~1)
skewness(subset_to_transform$Temperature_YSI^(1.2))
hist(subset_to_transform$Temperature_YSI^(1.2))
skewness((subset_to_transform$Temperature_YSI)-7)


## Dissoled Oxygen looks ok as is
hist(subset_to_transform$Dissolved_oxygen_percent)
skewness(subset_to_transform$Dissolved_oxygen_percent)

## pH  nothing works
hist(subset_to_transform$pH)
skewness(subset_to_transform$pH)
skewness(log(subset_to_transform$pH))
hist(log(subset_to_transform$pH))
skewness(sqrt(subset_to_transform$pH))
skewness(log10(subset_to_transform$pH))
skewness(1/(subset_to_transform$pH))
skewness(1/sqrt(subset_to_transform$pH))
boxcox(subset_to_transform$pH~1)
skewness(subset_to_transform$pH^(-0.3))
hist(subset_to_transform$pH^(-0.3))
skewness((subset_to_transform$pH)-7)


## Average viral abundance  log
hist(subset_to_transform$Average_viral_abundance)
skewness(subset_to_transform$Average_viral_abundance)
skewness(log(subset_to_transform$Average_viral_abundance))
hist(log(subset_to_transform$Average_viral_abundance))
skewness(sqrt(subset_to_transform$Average_viral_abundance))
skewness(log10(subset_to_transform$Average_viral_abundance))
skewness(1/(subset_to_transform$Average_viral_abundance))
skewness(1/sqrt(subset_to_transform$Average_viral_abundance))
boxcox(subset_to_transform$Average_viral_abundance~1)
skewness(subset_to_transform$Average_viral_abundance^(-0.6))
hist(subset_to_transform$Average_viral_abundance^(-0.6))
skewness((subset_to_transform$Average_viral_abundance)-7)

subset_to_transform$Average_viral_abundance <- log(subset_to_transform$Average_viral_abundance)

## Average bacterial abundance  box-cox transformation best
hist(subset_to_transform$Average_bacterial_abundance)
skewness(subset_to_transform$Average_bacterial_abundance)
skewness(log(subset_to_transform$Average_bacterial_abundance))
hist(log(subset_to_transform$Average_bacterial_abundance))
skewness(sqrt(subset_to_transform$Average_bacterial_abundance))
skewness(log10(subset_to_transform$Average_bacterial_abundance))
skewness(1/(subset_to_transform$Average_bacterial_abundance))
skewness(1/sqrt(subset_to_transform$Average_bacterial_abundance))
boxcox(subset_to_transform$Average_bacterial_abundance~1)
skewness(subset_to_transform$Average_bacterial_abundance^(0.2))
hist(subset_to_transform$Average_bacterial_abundance^(0.2))

subset_to_transform$Average_bacterial_abundance <- (subset_to_transform$Average_bacterial_abundance)^(0.2)

## Average chla  box-cox transformation best
hist(subset_to_transform$Average_chl_a)
skewness(subset_to_transform$Average_chl_a)
skewness(log(subset_to_transform$Average_chl_a))
hist(log(subset_to_transform$Average_chl_a))
skewness(sqrt(subset_to_transform$Average_chl_a))
skewness(log10(subset_to_transform$Average_chl_a))
skewness(1/(subset_to_transform$Average_chl_a))
skewness(1/sqrt(subset_to_transform$Average_chl_a))
boxcox(subset_to_transform$Average_chl_a~1)
skewness(subset_to_transform$Average_chl_a^(-0.2))
hist(subset_to_transform$Average_chl_a^(-0.2))

subset_to_transform$Average_chl_a <- (subset_to_transform$Average_chl_a)^(-0.2)

## Ave PO4 needs no transformation
hist(subset_to_transform$Average_PO4)
skewness(subset_to_transform$Average_PO4)
skewness(log(subset_to_transform$Average_PO4))
hist(log(subset_to_transform$Average_PO4))
skewness(sqrt(subset_to_transform$Average_PO4))
skewness(log10(subset_to_transform$Average_PO4))
skewness(1/(subset_to_transform$Average_PO4))
skewness(1/sqrt(subset_to_transform$Average_PO4))
boxcox(subset_to_transform$Average_PO4~1)
skewness(subset_to_transform$Average_PO4^(-0.2))
hist(subset_to_transform$Average_PO4^(-0.2))


## Average Silicate box-cox
hist(subset_to_transform$Average_SiO2)
skewness(subset_to_transform$Average_SiO2)
skewness(log(subset_to_transform$Average_SiO2))
hist(log(subset_to_transform$Average_SiO2))
skewness(sqrt(subset_to_transform$Average_SiO2))
skewness(log10(subset_to_transform$Average_SiO2))
skewness(1/(subset_to_transform$Average_SiO2))
skewness(1/sqrt(subset_to_transform$Average_SiO2))
boxcox(subset_to_transform$Average_SiO2~1)
skewness(subset_to_transform$Average_SiO2^(0.8))
hist(subset_to_transform$Average_SiO2^(0.8))

subset_to_transform$Average_SiO2 <- (subset_to_transform$Average_SiO2)^(0.8)

## Average NO3+NO2 again accoring to skewness it is nothing, but the distribution is inverted
hist(subset_to_transform$Average_NO3_NO2)
skewness(subset_to_transform$Average_NO3_NO2)
shapiro.test(subset_to_transform$Average_NO3_NO2)
skewness(log(subset_to_transform$Average_NO3_NO2))
hist(log(subset_to_transform$Average_NO3_NO2))
skewness(sqrt(subset_to_transform$Average_NO3_NO2))
skewness(log10(subset_to_transform$Average_NO3_NO2))
skewness(1/(subset_to_transform$Average_NO3_NO2))
skewness(1/sqrt(subset_to_transform$Average_NO3_NO2))
#boxcox(subset_to_transform$Average_NO3_NO2~1)
skewness(subset_to_transform$Average_NO3_NO2^(0.5))
hist(subset_to_transform$Average_NO3_NO2^(0.5))

summary(subset_to_transform)

write.csv(subset_to_transform, "../results/Jericho_env_data_mean_imputed_transformed.csv")

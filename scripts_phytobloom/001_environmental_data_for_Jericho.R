## Generate proper environmental table for Jericho

## Export only useful data, so that I don't have to keep using the same to import and reimport in different scripts. 
library(plotrix)
library(lubridate)


Jericho_Pier_all_data <- read.csv("../../JerichoTimeSeries/data/JerichoDataCollected.csv", header=TRUE)

# Convert the date to something more sensible and so that it might help with the panel chart. 
Jericho_Pier_all_data$Date <- as.Date(Jericho_Pier_all_data$Date, format='%d-%b-%y')

library(lubridate)
#adding season to table
Jericho_Pier_all_data$quarter <- quarter(Jericho_Pier_all_data$Date)
Jericho_Pier_all_data$season <- NA
Jericho_Pier_all_data <- within(Jericho_Pier_all_data, {season[quarter == 1] <- "spr"
                                      season[quarter == 2] <- "sum"
                                      season[quarter == 3] <- "fall"
                                      season[quarter == 4] <- "win"
})


#First want to do the averages of the viral abundance data, bacterial, chla, nutrients, etc and calculate the standard error of the mean. 
std <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

Jericho_Pier_all_data$Average_viral_abundance <-cbind(rowMeans(Jericho_Pier_all_data[,c("Raw_Viruses_rep_A", "Raw_Viruses_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_viral_abundance <- apply(Jericho_Pier_all_data[,c("Raw_Viruses_rep_A", "Raw_Viruses_rep_B")],
                                                              1,
                                                              std.error,
                                                              na.rm = TRUE) 
Jericho_Pier_all_data$Average_bacterial_abundance <-cbind(rowMeans(Jericho_Pier_all_data[,c("Raw_Bacteria_Rep_A", "Raw_Bacteria_Rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_bacterial_abundance <- apply(Jericho_Pier_all_data[,c("Raw_Bacteria_Rep_A", "Raw_Bacteria_Rep_B")],
                                                                  1,
                                                                  std.error,
                                                                  na.rm = TRUE) 


Jericho_Pier_all_data$Average_chl_a <-cbind(rowMeans(Jericho_Pier_all_data[,c("Chl_a_rep_A", "Chl_a_rep_B", "Chl_a_rep_C")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_chl_a <- apply(Jericho_Pier_all_data[,c("Chl_a_rep_A", "Chl_a_rep_B", "Chl_a_rep_C")],
                                                    1,
                                                    std.error,
                                                    na.rm = TRUE) 

Jericho_Pier_all_data$Average_PO4 <-cbind(rowMeans(Jericho_Pier_all_data[,c("PO4_rep_A", "PO4_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_PO4 <- apply(Jericho_Pier_all_data[,c("PO4_rep_A", "PO4_rep_B")],
                                                  1,
                                                  std.error,
                                                  na.rm = TRUE) 

Jericho_Pier_all_data$Average_SiO2 <-cbind(rowMeans(Jericho_Pier_all_data[,c("SiO2_rep_A", "SiO2_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_SiO2 <- apply(Jericho_Pier_all_data[,c("SiO2_rep_A", "SiO2_rep_B")],
                                                   1,
                                                   std.error,
                                                   na.rm = TRUE)

Jericho_Pier_all_data$Average_NO3_NO2 <-cbind(rowMeans(Jericho_Pier_all_data[,c("NO3_NO2_rep_A", "NO3_NO2_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_NO3_NO2 <- apply(Jericho_Pier_all_data[,c("NO3_NO2_rep_A", "NO3_NO2_rep_B")],
                                                      1,
                                                      std.error,
                                                      na.rm = TRUE)



## Fill in temp for YSI where it was missing from the thermometer data. 

find_na <- is.na(Jericho_Pier_all_data$Temperature_YSI)

Jericho_Pier_all_data$Temperature_YSI[is.na(Jericho_Pier_all_data$Temperature_YSI)] <- Jericho_Pier_all_data$Temperature_thermometer[find_na]


### Do with salinity
## what is the equation?
## use relationship to fill in missing values
test_equation <- line(Jericho_Pier_all_data$Salinity_ppt_refractometer,Jericho_Pier_all_data$Salinity_ppt_YSI)

plot(Jericho_Pier_all_data$Salinity_ppt_refractometer,Jericho_Pier_all_data$Salinity_ppt_YSI)
abline(coef(test_equation))


transformed_refractometer_data <- Jericho_Pier_all_data$Salinity_ppt_refractometer*coef(test_equation)[2] + coef(test_equation)[2]


find_na <- is.na(Jericho_Pier_all_data$Salinity_ppt_YSI)

Jericho_Pier_all_data$Salinity_ppt_YSI[find_na] <- transformed_refractometer_data[find_na]




## write out useful environmental data
## need to write out the secchi disk transparency measurements 

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

## write out talbe for use in analyses

filtered_Jericho_Pier <- Jericho_Pier_all_data[,columns_useful_for_analyses]

filtered_Jericho_Pier$pH <- as.character(filtered_Jericho_Pier$pH)
filtered_Jericho_Pier$pH[13]  <- "6.5"
filtered_Jericho_Pier$pH <- as.numeric(filtered_Jericho_Pier$pH)



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

write.csv(filtered_Jericho_Pier, "../results/Jericho_data_for_env_analysis_all_time_series.csv")




### only bloom data. ####

## Subset the  data so that it is only the bloom samples
Jericho_Pier_all_data <- Jericho_Pier_all_data[32:39,]

#First want to do the averages of the viral abundance data, bacterial, chla, nutrients, etc and calculate the standard error of the mean. 
std <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

Jericho_Pier_all_data$Average_viral_abundance <-cbind(rowMeans(Jericho_Pier_all_data[,c("Raw_Viruses_rep_A", "Raw_Viruses_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_viral_abundance <- apply(Jericho_Pier_all_data[,c("Raw_Viruses_rep_A", "Raw_Viruses_rep_B")],
                                                                  1,
                                                                  std.error,
                                                                  na.rm = TRUE) 
Jericho_Pier_all_data$Average_bacterial_abundance <-cbind(rowMeans(Jericho_Pier_all_data[,c("Raw_Bacteria_Rep_A", "Raw_Bacteria_Rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_bacterial_abundance <- apply(Jericho_Pier_all_data[,c("Raw_Bacteria_Rep_A", "Raw_Bacteria_Rep_B")],
                  1,
                  std.error,
                  na.rm = TRUE) 


Jericho_Pier_all_data$Average_chl_a <-cbind(rowMeans(Jericho_Pier_all_data[,c("Chl_a_rep_A", "Chl_a_rep_B", "Chl_a_rep_C")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_chl_a <- apply(Jericho_Pier_all_data[,c("Chl_a_rep_A", "Chl_a_rep_B", "Chl_a_rep_C")],
                                                                  1,
                                                                  std.error,
                                                                  na.rm = TRUE) 

Jericho_Pier_all_data$Average_PO4 <-cbind(rowMeans(Jericho_Pier_all_data[,c("PO4_rep_A", "PO4_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_PO4 <- apply(Jericho_Pier_all_data[,c("PO4_rep_A", "PO4_rep_B")],
                                                    1,
                                                    std.error,
                                                    na.rm = TRUE) 

Jericho_Pier_all_data$Average_SiO2 <-cbind(rowMeans(Jericho_Pier_all_data[,c("SiO2_rep_A", "SiO2_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_SiO2 <- apply(Jericho_Pier_all_data[,c("SiO2_rep_A", "SiO2_rep_B")],
                                                  1,
                                                  std.error,
                                                  na.rm = TRUE)

Jericho_Pier_all_data$Average_NO3_NO2 <-cbind(rowMeans(Jericho_Pier_all_data[,c("NO3_NO2_rep_A", "NO3_NO2_rep_B")], na.rm = TRUE))
Jericho_Pier_all_data$Standard_error_NO3_NO2 <- apply(Jericho_Pier_all_data[,c("NO3_NO2_rep_A", "NO3_NO2_rep_B")],
                                                   1,
                                                   std.error,
                                                   na.rm = TRUE)



## write out useful environmental data
## need to write out the secchi disk transparency measurements 

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

## write out talbe for use in analyses

filtered_Jericho_Pier <- Jericho_Pier_all_data[,columns_useful_for_analyses]

## correct values that look problematic

## viral abundance 
## how to throw out data. If standard error is greater than something replace the data point with x? 42450000 31050000 this is 70%
## 19165000 13435000 70%
## 31800000 14600000
## greater than 50% of the average then I reject it??




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


write.csv(filtered_Jericho_Pier, "../results/Jericho_data_for_env_analysis.csv")

 read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)

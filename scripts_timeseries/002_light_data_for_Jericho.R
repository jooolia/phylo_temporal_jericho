### pull in solar radiation data from UBC
## found at http://ibis.geog.ubc.ca/~achristn/infrastructure/climate-data/STSDA2-01-2011-to-12-2011.html

irradiation_data <- read.csv("../data/ibis_shortwave_irradiance_2010_2011.csv")


## want to get a good date so that I can pull in the data in to the Jericho Pier data
irradiation_data$Date <-as.Date(paste(irradiation_data$Year, irradiation_data$Month, irradiation_data$Day), "%Y %b %d")

write.csv(irradiation_data, file = "../data/ibis_shortwave_irradiance_2010_2011_date_fixed.csv")


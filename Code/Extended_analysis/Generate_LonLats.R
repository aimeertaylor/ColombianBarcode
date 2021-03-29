# Distances extended
rm(list = ls())
load("../../RData/metadata_extended.RData")
load("../../RData/geo_dist_info_cities.RData")
library(measurements) # For conv_unit

# Copied from Generate_geo_dist_info.R
# Import coordinate data
LonLat_city <- read.table('../../TxtData/Geo_coordinates.txt', skip = 21, header = TRUE, sep = ' ', row.names = )
rownames(LonLat_city) = LonLat_city$City

# Convert Longitude Latitude format from degrees, minutes and seconds to decimal degrees:
LonLat_city$Latitude = gsub('°', ' ', LonLat_city$Latitude)
LonLat_city$Longitude = gsub('°', ' ', LonLat_city$Longitude)
LonLat_city$Latitude = gsub('’', ' ', LonLat_city$Latitude)
LonLat_city$Longitude = gsub('’', ' ', LonLat_city$Longitude)
LonLat_city$Latitude = gsub('”N', '', LonLat_city$Latitude)
LonLat_city$Longitude = gsub('”W', '', LonLat_city$Longitude)
LonLat_city$Latitude = as.numeric(conv_unit(LonLat_city$Latitude, from = 'deg_min_sec', to = 'dec_deg'))
LonLat_city$Longitude = as.numeric(conv_unit(LonLat_city$Longitude, from = 'deg_min_sec', to = 'dec_deg'))

# Extract mean values from metadata where available
meanLonLat <- unique(metadata[, c("City", "LON", "LAT")]) %>%
  group_by(City) %>% 
  summarise(Latitude = mean(LAT, na.rm = T),
            Longitude = mean(LON, na.rm = T)) %>% 
  as.data.frame()
rownames(meanLonLat) <- meanLonLat$City

# Join 
LonLats <- rbind(meanLonLat[!meanLonLat$City %in% LonLat_city$City,], LonLat_city)

# Fill in missing values by hand
LonLats["Esmeraldas",c("Latitude","Longitude")]  = c(0.96817, -79.65172) 
LonLats["Orellana",c("Latitude","Longitude")]  = c(-0.71126, -77.15426)
LonLats["SanLorenzo",c("Latitude","Longitude")]  = c(1.26947, -78.84392)
LonLats["Sucumbios",c("Latitude","Longitude")]  = c(0.08892, -76.88975)
LonLats["TobarDonoso",c("Latitude","Longitude")]  = c(0.99959, -78.42971)

save(LonLats, file = "../../RData/LonLats.RData")
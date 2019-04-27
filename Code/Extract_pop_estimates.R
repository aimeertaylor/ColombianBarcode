###########################################################
# Script to extract world pop estimates 
# Last census in 2005 
# Tutorial following:
# http://www.guilles.website/2017/12/16/exploring-worldpop-gis-data-with-r/
###########################################################

library(rgdal)
library(raster)
library(sp)
library(colorRamps)
library(ggmap)
library(ggplot2)
PLOT = F

# Load Colombia data (data collected 1993 to 2007)
WorldPopCol <- raster("../PopData/Colombia_100m_Population/COL_pph_v2b_2010.tif")

# Shape file
gtmMunis = readOGR("../GISData/COL_adm_shp/COL_adm2.shp", encoding = "utf-8")

# Plot world pop data  (takes a few seconds)
if(PLOT){
  plot(log10(WorldPopCol), colNA="black", 
       col=colorRampPalette(c("#121212", "#123620", "#108650", "#80D6C0", "#DFDF3B"),0.5)(110))
}
 
# Convert adm2 name to character
gtmMunis$NAME_2 = as.character(gtmMunis$NAME_2)

#===============================================
# At the city level 
#===============================================
# Look for city names
gtmMunis$NAME_2[grepl('Guap', gtmMunis$NAME_2)] # "Guapí"
gtmMunis$NAME_2[grepl('Tumac', gtmMunis$NAME_2)] # "Tumaco"
gtmMunis$NAME_2[grepl('Buenaventura', gtmMunis$NAME_2)] # "Buenaventura"
gtmMunis$NAME_2[grepl('Quibd', gtmMunis$NAME_2)] # "Quibdó"
gtmMunis$NAME_2[grepl('Tad', gtmMunis$NAME_2)] # "Tadó"

# City names
sites = c("Guapí", "Tumaco", "Buenaventura", "Quibdó", "Tadó")

# Population per site  
populations = sapply(sites, function(site){
  sum(extract(WorldPopCol, gtmMunis[gtmMunis$NAME_2==site,], fun = sum, na.rm = TRUE))
})

# Name without accents
names(populations) = c("Guapi", "Tumaco", "Buenaventura", "Quibdo", "Tado")
save(populations, file = '~/Documents/ColombianBarcode/RData/populations_cities.RData')


#===============================================
# At the state level (takes a long time!)
#===============================================
# Look for city names
gtmMunis$NAME_1[grepl('Nari', gtmMunis$NAME_1)] # Nariño 
gtmMunis$NAME_1[grepl('Cauca', gtmMunis$NAME_1)] # Cauca and Valle del Cauca (latter is) 
gtmMunis$NAME_1[grepl('Valle', gtmMunis$NAME_1)] # Valle del Cauca
gtmMunis$NAME_1[grepl('Choc', gtmMunis$NAME_1)] # Chocó

# City names
sites = c("Nariño", "Cauca", "Valle del Cauca", "Chocó")

# Population per site  
populations = sapply(sites, function(site){
  sum(extract(WorldPopCol, gtmMunis[gtmMunis$NAME_1 ==site,], fun = sum, na.rm = TRUE))
})

# Name without accents
names(populations) = c("Narino", "Cauca", "Valle", "Choco")
save(populations, file = '~/Documents/ColombianBarcode/RData/populations_states.RData')

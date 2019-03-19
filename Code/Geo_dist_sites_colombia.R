############################################################
# Script to calculate distances between sites 
#
# To-do list
# Put gcd.hf.R into my own library
# Check lat and long on the map
############################################################

rm(list = ls())
source('/Users/aimeet/Documents/BroadLaptop/RFunctions/gcd.hf.R') 
require(gtools)
library(measurements)

# Import data
LonLat <- read.table('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/TxtData/Geo_coordinates.txt', 
                   skip = 21, header = TRUE, sep = ' ')

# Convert Longitude Latitude format from degrees, minutes and seconds to decimal degrees:
LonLat$Latitude = gsub('°', ' ', LonLat$Latitude)
LonLat$Longitude = gsub('°', ' ', LonLat$Longitude)
LonLat$Latitude = gsub('’', ' ', LonLat$Latitude)
LonLat$Longitude = gsub('’', ' ', LonLat$Longitude)
LonLat$Latitude = gsub('”N', '', LonLat$Latitude)
LonLat$Longitude = gsub('”W', '', LonLat$Longitude)
LonLat$Latitude = as.numeric(conv_unit(LonLat$Latitude, from = 'deg_min_sec', to = 'dec_deg'))
LonLat$Longitude = as.numeric(conv_unit(LonLat$Longitude, from = 'deg_min_sec', to = 'dec_deg'))

write.table(LonLat, file = '/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/TxtData/LonLat_dec_deg.txt', 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Calculate distances
sites <- LonLat$City
site_combinations <- combinations(n = length(sites), r = 2, as.character(sites))
pairwise_site_distance <- data.frame(site_combinations)
pairwise_site_distance$distance <- NA

for(i in 1:nrow(site_combinations)){
  site1 <- site_combinations[i,1]
  site2 <- site_combinations[i,2]
  pairwise_site_distance[i, 'distance'] <- gcd.hf(A = as.matrix(LonLat[LonLat$City == site1, c('Longitude', 'Latitude')]), 
                                                  B = as.matrix(LonLat[LonLat$City == site2, c('Longitude', 'Latitude')]))
}

# Sort pairwise distance
pairwise_site_distance <- pairwise_site_distance[sort(pairwise_site_distance$distance, 
                                                      index.return = TRUE)$ix,]
# Extract geo_order
geo_order <- apply(pairwise_site_distance[,c(1,2)], 1, function(x){paste(sort(x), sep = '_', collapse = '_')})

# Make matrix with all distances 
pairwise_site_distance_all <- rbind(pairwise_site_distance, 
                                    data.frame(X1 = pairwise_site_distance$X2,
                                               X2 = pairwise_site_distance$X1,
                                               distance = pairwise_site_distance$distance), 
                                    data.frame(X1 = sites, X2 = sites, distance = 0))

# Make a vector with all distances
names_pairwise_site_distance_all <- apply(pairwise_site_distance_all[,c(1,2)], 1, function(x){paste(x, sep = '_', collapse = '_')})
pairwise_site_distance_all <- pairwise_site_distance_all$distance
names(pairwise_site_distance_all) <- names_pairwise_site_distance_all

geo_dist_info <- list(geo_order = geo_order, 
                      pairwise_site_distance = pairwise_site_distance,
                      pairwise_site_distance_all = pairwise_site_distance_all)

save(geo_dist_info, file = '/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/geo_dist_info.RData')


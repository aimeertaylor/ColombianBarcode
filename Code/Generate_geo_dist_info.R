############################################################
# Script to calculate distances between sites and add gravity 
# estimates
############################################################
rm(list = ls())
source('./gcd.hf.R') # For calculating great circle distance
require(gtools) # For combinations
library(measurements) # For conv_unit

# Import coordinatesl data
LonLat_city <- read.table('../TxtData/Geo_coordinates.txt', skip = 21, header = TRUE, sep = ' ', row.names = )
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

# Make a state level version (where simply take midpoint between Tado and Quidbo)
# "Narino" "Cauca"  "Valle"  "Choco" 
LonLat_state = list("Narino" = LonLat_city[LonLat_city$City == 'Tumaco', c('Latitude', 'Longitude')], 
                   "Cauca" = LonLat_city[LonLat_city$City == 'Guapi', c('Latitude', 'Longitude')], 
                   "Valle" = LonLat_city[LonLat_city$City == 'Buenaventura', c('Latitude', 'Longitude')], 
                   "Choco" = colMeans(LonLat_city[LonLat_city$City == 'Quibdo' | LonLat_city$City == 'Tado', c('Latitude', 'Longitude')]))

# Match formats
LonLat_city = LonLat_city[,c('Latitude', 'Longitude')]
LonLat_state = do.call(rbind,LonLat_state)

# Needed for network plot 
save(LonLat_city, LonLat_state, file = '../RData/LonLat.RData')

# Calculate distances between cities and state
LonLat_list = list(cities = LonLat_city, states = LonLat_state)
for(l in names(LonLat_list)){
  
  LonLat = LonLat_list[[l]]
  sites <- rownames(LonLat)
  site_combinations <- combinations(n = length(sites), r = 2, as.character(sites))
  pairwise_site_distance <- data.frame(site_combinations)
  pairwise_site_distance$distance <- NA
  
  for(i in 1:nrow(site_combinations)){
    site1 <- site_combinations[i,1]
    site2 <- site_combinations[i,2]
    pairwise_site_distance[i, 'distance'] <- gcd.hf(A = as.matrix(LonLat[site1, c('Longitude', 'Latitude')]), 
                                                    B = as.matrix(LonLat[site2, c('Longitude', 'Latitude')]))
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
  

  # Add gravity estimates  
  load(sprintf('../RData/populations_%s.RData', l)) # Populations for gravity estimates 
  
  geo_dist_info$pairwise_site_distance$gravity_estimate = NA
  for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
    x1 = as.character(geo_dist_info$pairwise_site_distance[i,'X1'])
    x2 = as.character(geo_dist_info$pairwise_site_distance[i,'X2'])
    g = populations[x1] * populations[x2] / (geo_dist_info$pairwise_site_distance[i,'distance'])
    geo_dist_info$pairwise_site_distance[i, 'gravity_estimate'] = g
  }
  
  save(geo_dist_info, file = sprintf('../RData/geo_dist_info_%s.RData', l))
}


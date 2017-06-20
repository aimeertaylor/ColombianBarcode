# To-do list
# Update lon/lat with accurate lon/lat
# Remove Quidbo for the time being as we don't have site resolution for Choco in the data at present

# Script to calculate distances between sites 
rm(list = ls())
source('/Users/aimeet/Documents/BroadLaptop/RFunctions/gcd.hf.R')
require(gtools)

# Import data
Lonlat <- read.csv('/Users/aimeet/Documents/BroadLaptop/ColombianData/TempData/HMMDataFormat.csv')
Lonlat <- Lonlat[-2, ] # Remove Quidbo for the time being as we don't have site resolution for Choco in the data at present
sites <- Lonlat$State
site_combinations <- combinations(n = 4, r = 2, as.character(sites))
pairwise_site_distance <- data.frame(site_combinations)
pairwise_site_distance$distance <- NA

for(i in 1:nrow(site_combinations)){
  site1 <- site_combinations[i,1]
  site2 <- site_combinations[i,2]
  pairwise_site_distance[i, 'distance'] <- gcd.hf(A = as.matrix(Lonlat[Lonlat$State == site1, 4:3]), 
                                                  B = as.matrix(Lonlat[Lonlat$State == site2, 4:3]))
}

# Sort pairwise distance
pairwise_site_distance <- pairwise_site_distance[sort(pairwise_site_distance$distance, 
                                                      index.return = TRUE)$ix,]
# Extract geo_order
geo_order <- apply(pairwise_site_distance[,c(1,2)], 1, function(x){
  paste(sort(x), sep = '_', collapse = '_')})

# Make matrix with all distances 
pairwise_site_distance_all <- rbind(pairwise_site_distance, 
                                    data.frame(X1 = pairwise_site_distance$X2,
                                               X2 = pairwise_site_distance$X1,
                                               distance = pairwise_site_distance$distance), 
                                    data.frame(X1 = sites, X2 = sites, distance = 0))

# Make a vector with all distances
names_pairwise_site_distance_all <- apply(pairwise_site_distance_all[,c(1,2)], 1, paste, sep = '_', collapse = '_')
pairwise_site_distance_all <- pairwise_site_distance_all$distance
names(pairwise_site_distance_all) <- names_pairwise_site_distance_all

geo_dist_info <- list(geo_order = geo_order, 
                      pairwise_site_distance = pairwise_site_distance,
                      pairwise_site_distance_all = pairwise_site_distance_all)

save(geo_dist_info, file = '/Users/aimeet/Documents/BroadLaptop/ColombianData/RData/geo_dist_info.RData')


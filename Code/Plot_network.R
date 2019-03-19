############################################################################################
# This script is adapted from Maps_MORU.R
# 1) Map based on fraction of IBD proportions > 0.5 (2 plots, combined in manuscript)
# This map was created using QGIS (QGIS Development Team 2017) and the raster (Hijmans 2016),
# maptools (Bivand and Lewin-Koh 2017) and igraph (Csardi and Nepusz 2006) packages
############################################################################################
rm(list = ls())

library(maptools)
library(raster)
library(igraph)
PDF = F
Zoom = T

###################################################################################
# Colombian Pacific Coast
###################################################################################

# # Load shape files (downloaded from http://gadm.org/)
# Colombia <- getData("GADM",country="Colombia",level=0)
# Panama <- getData("GADM",country="Panama",level=0)
# Venezuela <- getData("GADM",country="Venezuela",level=0)
# Brazil <- getData("GADM",country="Brazil",level=0)
# Peru <- getData("GADM",country="Peru",level=0)
# Ecuador <- getData("GADM",country="Ecuador",level=0)
CP_coast <- readShapeLines("/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/QGIS/CP_coast.shp")

# Load data
load('~/Dropbox/IBD_IBS/PlasmodiumRelatedness/RData/Colombia_mles.RData') # load relatedness estimates
load('../RData/populations.RData') # Load population estimates 
load('../RData/geo_dist_info.RData') # Load distances

# Calculate gravity estimates and add
geo_dist_info$pairwise_site_distance$gravity_estimate = NA
for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
  x1 = as.character(geo_dist_info$pairwise_site_distance[i,'X1'])
  x2 = as.character(geo_dist_info$pairwise_site_distance[i,'X2'])
  g = populations[x1] * populations[x2] / (geo_dist_info$pairwise_site_distance[i,'distance'])^2
  geo_dist_info$pairwise_site_distance[i, 'gravity_estimate'] = g
}

# Vertices
site_comps <- geo_dist_info$geo_order

# Store in which to store edge weights
mylinks <- data.frame(array(dim = c(length(site_comps), 3),
                            dimnames = list(site_comps,c('site1', 'site2','weight'))))
mylinks$site1 <- as.character(geo_dist_info$pairwise_site_distance[,1])
mylinks$site2 <- as.character(geo_dist_info$pairwise_site_distance[,2])
X = geo_dist_info$pairwise_site_distance[,'gravity_estimate']
mylinks[, 'weight_barcode'] <- (X - min(X))/(max(X) - min(X))
  
# Generate net, remove loops
net <- graph.data.frame(mylinks[mylinks$site1 != mylinks$site2,], directed = FALSE)

# Coordinates
lonlat <- read.table('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/TxtData/LonLat_dec_deg.txt', 
                     header = TRUE)
rownames(lonlat) <- lonlat[,'City']
lonlat$lon <- (lonlat$Longitude)
lonlat$lat <- (lonlat$Latitude)
layout_lonlat <- as.matrix(lonlat[rownames(as.matrix(V(net))), c('lon','lat')])

# IBD barcode
E(net)$width <- E(net)$weight_barcode *10
E(net)$curved <- 0.69*c(1,-1.3,1.2,1.4,-1,-0.9,-1,0,0,0)
E(net)$loop.angle <- pi*c(1.68, 2, rep(1,8))
E(net)$color <- adjustcolor('blue', alpha.f = 0.75)

# Change labels, text colour etc.
print(V(net)) # Print order
V(net)$label <- c('Guapi','Buenaventura','Quibdó', 'Tadó','Tumaco')
V(net)$label.color <- "black"
V(net)$label.dist <- c(3, 7, 2, 3, 3)
V(net)$label.degree <- c(pi, 2*pi, 1.5*pi, 2*pi, 0.5*pi)
V(net)$label.cex <- 1.2
V(net)$color <- "black"
V(net)$size <- 8 
V(net)$frame.color <- "white"

# Plot Coast
E(net)$color <- adjustcolor('blue', alpha.f = 0)
if(PDF){pdf(file = './Colombia_coast.pdf', height = 7, width = 6)}
par(mfrow = c(1,1), family = "serif", mar = c(0,0,0,0))
plot(CP_coast, col = 'black', lwd = 2, xlim =  c(-78.5, -76.7))
plot(net, layout=layout_lonlat, add = TRUE, rescale = FALSE)
compassRose(x = -78.5, y = 5.74)
text(x = -78.52596, y = 3.72, labels = '50 km', cex = 1.2, pos = 3)
segments(y0 = 3.7, y1 = 3.7, x0 = -78.75928, x1 = -78.29265, lwd = 2, col = 'black')
if(PDF){dev.off()}

# Plot Net
E(net)$color <- adjustcolor('blue', alpha.f = 0.75)
if(PDF){pdf(file = './Colombia_Network.pdf', height = 7, width = 6)}
par(mfrow = c(1,1), family = "serif", mar = c(0,0,0,0))
plot(CP_coast, col = 'black', lwd = 2, xlim =  c(-78.5, -76.7))
plot(net, layout=layout_lonlat, add = TRUE, rescale = FALSE)
compassRose(x = -78.5, y = 5.74)
text(x = -78.52596, y = 3.72, labels = '50 km', cex = 1.2, pos = 3)
segments(y0 = 3.7, y1 = 3.7, x0 = -78.75928, x1 = -78.29265, lwd = 2, col = 'black')
if(PDF){dev.off()}

# # Plot map
# #png(file = './Colombia_Map.png', height = 7, width = 6, units = 'in', res = 400)
# par(mfrow = c(1,1), family = "serif", mar = c(0,0,0,0))
# plot(Colombia, col = 'gray', border = 'white', ylim = c(-1.5, 9.25), xlim = c(-82, -67))
# plot(Brazil, add = TRUE, col = 'gray', border = 'white')
# plot(Panama, add = TRUE, col = 'gray', border = 'white')
# plot(Peru, add = TRUE, col = 'gray', border = 'white')
# plot(Ecuador, add = TRUE, col = 'gray', border = 'white')
# plot(Venezuela, add = TRUE, col = 'gray', border = 'white')
# plot(CP_coast, col = 'black', lwd = 2, add = TRUE)
# text(x = c(-81.05882, -79.57647, -80.81176, -78.26588, -74.51765, -70.62118, -73.32000, -67.81059), 
#      y = c(8.401712, 10.388241,  4.041895, -1.252799, -2.392804,  8.416126, 3.531259, -2.230062), 
#      labels = c('Panama','Caribbean Sea', 'Pacific Ocean', 'Ecuador', 'Peru', 'Venezuela', 'Colombia', 'Brazil'), 
#      col = c('white', 'darkgray', 'darkgray', 'white', 'white', 'white', 'white', 'white'), cex = c(0.65,rep(1,7)))
# #dev.off()


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
library(sf) # For st_read instead of readShapeLines (returns warning: "use rgdal::readOGR or sf::st_read" )
PDF = F # Set to true to generate pdfs
DWN = T # Set to true to download country files
ZM = T # Set to true to zoomed out plot surrounding countries 
NET = T # Set to true to include network on coast

# Colombian Pacific Coast shape file
CP_coast <- readShapeLines("../GISData/CP_coast.shp")

# Load data to weight edges of the network
load('../RData/mles.RData') # load relatedness estimates
load('../RData/geo_dist_info.RData') # Load distances

# Creat a data.frame to covert into a graph
mylinks <- data.frame(site1 = as.character(geo_dist_info$pairwise_site_distance[,1]),
                      site2 = as.character(geo_dist_info$pairwise_site_distance[,2]), 
                      stringsAsFactors = F)

# ******* This bit will change depending on plot *******
X = geo_dist_info$pairwise_site_distance[,'gravity_estimate']
mylinks$edge_weight <- (X - min(X))/(max(X) - min(X))

# Generate net, remove loops
net <- graph.data.frame(mylinks[mylinks$site1 != mylinks$site2,], directed = FALSE)

# Coordinates for graph layout
lonlat <- read.table('../TxtData/LonLat_dec_deg.txt', header = TRUE)
rownames(lonlat) <- lonlat[,'City']
layout_lonlat <- as.matrix(lonlat[rownames(as.matrix(V(net))), c('Longitude','Latitude')])
colnames(layout_lonlat) = c('lon','lat')

# IBD barcode
E(net)$width <- E(net)$edge_weight * 10
E(net)$curved <- 0.69*c(1,-1.3,1.2,1.4,-1,-0.9,-1,0,0,0)
E(net)$loop.angle <- pi*c(1.68, 2, rep(1,8))
E(net)$color <- adjustcolor('blue', alpha.f = 0.75)

# Change labels, text colour etc.
print(V(net)) # Print order
V(net)$label <- c('Quibdó','Guapí','Buenaventura','Tadó','Tumaco')
V(net)$label.color <- "black"
V(net)$label.dist <- c(3, 3, 7, 3, 3)
V(net)$label.degree <- c(1.5*pi, pi, 2*pi, 2*pi, 0.5*pi)
V(net)$label.cex <- 1.2
V(net)$color <- "black"
V(net)$size <- 8 
V(net)$frame.color <- "white"

# Plot Net
if(PDF){pdf(file = '../Plots/Colombia_networ.pdf', height = 7, width = 6)}
if(!NET){E(net)$color <- adjustcolor('blue', alpha.f = 0)} # Plot coast only
par(mfrow = c(1,1), family = "serif", mar = c(0,0,0,0))
plot(CP_coast, col = 'black', lwd = 2, xlim =  c(-78.5, -76.7))
plot(net, layout=layout_lonlat, add = TRUE, rescale = FALSE)
compassRose(x = -78.5, y = 5.74)
text(x = -78.52596, y = 3.72, labels = '50 km', cex = 1.2, pos = 3)
segments(y0 = 3.7, y1 = 3.7, x0 = -78.75928, x1 = -78.29265, lwd = 2, col = 'black')
if(PDF){dev.off()}


#================================================
# Plot of surrounding countries
#================================================
if(ZM){
  
  if(PDF){png(file = '../Plots/Colombia_Map.png', height=7, width=6, units='in', res=300)}

  # Load or download country shape files from http://gadm.org/  
  PATH = '../GISData/'
  Colombia <- getData("GADM",download = DWN,country="Colombia",level=0,path = PATH)
  Panama <- getData("GADM",download = DWN,country="Panama",level=0,path = PATH)
  Venezuela <- getData("GADM",download = DWN,country="Venezuela",level=0,path = PATH)
  Brazil <- getData("GADM",download = DWN,country="Brazil",level=0,path = PATH)
  Peru <- getData("GADM",download = DWN,country="Peru",level=0,path = PATH)
  Ecuador <- getData("GADM",download = DWN,country="Ecuador",level=0,path = PATH)
  
  par(mfrow = c(1,1), family = "serif", mar = c(0,0,0,0))
  plot(Colombia, col = 'gray', border = 'white', ylim = c(-1.5, 9.25), xlim = c(-82, -67))
  plot(Brazil, add = TRUE, col = 'gray', border = 'white')
  plot(Panama, add = TRUE, col = 'gray', border = 'white')
  plot(Peru, add = TRUE, col = 'gray', border = 'white')
  plot(Ecuador, add = TRUE, col = 'gray', border = 'white')
  plot(Venezuela, add = TRUE, col = 'gray', border = 'white')
  plot(CP_coast, col = 'black', lwd = 2, add = TRUE)
  text(x = c(-81.05882, -79.57647, -80.81176, -78.26588, -74.51765, -70.62118, -73.32000, -67.81059),
       y = c(8.401712, 10.388241,  4.041895, -1.252799, -2.392804,  8.416126, 3.531259, -2.230062),
       labels = c('Panama','Caribbean Sea', 'Pacific Ocean', 'Ecuador', 'Peru', 'Venezuela', 'Colombia', 'Brazil'),
       col = c('white', 'darkgray', 'darkgray', 'white', 'white', 'white', 'white', 'white'), cex = c(0.65,rep(1,7)))
  if(PDF){dev.off()}
}



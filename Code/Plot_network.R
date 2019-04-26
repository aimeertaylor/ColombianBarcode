############################################################################################
# This script is adapted from Maps_MORU.R
# 1) Map based on fraction of IBD proportions > 0.5 (2 plots, combined in manuscript)
# This map was created using QGIS (QGIS Development Team 2017) and the raster (Hijmans 2016),
# maptools (Bivand and Lewin-Koh 2017) and igraph (Csardi and Nepusz 2006) packages
# Update for state possibility 
############################################################################################
rm(list = ls())

library(maptools)
library(raster)
library(igraph)
library(sf) # For st_read instead of readShapeLines (returns warning: "use rgdal::readOGR or sf::st_read" )
PDF = T # Set to true to generate pdfs
DWN = F # Set to true to download country files
ZM = F # Set to true to zoomed out plot surrounding countries 
NET = T # Set to true to include network on coast
Threshold = '0.25' # on edges
Filter = 'Unfiltered' # on graph

# Colombian Pacific Coast shape file
CP_coast <- readShapeLines("../GISData/CP_coast.shp")

# Load data to weight edges of the network
load('../RData/proportions_sensitivities.RData') # load relatedness estimates
load('../RData/LonLat.RData') # Coordinates for graph layout

if(PDF){png(file = '../Plots/Colombia_network%d.png', width = 1500, height = 2100, res = 300)}
for(Gravity in c(T,F)){ # Set to true to plot gravity vs genetic
  for(Cities in c(T,F)){ # Set to true to plot cities rather than states
    if(Cities){
      load('../RData/geo_dist_info_cities.RData')
      attach(geo_dist_info)
      layout_lonlat <- LonLat_city[,2:1]
      intra = sapply(rownames(layout_lonlat), function(x)paste(x,x,sep= "_"))
      proportions = proportions_cities['mean',geo_order,Threshold,Filter]
      Labels = array(c('Quibdó','Tadó','Buenaventura','Guapí','Tumaco'), 
                     dimnames = list(rownames(layout_lonlat)))
    } else {
      load('../RData/geo_dist_info_states.RData')  
      attach(geo_dist_info)
      layout_lonlat <- LonLat_state[,2:1]
      intra = sapply(rownames(layout_lonlat), function(x)paste(x,x,sep= "_"))
      proportions = proportions_states['mean',geo_order,Threshold,Filter]
      Labels = array(c('Nariño','Cauca','Valle del Cauca','Chocó'), 
                     dimnames = list(rownames(layout_lonlat)))
    } 
    
    if(Gravity){
      X = geo_dist_info$pairwise_site_distance[,'gravity_estimate']
    } else {
      X = proportions
    }
    
    # Creat a data.frame to covert into a graph
    mylinks <- data.frame(site1 = as.character(geo_dist_info$pairwise_site_distance[,1]),
                          site2 = as.character(geo_dist_info$pairwise_site_distance[,2]), 
                          stringsAsFactors = F)
    mylinks$edge_weight <- (X - min(X))/(max(X) - min(X))
    
    # Generate net, remove loops
    net <- graph.data.frame(mylinks[mylinks$site1 != mylinks$site2,], directed = FALSE)
    
    # IBD barcode
    E(net)$width <- E(net)$edge_weight * 10
    E(net)$curved <- 0.69*c(1,-1,-0.3,1,0.3,-1,-1,1,1,1) # Designed for cities
    E(net)$color <- adjustcolor('blue', alpha.f = 0.75)
    
    # Change labels, text colour etc.
    print(V(net)) # Print order
    V(net)$label <- Labels[attributes(V(net))$names]
    V(net)$label.color <- "black"
    V(net)$label.dist <- c(5,5,7,3,7)
    V(net)$label.degree <- rep(2*pi, length(Labels)) 
    V(net)$label.cex <- 1.2
    V(net)$color <- "black"
    V(net)$size <- 8 
    V(net)$frame.color <- "white"
    
    # Plot Net
    if(!NET){E(net)$color <- adjustcolor('blue', alpha.f = 0)} # Plot coast only
    par(mfrow = c(1,1), family = "serif", mar = c(0,0,1,0))
    plot(CP_coast, col = 'black', lwd = 2, xlim =  c(-78.5, -76.7), main = ifelse(Gravity, 'Gravity', 'Genetic'))
    plot(net, layout=layout_lonlat[attributes(V(net))$names,], add = TRUE, rescale = FALSE)
    compassRose(x = -78.5, y = 5.74)
    text(x = -78.52596, y = 3.72, labels = '50 km', cex = 1.2, pos = 3)
    segments(y0 = 3.7, y1 = 3.7, x0 = -78.75928, x1 = -78.29265, lwd = 2, col = 'black')
  }  
} 
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



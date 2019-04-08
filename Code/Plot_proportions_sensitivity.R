#############################################
# This script plots proportions of 
# highly related parasite sample comparisons
# using various thresholds for high relatedness
# with and without filtering 
#############################################

rm(list = ls())
load('../RData/proportions_sensitivities')
Thresholds = dimnames(proportions_cities)[[3]]
Filter_status = dimnames(proportions_cities)[[4]]
PDF = TRUE

if(PDF){pdf('../Plots/Proportions_sensitivity.pdf', height = 10, width = 10)}

# Cities (intra ordered by expectation based on transmission and isolation)
load('../RData/geo_dist_info_cities.RData')
inter = geo_dist_info$geo_order  
intra = c('Guapi_Guapi', 'Tado_Tado', 'Tumaco_Tumaco', 'Buenaventura_Buenaventura', 
          'Quibdo_Quibdo')
intra_inter = c(intra, inter)

par(mfrow = c(3,3), family = 'serif')
for(fs in Filter_status){
  for(threshold in Thresholds){
    
    X = proportions_cities[,intra_inter,threshold, fs]
    
    # Barplot
    Midpoints = barplot(X['mean',], names.arg = '', 
                        las = 2, ylim = c(0,max(X)), density = rep(c(100,25), c(5, 10)), 
                        main = sprintf('%s', fs), 
                        ylab = bquote('Proportion with'~hat(italic(r))>.(threshold)))
    title(xlab = 'Cities')
    
    # x labels rotate 60 degrees, srt=60
    text(x = Midpoints, y = -max(X)/20, srt = 40, adj= 1, xpd = TRUE, 
         labels = gsub('_', ' & ',colnames(X)), cex = 0.4)
    
    # CIs
    segments(x0 = Midpoints[,1], x1 = Midpoints[,1], 
             y0 = X['2.5%',], y1 = X['97.5%',], lty = 1)
  }}


# States (intra ordered by expectation based on transmission and isolation)
load('../RData/geo_dist_info_states.RData')
inter = geo_dist_info$geo_order  
intra = c('Cauca_Cauca', 'Narino_Narino', 'Valle_Valle', 'Choco_Choco')
intra_inter = c(intra, inter)

par(mfrow = c(3,3), family = 'serif')
for(fs in Filter_status){
  for(threshold in Thresholds){
    
    X = proportions_states[,intra_inter,threshold, fs]
    
    # Barplot
    Midpoints = barplot(X['mean',], names.arg = '', 
                        las = 2, ylim = c(0,max(X)), density = rep(c(100,25), c(4, 6)), 
                        main = sprintf('%s', fs), 
                        ylab = bquote('Proportion with'~hat(italic(r))>.(threshold)))
    title(xlab = 'States')
    
    # x labels rotate 60 degrees, srt=60
    text(x = Midpoints, y = -max(X)/20, srt = 40, adj= 1, xpd = TRUE, 
         labels = gsub('_', ' & ',colnames(X)), cex = 0.5)
    
    # CIs
    segments(x0 = Midpoints[,1], x1 = Midpoints[,1], 
             y0 = X['2.5%',], y1 = X['97.5%',], lty = 1)
  }}


# Times
par(mfrow = c(3,3), family = 'serif')
for(fs in Filter_status){
  for(threshold in Thresholds){
    
    # Extract result
    X = proportions_times[,,threshold, fs]
    times_sorted = as.character(sort(as.numeric(colnames(X))))
    
    # Barplot
    Midpoints = barplot(X['mean',times_sorted], names.arg = times_sorted, cex.names = 0.5,
                        las = 2, ylim = c(0,max(X)), 
                        main = sprintf('%s', fs), 
                        ylab = bquote('Proportion with'~hat(italic(r))>.(threshold)))
    title(xlab = expression(Delta~'Time (weeks)'))
    
    # CIs
    segments(x0 = Midpoints[,1], x1 = Midpoints[,1], 
             y0 = X['2.5%', times_sorted], y1 = X['97.5%', times_sorted], lty = 1)
  }}

if(PDF){dev.off()}

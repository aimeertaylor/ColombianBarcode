###############################################################################
#' Script to explore evidence of connectivity accross sampling dates and locations
###############################################################################
load("../../RData/fraction_highly_related_extended.RData")
library(RColorBrewer)
PNG <- FALSE
CEX.LAB = 1.15
CEX.AXIS = 1
CEX.LEGEND = 1.25
remove_small_samples <- TRUE

if(PNG) png(file = '../../Plots/connectivity%d.png', width = 1500, height = 1500, res = 300)
if(remove_small_samples){
  
  fraction_highly_related_year <- fraction_highly_related_year[fraction_highly_related_year$npairs >= 10, ]
  Error_bars_year <- Error_bars_year[as.character(fraction_highly_related_year$year_diff), ]
  
  fraction_highly_related_city <- fraction_highly_related_city[fraction_highly_related_city$min_n_pair >= 10, ]
  Error_bars_city <- Error_bars_city[fraction_highly_related_city$City12, ]
}


#============== By year =================
X <- barplot(fraction_highly_related_year$fHR, 
             col = brewer.pal(5, "GnBu")[3], 
             las = 2, xlab = 'Collection date yearly difference', 
             xaxt = 'n', 
             cex.lab = CEX.LAB,  
             cex.axis = CEX.AXIS,  
             ylab = expression('Fraction of highly-related'~italic('P. falciparum')~'sample pairs'), 
             cex.names = 1, 
             ylim = c(0,max(Error_bars_year)))

segments(y0 = Error_bars_year[,'2.5%'], y1 = Error_bars_year[,'97.5%'], x0 = X, x1 = X)
axis(side = 1, at = X[,1], labels = fraction_highly_related_year$year_diff)


#============== By City inter/intra =================
# For addition of distance line 
normalised_geo_dist = (fraction_highly_related_city$intercity_dist - min(fraction_highly_related_city$intercity_dist))/diff(range(fraction_highly_related_city$intercity_dist))

# Bar plot 
X <- barplot(fraction_highly_related_city$fHR, 
             las = 2, col = brewer.pal(5, "GnBu")[3], 
             density = c(35,NA)[(fraction_highly_related_city$intercity_dist == 0)+1], 
             xlab = '', xaxt = 'n', 
             cex.lab = CEX.LAB,  
             cex.axis = CEX.AXIS,  
             ylab = expression('Fraction of highly-related'~italic('P. falciparum')~'sample pairs'), 
             ylim = c(-0.1,max(Error_bars_city)+0.1))
segments(y0 = Error_bars_city[,'2.5%'], y1 = Error_bars_city[,'97.5%'],
         x0 = X, x1 = X)

# x labels rotate 60 degrees, srt=60
text(x = X, y = -max(Error_bars_city)/75, 
     srt = 35, adj= 1, xpd = TRUE, 
     cex = 0.25, 
     labels = gsub("_", " ",fraction_highly_related_city$City12))

# Add distance 
# note that par(new = T) resulted in expansion of plotting space for which I couldn't find any
# documentation. I'm therefore plotting distance directly onto the barplot
lines(x = X, y = normalised_geo_dist, 
      type = 'b', pch = 20, #panel.first = grid(nx = NA), 
      yaxt = 'n', xaxt = 'n', ylab = ' ', xlab = ' ', bty = 'n')
text(x = max(X), pos = 3, cex = CEX.AXIS, 
     y = max(normalised_geo_dist), 
     labels = bquote(.(round(max(fraction_highly_related_city$intercity_dist), 0))))

Esmeraldas_Quibdo_ind <- which(fraction_highly_related_city$City12 == "Esmeraldas_Quibdo")
text(x = X[Esmeraldas_Quibdo_ind,], pos = 4, cex = CEX.AXIS, 
     y = normalised_geo_dist[Esmeraldas_Quibdo_ind], 
     labels = bquote(.(round(fraction_highly_related_city$intercity_dist[Esmeraldas_Quibdo_ind], 0))))

# Legend
legend('top',density = c(100,35), fill = brewer.pal(5, "GnBu")[3], 
       bty = 'n', cex = CEX.LEGEND, inset = 0.10, 
       legend = c('within cities','across cities'))
legend('top', lty = 1, pch = 20, legend = 'Inter-city great-circle distance (km)', 
       bty = 'n', cex = 1)



#============== By City coast access =================
# Bar plot 
X <- barplot(fraction_highly_related_city$fHR, 
             las = 2, col = c(brewer.pal(5, "GnBu")[3],"#000000FF")[(!fraction_highly_related_city$InterCoast & fraction_highly_related_city$intercity_dist > 0)+1], 
             density = c(NA,35)[(fraction_highly_related_city$InterCoast & fraction_highly_related_city$intercity_dist > 0)+1], 
             xlab = '', xaxt = 'n', 
             cex.lab = CEX.LAB,  
             cex.axis = CEX.AXIS,  
             ylab = expression('Fraction of highly-related'~italic('P. falciparum')~'sample pairs'), 
             ylim = c(-0.1,max(Error_bars_city)+0.1))
segments(y0 = Error_bars_city[,'2.5%'], y1 = Error_bars_city[,'97.5%'],
         x0 = X, x1 = X)

# x labels rotate 60 degrees, srt=60
text(x = X, y = -max(Error_bars_city)/75, 
     srt = 35, adj= 1, xpd = TRUE, 
     cex = 0.25, 
     labels = gsub("_", " ",fraction_highly_related_city$City12))

# Add distance 
lines(x = X, y = normalised_geo_dist, 
      type = 'b', pch = 20, 
      yaxt = 'n', xaxt = 'n', ylab = ' ', xlab = ' ', bty = 'n')

if(!remove_small_samples){
text(x = max(X), pos = 3, cex = CEX.AXIS,
     y = max(normalised_geo_dist),
     labels = bquote(.(round(max(fraction_highly_related_city$intercity_dist), 0))))
}

Esmeraldas_Quibdo_ind <- which(fraction_highly_related_city$City12 == "Esmeraldas_Quibdo")
text(x = X[Esmeraldas_Quibdo_ind,], pos = 4, cex = CEX.AXIS, 
     y = normalised_geo_dist[Esmeraldas_Quibdo_ind], 
     labels = bquote(.(round(fraction_highly_related_city$intercity_dist[Esmeraldas_Quibdo_ind], 0))))

# Legend
legend('top',density = c(35,NA), fill = c(brewer.pal(5, "GnBu")[3], "black"), 
       bty = 'n', cex = CEX.LEGEND, inset = 0.10, 
       legend = c('InterCoastal','Not coastal'))
legend('top', lty = 1, pch = 20, legend = 'Inter-city great-circle distance (km)', 
       bty = 'n', cex = 1)

if(PNG) dev.off()
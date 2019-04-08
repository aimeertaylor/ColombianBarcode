#################################
# To-do: errors and aethestics
#
# Quidbo + Tado = 0.25
# Tumaco = 19
# Buenaventura = 0.14
# Guapi = 0.14
#################################
rm(list = ls())
library(plotrix) # For gap.barplot
library(RColorBrewer)
cols = brewer.pal(8, 'Dark2') # Color Scheme
r_threshold = 0.5 # This is the threshold used in Generate_proportions.R
PDF = T  

if(PDF){pdf('../Plots/Proportions_and_W_distance.pdf')}
par(family = 'serif')
load('../RData/proportions_time.RData')
load('../RData/proportions_geo.RData')
load('../RData/CIs_clonal_intra.RData')
load('../RData/mle_CIs.RData')

# Load geography
load('../RData/geo_dist_info_cities.RData')
attach(geo_dist_info)
inter = geo_order  
intra = c('Guapi_Guapi', 'Tado_Tado', 'Tumaco_Tumaco', 'Buenaventura_Buenaventura', 'Quibdo_Quibdo')
site_comps = c(intra, inter)
no_site_comps = length(site_comps)
no_time_bins = ncol(proportions_time)

#================================================================
# Function to automatically set gap in a gap plot given counts
#================================================================
gapfunction_counts <- function(x){
  Z <- sort(x, decreasing = TRUE)
  Eps = 0.2*sort(abs(diff(Z)), decreasing = T)[1]
  is <- c(0, sort(which(abs(diff(Z)) > 100))) + 1
  gap_temp <- sort(c(Z[is[2]]+Eps, Z[is[1]]-Eps))
  while(any((Z >= min(gap_temp)) & (Z <= max(gap_temp)))){
    gap_temp <- gap_temp-10
  }
  return(gap_temp)
}

#================================================================
# Histogram of mles 
#================================================================
par(family = 'serif', mfrow = c(1,1))
X = hist(mle_CIs$rhat, breaks = seq(0,1,0.025), plot = F)
Gap = gapfunction_counts(x = X$counts)
Yticks = round(c(max(X$counts), seq(0, floor(Gap[1]), length.out = 4)),0)
Y = gap.barplot(y = X$counts, 
                panel.first = grid(), 
                gap = Gap, xaxt = 'n',
                ytics = Yticks,
                main = '', 
                ylab = expression('Number of'~italic('P. falciparum')~'sample pairs'), 
                xlab = '',
                las = 2, col = rep(cols[1], length(X$mids)))
title(xlab = expression('Estimate of genetic relatedness,'~italic(widehat(r))), 
      line = 1.5, cex.lab = 1)
axis(side = 1, line = -1, at = seq(0,length(Y),length.out = 5), 
     labels = seq(0,1,length.out = 5), cex.axis = 1, tick = FALSE)


#===========================================================
# Regular plot of CIs 
#===========================================================
Ordered_r = sort.int(mle_CIs$rhat, index.return = T) # Order estimates

# NULL plot
plot(NULL, ylim = c(0,1), xlim = c(1,length(mle_CIs$rhat)), 
     ylab = expression('Relatedness estimate,'~hat(italic(r))), 
     xlab = expression("Comparison index ranked by"~hat(italic(r))), 
     bty = 'n', las = 1, panel.first = grid())

# CIs by segment
segments(x0 = 1:length(mle_CIs$rhat), x1 = 1:length(mle_CIs$rhat),
         y0 = mle_CIs$`2.5%`[Ordered_r$ix], y1 = mle_CIs$`97.5%`[Ordered_r$ix],
         col = adjustcolor(cols[1], alpha.f = 0.5),  
         lwd = 0.1)

lines(Ordered_r$x) # Add mles

#===========================================================
# Travel_time instead of distance
# Guapi Buenaventura: c(4, 12)/2 = 8
# Guapi Tumaco: (6, 3.5)/2 = 4.75
# Buenaventura Tumaco: (17.5, 7, 15)/3 = 13
# 1 = road; 2 = sea; 3 = both
#===========================================================
pairwise_site_distance$travel_time = c(1.5, 4.75, 8.5, 8, 10, 13, 8+8.5, 8+10, 19.5, 21)
pairwise_site_distance$travel_type = c(1,2,1,2,1,3,3,3,1,1)

plot(x = pairwise_site_distance$distance, 
     y = pairwise_site_distance$travel_time, 
     ylab = 'Travel time estimate (hr)', xlab = 'Distance (km)', 
     panel.first = grid(), col = cols[pairwise_site_distance$travel_type], 
     pch = 16, bty = 'n')
X = cor.test(pairwise_site_distance$distance, pairwise_site_distance$travel_time)
text(x = 400, y = 5, labels = sprintf('correlation = %s', round(X$estimate, 3)))
legend('topleft', legend = c('Road (Google maps)', 
                             '† Sea (averaged over cargo ship and boat)', 
                             'Combination of road and sea'), 
       bty = 'n', pch = c(16), col = cols[1:3])
text(y = 15, x = 30, labels = '† All info. is directly from the sites', pos = 4)


#===========================================================
# Histogram of site comparison partioned by clone
#===========================================================
par(mfrow = c(1,1), family = 'serif', mar = c(4,5,2,2))
X <- barplot(1-proportions_geo_cloned['#D3D3D3FF',intra,'all'], las = 2, xlab = '', 
             ylab = 'Proportion of clonal sample comparisons', 
             cex.names = 0.5, ylim = c(0, max(CIs_clonal_intra)),  
             col = cols[1], 
             xaxt = 'n')
text(x = X, y = -0.002, srt = 30, adj= 1, xpd = TRUE, 
     labels = do.call(rbind, strsplit(intra, split = '_'))[,1], cex=1)
segments(y0 = CIs_clonal_intra['2.5%',intra], y1 = CIs_clonal_intra['97.5%',intra],
         x0 = X, x1 = X, lwd = 1.5)


#===========================================================
# Histogram of time (not partioned by site comp)
#===========================================================
par(mfrow = c(2,1), family = 'serif', mar = c(4,5,1,1))
X <- barplot(proportions_time['mean',], las = 2, col = cols[1], 
             xaxt = 'n', xlab = '',
             ylab = bquote('Proportion with LCI of'~~italic(widehat(r))>.(r_threshold)), 
             cex.names = 0.75, ylim = c(0,max(proportions_time)))
segments(y0 = proportions_time['2.5%',], y1 = proportions_time['97.5%',], x0 = X, x1 = X)
axis(side = 1, at = X, labels = colnames(proportions_time), tick = F, line = -0.5, las = 2)
title(xlab = expression(Delta~'Time (weeks)'), line = 2)


#===========================================================
# Histogram of time (partioned by clone)
#===========================================================
X <- barplot(proportions_time_cloned[, ,'r_threshold'], 
             las = 2, xlab = expression(Delta~'Time (weeks)'), 
             xaxt = 'n',
             ylab = bquote('Proportion with LCI of'~italic(widehat(r))>.(r_threshold)), 
             cex.names = 1, ylim = c(0,max(proportions_time)), 
             col = rownames(proportions_time_cloned[, ,'all']))
axis(side = 1, at = X, labels = colnames(proportions_time), tick = F, line = -0.5, las = 2)
segments(y0 = proportions_time['2.5%',], y1 = proportions_time['97.5%',], x0 = X, x1 = X)


#===========================================================
# Histogram of time (partioned by site comp)
#===========================================================
par(mfrow = c(2,1), family = 'serif')
# Grouped version 
X <- barplot(proportions_time_grouped[, ,'all'], las = 2, xaxt = 'n', 
             xlab = expression(Delta~'Time (weeks)'), 
             ylab = 'Proportion of sample comparisons', 
             cex.names = 1, col = rainbow(no_site_comps))
axis(side = 1, at = X, labels = colnames(proportions_time), tick = F, line = -0.5)

X <- barplot(proportions_time_grouped[, ,'r_threshold'], 
             las = 2, xlab = expression(Delta~'Time (weeks)'), 
             xaxt = 'n',
             ylab = bquote('Proportion with LCI of'~italic(widehat(r))>.(r_threshold)), 
             cex.names = 1, ylim = c(0,max(proportions_time)), col = rainbow(no_site_comps))
axis(side = 1, at = X, labels = colnames(proportions_time), tick = F, line = -0.5)
segments(y0 = proportions_time['2.5%',], y1 = proportions_time['97.5%',], x0 = X, x1 = X)
legend('topright', fill = rainbow(no_site_comps), bty = 'n',legend = gsub('_', ' & ', site_comps), 
       cex = 0.5)


#===========================================================
# Histogram of site comparison
#===========================================================
# For addition of distance line in [0,0.25]
normalised_geo_dist = (pairwise_site_distance_all[site_comps]/max(pairwise_site_distance_all))/4

par(mfrow = c(2,1), family = 'serif', mar = c(8,5,1,1))

# Bar plot 
X <- barplot(proportions_geo['mean',site_comps], 
             las = 2, col = cols[1], 
             density = rep(c(100,25), c(length(intra), length(site_comps)-length(intra))),
             xlab = '', xaxt = 'n',
             ylab = bquote('Proportion of'~hat(italic(r))>.(r_threshold)), 
             cex.names = 0.5, ylim = c(0,max(proportions_geo)))
segments(y0 = proportions_geo['2.5%',site_comps], y1 = proportions_geo['97.5%',site_comps],
         x0 = X, x1 = X)

# x labels rotate 60 degrees, srt=60
text(x = X, y = -0.01, srt = 40, adj= 1, xpd = TRUE, labels = gsub('_', ' & ',site_comps), cex=0.5)

# Add distance 
# note that par(new = T) resulted in expansion of plotting space for which I couldn't find any
# documentation. I'm therefore plotting distance directly onto the barplot
lines(x = X, y = normalised_geo_dist, ylim = c(0,510),
      type = 'b', pch = 20, panel.first = grid(nx = NA),
      yaxt = 'n', xaxt = 'n', ylab = ' ', xlab = ' ', bty = 'n')
text(x = max(X), pos = 3, cex = 0.75, 
     y = max(normalised_geo_dist), 
     labels = bquote(.(round(max(pairwise_site_distance_all), 0))~'(km)'))

# Legend
legend('top',density = c(100,35), fill = cols[1], bty = 'n', 
       legend = c(expression('within sites ('~Delta~'distance = 0 km)'),
                  expression('across sites ('~Delta~'distance > 0 km)')))
legend('top', lty = 1, pch = 20, legend = expression(Delta~'distance'), 
       inset = 0.15, bty = 'n')


# Break down partitioned by clone
X <- barplot(proportions_geo_cloned[,site_comps,'r_threshold'], las = 2, xlab = '', 
             ylab = bquote('Proportion with LCI of'~italic(widehat(r))>.(r_threshold)), 
             cex.names = 0.5, ylim = c(0,max(proportions_geo)), 
             col = rownames(proportions_geo_cloned[,site_comps,'all']), 
             xaxt = 'n')
segments(y0 = proportions_geo['2.5%',site_comps], y1 = proportions_geo['97.5%',site_comps],
         x0 = X, x1 = X)
text(x = X, y = -0.01, srt = 30, adj= 1, xpd = TRUE, labels = gsub('_', ' & ',site_comps), cex=0.5)




#===========================================================
# Histogram of site comparison partioned by time
#===========================================================
par(mfrow = c(2,1), family = 'serif', mar = c(7,5,1,1))

# Grouped version
X <- barplot(proportions_geo_grouped[,site_comps,'all'], las = 2, xaxt = 'n',
             ylab = 'Proportion of sample comparisons',
             cex.names = 0.5, col = rainbow(no_site_comps), xlab = ' ')
text(x = X, y = -0.05, srt = 30, adj= 1, xpd = TRUE, labels = gsub('_', ' ',site_comps), cex=0.5)

# Break down it to site comparsion contributions 
X <- barplot(proportions_geo_grouped[,site_comps,'r_threshold'], las = 2, xlab = '', 
             ylab = bquote('Proportion with LCI of'~~italic(widehat(r))>.(r_threshold)), 
             cex.names = 0.5, ylim = c(0,max(proportions_geo)), col = rainbow(no_time_bins), 
             xaxt = 'n')
segments(y0 = proportions_geo['2.5%',site_comps], y1 = proportions_geo['97.5%',site_comps],
         x0 = X, x1 = X)
legend('topright', fill = rainbow(no_time_bins), bty = 'n', cex = 0.5, 
       legend = rownames(proportions_geo_grouped), title = expression(Delta~'Time (weeks)'))
# x labels 
text(x = X, y = -0.05, srt = 30, adj= 1, xpd = TRUE, labels = gsub('_', ' ',site_comps), cex=0.5)




#==================================================
# Plot transport distance results 
#==================================================
load('../RData/All_W_results.RData')
par(mfrow = c(2,2), family = 'serif', mar = c(6,4,3,1))
intra = site_comps[1:5]
for(j in 3:1){
  X = All_W_results[[j]]
  mpts = barplot(X['cost',], las = 2, xaxt = 'n', ylab = "1-Wasserstein distance", 
                 main = names(All_W_results)[j], ylim = c(0,1))
  text(x = mpts[,1], y = -0.02, srt = 40, adj= 1, xpd = TRUE,
       labels =  gsub('_', ' & ', colnames(X)), cex = 0.7)
  segments(x0 = mpts[,1], x1 = mpts[,1], y0 = X['2.5%',], y1 = X['97.5%',])
}

if(PDF){dev.off()}


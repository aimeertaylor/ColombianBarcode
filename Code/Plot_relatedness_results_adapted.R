
rm(list = ls())
library(plotrix) # For gap.barplot
library(RColorBrewer)
source('./simtests.R')
source('~/Dropbox/IBD_IBS/plotting_functions.R')
par(family = 'serif')

# Load geo_dist info and raw data 
load('../RData/geo_dist_info.RData')
attach(geo_dist_info)

# Load and summarise raw data 
load('../RData/SNPData.RData') 
SNPDataBinary <- SNPData[,6:255] # Extract the SNPData w/o meta data
numSNPs <- ncol(SNPDataBinary)
numSamples <- nrow(SNPDataBinary)
Frequencies <- colMeans(SNPDataBinary, na.rm = TRUE)

# Load IBD and IBS results
load('../RData/mle_CIs.RData') 

# Color Scheme
cols = brewer.pal(8, 'Dark2')

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
     xlab = expression("Sample comparison index ranked by"~hat(italic(r))), 
     bty = 'n', las = 1, panel.first = grid())

# CIs by segment
segments(x0 = 1:length(mle_CIs$rhat), x1 = 1:length(mle_CIs$rhat),
         y0 = mle_CIs$`2.5%`[Ordered_r$ix], y1 = mle_CIs$`97.5%`[Ordered_r$ix],
         col = 'gray', 
         lwd = 0.1)

lines(Ordered_r$x) # Add mles

#===========================================================
# Travel_time instead of distance
# Guapi Buenaventura: c(4, 12)/2 = 8
# Guapi Tumaco: (6, 3.5)/2 = 4.75
# Buenaventura Tumaco: (17.5, 7, 15)/3 = 13
#===========================================================
# 1 = road
# 2 = sea
# 3 = both
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


#+++++++++++++++++++++
# THIS IS WHERE I GOT TO

```{r}
#=================================================================
# Calculate proportions by time and site comp
#=================================================================
nrep <- 1000 # For bootstrap confidence intervals
interval <- 50 # The resolution of time bins in weeks
max_time_dist = max(as.numeric(Result$time_dist))
time_breaks <- seq(0, max_time_dist, interval) 
max_time_breaks = max(time_breaks)
if(max_time_breaks < max_time_dist){time_breaks = c(time_breaks, max_time_breaks + interval)}
time_bins = time_breaks[-length(time_breaks)] # The time bins value are the min for each bin
no_time_bins = length(time_bins)
site_comps = unique(Result$site_comp)
no_site_comps = length(site_comps)
set.seed(1) # For reproducibility

# Add time bins 
Result$time_bins = NA
for(i in 1:no_time_bins){
  ind <- (as.numeric(Result$time_dist) >= time_breaks[i]) & (as.numeric(Result$time_dist) < time_breaks[i+1])
  Result$time_bins[ind] = time_bins[i]
}

# Stores (inc. those where intervals are broken down into site_comps and time bins)
proportions_time <- array(dim = c(no_time_bins, 2), dimnames = list(time_bins, c('IBD', 'IBS')))
proportions_dist <- array(dim = c(no_site_comps, 2), dimnames = list(site_comps, c('IBD', 'IBS')))
proportions_time_deltas <- array(dim = c(no_time_bins, 2, nrep), dimnames = list(time_bins, c('IBD', 'IBS'), NULL))
proportions_dist_deltas <- array(dim = c(no_site_comps, 2, nrep), dimnames = list(site_comps,c('IBD', 'IBS'), NULL))
proportions_time_grouped = array(0, dim = c(no_site_comps, no_time_bins, 3), 
                                 dimnames = list(site_comps, time_bins, c('lowIBD','highIBD', 'all')))
proportions_dist_grouped = array(0, dim = c(no_time_bins, no_site_comps, 3), 
                                 dimnames = list(time_bins, site_comps, c('lowIBD','highIBD', 'all')))

# Function to extract proportion related based on bootstrap
X <- array(dim = c(2, nrep)) # For array in apply(fun_boot)
fun_boot <- function(x, denom, actual_IBD, actual_IBS){
  X <- cbind(x = mean(actual_IBD[sample(denom, size = denom, replace = TRUE)]),
             y = mean(actual_IBS[sample(denom, size = denom, replace = TRUE)]))
  return(X)
}

# Calculation proportions with time
for(i in time_bins){
  
  index_i = as.character(i)
  
  # Extract data 
  ind <- Result$time_bins == i
  actual_IBD <- Result$IBD_tail[ind]
  actual_IBS <- Result$IBS_tail[ind]
  denom = sum(ind)
  
  # Proportions of different site comps within
  site_comp_breakdown_all = table(Result$site_comp[ind])
  site_comp_breakdown_highIBD = table(Result$site_comp[ind][as.logical(Result$IBD_tail[ind])])
  site_comp_breakdown_lowIBD = table(Result$site_comp[ind][!as.logical(Result$IBD_tail[ind])])
  proportions_time_grouped[names(site_comp_breakdown_all), index_i, 'all'] = site_comp_breakdown_all/denom
  proportions_time_grouped[names(site_comp_breakdown_highIBD), index_i, 'highIBD'] = site_comp_breakdown_highIBD/denom
  proportions_time_grouped[names(site_comp_breakdown_lowIBD), index_i, 'lowIBD'] = site_comp_breakdown_lowIBD/denom
  
  # Observed proportion overall wrt time
  proportions_time[index_i, 'IBD'] <- mean(actual_IBD)
  proportions_time[index_i, 'IBS'] <- mean(actual_IBS)
  
  # Boostrap proportions wrt time
  bootstrap <- apply(X, 2, FUN = fun_boot, denom, actual_IBD, actual_IBS)
  deltas <- bootstrap - proportions_time[index_i, ]
  proportions_time_deltas[index_i, ,] <- deltas
}

# Check proptions wrt time the same if grouped or not: Yes 
# cbind(rowSums(t(proportions_time_grouped[,,'highIBD']), na.rm = TRUE), proportions_time[,'IBD'])

# Calculate proportions with site
for(i in site_comps){
  
  # Extract data
  ind <- Result$site_comp == i 
  actual_IBD <- Result$IBD_tail[ind]
  actual_IBS <- Result$IBS_tail[ind]
  denom = sum(ind)
  
  # Proportions of different times within
  time_bin_breakdown_all = table(Result$time_bins[ind])
  time_bin_breakdown_highIBD = table(Result$time_bins[ind][as.logical(Result$IBD_tail[ind])])
  time_bin_breakdown_lowIBD = table(Result$time_bins[ind][!as.logical(Result$IBD_tail[ind])])
  proportions_dist_grouped[names(time_bin_breakdown_all), i, 'all'] = time_bin_breakdown_all/denom
  proportions_dist_grouped[names(time_bin_breakdown_highIBD), i, 'highIBD'] = time_bin_breakdown_highIBD/denom
  proportions_dist_grouped[names(time_bin_breakdown_lowIBD), i, 'lowIBD'] = time_bin_breakdown_lowIBD/denom
  
  # Observed proportion 
  proportions_dist[i, 'IBD'] <- mean(actual_IBD)
  proportions_dist[i, 'IBS'] <- mean(actual_IBS)
  
  # Boostrap proportions
  bootstrap <- apply(X, 2, FUN = fun_boot, denom, actual_IBD, actual_IBS)
  deltas <- bootstrap - proportions_dist[i, ]
  proportions_dist_deltas[i, ,] <- deltas
}

# # Check proptions wrt time the same if grouped or not: Yes 
# cbind(rowSums(t(proportions_dist_grouped[,,'highIBD']), na.rm = TRUE), proportions_dist[,'IBD'])

# Extract 95% percentiles and save
deltaCIs_time <- apply(proportions_time_deltas, c(1,2), quantile, probs = c(0.025, 0.975), na.rm = TRUE)
deltaCIs_dist <- apply(proportions_dist_deltas, c(1,2), quantile, probs = c(0.025, 0.975), na.rm = TRUE)
# deltaCIs <- list(deltaCIs_time = deltaCIs_time, deltaCIs_dist = deltaCIs_dist)
# save(delataCIs, file = "/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/delataCIs.RData")
```


## Plot of proportions $\widehat{r}_m > 0.5$ and $\widehat{\text{IBS}}_m > 0.84$ against time 

```{r, include = TRUE}
par(mfrow = c(1,1), family = 'serif', mar = c(4,5,1,1))
CIs <- rbind(proportions_time[, 'IBD'], proportions_time[, 'IBD']) + deltaCIs_time[,,'IBD']
X <- barplot(proportions_time[, 'IBD'], las = 2,col = cols[1], xaxt = 'n', 
             xlab = '',
             ylab = bquote('Proportion of'~~italic(widehat(r))>.(round(threshold_IBD,2))), 
             cex.names = 0.75,ylim = c(0,max(CIs)))
segments(y0 = CIs[1,], y1 = CIs[2,], x0 = X, x1 = X)
axis(side = 1, at = X, labels = rownames(proportions_time), tick = F, line = -0.5)
title(xlab = expression(Delta~'Time (weeks)'), line = 2)
```

```{r, include = TRUE}
par(mfrow = c(1,1), family = 'serif', mar = c(4,5,1,1))
CIs <- rbind(proportions_time[, 'IBS'], proportions_time[, 'IBS']) + deltaCIs_time[,,'IBS']
X <- barplot(proportions_time[, 'IBS'], las = 2, col = cols[2], xaxt = 'n', 
             xlab = '', 
             ylab = bquote('Proportion of'~~widehat(IBS)>.(round(threshold_IBS, 1))), 
             cex.names = 0.75, ylim = c(0,max(CIs)))
segments(y0 = CIs[1,], y1 = CIs[2,], x0 = X, x1 = X)
axis(side = 1, at = X, labels = rownames(proportions_time), tick = F, line = -0.5)
title(xlab = expression(Delta~'Time (weeks)'), line = 2)
```

## Plot of proportions $\widehat{r}_m > 0.5$ in time broken down by site comparison  

```{r, include = TRUE, fig.height = 10}
par(mfrow = c(2,1), family = 'serif')
# Break down it to site comparsion contributions 
CIs <- rbind(proportions_time[, 'IBD'], proportions_time[, 'IBD']) + deltaCIs_time[,,'IBD']
X <- barplot(proportions_time_grouped[, ,'highIBD'], 
             las = 2, 
             xlab = expression(Delta~'Time (weeks)'), 
             xaxt = 'n',
             ylab = bquote('Proportion of'~~italic(widehat(r))>.(round(threshold_IBD, 1))), 
             cex.names = 1, ylim = c(0,max(CIs)), col = rainbow(no_site_comps))
axis(side = 1, at = X, labels = colnames(proportions_time_grouped))
segments(y0 = CIs[1,], y1 = CIs[2,], x0 = X, x1 = X)
legend('topright', fill = rainbow(no_site_comps), bty = 'n',
       legend = sapply(site_comps, function(x)gsub('_', ' & ',x)))

# Grouped version 
X <- barplot(proportions_time_grouped[, ,'all'], las = 2, xaxt = 'n', 
             xlab = expression(Delta~'Time (weeks)'), 
             ylab = expression('Proportion of'~~italic(widehat(r))), 
             cex.names = 1, col = rainbow(no_site_comps))
axis(side = 1, at = X, labels = colnames(proportions_time_grouped))
```

## Plot of proportions $\widehat{r}_m > 0.5$ and $\widehat{\text{IBS}}_m > 0.84$ against city comparisons 

```{r, include = TRUE}
# Prep
# This bit is a bit hacky 
inter = apply(pairwise_site_distance[,1:2], 1, function(x){paste(sort(x), collapse = '_')})
intra <- c('Guapi_Guapi',
           'Tumaco_Tumaco',
           'Tado_Tado', 
           'Buenaventura_Buenaventura', 
           'Quibdo_Quibdo')
inter_sort = sort.int(pairwise_site_distance$distance, index.return = T)
intra_inter = c(intra, inter[inter_sort$ix])

# For addition of distance line in [0,0.25]
normalised_geo_dist = (pairwise_site_distance_all[intra_inter]/max(pairwise_site_distance_all))/4
```

```{r, include = TRUE}
par(mfrow = c(1,1), family = 'serif', mar = c(7,5,1,1))
# Bar plot IBS
CIs <- rbind(proportions_dist[intra_inter, 'IBS'], proportions_dist[intra_inter, 'IBS']) + deltaCIs_dist[, intra_inter,'IBS']
X <- barplot(proportions_dist[intra_inter, 'IBS'], 
             las = 2, col = cols[2], 
             density = rep(c(100,25), c(length(intra), length(inter))),
             xlab = '', xaxt = 'n',
             ylab = bquote('Proportion of'~~widehat(IBS)>.(round(threshold_IBS, 1))), 
             cex.names = 0.5, ylim = c(0,max(CIs)))
segments(y0 = CIs[1,], y1 = CIs[2,], x0 = X, x1 = X)

# Legend
legend('top',density = c(100,35), fill = cols[2], bty = 'n', 
       legend = c(expression('within site ('~Delta~'distance = 0 km)'),
                  expression('across sites ('~Delta~'distance > 0 km)')))

# x labels rotate 60 degrees, srt=60
text(x = X, y = -0.01, srt = 40, adj= 1, xpd = TRUE,
     labels = c(do.call(rbind,strsplit(intra, split = '_'))[,1], 
                sapply(inter, function(x)gsub('_', ' ',x))), cex=1)


# Add distance 
# note that par(new = T) resulted in expansion of plotting space for which I couldn't find any
# documenetation. I'm therefore plotting distance directly onto the barplot
lines(x = X, y = normalised_geo_dist, ylim = c(0,510),
     type = 'b', pch = 20, panel.first = grid(nx = NA),
     yaxt = 'n', xaxt = 'n', ylab = ' ', xlab = ' ', bty = 'n')
text(x = max(X), pos = 3, cex = 0.75, 
     y = max(normalised_geo_dist), 
     labels = bquote(.(round(max(pairwise_site_distance_all), 0))~'(km)'))
text(x = 14.4,  y = 0.18, pos = 3, cex = 0.75, 
     labels = expression(Delta~'distance'), srt = 29)

```

```{r, include = TRUE}
par(mfrow = c(1,1), family = 'serif',mar = c(7,5,1,1))

# Bar plot IBD
CIs <- rbind(proportions_dist[intra_inter, 'IBD'], proportions_dist[intra_inter, 'IBD']) + deltaCIs_dist[, intra_inter,'IBD']
X = barplot(proportions_dist[intra_inter, 'IBD'],  
            las = 2,  col = cols[1],
            xlab = '', xaxt = 'n',
            ylab = bquote('Proportion of'~~italic(widehat(r))>.(round(threshold_IBD,1))), 
            density = rep(c(100,25), c(length(intra), length(inter))),
            cex.names = 0.5, ylim = c(0,max(CIs)))
segments(y0 = CIs[1,], y1 = CIs[2,], x0 = X, x1 = X)

# x labels rotate 60 degrees, srt=60
text(x = X, y = -0.01, srt = 40, adj= 1, xpd = TRUE,
     labels = c(do.call(rbind,strsplit(intra, split = '_'))[,1], 
                sapply(inter, function(x)gsub('_', ' ',x))), cex=1)

# Legend
legend('top',density = c(100,35), fill = cols[1], bty = 'n', y.intersp = 1.1, 
       legend = c(expression('within site,'~Delta~'distance = 0 km'),
                  expression('across sites,'~Delta~'distance > 0 km')))

# Add distance 
# note that par(new = T) resulted in expansion of plotting space for which I couldn't find any
# documenetation. I'm therefore plotting distance directly onto the barplot
lines(x = X, y = normalised_geo_dist, ylim = c(0,510),
     type = 'b', pch = 20, panel.first = grid(nx = NA),
     yaxt = 'n', xaxt = 'n', ylab = ' ', xlab = ' ', bty = 'n')
text(x = max(X), pos = 3, cex = 0.75, 
     y = max(normalised_geo_dist), 
     labels = bquote(.(round(max(pairwise_site_distance_all), 0))~'(km)'))
text(x = 14.4,  y = 0.18, pos = 3, cex = 0.75, 
     labels = expression(Delta~'distance'), srt = 39)
```

## Plot of proportions $\widehat{r}_m > 0.5$ in time broken down by time between collection dates (weeks)



```{r, include = TRUE}
par(mfrow = c(1,1), family = 'serif', mar = c(7,5,1,1))
# Break down it to site comparsion contributions 
CIs <- rbind(proportions_dist[inter,'IBD'],
             proportions_dist[inter, 'IBD']) + deltaCIs_dist[,inter,'IBD']
X <- barplot(proportions_dist_grouped[,inter ,'highIBD'], las = 2, xlab = '', 
             ylab = bquote('Proportion of'~~italic(widehat(r))>.(round(threshold_IBD, 1))), 
             cex.names = 0.5, ylim = c(0,max(CIs)), col = rainbow(no_time_bins), 
             xaxt = 'n')
segments(y0 = CIs[1,], y1 = CIs[2,], x0 = X, x1 = X)
legend('topleft', fill = rainbow(no_time_bins), bty = 'n',
       legend = time_bins, title = expression(Delta~'Time (weeks)'))
# x labels 
text(x = X, y = -0.002, 
     srt = 30, adj= 1, xpd = TRUE,
     labels = sapply(inter, function(x)gsub('_', ' ',x)), cex=1)

# # Grouped version 
# X <- barplot(proportions_dist_grouped[,inter,'all'], las = 2,xaxt = 'n',
#              ylab = expression('Proportion of'~italic(widehat(r))), 
#              cex.names = 0.5, col = rainbow(no_site_comps), xlab = ' ')
# #rotate 60 degrees, srt=60
# text(x = X, y = -0.01, 
#      srt = 30, adj= 1, xpd = TRUE,
#      labels = c(do.call(rbind,strsplit(intra, split = '_'))[,1], 
#                 sapply(inter, function(x)gsub('_', ' ',x))), cex=0.8)



```

```{r, include = TRUE}
par(mfrow = c(1,1), family = 'serif', mar = c(7,5,1,1))

# Break down it to site comparsion contributions 
CIs <- rbind(proportions_dist[intra_inter,'IBD'], proportions_dist[intra_inter, 'IBD']) + deltaCIs_dist[,intra_inter,'IBD']
X <- barplot(proportions_dist_grouped[,intra_inter ,'highIBD'], las = 2, xlab = '', xaxt = 'n',
             ylab = bquote('Proportion of'~~italic(widehat(r))>.(round(threshold_IBD, 1))), 
             cex.names = 0.5, ylim = c(0,max(CIs)), col = rainbow(no_time_bins))
segments(y0 = CIs[1,], y1 = CIs[2,], x0 = X, x1 = X)


# Max_time_diff_in_highly_related
max_time_diff = max(as.numeric(unlist(apply(proportions_dist_grouped[,intra_inter ,'highIBD'],2,function(x){names(x[x>0])}))))

legend('top', fill = rainbow(sum(time_bins <= max_time_diff)), bty = 'n',
       legend = time_bins[time_bins<=max_time_diff], title = expression(Delta~'Time (weeks)'))

# # Add distance 
# # note that par(new = T) resulted in expansion of plotting space for which I couldn't find any
# # documenetation. I'm therefore plotting distance directly onto the barplot
# lines(x = X, y = normalised_geo_dist, ylim = c(0,510),
#      type = 'b', pch = 20, panel.first = grid(nx = NA),
#      yaxt = 'n', xaxt = 'n', ylab = ' ', xlab = ' ', bty = 'n')
# text(x = max(X), pos = 3, cex = 0.75, 
#      y = max(normalised_geo_dist[-(1:5)]), 
#      labels = bquote(.(round(max(pairwise_site_distance_all), 0))~'(km)'))
# text(x = 14.4,  y = 0.18, pos = 3, cex = 0.75, 
#      labels = expression(Delta~'distance'), srt = 29)

# Site labels 
text(x = X, y = -0.01, srt = 40, adj= 1, xpd = TRUE,
     labels = c(do.call(rbind,strsplit(intra, split = '_'))[,1], 
                sapply(inter, function(x)gsub('_', ' ',x))), cex=1)

# # Grouped version 
# X <- barplot(proportions_dist_grouped[,intra_inter,'all'], las = 2, xlab = '', 
#              ylab = expression('Proportion of'~italic(widehat(r))), 
#              cex.names = 0.5, col = rainbow(no_site_comps))
```

## Plot of proportions $\widehat{r}_m > 0.5$ against inter-city distance

```{r, include = TRUE}
par(mfrow = c(1,1), family = 'serif', mar = c(4,5,1,1))
plot(y = proportions_dist[, 'IBD'], 
     x = pairwise_site_distance_all[rownames(proportions_dist)], 
     bty = 'n', pch = 21, bg = 'gray', cex = 1, 
     xlab = 'Inter-city distance (km)', 
     ylab = expression('Proportion of'~~italic(widehat(r))>0.5))
text(y = proportions_dist[, 'IBD'], # pos 1,2,3,4 = b, l, 
     x = pairwise_site_distance_all[rownames(proportions_dist)], 
     labels = rownames(proportions_dist), 
     cex = 0.3, 
     pos = 3)
segments(y0 = CIs[1,rownames(proportions_dist)], y1 = CIs[2,rownames(proportions_dist)], 
         x0 = pairwise_site_distance_all[rownames(proportions_dist)], 
         x1 = pairwise_site_distance_all[rownames(proportions_dist)])
```

## Plot of city level FST estimates against proportions $\widehat{r}_m > 0.5$ 

```{r, include = TRUE}
# pi_IBD versus FST (not clear pattern)
par(mfrow = c(1,1), family = 'serif', mar = c(4,5,1,1))
inter_city_comp <- names(Fst_barcode$Pair_wise_site_comparisons_Fst)
plot(y = Fst_barcode$Pair_wise_site_comparisons_Fst, 
     x = proportions_dist[names(Fst_barcode$Pair_wise_site_comparisons_Fst),'IBD'], 
     bty = 'n', pch = 21, bg = 'gray', 
     panel.first = grid(), 
     xlab = expression('Proportion of'~~italic(widehat(r))>0.5), 
     ylab = expression(italic(widehat(F)[ST]^Hudson)))
segments(x0 = CIs[1,names(Fst_barcode$Pair_wise_site_comparisons_Fst)], 
         x1 = CIs[2,names(Fst_barcode$Pair_wise_site_comparisons_Fst)], 
         y0 = Fst_barcode$Pair_wise_site_comparisons_Fst, 
         y1 = Fst_barcode$Pair_wise_site_comparisons_Fst)
segments(x0 = proportions_dist[names(Fst_barcode$Pair_wise_site_comparisons_Fst),'IBD'], 
         x1 = proportions_dist[names(Fst_barcode$Pair_wise_site_comparisons_Fst),'IBD'],
         y0 = Fst_barcode$Pair_wise_site_comparisons_Fst_CIs[comparison_names, '2.5%'], 
         y1 = Fst_barcode$Pair_wise_site_comparisons_Fst_CIs[comparison_names, '97.5%'])

```

### Are all Buenaventura Tomaco disproportionally from 2005? Was it a clonal outbreak? 

```{r, include = TRUE, fig.width=10, fig.height=5}
ind_all <- Result$site_comp == "Buenaventura_Tumaco"
ind_high <- Result$site_comp == "Buenaventura_Tumaco" & Result$IBD_tail == 1
par(mfrow = c(1,2), family = 'serif')

for(ind in list(ind_all, ind_high)){
  
  # Need to reorganise the samples for all in one site and other in other site 
  X <- array(dim = c(sum(ind), 2, 2), 
             dimnames = list(NULL, c('date', 'site'), c('sample1', 'sample2')))
  X[,'date','sample1'] <- as.character(SNPData[as.character(Result$sample1[ind]), 'COLLECTION.DATE'])
  X[,'date','sample2'] <- as.character(SNPData[as.character(Result$sample2[ind]), 'COLLECTION.DATE'])
  X[,'site','sample1'] <- SNPData[as.character(Result$sample1[ind]), 'City']
  X[,'site','sample2'] <- SNPData[as.character(Result$sample2[ind]), 'City']
  
  # Indeces to re-organise by site
  ind_B1 <- which(X[,'site','sample1'] == 'Buenaventura')
  ind_B2 <- which(X[,'site','sample2'] == 'Buenaventura')
  
  # All those collected 
  plot(x = as.Date(c(X[ind_B1,'date','sample1'], X[ind_B2,'date','sample2'])), 
       y = as.Date(c(X[ind_B1,'date','sample2'], X[ind_B2,'date','sample1'])), 
       bty = 'n', 
       xlab = 'Buenaventura date of collection', ylab = 'Tumaco date of collection', 
       ylim = as.Date(c("1993-01-01", "2007-12-31")), 
       xlim = as.Date(c("1993-01-01", "2007-12-31")))
  
  # Add source sink some how and intra relatedness
  polygon(x = as.Date(c("1993-01-01", "2007-12-31", "2007-12-31")), 
          y = as.Date(c("1993-01-01", "2007-12-31", "1993-01-01")), col = 'blue', border = NA)
  polygon(x = as.Date(c("1993-01-01", "1993-01-01", "2007-12-31")), 
          y = as.Date(c("1993-01-01", "2007-12-31", "2007-12-31")), col = 'pink', border = NA)
  
  # Points
  points(x = as.Date(c(X[ind_B1,'date','sample1'], X[ind_B2,'date','sample2'])) , 
         y = as.Date(c(X[ind_B1,'date','sample2'], X[ind_B2,'date','sample1'])) , 
         pch = 20, col = adjustcolor('black', alpha = 0.5))
}
```




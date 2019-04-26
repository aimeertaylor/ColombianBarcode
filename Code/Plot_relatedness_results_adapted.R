#################################
# To-do
# Decide on threshold to show (0.25 or 0.5)
# Call this say plots of histrograms

# Quidbo + Tado = 0.25
# Tumaco = 19
# Buenaventura = 0.14
# Guapi = 0.14
#################################
rm(list = ls())
library(plotrix) # For gap.barplot
library(RColorBrewer)
par(family = 'serif')

# Load geo_dist info and raw data 
load('../RData/geo_dist_info.RData')
attach(geo_dist_info)
intra = c('Guapi_Guapi','Buenaventura_Buenaventura','Tumaco_Tumaco','Quibdo_Quibdo', 'Tado_Tado')
intra_inter = c(intra, geo_order[!grepl('Tado', geo_order)])

# Load and summarise raw data 
load('../RData/SNPData.RData') 
SNPDataBinary <- SNPData[,6:255] # Extract the SNPData w/o meta data
numSNPs <- ncol(SNPDataBinary)
numSamples <- nrow(SNPDataBinary)
Frequencies <- colMeans(SNPDataBinary, na.rm = TRUE)

# Load IBD and IBS results
load('../RData/mle_CIs.RData') 

mle_CIs[mle_CIs$City1 == "Tado", 'City1'] = 'Quibdo' 
mle_CIs[mle_CIs$City2 == "Tado", 'City2'] = 'Quibdo' 
mle_CIs$City12 = apply(mle_CIs[,c('City1','City2')], 1, function(x)paste(sort(x),collapse="_"))

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


#=================================================================
# Calculate proportions by time and site comp
# This is where I got to 
#=================================================================
nrep <- 100 # For bootstrap confidence intervals
interval <- 50 # The resolution of time bins in weeks
max_time_dist = max(as.numeric(mle_CIs$time_dist))
time_breaks <- seq(0, max_time_dist, interval) 
max_time_breaks = max(time_breaks)
if(max_time_breaks < max_time_dist){time_breaks = c(time_breaks, max_time_breaks + interval)}
time_bins = time_breaks[-length(time_breaks)] # The time bins value are the min for each bin
no_time_bins = length(time_bins)
site_comps = unique(mle_CIs$City12)
no_site_comps = length(site_comps)
set.seed(1) # For reproducibility

# Add time bins 
mle_CIs$time_bins = NA
for(i in 1:no_time_bins){
  ind <- (as.numeric(mle_CIs$time_dist) >= time_breaks[i]) & (as.numeric(mle_CIs$time_dist) < time_breaks[i+1])
  mle_CIs$time_bins[ind] = time_bins[i]
}

# Add membership by colour
load('../RData/edge_cols.RData')
mle_CIs$sample_comp = apply(mle_CIs[, c("individual1", "individual2")], 1, function(x)paste(sort(x), collapse = "_"))
mle_CIs$comp_col = "#D3D3D3FF" # gray
if(all(names(edge_cols) %in% mle_CIs$sample_comp)){
  ind = mle_CIs$sample_comp %in% names(edge_cols) # reduce df to relevant comp of 
  for(comp in names(edge_cols)){
    ind_comp = which(mle_CIs$sample_comp[ind] == comp)
    mle_CIs$comp_col[ind][ind_comp] = edge_cols[comp]
  }} else {stop('Not all edges present')}

comp_cols = unique(mle_CIs$comp_col)

# Stores (inc. those where intervals are partioned into site_comps and time bins)
proportions_time = array(0, dim = c(3, no_time_bins), dimnames = list(c('mean', '2.5%', '97.5%'), time_bins))
proportions_geo = array(0, dim = c(3, no_site_comps), dimnames = list(c('mean', '2.5%', '97.5%'), site_comps))
proportions_time_grouped = array(0, dim = c(no_site_comps, no_time_bins, 2), dimnames = list(site_comps, time_bins, c('all', 'r_threshold')))
proportions_geo_grouped = array(0, dim = c(no_time_bins, no_site_comps, 2), dimnames = list(time_bins, site_comps, c('all', 'r_threshold')))
proportions_time_cloned = array(0, dim = c(length(comp_cols), no_time_bins, 2), dimnames = list(comp_cols, time_bins, c('all', 'r_threshold')))
proportions_geo_cloned = array(0, dim = c(length(comp_cols), no_site_comps, 2), dimnames = list(comp_cols, site_comps, c('all', 'r_threshold')))

r_threshold = 0.01

# Calculation proportions with time
for(i in time_bins){
  
  # Extract indices, data, etc. 
  index_i = as.character(i)
  ind <- mle_CIs$time_bins == i
  n_time_bin <- sum(ind)
  rhats <- mle_CIs$`2.5%`[ind]
  
  # Proportions of different site comps within
  site_comp_breakdown_all = table(mle_CIs$City12[ind])
  site_comp_breakdown_r_threshold = table(mle_CIs$City12[ind][rhats > r_threshold])
  proportions_time_grouped[names(site_comp_breakdown_all), index_i, 'all'] = site_comp_breakdown_all/n_time_bin
  proportions_time_grouped[names(site_comp_breakdown_r_threshold), index_i, 'r_threshold'] = site_comp_breakdown_r_threshold/n_time_bin
  
  # Proportions of different clonal components within
  clone_comp_breakdown_all = table(mle_CIs$comp_col[ind])
  clone_comp_breakdown_r_threshold = table(mle_CIs$comp_col[ind][rhats > r_threshold])
  proportions_time_cloned[names(clone_comp_breakdown_all), index_i, 'all'] = clone_comp_breakdown_all/n_time_bin
  proportions_time_cloned[names(clone_comp_breakdown_r_threshold), index_i, 'r_threshold'] = clone_comp_breakdown_r_threshold/n_time_bin
  
  # Boostrap proportions wrt time
  prob_b <- sapply(1:nrep, function(b)mean(sample(rhats, size = n_time_bin, replace = TRUE) > r_threshold))
  
  # Observed proportion overall wrt time and quantiles
  proportions_time[,index_i] = c(mean(rhats > r_threshold), quantile(prob_b, probs = c(0.025, 0.975)))
}

# # Check proptions wrt time the same if grouped or not: Yes 
# cbind(colSums(proportions_time_grouped[,,'r_threshold'], na.rm = TRUE),proportions_time['mean',])

# Calculate proportions with site
for(i in site_comps){
  
  # Extract data
  ind <- mle_CIs$City12 == i 
  n_city12 = sum(ind)
  rhats <- mle_CIs$`2.5%`[ind]
  
  # Proportions of different times within
  time_bin_breakdown_all = table(mle_CIs$time_bins[ind])
  time_bin_breakdown_r_threshold = table(mle_CIs$time_bins[ind][rhats > r_threshold])
  proportions_geo_grouped[names(time_bin_breakdown_all), i, 'all'] = time_bin_breakdown_all/n_city12
  proportions_geo_grouped[names(time_bin_breakdown_r_threshold), i, 'r_threshold'] = time_bin_breakdown_r_threshold/n_city12
  
  # Proportions of different clonal components within
  clone_comp_breakdown_all = table(mle_CIs$comp_col[ind])
  clone_comp_breakdown_r_threshold = table(mle_CIs$comp_col[ind][rhats > r_threshold])
  proportions_geo_cloned[names(clone_comp_breakdown_all), i, 'all'] = clone_comp_breakdown_all/n_city12
  proportions_geo_cloned[names(clone_comp_breakdown_r_threshold), i, 'r_threshold'] = clone_comp_breakdown_r_threshold/n_city12
  
  # Boostrap proportions 
  prob_b <- sapply(1:nrep, function(b)mean(sample(rhats, size = n_city12, replace = TRUE) > r_threshold))
  
  # Observed proportion 
  proportions_geo[,i] <- c('mean' =  mean(rhats > r_threshold), quantile(prob_b, probs = c(0.025, 0.975)))
}

# # Check proptions wrt time the same if grouped or not: Yes 
# cbind(colSums(proportions_geo_grouped[,,'r_threshold'], na.rm = TRUE), proportions_geo['mean',])

# CIs for intra clonal proportion
CIs_clonal_intra = sapply(intra, function(i){
  # Extract data
  ind <- mle_CIs$City12 == i 
  n_city12 = sum(ind)
  clones <- mle_CIs$comp_col[ind]
  
  # Boostrap proportions not "#D3D3D3FF" (i.e. not singlton)
  prob_b <- sapply(1:nrep, function(b)mean(sample(clones, size = n_city12, replace = TRUE) != "#D3D3D3FF"))
  quantile(prob_b, probs = c(0.025, 0.975))
})

#----------------------------------------------------
# Histogram of site comparison partioned by clone
#----------------------------------------------------
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


#----------------------------------------------------
# Histogram of time (not partioned by site comp)
#----------------------------------------------------
par(mfrow = c(2,1), family = 'serif', mar = c(4,5,1,1))
X <- barplot(proportions_time['mean',], las = 2, col = cols[1], 
             xaxt = 'n', xlab = '',
             ylab = bquote('Proportion of'~~italic(widehat(r))>.(r_threshold)), 
             cex.names = 0.75, ylim = c(0,max(proportions_time)))
segments(y0 = proportions_time['2.5%',], y1 = proportions_time['97.5%',], x0 = X, x1 = X)
axis(side = 1, at = X, labels = colnames(proportions_time), tick = F, line = -0.5)
title(xlab = expression(Delta~'Time (weeks)'), line = 2)

#----------------------------------------------------
# Histogram of time (partioned by clone)
#----------------------------------------------------
X <- barplot(proportions_time_cloned[, ,'r_threshold'], 
             las = 2, xlab = expression(Delta~'Time (weeks)'), 
             xaxt = 'n',
             ylab = bquote('Proportion of sample comparisons with'~~italic(widehat(r))>.(r_threshold)), 
             cex.names = 1, ylim = c(0,max(proportions_time)), 
             col = rownames(proportions_time_cloned[, ,'all']))
axis(side = 1, at = X, labels = colnames(proportions_time), tick = F, line = -0.5)
segments(y0 = proportions_time['2.5%',], y1 = proportions_time['97.5%',], x0 = X, x1 = X)


#----------------------------------------------------
# Histogram of time (partioned by site comp)
#----------------------------------------------------
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
             ylab = bquote('Proportion of sample comparisons with'~~italic(widehat(r))>.(r_threshold)), 
             cex.names = 1, ylim = c(0,max(proportions_time)), col = rainbow(no_site_comps))
axis(side = 1, at = X, labels = colnames(proportions_time), tick = F, line = -0.5)
segments(y0 = proportions_time['2.5%',], y1 = proportions_time['97.5%',], x0 = X, x1 = X)
legend('topright', fill = rainbow(no_site_comps), bty = 'n',legend = gsub('_', ' & ', site_comps))








#----------------------------------------------------
# Histogram of site comparison
#----------------------------------------------------
# For addition of distance line in [0,0.25]
normalised_geo_dist = (pairwise_site_distance_all[intra_inter]/max(pairwise_site_distance_all))/4

par(mfrow = c(2,1), family = 'serif', mar = c(8,5,1,1))

# Bar plot 
X <- barplot(proportions_geo['mean',intra_inter], 
             las = 2, col = cols[1], 
             density = rep(c(100,25), c(length(intra), length(intra_inter)-length(intra))),
             xlab = '', xaxt = 'n',
             ylab = bquote('Proportion of'~hat(italic(r))>.(r_threshold)), 
             cex.names = 0.5, ylim = c(0,max(proportions_geo)))
segments(y0 = proportions_geo['2.5%',intra_inter], y1 = proportions_geo['97.5%',intra_inter],
         x0 = X, x1 = X)

# x labels rotate 60 degrees, srt=60
text(x = X, y = -0.01, srt = 40, adj= 1, xpd = TRUE, labels = gsub('_', ' & ',intra_inter), cex=1)

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
X <- barplot(proportions_geo_cloned[,intra_inter,'r_threshold'], las = 2, xlab = '', 
             ylab = bquote('Proportion of sample comparisons with'~~italic(widehat(r))>.(r_threshold)), 
             cex.names = 0.5, ylim = c(0,max(proportions_geo)), 
             col = rownames(proportions_geo_cloned[,intra_inter,'all']), 
             xaxt = 'n')
segments(y0 = proportions_geo['2.5%',intra_inter], y1 = proportions_geo['97.5%',intra_inter],
         x0 = X, x1 = X)
text(x = X, y = -0.01, srt = 30, adj= 1, xpd = TRUE, labels = gsub('_', ' & ',intra_inter), cex=1)



#----------------------------------------------------
# Histogram of site comparison partioned by time
#----------------------------------------------------
par(mfrow = c(2,1), family = 'serif', mar = c(7,5,1,1))

# Grouped version
X <- barplot(proportions_geo_grouped[,intra_inter,'all'], las = 2, xaxt = 'n',
             ylab = 'Proportion of sample comparisons',
             cex.names = 0.5, col = rainbow(no_site_comps), xlab = ' ')
text(x = X, y = -0.05, srt = 30, adj= 1, xpd = TRUE, labels = gsub('_', ' ',intra_inter), cex=1)

# Break down it to site comparsion contributions 
X <- barplot(proportions_geo_grouped[,intra_inter,'r_threshold'], las = 2, xlab = '', 
             ylab = bquote('Proportion of sample comparisons with'~~italic(widehat(r))>.(r_threshold)), 
             cex.names = 0.5, ylim = c(0,max(proportions_geo)), col = rainbow(no_time_bins), 
             xaxt = 'n')
segments(y0 = proportions_geo['2.5%',intra_inter], y1 = proportions_geo['97.5%',intra_inter],
         x0 = X, x1 = X)
legend('topright', fill = rainbow(no_time_bins), bty = 'n', cex = 0.75, 
       legend = time_bins, title = expression(Delta~'Time (weeks)'))
# x labels 
text(x = X, y = -0.05, srt = 30, adj= 1, xpd = TRUE, labels = gsub('_', ' ',intra_inter), cex=1)




#==================================================
# Plot transport distance results 
#==================================================
load('../RData/All_W_results.RData')
par(mfrow = c(2,2), family = 'serif', mar = c(6,4,3,1))
for(j in 3:1){
  X = All_W_results[[j]]
  mpts = barplot(X['cost',], las = 2, xaxt = 'n', ylab = "1-Wasserstein distance", 
                 main = names(All_W_results)[j], ylim = c(0,1), col = cols[1])
  text(x = mpts[,1], y = -0.02, srt = 40, adj= 1, xpd = TRUE,
       labels =  gsub('_', ' & ', colnames(X)), cex = 0.7)
  segments(x0 = mpts[,1], x1 = mpts[,1], y0 = X['2.5%',], y1 = X['97.5%',])
}

# Compare transport and genetic
plot(proportions_geo['mean',intra_inter[-(1:length(intra))]], All_W_results[[3]]['cost',intra_inter[-(1:length(intra))]]) 



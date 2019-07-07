###########################################
# This script generates plots generated for 
# the Gordon Research Conference 2019
###########################################
rm(list = ls())
library(plotrix) # For gap.barplot
library(RColorBrewer)
cols = brewer.pal(8, 'Dark2') # Codlor Scheme
r_threshold = 0.25 # This is the threshold used in Generate_proportions.R
PDF = T

if(PDF){png('../CIS_GRC_%d.png', res = 500, height = 6, width = 6, units = 'in')}
par(family = 'serif')
load('../RData/mle_CIs.RData')

#===========================================================
# Regular plot of CIs 
#===========================================================
Ordered_r = sort.int(mle_CIs$rhat, index.return = T) # Order estimates

# NULL plot
plot(NULL, ylim = c(0,1), xlim = c(1,length(mle_CIs$rhat)), 
     ylab = 'Relatedness estimate', cex.lab = 1.25, 
     xlab = "Sample comparison ranked by relatedness estimate", 
     las = 1, panel.first = grid())

segments(x0 = 1:length(mle_CIs$rhat), x1 = 1:length(mle_CIs$rhat),
         y0 = mle_CIs$`2.5%`[Ordered_r$ix], y1 = mle_CIs$`97.5%`[Ordered_r$ix],
         col = 'lightgray', lwd = 0.1)

lines(Ordered_r$x, lwd = 2) # Add mles

#===========================================================
# Plot of CIs by r_threshold 
#===========================================================
Ordered_r = sort.int(mle_CIs$rhat, index.return = T) # Order estimates

# NULL plot
plot(NULL, ylim = c(0,1), xlim = c(1,length(mle_CIs$rhat)), 
     ylab = 'Relatedness estimate', #expression('Relatedness estimate'~hat(italic(r))), 
     xlab = "Sample comparison ranked by relatedness estimate", #expression("Sample comparison ranked by"~hat(italic(r))), 
     las = 1, panel.first = grid(), cex.axis = 1.2, cex.lab = 1.25)

# Make transparency a function of lower CI
Thresholds = c(0.25)
cols_inds = apply(sapply(1:length(Thresholds), function(i){
  a = rep(1, nrow(mle_CIs))
  
  # When theshold is 0.01, 0.25 or 0.5
  ind = mle_CIs$`2.5%` > Thresholds[i] 

  if(Thresholds[i] == 0.99){ # Overwrite
    ind = mle_CIs$`97.5%` > Thresholds[i] & mle_CIs$`2.5%` > 0.01}
  
  a[ind] <- (i+1)
  return(a)}), 1, max)

# Calculate "colour"
cols_CIs = c(adjustcolor('lightgray'),cols[1])[cols_inds]

# CIs by segment
segments(x0 = 1:length(mle_CIs$rhat), x1 = 1:length(mle_CIs$rhat),
         y0 = mle_CIs$`2.5%`[Ordered_r$ix], y1 = mle_CIs$`97.5%`[Ordered_r$ix],
         col = cols_CIs[Ordered_r$ix])

lines(Ordered_r$x, lwd = 2) # Add mles

# Add legend
legend('topleft', fill = c(adjustcolor('lightgray'),cols[1]), 
       inset = 0.01, 
       legend = c('Statistically indistinguishable from 0', 
                  'Statistically distinguishable from 0.25'))



if(PDF){dev.off()}




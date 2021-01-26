#################################################################
# README!!!!!!!!!!!!!!!!
# This script adds metadata to the dataframe of relatedness estimates
# We also remove overly precise estimates close to zero and one
# Must do this part by hand!!!
#################################################################
rm(list = ls())
library(stringr)
eps <- 0.01 # Used to define clones in Taylor et al. 2020

# Function to summarise mles to track formatting
summarise_mles <- function(x, snp_count_threshold = 10, 
                           zoom = F, eps = 0.01){
  
  x <- dplyr::arrange(x, rhat, CI_width)
  
  # Plot relatedness values
  plot(NULL, main = "", 
       ylim = c(0,1), 
       xlim = c(1,length(x$rhat)), 
       ylab = 'Relatedness estimate', 
       xlab = "Sample pair index", 
       las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
  axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
  segments(x0 = 1:nrow(x), x1 = 1:nrow(x), # Add CIs
           y0 = x$r2.5., y1 = x$r97.5.,
           lwd = 0.1, col = 'gray')
  lines(x$rhat, lwd = 1, col = "black") # Add mles
  abline(h = 1-eps, lty = "dashed", col = 'blue', lwd = 0.5)
  
  if(zoom) {
    # Plot relatedness values
    plot(NULL, main = "", 
         ylim = c(0,1), 
         xlim = c(sum(x$rhat < 0.8, na.rm = T),length(x$rhat)), 
         ylab = 'Relatedness estimate', 
         xlab = "Sample pair index", 
         las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
    axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
    segments(x0 = 1:nrow(x), x1 = 1:nrow(x), # Add CIs
             y0 = x$r2.5., y1 = x$r97.5.,
             lwd = 0.1, col = 'gray')
    lines(x$rhat, lwd = 1, col = "black") # Add mles
    abline(h = 1-eps, lty = "dashed", col = 'blue', lwd = 0.5)
  }
  
  # Summarise data paucity
  snp_counts <- array(0, dim = c(2,snp_count_threshold), 
                      dimnames = list(c("per_sample", "per_sample_pair"), 1:snp_count_threshold))
  snp_counts[1, ] <- table(metadata[unique(c(x$individual1, x$individual2)), "snp_count"])[colnames(snp_counts)]
  snp_counts[2, ] <- table(x$snp_count)[colnames(snp_counts)]
  
  writeLines(paste(sprintf("No. sample pairs: %s", nrow(x)), 
                   sprintf("No. samples: %s", length(unique(c(x$individual1, x$individual2)))), 
                   "Numbers of samples and sample pairs with sparse data:", sep = "\n"))
  print(snp_counts)
}

# Load metadata
load(file = "../../RData/metadata_extended.RData")

# Load relatedness estimates
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s.RData', freqs_used)) # Load All_results
mle_CIs$CI_width <- (mle_CIs$r97.5. - mle_CIs$r2.5.) # Add CI width to relatedness estimates

summarise_mles(mle_CIs)
mle_CIs <- mle_CIs[!is.na(mle_CIs$rhat),] # Remove NAs
summarise_mles(mle_CIs)



# ========== Check missing loci count ==========
# Load the extended data set
load(file = "../../RData/snpdata_extended.RData")
missing_loci_count <- sapply(1:nrow(mle_CIs), function(i){
  ys <- snpdata[,c(mle_CIs$individual1[i], mle_CIs$individual2[i])]
  sum(rowSums(is.na(ys)) > 0)
})
all((250-mle_CIs$missing_loci_count) == missing_loci_count)

# ========== Add meta data to the results ==========

# Factors to characters
if(class(mle_CIs$individual1) == 'factor'){mle_CIs$individual1 = as.character(mle_CIs$individual1)}
if(class(mle_CIs$individual2) == 'factor'){mle_CIs$individual2 = as.character(mle_CIs$individual2)}

# Add cities
# Check NA returned for missing sids: 
is.na(metadata['agsdhfsdhf', 'City'])
mle_CIs$City1 = metadata[mle_CIs$individual1, 'City']
mle_CIs$City2 = metadata[mle_CIs$individual2, 'City']
mle_CIs$State1 = metadata[mle_CIs$individual1, 'STATE']
mle_CIs$State2 = metadata[mle_CIs$individual2, 'STATE']

# Add city comparisons
mle_CIs$City12 = apply(mle_CIs[,c('City1','City2')], 1, function(x)paste(sort(x),collapse="_"))

# Check NA returned for missing sids: yes
metadata['agsdhfsdhf', 'COLLECTION.DATE']
any(is.na(metadata$Year)) # Check for missing years

# All the samples with missing dates are recent
range(metadata$Year[is.na(metadata$COLLECTION.DATE)])

# Add dates and years
mle_CIs$date1 = metadata[mle_CIs$individual1, 'COLLECTION.DATE']
mle_CIs$date2 = metadata[mle_CIs$individual2, 'COLLECTION.DATE']
mle_CIs$year1 = metadata[mle_CIs$individual1, 'Year']
mle_CIs$year2 = metadata[mle_CIs$individual2, 'Year']

# Add time dist based on date
mle_CIs$time_dist = abs(difftime(as.Date(mle_CIs$date1), 
                                 as.Date(mle_CIs$date2), units = 'weeks'))

# Add sample comp
mle_CIs$sample_comp = apply(mle_CIs[, c("individual1", "individual2")], 1, 
                            function(x) paste(sort(x), collapse = "_"))


##############################################################
# Removal of overly precise estimates close to zero and one
# Must do this part by hand!!!
#
# Explaination: intuitively, one expects a smaller CI width with higher SNP count
# However, due to limitations of the parametric bootstrap given an mle plug-in,
# there can be some anomalies where rhat is equal to zero or one
# we remove these anomalies by bespoke visual inspection. 
##############################################################

#===== Boundry close to zero =====
unrelated_ind <- mle_CIs$r2.5. < eps & mle_CIs$r97.5. < 0.5

# Extract point estimates and map to zero one to project to grayscale
rhats_unrelated <- mle_CIs$rhat[unrelated_ind] 
rhats_unrelated_01 <- (rhats_unrelated  - min(rhats_unrelated )) / diff(range(rhats_unrelated))

# Plot: no outliers
plot(y = mle_CIs$CI_width[unrelated_ind], 
     x = mle_CIs$snp_count[unrelated_ind], 
     bg = gray(rhats_unrelated_01), 
     pch = 21, cex = 0.5, bty = "n", 
     ylab = "CI width", xlab = "SNP count") 



#===== Boundry close to one =====
clone_ind <- mle_CIs$r2.5. > 0.5 & mle_CIs$r97.5. > (1-eps)

# Extract point estimates and map to zero one to project to grayscale
rhats_clones <- mle_CIs$rhat[clone_ind] 
rhats_clones_01 <- (rhats_clones - min(rhats_clones)) / diff(range(rhats_clones))

# All: problem seems to dissappear above 50 SNPs
plot(y = mle_CIs$CI_width[clone_ind], 
     x = mle_CIs$snp_count[clone_ind], 
     bg = gray(1-rhats_clones_01), 
     pch = 21, cex = 0.5, bty = "n", 
     ylab = "CI width", xlab = "SNP count") 
abline(v = 30, lty = "dashed")

# Zoom in: problem up to 30 SNPs 
plot(y = mle_CIs$CI_width[clone_ind], 
     x = mle_CIs$snp_count[clone_ind], 
     bg = gray(1-rhats_clones_01), 
     pch = 21, cex = 0.5, bty = "n", 
     ylab = "CI width", xlab = "SNP count", 
     xlim = c(0,50), ylim = c(0,0.1)) 
abline(v = 30, lty = "dashed")

# Lets manually replace point estimates with CIs < 0.02 and SNP count < 30 with NAs
to_filter <- mle_CIs$CI_width < 0.02 & # Too tight CI
  mle_CIs$snp_count <= 30 & # Too few data
  mle_CIs$r2.5. > 0.5 & mle_CIs$r97.5. > (1-eps) # Clone

#===============================
# Filtering mles
#===============================
par(mfrow = c(1,2))
summarise_mles(mle_CIs, zoom = T)
mle_CIs <- mle_CIs[!to_filter,] # Remove boundary problems
summarise_mles(mle_CIs, zoom = T)
mle_CIs <- mle_CIs[!(mle_CIs$r2.5. < eps & mle_CIs$r97.5. > (1-eps)), ] # Remove uninformative
summarise_mles(mle_CIs, zoom = T)

# Save data frame
save(mle_CIs, file = sprintf("../../RData/mles_CIs_extended_freqs%s_meta.RData", freqs_used))



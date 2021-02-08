#################################################################
# README!!!!!!!!!!!!!!!!
# This script adds metadata to the dataframe of relatedness estimates
# We also remove overly precise estimates close to zero and one
# Must do this part by hand!!!
#################################################################
rm(list = ls())
source("../../Code/Extended_analysis/summarise_mles.R")
eps <- 0.01 # Used to define clones in Taylor et al. 2020

# Load metadata
load(file = "../../RData/metadata_extended.RData")

# Load relatedness estimates
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s.RData', freqs_used)) # Load All_results
mle_CIs$CI_width <- (mle_CIs$r97.5. - mle_CIs$r2.5.) # Add CI width to relatedness estimates

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
rhats_unrelated_01 <- (rhats_unrelated  - min(rhats_unrelated, na.rm = T)) / diff(range(rhats_unrelated, na.rm = T))

# Plot: no outliers
plot(y = mle_CIs$CI_width[unrelated_ind], 
     x = mle_CIs$snp_count[unrelated_ind], 
     bg = sapply(1-rhats_unrelated_01, function(x) ifelse(is.na(x), NA, gray(x))), 
     pch = 21, cex = 0.5, bty = "n", 
     ylab = "CI width", xlab = "SNP count") 



#===== Boundry close to one =====
clone_ind <- mle_CIs$r2.5. > 0.5 & mle_CIs$r97.5. > (1-eps)

# Extract point estimates and map to zero one to project to grayscale
rhats_clones <- mle_CIs$rhat[clone_ind] 
rhats_clones_01 <- (rhats_clones - min(rhats_clones, na.rm = T)) / diff(range(rhats_clones, na.rm = T))

# All: problem seems to dissappear above 50 SNPs
plot(y = mle_CIs$CI_width[clone_ind], 
     x = mle_CIs$snp_count[clone_ind], 
     bg = sapply(1-rhats_clones_01, function(x) ifelse(is.na(x), NA, gray(x))), 
     pch = 21, cex = 0.5, bty = "n", 
     ylab = "CI width", xlab = "SNP count") 
abline(v = 30, lty = "dashed")

# Zoom in: problem up to 30 SNPs 
plot(y = mle_CIs$CI_width[clone_ind], 
     x = mle_CIs$snp_count[clone_ind], 
     bg = sapply(1-rhats_clones_01, function(x) ifelse(is.na(x), NA, gray(x))), 
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
mle_CIs <- mle_CIs[!is.na(mle_CIs$rhat),] # Remove NAs
summarise_mles(mle_CIs, zoom = T)

# Save data frame
save(mle_CIs, file = sprintf("../../RData/mles_CIs_extended_freqs%s_meta.RData", freqs_used))


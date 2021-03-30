###############################################################################
#' READ ME (Must do this part by hand!) 
#' Removal of NAs and overly precise estimates close to zero and one.
#'
#' Explanation: intuitively, one expects a smaller CI width with higher SNP
#' count. However, due to limitations of the parametric bootstrap given an mle
#' plug-in, there can be some anomalies where rhat is equal to zero or one (but
#' less so for zero because IBS is compatible with "not IBD" but "not IBS" is
#' not compatible with IBD, besides when there are errors). We remove these
#' anomalies by visual inspection.
###############################################################################
rm(list = ls())
source("../../Code/Extended_analysis/summarise_mles.R")
eps <- 0.01 # Used to define clones in Taylor et al. 2020
PDF <- TRUE

if(PDF) pdf("../../Plots/filter_mles_CIs_extended.pdf")

# Load metadata (to pass snp counts to summarise_mles function)
load(file = "../../RData/metadata_extended.RData")

# Load relatedness estimates
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s.RData', freqs_used)) # Load All_results
mle_CIs$CI_width <- (mle_CIs$r97.5. - mle_CIs$r2.5.) # Add CI width to relatedness estimates

#===== Boundry close to zero =====
unrelated_ind <- mle_CIs$r2.5. < eps & mle_CIs$r97.5. < 0.5

# Extract point estimates and map to range zero to one to project to grayscale
rhats_unrelated <- mle_CIs$rhat[unrelated_ind] 
rhats_unrelated_01 <- (rhats_unrelated  - min(rhats_unrelated, na.rm = T)) / diff(range(rhats_unrelated, na.rm = T))

# Visually identify outliers: no outliers
# ?gray: zero indicates "black" and one indicates "white"
par(mfrow = c(1,1))
plot(y = mle_CIs$CI_width[unrelated_ind], 
     x = mle_CIs$snp_count[unrelated_ind], 
     bg = sapply(rhats_unrelated_01, function(x) ifelse(is.na(x), NA, gray(x))), 
     pch = 21, cex = 0.5, bty = "n", cex.main = 1, 
     main = sprintf("Lower CI limit < %s and upper CI limit < 0.5",eps), 
     ylab = "95% confidence interval (CI) width", xlab = "SNP count") 
legend("bottomleft", pch = 21, pt.bg = c("black", "white"), 
       bty = "n", title = "estimated \n relatedness",  
       legend = c(format(min(rhats_unrelated), digits = 2), 
                  format(max(rhats_unrelated), digits = 2)))

writeLines(sprintf("%s estimates are classified unrelated here", length(rhats_unrelated)))


#===== Boundry close to one =====
clone_ind <- mle_CIs$r2.5. > 0.5 & mle_CIs$r97.5. > (1-eps)

# Extract point estimates and map to zero one to project to grayscale
rhats_clones <- mle_CIs$rhat[clone_ind] 
rhats_clones_01 <- (rhats_clones - min(rhats_clones, na.rm = T)) / diff(range(rhats_clones, na.rm = T))

# All: problem seems to dissappear above 30 SNPs
plot(y = mle_CIs$CI_width[clone_ind], 
     x = mle_CIs$snp_count[clone_ind], 
     bg = sapply(1-rhats_clones_01, function(x) ifelse(is.na(x), NA, gray(x))), 
     pch = 21, cex = 0.5, bty = "n", cex.main = 1, 
     main = sprintf("Lower CI limit > 0.5 and upper CI limit > %s",1-eps), 
     ylab = "95% confidence interval (CI) width", xlab = "SNP count") 
legend("topright", pch = 21, pt.bg = c("black", "white"), inset = 0.05, 
       bty = "n", title = "estimated \n relatedness", 
       legend = c(format(max(rhats_clones), digits = 2, drop0trailing = F), 
                  format(min(rhats_clones), digits = 2)))

# Zoom in: 
plot(y = mle_CIs$CI_width[clone_ind], 
     x = mle_CIs$snp_count[clone_ind], 
     bg = sapply(1-rhats_clones_01, function(x) ifelse(is.na(x), NA, gray(x))), 
     pch = 21, cex = 0.5, bty = "n", 
     ylab = "95% confidence interval (CI) width", xlab = "SNP count", 
     xlim = c(0,50), ylim = c(0,0.15)) 
abline(v = 30, lty = "dashed")
abline(h = 0.07, lty = "dashed")


# Lets manually filter point estimates with CIs < 0.02 and SNP count < 30 with NAs
to_filter <- mle_CIs$CI_width < 0.07 & # Too tight CI
  mle_CIs$snp_count <= 30 & # Too few data
  mle_CIs$r2.5. > 0.5 & mle_CIs$r97.5. > (1-eps) # Clone

writeLines(sprintf("Among %s estimates clonally-classified estimates, %s are outliers", 
                   length(rhats_clones), sum(to_filter)))



#===============================
# Filtering mles
# Sample comparisons with only one
# marker remain but these are not 
# among the clones. 
#===============================
if(!PDF) par(mfrow = c(1,2))
summarise_mles(mle_CIs, metadata_ = metadata, zoom = T)
mle_CIs <- mle_CIs[!to_filter,] # Remove boundary problems
summarise_mles(mle_CIs, metadata_ = metadata,  zoom = T)

# How many still uninformative? 
par(mfrow = c(1,1))
hist(mle_CIs$CI_width, 
     xlab = "95% Confidence interval width", main = "")
mean(mle_CIs$CI_width > 0.5)
mean(mle_CIs$CI_width > 0.95)
abline(v = 0.95, lty = 'dotted')

# Replot all to check all problematic removed
# Looks like the negative trends is somehow stratified
par(mfrow = c(1,1))
clone_ind <- mle_CIs$r2.5. > 0.5 & mle_CIs$r97.5. > (1-eps)
rhats_clones <- mle_CIs$rhat[clone_ind] 
rhats_clones_01 <- (rhats_clones - min(rhats_clones, na.rm = T)) / diff(range(rhats_clones, na.rm = T))
plot(y = mle_CIs$CI_width[clone_ind], 
     x = mle_CIs$snp_count[clone_ind], 
     bg = sapply(1-rhats_clones_01, function(x) ifelse(is.na(x), NA, gray(x))), 
     pch = 21, cex = 0.5, bty = "n", 
     ylab = "CI width", xlab = "SNP count") 
abline(v = 50, lty = "dashed")

# Save data frame
save(mle_CIs, file = sprintf("../../RData/mles_CIs_extended_freqs%s_filtered.RData", freqs_used))

if(PDF) dev.off()

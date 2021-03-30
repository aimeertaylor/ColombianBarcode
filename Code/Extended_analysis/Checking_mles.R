################################################################################
#' This script compares estimates based on marker-reordered samples from Taylor
#' et al. 2020 and those for the extended analysis. We then looks for CI limit
#' inconsistencies, of which there are a few on the LCI, due to underflow
#' problems
################################################################################
rm(list = ls())

# Original mles (NAs were not considered in CIs
# and there was a very minor NA problem in hmmloglikelihood.cpp)
load('../../RData/mles_CIs.RData')
mles_reordered <- mle_CIs

# Extended mles (NAs were considered in CIs)
load("../../RData/mles_CIs_extended_freqsTaylor2020.RData")
mles_extended <- mle_CIs

# Match names 
load('../../RData/SNPData.RData')
key <- as.character(SNPData$SAMPLE.CODE)
names(key) <- rownames(SNPData)

nrow(mles_reordered)
nrow(mles_extended)

rownames(mles_reordered) <- paste(key[as.character(mles_reordered$individual1)], 
                                 key[as.character(mles_reordered$individual2)], sep = "_")

rownames(mles_extended) <- paste(mles_extended$individual1, 
                                 mles_extended$individual2, sep = "_")

# Check all there
all(rownames(mles_reordered) %in% rownames(mles_extended))

# Compare estimates: don't expect exact match (248 vs 250 SNPs)
plot(x = mles_reordered$rhat, 
     y = mles_extended[rownames(mles_reordered), "rhat"])

# Compare CIs: don't expect exact match (248 vs 250 SNPs)
plot(x = mles_reordered$r97.5., 
     y = mles_extended[rownames(mles_reordered), "r97.5."])

# Compare exactly 0.5 
any(mles_reordered$rhat == 0.5) # No 
any(mles_extended[rownames(mles_reordered), "rhat"] == 0.5) # No 

# Compare uninformative among common samples
eps <- 0.01
sum((mles_reordered$r2.5. < eps) & (mles_reordered$r97.5. > (1-eps)), na.rm = T) 
sum((mles_extended[rownames(mles_reordered), "r2.5."] < eps) & 
      (mles_extended[rownames(mles_reordered), "r97.5."] > (1-eps)), na.rm = T) 

# Compare both related and clonal
sum((mles_reordered$r2.5. < 0.5) & (mles_reordered$r97.5. > (1-eps)), na.rm = T) 
sum((mles_extended[rownames(mles_reordered), "r2.5."] < 0.5) & 
      (mles_extended[rownames(mles_reordered), "r97.5."] > (1-eps)), na.rm = T) 

# Are clones the same? Not exactly, but only four differences (with 3069D0R in 3 of 4)
clones_original <- rownames(mles_reordered)[which(mles_reordered$r97.5. > (1-eps))]
clones_extended <- rownames(mles_reordered)[mles_extended[rownames(mles_reordered), "r97.5."] > (1-eps)]
clones_extended_only <- setdiff(clones_extended, clones_original)
clones_original_only <- setdiff(clones_original, clones_extended)
c(clones_extended_only, clones_original_only) # Only four differences

# Compare CI widths 
mles_reordered$CI_width <- mles_reordered$r97.5. - mles_reordered$r2.5.
mles_extended$CI_width <- mles_extended$r97.5. - mles_extended$r2.5.
summary(mles_reordered$CI_width)
summary(mles_extended[rownames(mles_reordered), "CI_width"])
plot(x = mles_reordered$CI_width, 
     y = mles_extended[rownames(mles_reordered), "CI_width"])

#===== Check for inconsistent CI limits =====#

# Original estimates: 
ind_rhat_less_thn_LCI_original <- mles_reordered$rhat < mles_reordered$r2.5. 
ind_rhat_bigg_thn_UCI_original <- mles_reordered$rhat > mles_reordered$r97.5.
# Seem to be problems with LCI but not UCI
sum(ind_rhat_less_thn_LCI_original)
sum(ind_rhat_bigg_thn_UCI_original)
# Problem with original estimates: seems to be an underflow issue
mles_rhat_less_thn_LCI_original <- mles_reordered[ind_rhat_less_thn_LCI_original, ]
unique(mles_rhat_less_thn_LCI_original[,c("r2.5.", "rhat")])

# Check go away when round
UCI <- round(mles_reordered$r97.5., 5)
LCI <- round(mles_reordered$r2.5., 5)
est <- round(mles_reordered$rhat, 5)
mles_reordered[which((LCI > est) | (UCI < est)), ]

# Across all matched extended estimates: also seems to be an underflow issue
ind_rhat_less_thn_LCI_extended <- mles_extended[rownames(mles_reordered), 'rhat'] < 
  mles_extended[rownames(mles_reordered), "r2.5."] 
ind_rhat_bigg_thn_UCI_extended <- mles_extended[rownames(mles_reordered), 'rhat'] > 
  mles_extended[rownames(mles_reordered), 'r97.5.']
# Seem to be problems with LCI but not UCI
sum(ind_rhat_less_thn_LCI_extended, na.rm = T)
sum(ind_rhat_bigg_thn_UCI_extended, na.rm = T)
# Problem with original estimates: seems to be an underflow issue
mles_rhat_less_thn_LCI_extended <- mles_extended[rownames(mles_reordered),][ind_rhat_less_thn_LCI_extended, ]
unique(mles_rhat_less_thn_LCI_extended[,c("r2.5.", "rhat")])

# Across all extended estimates: also seems to be an underflow issue
ind_rhat_less_thn_LCI_extended <- mles_extended$rhat < mles_extended$r2.5.
ind_rhat_bigg_thn_UCI_extended <- mles_extended$rhat > mles_extended$r97.5.
# problem on LCI only
sum(ind_rhat_less_thn_LCI_extended, na.rm = T)
sum(ind_rhat_bigg_thn_UCI_extended, na.rm = T)
# Problem on inspection appears to be under/over flow for rhat < r2.5.
mles_rhat_less_thn_LCI_extended <- mles_extended[ind_rhat_less_thn_LCI_extended, ]
unique(mles_rhat_less_thn_LCI_extended[,c("rhat","r2.5.")])

# Check it goes away with rounding: yes, 
# goes away if we round to two decimal places
# can check by manually changing the precision to which we round 
UCI <- round(mles_extended$r97.5., 2) 
LCI <- round(mles_extended$r2.5., 2)
est <- round(mles_extended$rhat, 2)
mles_extended[which((LCI > est) | (UCI < est)), ]

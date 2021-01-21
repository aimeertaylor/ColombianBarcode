rm(list = ls())

# Original mles (NAs were not considered in CIs)
load('../../RData/mles_CIs.RData')
mles_original <- mle_CIs

# Extended mles (NAs were considered in CIs)
load("../../RData/mles_CIs_extended_freqsTaylor2020.RData")
mles_extended <- mle_CIs

# Match names 
load('../../RData/SNPData.RData')
key <- as.character(SNPData$SAMPLE.CODE)
names(key) <- rownames(SNPData)

nrow(mles_original)
nrow(mles_extended)

rownames(mles_original) <- paste(key[as.character(mles_original$individual1)], 
                                 key[as.character(mles_original$individual2)], sep = "_")

rownames(mles_extended) <- paste(mles_extended$individual1, 
                                 mles_extended$individual2, sep = "_")

# Check all there
all(rownames(mles_original) %in% rownames(mles_extended))

# Compare estimates: yes, match. 
plot(x = mles_original$rhat, 
     y = mles_extended[rownames(mles_original), "rhat"])

# Compare CIs: no, not exact match (random numbers different)
plot(x = mles_original$r97.5., 
     y = mles_extended[rownames(mles_original), "r97.5."])

# Compare uninformative
eps <- 0.01
sum((mles_original$r2.5. < eps) & (mles_original$r97.5. > (1-eps)), na.rm = T) 
sum((mles_extended[rownames(mles_original), "r2.5."] < eps) & 
      (mles_extended[rownames(mles_original), "r97.5."] > (1-eps)), na.rm = T) 

# Compare both related and clonal
sum((mles_original$r2.5. < 0.5) & (mles_original$r97.5. > (1-eps)), na.rm = T) 
sum((mles_extended[rownames(mles_original), "r2.5."] < 0.5) & 
      (mles_extended[rownames(mles_original), "r97.5."] > (1-eps)), na.rm = T) 

# Are clones the same? Not exactly, but only four differences
clones_original <- rownames(mles_original)[which(mles_original$r97.5. > (1-eps))]
clones_extended <- rownames(mles_original)[mles_extended[rownames(mles_original), "r97.5."] > (1-eps)]
clones_extended_only <- setdiff(clones_extended, clones_original)
clones_original_only <- setdiff(clones_original, clones_extended)
c(clones_extended_only, clones_original_only) # Only four differences

# Compare CI widths (expect extended to be wider but not)
mles_original$CI_width <- mles_original$r97.5. - mles_original$r2.5.
mles_extended$CI_width <- mles_extended$r97.5. - mles_extended$r2.5.
summary(mles_original$CI_width)
summary(mles_extended[rownames(mles_original), "CI_width"])
plot(x = mles_original$CI_width, 
     y = mles_extended[rownames(mles_original), "CI_width"])


#===== Check for inconsistent CI limits =====#

# Original estimates: 
ind_rhat_less_thn_LCI_original <- mles_original$rhat < mles_original$r2.5. 
ind_rhat_bigg_thn_UCI_original <- mles_original$rhat > mles_original$r97.5.
# Seem to be problems with LCI but not UCI
sum(ind_rhat_less_thn_LCI_original)
sum(ind_rhat_bigg_thn_UCI_original)
# Problem with original estimates: seems to be an underflow issue
mles_rhat_less_thn_LCI_original <- mles_original[ind_rhat_less_thn_LCI_original, ]
unique(mles_rhat_less_thn_LCI_original[,c("r2.5.", "rhat")])

# Check its really an under/overflow problem by rounding: yes
UCI <- round(mles_original$r97.5., 5)
LCI <- round(mles_original$r2.5., 5)
est <- round(mles_original$rhat, 5)
mles_original[which((LCI > est) | (UCI < est)), ]


# Across all matched extended estimates: also seems to be an underflow issue
ind_rhat_less_thn_LCI_extended <- mles_extended[rownames(mles_original), 'rhat'] < 
  mles_extended[rownames(mles_original), "r2.5."] 
ind_rhat_bigg_thn_UCI_extended <- mles_extended[rownames(mles_original), 'rhat'] > 
  mles_extended[rownames(mles_original), 'r97.5.']
# Seem to be problems with LCI but not UCI
sum(ind_rhat_less_thn_LCI_extended, na.rm = T)
sum(ind_rhat_bigg_thn_UCI_extended, na.rm = T)
# Problem with original estimates: seems to be an underflow issue
mles_rhat_less_thn_LCI_extended <- mles_extended[rownames(mles_original),][ind_rhat_less_thn_LCI_extended, ]
unique(mles_rhat_less_thn_LCI_extended[,c("r2.5.", "rhat")])


# Across all matched extended estimates: also seems to be an underflow issue
ind_rhat_less_thn_LCI_extended <- mles_extended$rhat < mles_extended$r2.5.
ind_rhat_bigg_thn_UCI_extended <- mles_extended$rhat > mles_extended$r97.5.
# problem on both ends
sum(ind_rhat_less_thn_LCI_extended, na.rm = T)
sum(ind_rhat_bigg_thn_UCI_extended, na.rm = T)
# Problem on inspection appears to be under/over flow for rhat < r2.5.
mles_rhat_less_thn_LCI_extended <- mles_extended[ind_rhat_less_thn_LCI_extended, ]
unique(mles_rhat_less_thn_LCI_extended[,c("rhat","r2.5.")])
# Problem on inspection appears to be under/over flow for rhat > r97.5.
mles_rhat_bigg_thn_UCI_extended <- mles_extended[ind_rhat_bigg_thn_UCI_extended, ]
unique(mles_rhat_bigg_thn_UCI_extended[,c("rhat","r97.5.")])

# Is it always associated with a low snpcount? No!
sort(unique(mles_rhat_bigg_thn_UCI_extended$snp_count))
example_pairs <- rownames(unique(mles_rhat_bigg_thn_UCI_extended[,c("rhat","r97.5.")]))
mles_rhat_bigg_thn_UCI_extended[example_pairs, c("rhat","r97.5.","snp_count")]

# Check its really an under/overflow problem by rounding: yes, 
# goes away if we round to two decimal places
# can check by manually changing the precision to which we round 
UCI <- round(mles_extended$r97.5., 2) 
LCI <- round(mles_extended$r2.5., 2)
est <- round(mles_extended$rhat, 2)
mles_extended[which((LCI > est) | (UCI < est)), ]

# Example problems: any of example_pairs
# (check by hand computation returned by Generate_mles_CIs_extended):
which(rownames(mles_extended) == "PW0082-C_Pf007") # 124402 # snp_count = 8
which(rownames(mles_extended) == "PW0104-C_SPT26297") # 127456 # snp_count = 165!

##############################################################
# Script to generate matrices of relatedness between 
# "new" samples/CCs and "old" CCs (CCs of Taylor et al.)
# 
# Notes: 
# - Takes a few moments to run
# - Could extend to samples of Tayor et al. 
##############################################################
rm(list = ls())
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019

# Load relatedness results
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s_meta.RData', freqs_used)) 
mle_CIs <- mle_CIs[!is.na(mle_CIs$rhat), ] # remove NAs 
mle_CIs <- mle_CIs[!(mle_CIs$r2.5. < eps & mle_CIs$r97.5. > (1-eps)), ] # Remove uninformative

# Load Clonal components 
load("../../RData/Clonal_components.RData")
CC_original <- Clonal_components
load("../../RData/Clonal_components_extended_FSVC_sid.RData")
CC_extended <- Clonal_components

# Load metadata for sorting sids below
load("../../RData/metadata_extended.RData")

# Get IDs of samples that remain after removing 
# overly precise estimates close to zero and one; see Format_mles_CIs_extended.R
# NAs and uninformative (above)
FSVC_sid_all <- metadata$SAMPLE.CODE[!metadata$PloSGen2020]
FSVC_sid_ind <- FSVC_sid_all %in% c(mle_CIs$individual1, mle_CIs$individual2) 
FSVC_sid <- FSVC_sid_all[FSVC_sid_ind]
names(FSVC_sid) = FSVC_sid

#======================================================
# Relatedness estimates between 
# Clonal components of FSVC and Taylor et al 2020
#======================================================
cc_extended_to_cc_original <- lapply(CC_extended, function(cc_e){
  sapply(CC_original, function(cc_o){
    inds <- mle_CIs$individual1 %in% cc_o & 
      mle_CIs$individual2 %in% cc_e |
      mle_CIs$individual1 %in% cc_e & 
      mle_CIs$individual2 %in% cc_o
    if(any(inds)){
      c(av_rhat = mean(mle_CIs[inds,"rhat"], na.rm = T), 
        av_r2.5 = mean(mle_CIs[inds,"r2.5."], na.rm = T), 
        av_r97.5 = mean(mle_CIs[inds,"r97.5."], na.rm = T), 
        mn_r2.5 = min(mle_CIs[inds,"r2.5."], na.rm = T), 
        mx_r97.5 = max(mle_CIs[inds,"r97.5."], na.rm = T))
    } else {
      c(av_rhat = NA, av_r2.5 = NA, av_r97.5 = NA, mn_r2.5 = NA, mx_r97.5 = NA)
    }
  })
}) 

#======================================================
# Relatedness estimates between FSVC samples and 
# Clonal components of Taylor et al 2020
#======================================================
sid_extended_to_cc_original <- lapply(FSVC_sid, function(sid){
  sapply(CC_original, function(cc){
    inds <- mle_CIs$individual1 == sid & 
      mle_CIs$individual2 %in% cc |
      mle_CIs$individual1 %in% cc & 
      mle_CIs$individual2 == sid
    if(any(inds)){
      c(av_rhat = mean(mle_CIs[inds,"rhat"], na.rm = T), 
        av_r2.5 = mean(mle_CIs[inds,"r2.5."], na.rm = T), 
        av_r97.5 = mean(mle_CIs[inds,"r97.5."], na.rm = T), 
        mn_r2.5 = min(mle_CIs[inds,"r2.5."], na.rm = T), 
        mx_r97.5 = max(mle_CIs[inds,"r97.5."], na.rm = T))
    } else {
      c(av_rhat = NA, av_r2.5 = NA, av_r97.5 = NA, mn_r2.5 = NA, mx_r97.5 = NA)
    }
  })
}) 

# Remove any comparisons that are entirely NA because a given sample
# had no comparisons for a given CC
to_keep <- !sapply(sid_extended_to_cc_original, function(x) all(is.na(x)))
sid_extended_to_cc_original <- sid_extended_to_cc_original[to_keep]


# Save results
save(sid_extended_to_cc_original, cc_extended_to_cc_original, 
     file = "../../RData/relatedness_to_CCs.RData")


# NB: although uninformative relatedness values are removed,
# the averaged values range from zero to one
range(cc_extended_to_cc_original)
range(sid_extended_to_cc_original, na.rm = T)




##############################################################################
#' This script was written to find, among FSVC samples, those to filter before
#' generating graphs. We do this by computing the number of NA relatedness
#' estimates per sample, removing the sample with the highest count and
#' iterating. NB, we allow any CI width and see what CI widths we're left with
#' after removing samples. NB, although there is a clear association between
#' per-sample NA relatedness counts and per-sample marker data count, the
#' removed sample doesn't always have the minimum snp count.
#'
#' Could we use the transitivity property of clones etc. to impute NAs? This
#' third point is a project in itself.
##############################################################################
rm(list = ls())
library(igraph) 
source('../igraph_functions.R') # construct_adj_matrix

# Load relatedness results and metadata
load('../../RData/mles_CIs_extended_freqsTaylor2020_meta.RData') 
load('../../RData/metadata_extended.RData')

# Extract sids ordered by date
FSVC_sid <- dplyr::arrange(metadata[!metadata$PloSGen2020, ], COLLECTION.DATE)$SAMPLE.CODE
FSVC_pair_ind <- (mle_CIs$individual1 %in% FSVC_sid) & (mle_CIs$individual2 %in% FSVC_sid)

# Extract estimates into an upper adjacency matrix and make symmetric
A_est <- construct_adj_matrix(mle_CIs[FSVC_pair_ind, ], Entry = 'rhat')
A_est_full <- A_est
A_est_full[lower.tri(A_est_full)] <- t(A_est_full)[lower.tri(A_est_full)]
diag(A_est_full) <- 1
fields::image.plot(A_est)
fields::image.plot(A_est_full) 

# Compute number of NA comparisons per sample
na_count_per_sample <- rowSums(is.na(A_est_full))
snp_count_per_sample <- metadata[names(na_count_per_sample), "snp_count"]

# Initiate progressive removal of samples
A_est_full_new <- A_est_full
na_count_per_sample_new <- na_count_per_sample
snp_count_per_sample_new <- snp_count_per_sample
A_est_fulls <- list(A_est_full_new); i <- 2

# Progressively remove samples
while(max(na_count_per_sample_new) > 0){
  to_remove <- names(sort(na_count_per_sample_new, decreasing = T)[1])
  A_est_full_new <- A_est_full_new[setdiff(rownames(A_est_full_new), to_remove), 
                                   setdiff(colnames(A_est_full_new), to_remove)]
  na_count_per_sample_new <- rowSums(is.na(A_est_full_new))
  A_est_fulls[[i]] <- A_est_full_new
  i <- i + 1
}

# Visualise the progressive removal of samples
for(i in 1:length(A_est_fulls)){
  
  A_est_full_new <- A_est_fulls[[i]]
  na_count_per_sample_new <- rowSums(is.na(A_est_full_new))
  snp_count_per_sample_new <- metadata[names(na_count_per_sample_new), "snp_count"]
  names(snp_count_per_sample_new) <- names(na_count_per_sample_new)
  
  plot(x = snp_count_per_sample_new, 
       y = na_count_per_sample_new, 
       ylim = c(0,max(na_count_per_sample)), 
       xlim = c(1,max(snp_count_per_sample)), 
       xlab = "Per-sample marker data count", 
       ylab = "Per-sample NA relatedness count",
       bty = "n", pch = 20, cex.main = 1, cex = 0.5)
  abline(v = min(snp_count_per_sample_new), lty = "dashed")

  if(i < length(A_est_fulls)) {
    sid_removed <- setdiff(rownames(A_est_full_new), rownames(A_est_fulls[[i+1]])) 
    points(x = snp_count_per_sample_new[sid_removed], 
           y = na_count_per_sample_new[sid_removed])
  }
}

# Summarise samples kept and removed
sids_keep <- unique(unlist(dimnames(A_est_full_new)))

# Samples excluded in Filter_mles_CIs_extended.R: 
setdiff(FSVC_sid, c(mle_CIs$individual1[FSVC_pair_ind], mle_CIs$individual2[FSVC_pair_ind]))

# Samples removed here only: 
setdiff(c(mle_CIs$individual1[FSVC_pair_ind], mle_CIs$individual2[FSVC_pair_ind]), sids_keep)

# Samples removed in Filter_mles_CIs_extended.R and here: 
sids_remv <- setdiff(FSVC_sid, sids_keep)

# Metadata of removed SNPs
table(metadata[FSVC_sid, "City"])
table(metadata[sids_keep, "City"])
table(metadata[sids_remv, "City"])

# If we use this threshold, do we remove all the uninformative estimates? No
inds <- mle_CIs$individual1 %in% sids_keep & mle_CIs$individual2 %in% sids_keep
range(mle_CIs$CI_width[inds])

save(sids_remv, file = "../../RData/sids_to_remove_from_graphs.RData")

source("../../Code/Extended_analysis/summarise_mles.R")
summarise_mles(mle_CIs[FSVC_pair_ind,], metadata_ = metadata)
summarise_mles(mle_CIs[inds, ], metadata_ = metadata)







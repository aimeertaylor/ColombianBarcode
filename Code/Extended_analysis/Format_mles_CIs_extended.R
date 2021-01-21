#################################################################
# This script adds metadata to the dataframe of relatedness estimates
#################################################################
rm(list = ls())
library(stringr)

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
# Load metadata
load(file = "../../RData/metadata_extended.RData")
rownames(metadata) <- metadata$SAMPLE.CODE

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

# Save data frame
save(mle_CIs, file = sprintf("../../RData/mles_CIs_extended_freqs%s_meta.RData", freqs_used))



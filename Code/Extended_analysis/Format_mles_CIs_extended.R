#################################################################
#################################################################

rm(list = ls())
library(stringr)

# ======================= Load results/data =======================
# Relatedness estimates
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s.RData', freqs_used)) # Load All_results
# nrow(mle_CIs) == choose(577,2)

# Load metadata 
metadata <- read.csv("../../fwdgoldengaterecodeddata/Colombia-Ecuador_Barcode_Clusters_Missing40_Master.csv")
rownames(metadata) <- as.character(metadata$Sample)

# ========== Plot relatedness estimates and CIs ==========
Ordered_r = sort.int(mle_CIs$rhat, index.return = T) # Order estimates

# NULL plot
plot(NULL, ylim = c(0,1), xlim = c(1,length(mle_CIs$rhat)), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 1, panel.first = grid())

# Add CIs
segments(x0 = 1:length(mle_CIs$rhat), x1 = 1:length(mle_CIs$rhat),
         y0 = mle_CIs$`r2.5.`[Ordered_r$ix], y1 = mle_CIs$`r97.5.`[Ordered_r$ix],
         col = adjustcolor('gray', alpha.f = 0.5), lwd = 0.1)

# Add mles
lines(Ordered_r$x, lwd = 2) 
abline(h = 1-0.01, lty = 'dashed')


# ========== Add meta data to the results ==========
#
# # =========== Aside: check meta-data agrees with SNPData =========== 
# load("../../RData/SNPData.RData")
# rownames(SNPData) <- as.character(SNPData$SAMPLE.CODE)
# intersect_sids <- intersect(metadata$Sample, rownames(SNPData))
# length(intersect_sids)
# # Check all same: yes
# all(SNPData[intersect_sids, "City"] == metadata[intersect_sids, "Location"])
# all(SNPData[intersect_sids, "Year"] == metadata[intersect_sids, "Year"])
# #===================================================================

# Factors to characters
if(class(mle_CIs$individual1) == 'factor'){mle_CIs$individual1 = as.character(mle_CIs$individual1)}
if(class(mle_CIs$individual2) == 'factor'){mle_CIs$individual2 = as.character(mle_CIs$individual2)}

# Add cities
# Check NA returned for missing sids: metadata['agsdhfsdhf', 'Location']
mle_CIs$City1 = metadata[mle_CIs$individual1, 'Location']
mle_CIs$City2 = metadata[mle_CIs$individual2, 'Location']

# Replace NDs and NAs with "Unknown"
mle_CIs$City1[grepl("ND", mle_CIs$City1)] <- "Unknown"
mle_CIs$City2[grepl("ND", mle_CIs$City2)] <- "Unknown"
mle_CIs$City1[is.na(mle_CIs$City1)] <- "Unknown"
mle_CIs$City2[is.na(mle_CIs$City2)] <- "Unknown"

# Add city comparisons
mle_CIs$City12 = apply(mle_CIs[,c('City1','City2')], 1, function(x)paste(sort(x),collapse="_"))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# HAC until metadata discussed and updated
# Add dates and distance in time 
load("../../RData/SNPData.RData") 
rownames(SNPData) <- as.character(SNPData$SAMPLE.CODE)

# Check NA returned for missing sids: yes
# as.character(SNPData['agsdhfsdhf', 'COLLECTION.DATE'])
mle_CIs$date1 = as.character(SNPData[mle_CIs$individual1, 'COLLECTION.DATE'])
mle_CIs$date2 = as.character(SNPData[mle_CIs$individual2, 'COLLECTION.DATE'])

# For those missing a date, try to add the year
# Logical vector for sample names who only have year info (some not even that "ND")
year_only1 <- !mle_CIs$individual1 %in% rownames(SNPData) 
year_only2 <- !mle_CIs$individual2 %in% rownames(SNPData) 
mle_CIs$date1[year_only1] = paste0(metadata[mle_CIs$individual1[year_only1], "Year"], "-01-", "01")
mle_CIs$date2[year_only2] = paste0(metadata[mle_CIs$individual2[year_only2], "Year"], "-01-", "01")

# Otherwise add today's date
today <- Sys.Date() # format = "%Y-%b-%d"
mle_CIs$date1[year_only1][grepl("ND", mle_CIs$date1[year_only1])] <- today
mle_CIs$date2[year_only2][grepl("ND", mle_CIs$date2[year_only2])] <- today
mle_CIs$date1[year_only1][grepl("NA", mle_CIs$date1[year_only1])] <- today
mle_CIs$date2[year_only2][grepl("NA", mle_CIs$date2[year_only2])] <- today

mle_CIs$time_dist = abs(difftime(as.Date(mle_CIs$date1), 
                                 as.Date(mle_CIs$date2), units = 'weeks'))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Add sample comp
mle_CIs$sample_comp = apply(mle_CIs[, c("individual1", "individual2")], 1, 
                            function(x) paste(sort(x), collapse = "_"))

# Save data frame
save(mle_CIs, file = sprintf("../../RData/mles_CIs_meta_extended_freqs%s.RData", freqs_used))




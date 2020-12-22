rm(list = ls())

load('../../RData/SNPData.RData')
load('../../RData/mles_CIs.RData')
mles_original <- mle_CIs
load("../../RData/mles_CIs_extended_freqsTaylor2020.RData")
mles_extended <- mle_CIs

# Match names 
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
     
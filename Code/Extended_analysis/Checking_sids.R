rm(list = ls())

# Meta data used in Taylor et al. 2020 (Plos Genetics)
load('../../RData/SNPData.RData')
SNPData$SAMPLE.CODE <- as.character(SNPData$SAMPLE.CODE)
length(unique(SNPData$SAMPLE.CODE)) == nrow(SNPData) # No duplicate sample names
nrow(SNPData) # 325 samples in snpdata used in Taylor et al. 2020

# SNP data provided by Angela
snpdata <- read.delim("../../recodedgoldengatedata/Diego-Vladimir-Fabian_GG3D7_Recode_21Dec2020.txt", 
                      stringsAsFactors = F, check.names = F)


# 1) Check that all the 325 samples from Taylor et al. 2020 are present
all(SNPData$SAMPLE.CODE %in% colnames(snpdata))

# 2) Check that the freequencies estimated using these samples are the same
sid_inc <- colnames(snpdata) %in% SNPData$SAMPLE.CODE # Pull out 325 samples
freqs_recoded <- rowMeans(snpdata[, sid_inc], na.rm = T)
freqs_original <- colMeans(SNPData[, -c(1:5, ncol(SNPData)-1, ncol(SNPData))], na.rm = T)

# Yes, the two match exactly
plot(freqs_recoded, freqs_original)
abline(a = 0, b = 1)
abline(a = 1, b = -1)


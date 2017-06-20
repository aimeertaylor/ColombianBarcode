##############################################################################################
# Prelim-analysis of Colombian data: formatting the data s.t. it can be analysed under the HMM
# Format for HMM: chrom pos sample (0,1,-1)
#
# To-do list:
# Confirm missing ('' or '--')
# Decide how to deal with mixed (185 samples with at least one mixed and 154 without)
# Use proper positions, chrom and pos
##############################################################################################
rm(list = ls())

# Load data 
SNPMeta <- read.csv('/Users/aimeet/Documents/BroadLaptop/ColombianData/OriginalData/12863_2012_1071_MOESM1_ESM.csv', skip = 1)
SNPData <- read.csv('/Users/aimeet/Documents/BroadLaptop/ColombianData/OriginalData/ColombianDataset_from_DiegoEcheverry.csv',nrow = 325)
nsamp <- nrow(SNPData)

# Remove surplus columns and convert SNPData to a matrix
SNPData <- SNPData[,1:255]
rownames(SNPData) <- paste('sid', 1:nsamp, sep = '')
X <- as.matrix(SNPData[,6:255])
X[X == 'AA'] <- '0'
X[X == 'BB'] <- '1'
X[X == 'AB'] <- '-1' # Treat mixed as missing 
X[X == '--'] <- '-1' # Treat '--' as missing 
SNPData <- as.data.frame(SNPData)
SNPData[,6:255] <- apply(X, 1, as.numeric)
SNPData[,6:255][SNPData[,6:255] == '-1'] <- NA # Encode missing as NA for RData

# Assume SNP names are given by SNP ID in the same order
SNPnames <- as.character(SNPMeta$X.SNP.ID.[as.character(SNPMeta$X250.most.informative.SNPs) == 'x'])
pos <- sapply(strsplit(SNPnames, '-'), function(x){x[2]})
chrom <- as.numeric(sapply(sapply(strsplit(SNPnames, '-'), function(x){y <- strsplit(x[1], 'L')}), function(x){x[2]}))

# Combine chrom, pos and SNPdata into an array 
HMMData_temp <- cbind(chrom, pos, t(X))
HMMData <- NULL

# Sort order of SNPs within chromosomes
for(crom in unique(chrom)){
  ind <- chrom == crom
  X <- sort(as.numeric(pos[ind]), decreasing = FALSE, ind = TRUE)
  HMMData <- rbind(HMMData, HMMData_temp[ind, ][X$ix, ])
}

# Save data in format for HMM and RData
colnames(HMMData) <- c('chrom', 'pos', paste('sid', 1:nsamp, sep = ''))
write.table(HMMData, file = '/Users/aimeet/Documents/BroadLaptop/ColombianData/HMMInput/HMMData.txt', 
            quote = FALSE, row.names = FALSE, sep = '\t')
save(SNPData, file = '/Users/aimeet/Documents/BroadLaptop/ColombianData/RData/SNPData.RData')

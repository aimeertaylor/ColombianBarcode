################################################
# Generate extended data set
################################################
rm(list = ls())

#=====================================
# Load and process data 
#=====================================
data_set0 = read.delim("../../TxtData/hmmInputRecode.txt", stringsAsFactors = F) # Re-coded data set from Echeverry et al. 
data_set1 = read.csv("../../OriginalData/Guapi_WGStoBarcode.csv") # Data set provided by Vladimir courtesy of Manuela's email dated 24th Jan 2020 
data_set2 = read.delim("../../OriginalData/Ecuador_Colombia_2020-02-24.DiegoBarcode.txt", 
                       stringsAsFactors = F) # Data set provided by Angela (renders above obsolete) 

# Check contents
unique(unlist(data_set0[,-(1:2)]))
unique(unlist(data_set1[,-(1:2)]))
unique(unlist(data_set2[,-(1:2)])) # Why does dataset2 contain 2? How many? 

# Why does dataset2 contain 2? How many? 2
sum(data_set2[,-(1:2)] == "2", na.rm = T)

# Replace missing (-1) and het ("H") calls with NA and convert to numeric df
data_set2[,-(1:2)][data_set2[,-(1:2)] == "2"] <- NA # ****** This is a hack!!!! ******
data_set2[,-(1:2)][data_set2[,-(1:2)] == "-1"] <- NA
data_set2[,-(1:2)][data_set2[,-(1:2)] == "H"] <- NA
data_set2 <- as.data.frame(data.matrix(data_set2))
unique(unlist(data_set2[,-(1:2)])) # Check now 0,1,NA only

# Reformat colnames
colnames(data_set1)[1:2] = colnames(data_set0)[1:2]
colnames(data_set2)[1:2] = colnames(data_set0)[1:2]

# Add rownames 
rownames(data_set0) = apply(data_set0[,c("chrom", "pos")], 1, paste, collapse = "_")
rownames(data_set1) = apply(data_set1[,c("chrom", "pos")], 1, paste, collapse = "_")
rownames(data_set2) = apply(data_set2[,c("chrom", "pos")], 1, paste, collapse = "_")

# Check that data_set1 in contained within data_set2 as mentioned in Angela's email: yes
all(colnames(data_set1) %in% colnames(data_set2)) 

# Check that the data contained are the same besides NA: yes
all(data_set2[rownames(data_set1), colnames(data_set1)] == 
      data_set1[rownames(data_set1), ], na.rm = T)

# Check there are the same number of NAs: more in Angela's
sum(is.na(data_set2[,colnames(data_set1)]))
sum(is.na(data_set1[,colnames(data_set1)]))

# Plot the Colombian data from Manuela and Angela (code copied from Data_summary.Rmd)
SNPDataBinary1 <- data_set1[,-(1:2)]
SNPDataBinary2 <- data_set2[,colnames(data_set1)][,-(1:2)]
par(mfrow = c(2,1), family = 'serif', mar = c(2,2,1,1))
cols <- RColorBrewer::brewer.pal(3,'Dark2')
image(as.matrix(SNPDataBinary1[names(sort(rowSums(SNPDataBinary1, na.rm = T))),]), 
      ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)
title(ylab = sprintf('Sample ID (%s samples)', ncol(SNPDataBinary1)-2), 
      xlab = sprintf('SNP ID (%s SNPs)', nrow(SNPDataBinary1)), line = 1)
image(as.matrix(SNPDataBinary2[names(sort(rowSums(SNPDataBinary1, na.rm = T))),]), 
      ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)
title(ylab = sprintf('Sample ID (%s samples)', ncol(SNPDataBinary2)-2), 
      xlab = sprintf('SNP ID (%s SNPs)', nrow(SNPDataBinary2)), line = 1)

# Work with data_set2 instead but order samples by data_set1
columbia_samp <- which(colnames(data_set2) %in% colnames(data_set1))
data_set2 <- cbind(data_set2[, columbia_samp], data_set2[, -columbia_samp])

# Extract the missing SNPs
missing_SNPs = which(!rownames(data_set0) %in% rownames(data_set2))

# Add rows of missing SNPs to data_set2
X = cbind(data_set0$chrom[missing_SNPs], 
          data_set0$pos[missing_SNPs],  
          matrix(NA, nrow = length(missing_SNPs), ncol = ncol(data_set2)-2))
rownames(X) = rownames(data_set0)[missing_SNPs]
colnames(X) = colnames(data_set2)
data_set2 = rbind(data_set2, X)

# Concatenate data and remove rownames
data_set2 <- data_set2[rownames(data_set0),]
if(all(rownames(data_set2) == rownames(data_set0))){
  data_set = cbind(data_set0, data_set2[,-c(1,2)])
} else { stop('data sets mismatched') }

# Save the extended data_set 
save(data_set, file = "../../RData/data_set_extended.RData")

# Plot the extended data (code copied from Data_summary.Rmd)
SNPDataBinary <- data_set[,-(1:2)]
par(mfrow = c(1,1), family = 'serif', mar = c(2,2,1,1))
cols <- RColorBrewer::brewer.pal(3,'Dark2')
image(as.matrix(SNPDataBinary[names(sort(rowSums(SNPDataBinary, na.rm = T))),]), 
      ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)
title(ylab = sprintf('Sample ID (%s samples ordered by data set: Diego, Vladimir, Fabian)', ncol(data_set)-2), 
      xlab = sprintf('SNP ID (%s SNPs) (ordered by increasing alt proportion)', nrow(data_set)), line = 1)




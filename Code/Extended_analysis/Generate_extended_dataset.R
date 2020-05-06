################################################
# Generate extended data set
# For conversation, see email thread entitled 
# WGS-to-barcode data for extended Colombia-Ecuador analysis  
################################################
rm(list = ls())

#=====================================
# Load and process data 
#=====================================
data_set0 = read.delim("../../TxtData/hmmInputRecode.txt", stringsAsFactors = F) # Re-coded data set from Echeverry et al. 
data_set1 = read.csv("../../OriginalData/Guapi_WGStoBarcode.csv") # Data set provided by Vladimir courtesy of Manuela's email dated 24th Jan 2020 
data_set2 = read.delim("../../OriginalData/Ecuador_Colombia_2020-02-24.DiegoBarcode.txt", 
                       stringsAsFactors = F) # Data set provided by Angela (renders above obsolete) 

#' Data set provided by Angela Tue 6th May (should replace above) Within it, het
#' calls have been replaced by missing (since all samples are deemed monoclonal)
#' There is an additional column "GG" of SNP Names. Angela added an additional
#' read depth filter to minimize reference bias. A.s. more missing, but more
#' reliable. Also seems to have resulted in the removal of a multiallelic call
#' in sample Pf024 Previously missing SNPs (chrom = 4, pos = 156685; and chrom =
#' 13, pos = 2675856) are included. While the addition of the previously missing
#' SNP (chrom = 13, pos = 2675856, old name = MAL13-2672190) has introduced two
#' multiallelic second alt calls (one in sample Pf050, another in sample
#' SPT26229)
data_set3 = read.delim("../../OriginalData/Ecuador_Colombia.DiegoBarcode.Complete.RD5.txt", 
                       stringsAsFactors = F)

# Check contents
unique(unlist(data_set0[,-(1:2)]))
unique(unlist(data_set1[,-(1:2)]))
unique(unlist(data_set2[,-(1:2)])) # Why does dataset2 contain 2? How many? 
unique(unlist(data_set3[,-(1:3)])) # Dataset3 is numeric because doesn't contain "H"

# Why does dataset2 contain 2? How many? These are multiallelic calls confirmed by Angela over email
sum(data_set2[,-(1:2)] == "2", na.rm = T)

#' Why does dataset3 contain an additional multialleleic calls relative to
#' dataset2? The addition of the previously missing SNP (chrom = 13, pos =
#' 2675856, old name = MAL13-2672190) has introduced two multiallelic second alt
#' calls (one in sample Pf050, another in sample SPT26229). There has also been
#' the removal of a multiallelic call associated with existing SNP (chrom = 4,
#' pos = 336020, name = MAL04-342497) in sample Pf024
sum(data_set3[,-(1:3)] == 2, na.rm = T) 
j_multi2 <- which(colSums(data_set2[,-(1:2)] == 2, na.rm = T) > 0)
i_multi2 <- which(rowSums(data_set2[,-(1:2)] == 2, na.rm = T) > 0)
j_multi3 <- which(colSums(data_set3[,-(1:3)] == 2, na.rm = T) > 0)
i_multi3 <- which(rowSums(data_set3[,-(1:3)] == 2, na.rm = T) > 0)
data_set2[i_multi2, c(1:2, j_multi2+2)]
data_set3[i_multi3, c(1:3, j_multi3+3)]

# Replace multiallelic alt (2), missing (-1) and het ("H") calls with NA and convert to numeric df
data_set2[,-(1:2)][data_set2[,-(1:2)] == "2"] <- NA 
data_set2[,-(1:2)][data_set2[,-(1:2)] == "-1"] <- NA
data_set2[,-(1:2)][data_set2[,-(1:2)] == "H"] <- NA
data_set2 <- as.data.frame(data.matrix(data_set2))
unique(unlist(data_set2[,-(1:2)])) # Check now 0,1,NA only

# Replace multiallelic alt (2) in dataset3 and check 
# After discussion with, Angela this choice seems reasonable,
# "especially because Goldengate wouldn't have been able to detect these anyway"
data_set3[,-(1:3)][data_set3[,-(1:3)] == 2] <- NA 
unique(unlist(data_set3[,-(1:3)]))

# Remove additional column 
data_set3 <- data_set3[,-3]

# Reformat colnames
colnames(data_set1)[1:2] = colnames(data_set0)[1:2]
colnames(data_set2)[1:2] = colnames(data_set0)[1:2]
colnames(data_set3)[1:2] = colnames(data_set0)[1:2]

# Add rownames 
rownames(data_set0) = apply(data_set0[,c("chrom", "pos")], 1, paste, collapse = "_")
rownames(data_set1) = apply(data_set1[,c("chrom", "pos")], 1, paste, collapse = "_")
rownames(data_set2) = apply(data_set2[,c("chrom", "pos")], 1, paste, collapse = "_")
rownames(data_set3) = apply(data_set3[,c("chrom", "pos")], 1, paste, collapse = "_")

# Check that data_set1 in contained within data_set2 and data_set3 as mentioned in Angela's email: yes
all(colnames(data_set1) %in% colnames(data_set2)) 
all(colnames(data_set1) %in% colnames(data_set3)) 

# Check that the data contained are the same besides NA: yes
all(data_set2[rownames(data_set1), colnames(data_set1)] == 
      data_set1[rownames(data_set1), ], na.rm = T)
# Check that the data contained are the same besides NA: yes
all(data_set3[rownames(data_set1), colnames(data_set1)] == 
      data_set1[rownames(data_set1), ], na.rm = T)

# Check there are the same number of NAs: more in Angela's most recent
sum(is.na(data_set3[,colnames(data_set1)]))
sum(is.na(data_set2[,colnames(data_set1)]))
sum(is.na(data_set1[,colnames(data_set1)]))

# Plot the Colombian data from Manuela and Angela (code copied from Data_summary.Rmd)
SNPDataBinary1 <- data_set1[,-(1:2)]
SNPDataBinary2 <- data_set2[,colnames(data_set1)][,-(1:2)]
SNPDataBinary3 <- data_set3[,colnames(data_set1)][,-(1:2)]
par(mfrow = c(3,1), family = 'serif', mar = c(2,2,1,1))
cols <- RColorBrewer::brewer.pal(3,'Dark2')
image(as.matrix(SNPDataBinary1[names(sort(rowSums(SNPDataBinary1, na.rm = T))),]), 
      ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)
title(ylab = sprintf('Sample ID (%s samples)', ncol(SNPDataBinary1)-2), 
      xlab = sprintf('SNP ID (%s SNPs)', nrow(SNPDataBinary1)), line = 1)
image(as.matrix(SNPDataBinary2[names(sort(rowSums(SNPDataBinary1, na.rm = T))),]), 
      ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)
title(ylab = sprintf('Sample ID (%s samples)', ncol(SNPDataBinary2)-2), 
      xlab = sprintf('SNP ID (%s SNPs)', nrow(SNPDataBinary2)), line = 1)
image(as.matrix(SNPDataBinary3[names(sort(rowSums(SNPDataBinary1, na.rm = T))),]), 
      ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)
title(ylab = sprintf('Sample ID (%s samples)', ncol(SNPDataBinary3)-2), 
      xlab = sprintf('SNP ID (%s SNPs)', nrow(SNPDataBinary3)), line = 1)

# Work with data_set3 instead but order samples by data_set1
columbia_samp <- which(colnames(data_set3) %in% colnames(data_set1))
data_set3 <- cbind(data_set3[, columbia_samp], data_set3[, -columbia_samp])

# Extract the missing SNPs: no longer any missing
missing_SNPs = which(!rownames(data_set0) %in% rownames(data_set3))

# # Add rows of missing SNPs to data_set3
# X = cbind(data_set0$chrom[missing_SNPs], 
#           data_set0$pos[missing_SNPs],  
#           matrix(NA, nrow = length(missing_SNPs), ncol = ncol(data_set3)-2))
# rownames(X) = rownames(data_set0)[missing_SNPs]
# colnames(X) = colnames(data_set3)
# data_set3 = rbind(data_set3, X)

# Concatenate data and remove rownames
data_set3 <- data_set3[rownames(data_set0),]
if(all(rownames(data_set3) == rownames(data_set0))){
  data_set = cbind(data_set0, data_set3[,-c(1,2)])
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




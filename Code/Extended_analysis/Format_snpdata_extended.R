#################################################################
# Script to format the extended data set
#################################################################
rm(list = ls())

# =============== Load meta data ===============
# To check against snp data
load(file = "../../RData/metadata_extended.RData")
sids <- metadata$SAMPLE.CODE

# =============== Load data provided via Angela ===============
# email sent 21st Dec 2020 in thread entitled "Golden Gate Recoded Data"
dataset <- read.delim("../../TxtData/Diego-Vladimir-Fabian_GG3D7_Recode_V2_07March2021.txt", 
                      check.names = F, # Stops Xs being added to names beginning with a digit and conversion of "-" to "."
                      stringsAsFactors = F) 
head(dataset[,1:10]) 
colnames(dataset)[1:2] <- c("chrom", "pos") # Rename first two cols

# Overlap: exact
setequal(sids, colnames(dataset)[-(1:3)])

# Remove marker name
snpdata <- dataset[,-3]  

# Order markers 
snpdata <- dplyr::arrange(snpdata, "chrom", "pos")

# =============== Replace three 2s with NA ===============
unique(unlist(snpdata[,-(1:2)])) # "0"   "1"   NA  "2"  
table(unlist(snpdata[,-(1:2)])) # "0"   "1"  "2"  
snpdata[,-(1:2)][snpdata[,-(1:2)] == "2"] <- NA

# Check: 
unique(unlist(snpdata[,-(1:2)]))
table(unlist(snpdata[,-(1:2)]))

# =============== Plot the extended data =============== 
par(mfrow = c(1,1), family = 'serif', mar = c(2,2,1,1))
cols <- RColorBrewer::brewer.pal(3,'Dark2')
image(as.matrix(snpdata[,-(1:2)]), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)

# Add snp_count to metadata
metadata$snp_count <- apply(snpdata[, metadata$SAMPLE.CODE], 2, function(x) sum(!is.na(x)))

missing_snps <- which(apply(snpdata[,-(1:2)], 1, function(x) all(is.na(x)))) # Any markers with zero data: yes, 2
missing_sids <- which(apply(snpdata, 2, function(x) all(is.na(x)))) # Any samples with zero data: yes, 9

dim(snpdata)
snpdata <- snpdata[-missing_snps, -missing_sids]
dim(snpdata)

# =============== Save the extended data_set  ===============
save(snpdata, file = "../../RData/snpdata_extended.RData")
save(metadata, file = "../../RData/metadata_extended.RData")


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
dataset <- read.delim("../../TxtData/Diego-Vladimir-Fabian_GG3D7_Recode_V2_05March2021.txt", 
                      check.names = F, # Stops Xs being added to names beginning with a digit and conversion of "-" to "."
                      stringsAsFactors = F) 
head(dataset[,1:10]) 
colnames(dataset)[1:3] <- c("chrom", "pos", "marker_name") # Rename first three cols

# Overlap: exact
setequal(sids, colnames(dataset)[-(1:3)])

# # Redundant: samples in snpdata but not in metadata and vice versa 
# snpdata_additional_sids <- colnames(dataset)[!colnames(dataset) %in% sids][-(1:3)] 
# metadata_additional_sids <- sids[!sids %in% colnames(dataset)[-(1:2)]]



# 112 additional sids in snpdata
length(snpdata_additional_sids)
# These include 8 ref (3D7 missing) and 104 primary infection
# Echeverry et al. 2013: 384 SNPs and 447 primary infections + 9 reference
# strains attempted; 32 sids excluded due to incomplete data, 15 sids excluded due to
# identification, leaves 400 primary infections (75 polyclonal, 325 monoclonal)
# 325 monolconal are in metadata: 
metadata$SAMPLE.CODE[metadata$PloSGen2020]

# 61 of the 325 samples of Echeverry et al. 2013 are missing in the new data set 
length(metadata_additional_sids)

# Remove samples not in metadata (remove Echeverry 2013 polyclonal and excluded)
snpdata <- dataset[,c("chrom", "pos", sids[sids %in% colnames(dataset)])]  

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

which(rowSums(snpdata[,-(1:2)], na.rm = T) == 0) # Any markers with zero data: yes, one
which(colSums(snpdata, na.rm = T) == 0) # Any samples with zero data: yes

# =============== Save the extended data_set  ===============
save(snpdata, file = "../../RData/snpdata_extended.RData")



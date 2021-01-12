library(stringr) # for str_detect
rm(list = ls())

# =============== Load data provided by Angela ===============
# email sent 21st Dec 2020 in thread entitled "Golden Gate Recoded Data"
dataset <- read.delim("../../recodedgoldengatedata/Diego-Vladimir-Fabian_GG3D7_Recode_21Dec2020.txt", 
                      check.names = F, # Stops Xs being added to names beginning with a digit and conversion of "-" to "."
                      stringsAsFactors = F) 
head(dataset[,1:10]) #
colnames(dataset)[1:3] <- c("chrom", "pos", "marker_name") # Rename first three cols
snpdata <- dataset[,!colnames(dataset) %in% c("marker_name", "GGSite")]  # Remove marker name (also seems to be snpdata$GGSite)

# =============== Replace "2" with NA ===============
# SNP calls in the extended data set: 
unique(unlist(snpdata[,-(1:2)])) # "0"   "1"   NA  "2"  
table(unlist(snpdata[,-(1:2)])) # "0"   "1"  "Het" "2"  
snpdata[,-(1:2)][snpdata[,-(1:2)] == "2"] <- NA

# Check: 
unique(unlist(snpdata[,-(1:2)]))
table(unlist(snpdata[,-(1:2)]))

# =============== Plot the extended data =============== 
par(mfrow = c(1,1), family = 'serif', mar = c(2,2,1,1))
cols <- RColorBrewer::brewer.pal(3,'Dark2')
image(as.matrix(snpdata[,-(1:2)]), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)

# =============== Check if there any markers with zero data: no
any(rowSums(snpdata[,-(1:2)], na.rm = T) == 0)
which(colSums(snpdata, na.rm = T) == 0)

# =============== Save the extended data_set  ===============
# Replace previous dataset
save(snpdata, file = "../../RData/snpdata_extended.RData")



#' As of Dec 10th, waiting for complete data from Angela (61 extra samples and
#' remove low quality)
library(stringr) # for str_detect
rm(list = ls())

# =============== Load data provided by Angela ===============
# Data are stored in fwdgoldengaterecodeddata
# email sent 2nd Dec 2020 in thread entitled "Golden Gate Recoded Data"
dataset <- read.delim("../../fwdgoldengaterecodeddata/GGCalls_based_on_3d7_w_Guapi.txt", 
                      check.names = F, # Stops Xs being added to names beginning with a digit and conversion of "-" to "."
                      stringsAsFactors = F) 
head(dataset[,1:10]) # First six columns are not SNP data



# =============== Replace "het" and "2" with NA ===============
# SNP calls in the extended data set: 
unique(unlist(dataset[,-(1:6)])) # "0"   "1"   NA    "Het" "2"  
table(unlist(dataset[,-(1:6)])) # "0"   "1"  "Het" "2"  

# Create a new data frame and replace "het" and "2"
snpdata <- dataset[,-(3:6)]
snpdata[,-(1:2)][dataset[,-(1:6)] == "2" | dataset[,-(1:6)] == "Het"] <- NA

# Check: 
unique(unlist(snpdata[,-(1:2)]))
table(unlist(snpdata[,-(1:2)]))


# =============== Rename columns ===============
colnames(snpdata)[1] <- "chrom"
colnames(snpdata)[2] <- "pos"



# =============== Convert into numeric =============== 
apply(snpdata, 2, class) 
apply(dataset, 2, class) 
snpdata <- apply(snpdata, 2, as.numeric) 

# Check all the same after conversion 
all(sapply(1:ncol(dataset[,-(1:6)]), function(j){
  all(dataset[,-(1:6)][,j] == snpdata[,-(1:2)][,j], na.rm = T)
}))



# =============== Plot the extended data =============== 
par(mfrow = c(1,1), family = 'serif', mar = c(2,2,1,1))
cols <- RColorBrewer::brewer.pal(3,'Dark2')
image(snpdata[,-(1:2)], ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)


# =============== Remove row (marker) with no data =============== 
# Also rename and convert to data frame for consistency with Generate_mles_CIs_extended
snpdata <- as.data.frame(snpdata[rowSums(snpdata[,-(1:2)], na.rm = T) > 0, ]) 


# =============== Save the extended data_set  ===============
# Replace previous dataset
save(snpdata, file = "../../RData/snpdata_extended.RData")



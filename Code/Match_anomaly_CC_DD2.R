rm(list = ls())

# Load data
load("../RData/anomaly_CC_SNPData.RData")
ControlDataRaw <- read.csv("../OriginalData/12863_2012_1071_MOESM2_ESM.csv", 
                          header = T, skip = 2, stringsAsFactors = F)

# Process data
head(ControlDataRaw[,1:10])
tail(ControlDataRaw[,1:10], n = 50) # rows 385 onwards are blank so remove in next step
DD2Data = ControlDataRaw[-c(385:nrow(ControlDataRaw)),c("DD2.plate1.", "DD2.plate4.")]
rownames(DD2Data) = ControlDataRaw$Chromosome.SNP[-c(385:nrow(ControlDataRaw))]
anomaly_CC_SNPOnlyData = t(anomaly_CC_SNPData[,grepl("MAL", colnames(anomaly_CC_SNPData))])


# Colate DD2 and anomaly data 
DD2_plus_anomaly_CC <- cbind(DD2Data[rownames(anomaly_CC_SNPOnlyData),], anomaly_CC_SNPOnlyData)

# Inspect by eye (could re-ecode DD2 to major/minor based on full sample data set but unnecessary)
DD2_plus_anomaly_CC



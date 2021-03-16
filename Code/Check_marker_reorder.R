load('../RData/BackupForReAnalysis/SNPData.RData')
load('../RData/SNPData.RData')

# Check fix of alphabetically miss-ordering of markers
datarecoded <- read.delim("../TxtData/Diego-Vladimir-Fabian_GG3D7_Recode_V2_07March2021.txt", 
                          check.names = F, # Stops Xs being added to names beginning with a digit and conversion of "-" to "."
                          stringsAsFactors = F) 
datarecoded <- dplyr::arrange(datarecoded,chrom, pos)

all(datarecoded$GGMarker == colnames(SNPData[,6:255]))

# Recode the SNPData
ind_pivot <- which(apply(SNPData[, 6:255], 1, function(x) !any(is.na(x))))[1]
pivot_data <- cbind(old = as.numeric(SNPData[ind_pivot,6:255]), 
                    new = datarecoded[, as.character(SNPData$SAMPLE.CODE[ind_pivot])])
SNPData_recoded <- SNPData[, 6:255]

for(j in 1:250){
  if(any(is.na(pivot_data[j,]))){
    SNPData_recoded[,j] <- NA
  } else if(pivot_data[j,1] != pivot_data[j,2]) {
    SNPData_recoded[,j] <- abs(SNPData_recoded[,j]-1) 
  } else {
    next()
  }
}
all((as.matrix(t(SNPData_recoded)) == as.matrix(datarecoded[,SNPData$SAMPLE.CODE])), na.rm = T)
image(as.matrix(t(SNPData_recoded)))
image(as.matrix(datarecoded[,as.character(SNPData$SAMPLE.CODE)]))


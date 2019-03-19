##############################################################################################
# Script to generate Results using Colombian data
# 407 seconds
##############################################################################################
rm(list = ls())
library(bbmle) # for mle2()
library(proxy) # for dist() with modifiable function 
library(tictoc)
tic()
source('/Users/aimeet/Documents/BroadLaptop/TM_border/Rscripts/Archived_for_future_ref/distance_functions.R')
source('/Users/aimeet/Documents/BroadLaptop/TM_border/QuantLocalPfConnIBD/FunctionFiles/simtests.R')
source('/Users/aimeet/Documents/BroadLaptop/IBD_v_IBS_across_the_globe/functions.R') # For mle_estimate etc.
load('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/geo_dist_info.RData') # Load geo_dist info and raw data 
load('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/SNPData.RData') # Load and summarise raw data 
SNPDataBinary <- SNPData[,6:255] # Extract the SNPData w/o meta data
attach(geo_dist_info)

# Run HMM (12 seconds) then reload
system.time(system('~/Documents/BroadLaptop/hmmIBD-2.0.2/hmmIBD -i ../TxtData/hmmIBD_input.txt -o ../TxtData/hmmIBD_output')) # 66.905 sec
Result <- read.delim('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/TxtData/hmmIBD_output.hmm_fract.txt')

# Order of as.vector(dist)
nSamples <- nrow(SNPData)
ind_j <- rep(1:nSamples, (nSamples:1)-1)
ind_i <- c()
for(i in 2:nSamples){
  ind_i <- c(ind_i, i:nSamples)
}

# Generate IBS 
IBS <- as.vector(dist(SNPDataBinary, method = bs_distance_function))
names(IBS) <- apply(cbind(rownames(SNPData)[ind_i], rownames(SNPData)[ind_j]), 1, 
                    function(x){paste(sort(x), collapse = '', sep = '')})

# Add sample columns
Result$sample_comp <- apply(Result[,c('sample1', 'sample2')], 1, 
                            function(x){paste(sort(x), collapse = '', sep = '')})

# Add IBS, IBD and tails 
Result$IBD <- Result$fract_sites_IBD
Result$IBS <- IBS[Result$sample_comp]
threshold_IBD <- 0.5 # Define magic numbers
Result$IBD_tail <- 1*(Result$IBD > threshold_IBD) 
IBD_cdf <- ecdf(Result$IBD) # What is the empirical Pr(IBD < 0.5)
threshold_IBS <- quantile(Result$IBS, probs = IBD_cdf(threshold_IBD)) # IBS sample quantile
Result$IBS_tail <- 1*(Result$IBS > threshold_IBS) 
Thresholds = c('threshold_IBD' = threshold_IBD, 'threshold_IBS' = as.numeric(threshold_IBS))

# Add additional columns/covariates 
Result$site_comp <- apply(cbind(SNPData[as.character(Result$sample1), 'City'], 
                                SNPData[as.character(Result$sample2), 'City']), 1, 
                          function(x){paste(sort(x), sep = '', collapse = '_')})
Result$geo_dist <- pairwise_site_distance_all[Result$site_comp]
Result$Tumaco <- Result$site_comp == "Tumaco_Tumaco"
Result$Guapi <- Result$site_comp == "Guapi_Guapi"
Result$Buenaventura <- Result$site_comp == "Buenaventura_Buenaventura"
Result$Quibdo <- Result$site_comp == "Quibdo_Quibdo"
Result$Tado <- Result$site_comp == "Tado_Tado"
Result$Within <- apply(Result[,c('Tumaco', 'Guapi', 'Buenaventura', 'Quibdo', 'Tado')], 1, any)
Result$time_dist <- abs(difftime(SNPData[as.character(Result$sample1), 'COLLECTION.DATE'], 
                                 SNPData[as.character(Result$sample2), 'COLLECTION.DATE'],
                                 units = 'weeks'))

# Add mle IBD estimated under indpendence (added 5th June 2018)
Frequencies <- colMeans(SNPDataBinary, na.rm = TRUE)
system.time(mle_distance <- dist(SNPDataBinary, method = mle_function, f = Frequencies)) # Takes 348 sec
Result$IBD_indep <- as.vector(mle_distance)

# Add corrected IBS estimate (added 19th June 2018)
h_constant <- mean(2*Frequencies*(1-Frequencies), na.rm = TRUE)
Result$IBSc <- 1 + ((Result$IBS - 1)/h_constant) 

# Save Result
save(Result, file = '/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/Result.RData')
save(Thresholds, file = '/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/Thresholds.RData')
toc()
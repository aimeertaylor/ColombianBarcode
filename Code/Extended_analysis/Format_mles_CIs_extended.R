###############################################################################
#' This script adds metadata to the filtered dataframe of relatedness estimates
###############################################################################
rm(list = ls())

# Load metadata
load(file = "../../RData/metadata_extended.RData")
load(file = "../../RData/LonLats.RData")
source("../../Code/gcd.hf.R")

# Load filtered relatedness estimates
freqs_used <- "Taylor2020"
load(sprintf("../../RData/mles_CIs_extended_freqs%s_filtered.RData", freqs_used)) 

# Factors to characters
if(class(mle_CIs$individual1) == 'factor'){mle_CIs$individual1 = as.character(mle_CIs$individual1)}
if(class(mle_CIs$individual2) == 'factor'){mle_CIs$individual2 = as.character(mle_CIs$individual2)}

# Compute and add the fraction highly related using 0.25 (threshold used in main analyses of Taylor et al. 2020)
high_relatedness_threshold <- 0.25
mle_CIs$highly_related <- as.numeric(mle_CIs$r2.5. > high_relatedness_threshold)

# Add State and country, noting that metadata STATE contains a mix of both
mle_CIs$State1 = metadata[mle_CIs$individual1, 'STATE']
mle_CIs$State2 = metadata[mle_CIs$individual2, 'STATE']
mle_CIs$Country1 = mle_CIs$State1
mle_CIs$Country2 = mle_CIs$State2
ind_EV1 <- (mle_CIs$Country1 == "Ecuador" | mle_CIs$Country1 == "Venezuela")
ind_EV2 <- (mle_CIs$Country2 == "Ecuador" | mle_CIs$Country2 == "Venezuela")
mle_CIs$Country1[!ind_EV1] <- "Colombia"
mle_CIs$Country2[!ind_EV2] <- "Colombia"
mle_CIs$State1[ind_EV1 | mle_CIs$State1 == "Colombia"] <- NA
mle_CIs$State2[ind_EV2 | mle_CIs$State2 == "Colombia"] <- NA
unique(c(mle_CIs$Country1, mle_CIs$Country2)) # Check countries
unique(c(mle_CIs$State1, mle_CIs$State2)) # Check states

# Add cities
# Check NA returned for missing sids: 
is.na(metadata['agsdhfsdhf', 'City'])
mle_CIs$City1 = metadata[mle_CIs$individual1, 'City']
mle_CIs$City2 = metadata[mle_CIs$individual2, 'City']

# Add city comparisons
mle_CIs$City12 = apply(mle_CIs[,c('City1','City2')], 1, function(x)paste(sort(x),collapse="_"))

# Add city distance 
city_comb <- unique(mle_CIs[, c("City1", "City2", "City12")])
intercity_dists <- sapply(1:nrow(city_comb), function(i){
  gcd.hf(A = as.matrix(LonLats[city_comb[i, "City1"], c('Longitude', 'Latitude')]), 
         B = as.matrix(LonLats[city_comb[i, "City2"], c('Longitude', 'Latitude')]))
})
names(intercity_dists) <- city_comb$City12
mle_CIs$geo_dist <- intercity_dists[mle_CIs$City12]

# Add InterCoast and InterPort
mle_CIs$InterCoast <- (mle_CIs$City1 != mle_CIs$City2) & (LonLats[mle_CIs$City1, "CoastAccess"] == 1) & (LonLats[mle_CIs$City2, "CoastAccess"] == 1) 
mle_CIs$InterPort <- (mle_CIs$City1 != mle_CIs$City2) & (LonLats[mle_CIs$City1, "Port"] == 1) & (LonLats[mle_CIs$City2, "Port"] == 1)  

# Check NA returned for missing sids: yes
metadata['agsdhfsdhf', 'COLLECTION.DATE']
any(is.na(metadata$Year)) # Check for missing years

# All the samples with missing dates are recent
range(metadata$Year[is.na(metadata$COLLECTION.DATE)])

# Add dates and years
mle_CIs$date1 = metadata[mle_CIs$individual1, 'COLLECTION.DATE']
mle_CIs$date2 = metadata[mle_CIs$individual2, 'COLLECTION.DATE']
mle_CIs$year1 = metadata[mle_CIs$individual1, 'Year']
mle_CIs$year2 = metadata[mle_CIs$individual2, 'Year']

# Add time dist based on date
mle_CIs$time_dist = abs(difftime(as.Date(mle_CIs$date1), 
                                 as.Date(mle_CIs$date2), units = 'weeks'))

# Add year dist based on yearly label (not rounded to nearest year)
mle_CIs$year_diff <- abs(as.numeric(mle_CIs$year1) - as.numeric(mle_CIs$year2))

# Add sample comp
mle_CIs$sample_comp = apply(mle_CIs[, c("individual1", "individual2")], 1, 
                            function(x) paste(sort(x), collapse = "_"))

# Save data frame
save(mle_CIs, file = sprintf("../../RData/mles_CIs_extended_freqs%s_meta.RData", freqs_used))




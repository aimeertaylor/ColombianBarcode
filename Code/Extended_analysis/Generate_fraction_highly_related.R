###############################################################################
#' Script to compute fraction highly related partioned by space and time
###############################################################################
library(tidymodels)
load("../../RData/mles_CIs_extended_freqsTaylor2020_meta.RData")
load("../../RData/metadata_extended.RData")
source("./summarise_mles.R")
nboot <- 100
eps <- 0.01 # Used to define clones in Taylor et al. 2020

# Filter uninformative
uninformative_ind <- mle_CIs$r2.5. < eps & mle_CIs$r97.5. > 1-eps
summarise_mles(mle_CIs, metadata_ = metadata)
mle_CIs <- mle_CIs[!uninformative_ind,] 
summarise_mles(mle_CIs, metadata_ = metadata)
writeLines(sprintf("Minimum SNP count among remaining sample pairs: %s", min(mle_CIs[, "snp_count"])))
writeLines(sprintf("Minimum SNP count among remaining samples: %s", 
                   min(metadata[unique(c(mle_CIs$individual1, mle_CIs$individual2)), "snp_count"])))


# ============== By year =================
fraction_highly_related_year <- mle_CIs %>%
  group_by(year_diff) %>%
  summarise(fHR = mean(highly_related), 
            npairs = length(highly_related), 
            propCoastTRUE = mean(Coast), 
            propPortTRUE = mean(Port)) %>%
  arrange(year_diff)

# Compute CIs using the bootstrap 
# (use group by instead of argument strata because "strata below 10% of the total are pooled together")
boots <- mle_CIs %>%
  group_by(year_diff) %>%
  bootstraps(times = nboot) 

# Matrix to collect fractions highly related
boots_fHR_matrix = array(dim = c(nrow(fraction_highly_related_year), nboot), 
                         dimnames = list(fraction_highly_related_year$year_diff, NULL))

# Extract fractions highly related into the matrix  
for(i in 1:length(boots$splits)){ # using forloop because boots$splits not all same length
  x <- analysis(boots$splits[[i]])
  booti <- x %>%
    group_by(year_diff) %>%
    summarise(fHR = mean(highly_related))  
  boots_fHR_matrix[as.character(booti$year_diff), i] <- booti$fHR
}

# Compute error bars
Error_bars_year <- t(apply(boots_fHR_matrix, 1, quantile, probs = c(0.025, 0.975), na.rm = T))


#============== By city =================
# Extract the fraction related and number of pairs
fraction_highly_related_city <- mle_CIs %>%
  group_by(City12) %>%
  summarise(fHR = mean(highly_related), 
            intercity_dist = mean(geo_dist), 
            Coast = unique(Coast)) %>%
  arrange(intercity_dist)

# Add minmum per-city sample size among informative relatedness estimates
fraction_highly_related_city$min_n_pair <- sapply(fraction_highly_related_city$City12, function(City12){
  x <- mle_CIs[mle_CIs$City12 == City12,]
  min(table(unique(rbind(as.matrix(x[,c("individual1", "City1")]), 
                         as.matrix(x[,c("individual2", "City2")])))[,"City1"]))})

# Compute CIs using the bootstrap 
boots <- mle_CIs %>%
  group_by(City12) %>%
  bootstraps(times = nboot) 

# Matrix to collect fractions highly related
boots_fHR_matrix = array(dim = c(nrow(fraction_highly_related_city), nboot), 
                         dimnames = list(fraction_highly_related_city$City12, NULL))

# Extract fractions highly related into the matrix  
for(i in 1:length(boots$splits)){ # using forloop because boots$splits not all same length
  x <- analysis(boots$splits[[i]])
  booti <- x %>%
    group_by(City12) %>%
    summarise(fHR = mean(highly_related))  
  boots_fHR_matrix[booti$City12, i] <- booti$fHR
}

# Compute CIs
Error_bars_city <- t(apply(boots_fHR_matrix, 1, quantile, probs = c(0.025, 0.975), na.rm = T))

# Save
save(fraction_highly_related_year, Error_bars_year, 
     fraction_highly_related_city, Error_bars_city, 
     file = "../../RData/fraction_highly_related_extended.RData")
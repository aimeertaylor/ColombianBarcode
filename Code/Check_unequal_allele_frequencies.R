###################################################################
# Check that departure form isolation-by-distance remains (yes) when
# the samples from different cities are treated as being from diff. 
# populations (as FST results suggest) using hmmIBD 
###################################################################
rm(list = ls())
load('../RData/SNPData.RData') 
load('../RData/geo_dist_info_cities.RData')
require(gtools)
RUN_hmmIBD <- FALSE
high_relatedness_threshold <- 0.25

# Import SNP data
data_set = read.delim("../TxtData/hmmInput.txt") 

# Reformat for hmmIBD
data_set[is.na(data_set)] <- -1

# Break into different files per city
cities <- unique(SNPData$City)
for(city in cities){
  samples <- rownames(SNPData)[SNPData$City == city]
  data_city <- data_set[,c('chrom', 'pos', samples)]
  write.table(data_city, file = sprintf("../TxtData/hmmIBD_%s.txt", city), 
              quote = FALSE, row.names = FALSE, sep = "\t")
}

# Create of pairwise city comparisons
city_comparisons <- combinations(n = length(cities), r = 2, v = cities, repeats.allowed = T)

# Run hmmIBD for each city comparison
if(RUN_hmmIBD){
  for(i in 1:nrow(city_comparisons)){
    city_pair <- city_comparisons[i,]
    writeLines(sprintf('\n %s %s',  city_pair[1], city_pair[2]))
    system(sprintf('~/Dropbox/hmmIBD-v2.0.4/hmmIBD -i ../TxtData/hmmIBD_%s.txt -I ../TxtData/hmmIBD_%s.txt -o ../TxtData/%s_%s', 
                   city_pair[1], city_pair[2], city_pair[1], city_pair[2]))
    
    system(sprintf('rm ../TxtData/%s_%s.hmm.txt', city_pair[1], city_pair[2])) # Delete unwanted file
  }
}

# Extract results
fraction_highly_related <- array(NA, dim = nrow(city_comparisons), 
                                 dimnames = list(apply(city_comparisons, 1, 
                                                       function(x){paste(sort(x), collapse = "_")})))
for(i in 1:nrow(city_comparisons)){
  city_pair <- city_comparisons[i,]
  results <- read.delim(sprintf('../TxtData/%s_%s.hmm_fract.txt', city_pair[1], city_pair[2]))
  fraction_highly_related[i] <- mean(results$fract_sites_IBD > high_relatedness_threshold)
}

# Plot results to check departure from isolation by distance still exists: yes 
plot(y = fraction_highly_related, 
     x = geo_dist_info$pairwise_site_distance_all[names(fraction_highly_related)], 
     pch = 20, ylab = '', xlab = 'Distance between cities (km)')
text(y = fraction_highly_related, 
     x = geo_dist_info$pairwise_site_distance_all[names(fraction_highly_related)], 
     labels = names(fraction_highly_related), cex = 0.5, pos = 4)
title(ylab = sprintf('Fraction highly related with r estimate > %s \n not accounting for uncertainty but using site specific frquencies',
                     high_relatedness_threshold), line = 2)
           
# Correlate with results that account for uncertainity 
load('../RData/proportions_geo.RData')
plot(y = proportions_geo['mean', names(fraction_highly_related)], 
     x = fraction_highly_related, pty = 's', 
     xlim = c(0,0.5), ylim = c(0,0.5), pch = 20, 
     ylab = 'With uncertainty using site averaged frquencies', 
     xlab = 'Without uncertainty using site specific frquencies', 
     main = sprintf('Fraction highly related with r estimate > %s',
                    high_relatedness_threshold))
abline(a = 0, b = 1)
text(y = proportions_geo['mean', names(fraction_highly_related)], 
     x = fraction_highly_related, 
     labels = names(fraction_highly_related), cex = 0.5, pos = 4)


#============================================================
# Modify error bars on binomial proportions (use bootstrap given not independent)
# Could use Lyndals test for signal of selection (or at least cite)
#============================================================
rm(list = ls())
library(RColorBrewer)
load('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/ViterbiPaths.RData')
load('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/Result.RData') 
chrom_lengths <- read.csv(file = '~/Documents/BroadLaptop/hmmIBD_benchmark/appendix3/pf3k_data/Pf_v3_chrom_length.csv', skip = 2, header = FALSE)
site_comps = unique(Result$site_comp) 
names_snps = colnames(ViterbiPaths)
no_samp_comps = nrow(ViterbiPaths)
no_snps = ncol(ViterbiPaths)
chr = as.numeric(do.call(rbind, strsplit(names_snps, split = '_'))[,1]) # Extract chr info
pos = as.numeric(do.call(rbind, strsplit(names_snps, split = '_'))[,2]) # Extract pos info
VP = abs(ViterbiPaths - 1) # Recode IBD as 1 and nIBD as 0 
VP_col = VP * rep(chr, each = no_samp_comps) # Recode as IBD 1:14 (depending on chr) and nIBD as 0
cumsum_chr = cumsum(c(0,chrom_lengths$V3))
x_coord = cumsum_chr[chr] + pos # Calculate X coordinates for plots along genome
x_chr_mid = colMeans(rbind(cumsum_chr[-length(cumsum_chr)], cumsum_chr[-1]))
  
#============================================================
# Signatures of clonal expansion: proportion of viterbi IBD across the 250 SNPs
#============================================================
par(mfrow = c(4,4), family = 'serif')
# Across all comparisons
plot(y = colMeans(VP), x = x_coord, type = 'l', cex.lab = 0.7, ylim = c(0,0.6),  
     main = 'All pairwise comparisons', xlab = 'SNP position on chromosome', xaxt = 'n', 
     ylab = sprintf('Propotion Viterbi IBD (%s comparisons)', no_samp_comps)) 
abline(v = cumsum_chr, lty = 'dotted')
axis(side = 1, at = x_chr_mid, labels = 1:14, cex.axis = 0.5, las = 2)

# Add CIs
half_CI = sqrt((-log(0.025))/(2*no_samp_comps)) # Check
segments(y0 = colMeans(VP)-half_CI, 
         y1 = colMeans(VP)+half_CI, 
         x0 = x_coord,  x1 = x_coord, lwd = 0.2)

for(site_comp in site_comps){
  sample_comps = as.character(Result$sample_comp[Result$site_comp == site_comp])
  n = length(sample_comps) 
  plot(y = colMeans(VP[sample_comps,]), x = x_coord, type = 'l', cex.lab = 0.7, ylim = c(0,0.6),  
       main = site_comp, xlab = 'SNP position on chromosome', xaxt = 'n', 
       ylab = sprintf('Propotion Viterbi IBD (%s comparisons)', n)) 
  
  # Add 95% CIs
  half_CI = sqrt((-log(0.025))/(2*n)) # Check
  segments(y0 = colMeans(VP[sample_comps,])-half_CI, 
           y1 = colMeans(VP[sample_comps,])+half_CI, 
           x0 = x_coord,  x1 = x_coord, lwd = 0.2)
  
  abline(v = cumsum_chr, lty = 'dotted')
  axis(side = 1, at = x_chr_mid, labels = 1:14, cex.axis = 0.5, las = 2)
}


#============================================================
# Heatmaps of Viterbi paths 
# Takes a long time (just plot Buenaventua and Tamaco with some IBD)
#============================================================
par(mfrow = c(2,1), family = 'serif', mar = c(4,4,1,1))
cols_chr = c('#FFFFFF', rainbow(14)) # To colour by chromosome

for(site_comp in "Buenaventura_Tumaco"){ # "Buenaventura_Buenaventura", "Tumaco_Tumaco"
  sample_comps = as.character(Result$sample_comp[Result$site_comp == site_comp])
  n = length(sample_comps) 
  plot(y = colMeans(VP[sample_comps,]), x = x_coord, type = 'l', cex.lab = 0.7, ylim = c(0,0.2),  
       xlab = 'SNP position on chromosome', xaxt = 'n', 
       ylab = sprintf('Propotion Viterbi IBD (%s comparisons)', n)) 
  
  # Add 95% CIs
  half_CI = sqrt((-log(0.025))/(2*n)) # Check
  segments(y0 = colMeans(VP[sample_comps,])-half_CI, 
           y1 = colMeans(VP[sample_comps,])+half_CI, 
           x0 = x_coord,  x1 = x_coord, lwd = 0.3, 
           col = cols_chr[chr+1])
  
  abline(v = cumsum_chr, lty = 'dotted')
  for(i in 1:14){
    axis(side = 1, at = x_chr_mid[i], labels = i, cex.axis = 1, col.axis = cols_chr[i+1])
  }
   
  # Restrict to comparisons with one or more IBD segments
  sample_comps_IBD = sample_comps[rowSums(VP[sample_comps,]) > 0]
  
  # # Uncomment to order by time dist
  # time_dists = Result$time_dist[Result$sample_comp %in% sample_comps_IBD]
  # new_order = sort.int(time_dists, index.return = T)$ix
  
  # Uncomment to order by genomic fraction IBD
  IBD_dists = Result$IBD[Result$sample_comp %in% sample_comps_IBD]
  new_order = sort.int(IBD_dists, index.return = T)$ix
  
  # Colour by chr
  image(t(VP_col[sample_comps_IBD,][new_order,]), xlab = '', ylab = '',
        xaxt = 'n', yaxt = 'n', col = cols_chr)
  
  title(line = 1, xlab = '250 SNPs', 
        ylab = sprintf('%s pairwise comparisons with one or more IBD segments', 
                                          length(sample_comps_IBD)))
}


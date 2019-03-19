############################################################################
# In this script I intend to extract the Viterbi paths for all comparisons
# Takes 1306 sec to run 
# Seems to a problem with X$Nsnp[i] - tell steve it was solved with hmmIBD-2.0.2
############################################################################

#============================================
# Set up
#============================================
library(dplyr) # For filter 
library(tictoc)

tic()

# Load IBD and IBS results
load('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/Result.RData') 
# Import segment data
Segments <- read.delim('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/TxtData/hmmIBD_output.hmm.txt')

#============================================
# Data re-formating and checks
#============================================
# Prefix chr to avoid potentially problems with non-unique pos and add comparison
Segments$SNPs_str = paste(Segments$chr, Segments$start, sep = "_")
Segments$SNPs_end = paste(Segments$chr, Segments$end, sep = "_")
Segments$sample_comp = apply(Segments[,c('sample1', 'sample2')], 1, function(x){paste(sort(x), collapse = '', sep = '')})

# Check naming of segments covers all sample comparisons
all(Segments$sample_comp %in% Result$sample_comp) & all(Result$sample_comp %in% Segments$sample_comp) 

# Check start and end positions are only those provided (seems so - obvious)
HMMdata = read.delim('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/TxtData/hmmIBD_input.txt')
all(HMMdata$pos %in% X) # FALSE
all(X  %in% HMMdata$pos) # TRUE
length(HMMdata$pos) # 250

#============================================
# Convert into long format
#============================================
comps = unique(Segments$sample_comp)
SNPs = paste(HMMdata$chr, HMMdata$pos, sep = '_') 
no_comp = length(comps)
no_SNPs = length(SNPs)
ViterbiPaths = array(dim = c(no_comp, no_SNPs), dimnames = list(comps, SNPs))

# Relys on correct ordering of SNPs in Segments and HMM
pb <- txtProgressBar(min = 0, max = no_comp, style = 3)

for(comp in comps){
  setTxtProgressBar(pb, which(comps == comp))
  comp_ind = Segments$sample_comp == comp
  X = Segments[comp_ind, ]
  
  # Check 1: SNPs are in order
  Z = sapply(X$chr, function(x){
    chr_ind = (X$chr == x)
    str_inc = all(X[chr_ind, 'start'] == cummax(X[chr_ind, 'start']))
    end_inc = all(X[chr_ind, 'end'] == cummax(X[chr_ind, 'end']))
    if(!all(str_inc, end_inc)){stop('SNPs not in monotonically increasing order')}
    })
  
  for(i in 1:nrow(X)){
    COLs = which(SNPs == X$SNPs_str[i]):which(SNPs == X$SNPs_end[i])
    # Check 2: expected number of SNPs
    if(length(COLs) != X$Nsnp[i])stop('Number of SNPs do not agree')
    
    ViterbiPaths[comp,COLs] = X$different[i] # 1 = difference = nIBD; 0 = IBD  
  }
}

# Save (1 = difference = nIBD; 0 = IBD)
save(ViterbiPaths, file = '/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/ViterbiPaths.RData')
toc()


#==============================================================================================
# Function to calculate Wier and Cockerham's theta per SNP and theta averaged over multiple SNPs
# Assumes biallelic and pairwise comparisons
# Note that when n_bar >> n1, n_c goes to zero and theta goes to +/- inf
#==============================================================================================
calculate_pairwise_WandC <- function(P1, P2, SNPDataBinary){ 
  
  r <- 2 # Assumes pairwise comparisons 
  numSNPs <- ncol(SNPDataBinary) 
  
  # Sample size per population (accounting for different numbers of missing data per SNP, assumes missing N)
  sample_sizes_per_subpopulation <- matrix(nrow = 2, ncol = numSNPs, dimnames = list(c('P1', 'P2'), NULL))
  sample_sizes_per_subpopulation['P1',] <- apply(SNPDataBinary[P1,], 2, FUN = function(x){sum(!is.na(x))})
  sample_sizes_per_subpopulation['P2',] <- apply(SNPDataBinary[P2,], 2, FUN = function(x){sum(!is.na(x))})
  
  # Are any sample sizes < 4? 
  snps_with_too_few_data <- as.logical(colSums(sample_sizes_per_subpopulation < 4))
  
  # Average sample size and n_c
  n_bar <- apply(sample_sizes_per_subpopulation, 2, mean)
  n_c <- (r*n_bar - apply((sample_sizes_per_subpopulation^2)/(r*n_bar), 2, sum))/(r-1)
  
  # Reference allele frequencies
  X <- rbind(colMeans(SNPDataBinary[P1,], na.rm = TRUE), colMeans(SNPDataBinary[P2,], na.rm = TRUE))
  
  # Are there any non-variant X? snps_with_too_few_data return NA
  invariant_snps <- (colSums(X) == 0 | colSums(X) == 1 | colSums(X) == 2)
  
  # Weighted mean of each allele frequency over subpopulations 
  # Allow SNPs with no data at one or both sites to return NaN values
  p_bar <- colSums(X*sample_sizes_per_subpopulation)/(r*n_bar)
  
  # Weighted sample variance of each allele frequency over subpopulations
  sq_diff <-  apply(X, 1, FUN = function(x){(x-p_bar)^2})
  s2 <- colSums(t(sq_diff)*sample_sizes_per_subpopulation)/(n_bar*(r-1))
  
  # Average heterozygosity ()
  hi <- (1- X^2 - (1-X)^2)
  h_bar <- colSums(hi * sample_sizes_per_subpopulation)/(r*n_bar)
  
  # Calculate a, b, c (Weir and Cockerham 1984)
  a = (n_bar/n_c) * (s2 - (1/(n_bar -1))*((p_bar*(1-p_bar)) - (((r-1)/r)*s2) - (h_bar/4)))
  b = (n_bar/(n_bar -1))*((p_bar*(1 - p_bar)) - (((r-1)/r) * s2) - (((2 * n_bar) - 1)/(4*n_bar))*h_bar)
  c = 0.5*h_bar
  
  NAsnps <- (invariant_snps | snps_with_too_few_data | n_c == 0) 
  
  # Calculate theta = Fst (Weir and Cockerham 1984)
  # Remove SNPs that have no data in one or both populations
  theta_hat = a/(a + b + c)
  theta_hat[NAsnps] <- NA
  theta_loci_hat_combined = sum(a[!NAsnps])/sum(a[!NAsnps] + b[!NAsnps] + c[!NAsnps])
  
  # ---------------------- End of function ----------------------
  return(list(theta_hat = theta_hat,
              theta_loci_hat_combined = theta_loci_hat_combined,
              #variance_theta_loci_hat_combined = variance_theta_loci_hat_combined,
              #Fst_miotto2015 = Fst_miotto2015,
              n_bar = n_bar, 
              n_c = n_c))
}


#==============================================================================================
# Function to calculate Fst using Hudsons estimator
# Assumes biallelic and pairwise comparisons
#==============================================================================================
calculate_pairwise_reich <- function(P1, P2, SNPDataBinary){ 
  
  # Sample size per population (accounting for different numbers of missing data per SNP, assumes missing N)
  n1 <- apply(SNPDataBinary[P1,], 2, FUN = function(x){sum(!is.na(x))})
  n2 <- apply(SNPDataBinary[P2,], 2, FUN = function(x){sum(!is.na(x))})
  
  sqrd_term <- (colMeans(SNPDataBinary[P1,], na.rm = TRUE) - colMeans(SNPDataBinary[P2,], na.rm = TRUE))^2
  a1 <- colSums(SNPDataBinary[P1,], na.rm = TRUE)
  a2 <- colSums(SNPDataBinary[P2,], na.rm = TRUE)
  h1 <- (a1*(n1-a1))/(n1*(n1-1))
  h2 <- (a2*(n2-a2))/(n2*(n2-1))
  N <- sqrd_term - (h1/n1) - (h2/n2)
  D <- N + h1 + h2
  
  # Snps to remove
  snps_with_too_few_data <- as.logical(colSums(rbind(n1,n2) < 4)) # Are any sample sizes < 4? 
  invariant_snps <- (a1 == n1 | a2 == n2) # Are there any non-variant X? snps_with_too_few_data return NA
  NAsnps <- (invariant_snps | snps_with_too_few_data) 
  
  F_st_snp <- N/D 
  F_st <- sum(N[!NAsnps])/sum(D[!NAsnps])
  
  # ---------------------- End of function ----------------------
  return(list(F_st_snp = F_st_snp,
              F_st = F_st))
}



#==============================================================================================
# Function to calculate Fst based on inferred frequncies 
# Using Wright's definition of Fst as written in appendix of Reich et al. 2009
#==============================================================================================
calculate_pairwise_inferred <- function(P1_name, P2_name, bootstrap_snps){ 
  
  # Ad-hoc fix for distingusing between site and year ie, MKK:06 and just site ie 'MKK
  if(grepl(':',P1_name)){
    # Load inferred frequencies
    load('/Users/aimeet/Documents/BroadLaptop/TM_border/RData/COIext_f_extracted_site_year.RData')  
  } else {
    load('/Users/aimeet/Documents/BroadLaptop/TM_border/RData/COIext_f_extracted_site.RData')  
  }
  
    # Ad-hoc fix due to changing PLU to MKK
    dimnames(frequency_quants)[[3]] <- sub('PLU', 'MKK', dimnames(frequency_quants)[[3]],)
    
    # Calculate reference allele frequencies
    p1 <- frequency_quants[bootstrap_snps,"50%", P1_name, drop = FALSE]
    p2 <- frequency_quants[bootstrap_snps,"50%", P2_name, drop = FALSE]
    q1 <- 1 - p1
    q2 <- 1 - p2
    
    N <- p1*(q2-q1) + p2*(q1 - q2)
    D <- p1*q2 + q1 * p2
    Fst_snp = N/D
    Fst <- mean(Fst_snp)
    return(list(Fst_snp = Fst_snp, Fst = Fst))
}




#==============================================================================================
# Function to calculate Weir and Hill 
#==============================================================================================
calculate_pairwise_WandH <- function(P1, P2, SNPDataBinary){

  # Sample size per population (accounting for different numbers of missing data per SNP, assumes missing N)
  ni <- apply(SNPDataBinary[P1,], 2, FUN = function(x){sum(!is.na(x))})
  nj <- apply(SNPDataBinary[P2,], 2, FUN = function(x){sum(!is.na(x))})

  nai_1 <- colSums(SNPDataBinary[P1,], na.rm = TRUE)
  naj_1 <- colSums(SNPDataBinary[P2,], na.rm = TRUE)

  nai_0 <- ni - nai_1
  naj_0 <- nj - naj_1
  p_0 = nai_0 / ni
  p_1 = naj_0 / nj
  pmean = (ni * p_0 + nj * p_1) / (ni + nj)
  nic = ni - ni * ni / (ni + nj)
  njc = nj - nj * nj / (ni + nj)
  nc = nic + njc
  msp = ni * (p_0 - pmean) * (p_0 - pmean) + nj * (p_1 - pmean) * (p_1 - pmean)
  msg = (ni * p_0 * (1 - p_0) + nj * p_1 * (1 - p_1)) / (ni - 1 + nj - 1)
  num = msp - msg;
  denom = msp + (nc - 1) * msg

  fstnum = sum(num, na.rm = TRUE)
  fstdenom = sum(denom, na.rm = TRUE);
  fst = fstnum / fstdenom

  # ---------------------- End of function ----------------------
  return(fst)
}
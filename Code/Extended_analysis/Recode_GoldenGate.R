##################################################################################
#' Formatting the Golden Gate Colombian data for analyses inc. more contemporary
#' WGS data from Guapi (Vladimir and co.) and Ecuador (Fabian and co.). This
#' script is not automated: go through each line manually
#'
#' Coding used in original data formating script; see Code/DataFormat.R:
#' AA = 1, BB = 0, AB = NA, -- = NA.
#'
#' Where AA, BB etc. is the GoldenGate code whose mapping onto 3D7 ref/alt (as
#' used by Manuela to process the WGS Guapi data; see email dated Feb 10th) is
#' non trivial. Nonetheless, Diego figured out the mapping (see notes below
#' copied from Diego's email dated via email dated April 9th 20 within the
#' thread whose subject is Quick question about funders and Colombia preprint)
#' and shared it via google drive (see email dated April 9th 2020 whose subject
#' is "revealing SNPs Goldengate April2020.xlsx")
#'
#' Please see in your email an access via google-drive to get an excel doc with
#' some sheets. I figured out the nucleotides for each SNP, in other words, we
#' now know what does mean AA or BB in each allele. to do so I combined
#' information such as:
#' 1) The order that Tim sent to Illumina in order to build the chip-array (it
#' contain the flanking regions and SNPs)
#' 2) The output information of a whole plate after the genotyping. With this, I
#' was able to know what is AA or BB and the equivalent SNPs
#' 3) As Tim mentioned to me, some SNPs were identified in the 5'-->3 or 3'-->5
#' strains (that was the very difficult part)  but as we had the reference
#' strains (specially 3D7) I figured out using PlasmoDB and the supplementary
#' figure 2 of my paper
#' 4) A lot of patients, I checked the files 3X, and in different ways in order to
#' avoid painful mistakes in further analyses. to the very right of the last
#' sheet "final table," you will see the most important table, which contains
#' the meaning of each AA or BB per allele and the respective nucleotide. You
#' may find previous files kind of busy and probably difficult to follow, I will
#' be happy to try to explain what I did.
##################################################################################
rm(list = ls())

# ===================================================
# Load mapping (see comments above for provenance)
# ===================================================
X <- read.csv(file = '../../OriginalData/OriginalDataGoldenGate/revealing SNPs Goldengate April2020_final_table.csv', 
              skip = 7, nrows = 250) 

# Check all rows imported and no more
# Last SNP imported should be MAL14-3126540
tail(X)[1:5]

# Inspect original excel sheet and colnames by hand and extract the relevant columns
colnames(X)

# [1] "chrom"                                         
# [2] "pos"  
# [3] "ALL.SNPs.from.Goldengate..previous.position."
# [16] "X3D7.ref.Alele.according.plasmoDB"  
# [20] "chrom.1"                                       
# [21] "pos.1"                                         
# [22] "Goldengate.Code.for.the.SNP"                   
# [23] "nucleotide..based.on.3D7."                     
# [24] "Goldengate.Code.for.the.SNP.1"                 
# [25] "nucleotide..based.on.3D7.and.other.ref.samples"

# Check chrom and pos match: yes
all(X[,"chrom"] == X[,"chrom.1"], na.rm = T)
all(X[,"pos"] == X[,"pos.1"], na.rm = T)

# Cherry pick columns: first check
head(X[, c(3,1,2,22:25,16)])

GoldenGate_to_3D7_map <- data.frame('GoldenGateSNP' = as.character(X[,3]), 
                                    'chrom' = X[,1], 
                                    'pos' = X[,2], 
                                    'Code_GoldenGate_1' = as.character(X[,22]), 
                                    'Nucleotide_3D7_1' = as.character(X[,23]), 
                                    'Code_GoldenGate_2' = as.character(X[,24]), 
                                    'Nucleotide_3D7_2' = as.character(X[,25]),
                                    'Ref_nucleotide_3D7' = as.character(X[,16]), 
                                    stringsAsFactors = FALSE)

GoldenGate_to_3D7_map$Ref_code_3D7 <- 
  apply(GoldenGate_to_3D7_map, 1, function(x){
    if (x['Ref_nucleotide_3D7'] == x['Nucleotide_3D7_1']) {
      return(x['Code_GoldenGate_1'])
    } else if (x['Ref_nucleotide_3D7'] == x['Nucleotide_3D7_2']) {
      return(x['Code_GoldenGate_2']) 
    } else {
      stop("3D7 reference nucleotide matches neither")
    }
  })


# ===================================================
# Load data used in my analysis that features in bioRxiv preprint entitled
# "Identity-by-descent relatedness estimates with uncertainty characterise
# departure from isolation-by-distance between Plasmodium falciparum
# populations on the Colombian-Pacific coast" and convert back to AA and BB
# ===================================================
hmmInput <- read.delim(file = '../../TxtData/hmmInput.txt')
load(file = '../../RData/SNPData.RData')

# Calculate frequencies before converting 
SNP_inds <- grepl("MAL", colnames(SNPData))
colnames(SNPData)[SNP_inds] # Check
SNPData_freq <- colMeans(SNPData[, SNP_inds], na.rm = TRUE)
hmmInput_freq <- rowMeans(hmmInput[,-(1:2)], na.rm = TRUE)

# Convert SNPData 
SNPData[, SNP_inds][SNPData[, SNP_inds] == 1] <- "AA" 
SNPData[, SNP_inds][SNPData[, SNP_inds] == 0] <- "BB" 
rbind(head(SNPData), tail(SNPData))

# Convert hmmInput
hmmInput[,-(1:2)][hmmInput[,-(1:2)] == 1] <- "AA"
hmmInput[,-(1:2)][hmmInput[,-(1:2)] == 0] <- "BB"
rbind(head(hmmInput), tail(hmmInput))

# Recode SNPData to 0 = ref and 1 = alt
SNPDataRecode <- SNPData
SNPDataRecode[,SNP_inds] <- NA 
for(SNP in colnames(SNPData)[SNP_inds]){
  SNP_ind <- which(GoldenGate_to_3D7_map$GoldenGateSNP == SNP)
  SNP_ref_code <- GoldenGate_to_3D7_map$Ref_code_3D7[SNP_ind]
  SNP_ref_ind <- which(SNPData[, SNP] == SNP_ref_code)
  SNP_alt_ind <- which(SNPData[, SNP] != SNP_ref_code)
  SNPDataRecode[SNP_ref_ind, SNP] <- 0
  SNPDataRecode[SNP_alt_ind, SNP] <- 1
}

# Recode hmmInput to 0 = ref and 1 = alt
hmmInputRecode <- hmmInput
hmmInputRecode[,-(1:2)] <- NA 
for(i in 1:nrow(hmmInput)){
  SNP_ind <- which(GoldenGate_to_3D7_map$chrom == hmmInput$chrom[i] & 
                     GoldenGate_to_3D7_map$pos == hmmInput$pos[i])
  
  if (length(SNP_ind) != 1) {
    stop("SNP ind is not unique") } else {
      SNP_ref_code <- GoldenGate_to_3D7_map$Ref_code_3D7[SNP_ind]
      SNP_ref_ind <- which(hmmInput[i,-(1:2)] == SNP_ref_code)
      SNP_alt_ind <- which(hmmInput[i,-(1:2)] != SNP_ref_code)
      hmmInputRecode[i,-(1:2)][SNP_ref_ind] <- 0
      hmmInputRecode[i,-(1:2)][SNP_alt_ind] <- 1
    }
}

# Check numeric
str(SNPDataRecode)
str(hmmInputRecode)

# Check frequencies are unchanged (f) or 1-f
SNPDataRecode_freq <- colMeans(SNPDataRecode[, SNP_inds], na.rm = TRUE)
hmmInputRecode_freq <- rowMeans(hmmInputRecode[,-(1:2)], na.rm = TRUE)

plot(SNPData_freq, SNPDataRecode_freq)
abline(a = 0, b = 1, col = 'blue'); abline(a = 1, b = -1, col = 'blue')
plot(hmmInput_freq, hmmInputRecode_freq)
abline(a = 0, b = 1, col = 'blue'); abline(a = 1, b = -1, col = 'blue')

# Save Recoded data in format for HMM and RData
write.table(hmmInputRecode, file = '../../TxtData/hmmInputRecode.txt', 
            quote = FALSE, row.names = FALSE, sep = '\t')
save(SNPDataRecode, file = '../../RData/SNPDataRecode.RData')



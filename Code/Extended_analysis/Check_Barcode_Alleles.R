##################################################################################
#' Mon 7th Dec 2020: I belive this script is obsolete since Angela re-mapped
#' data while I was on maternity leave 
#' (see email thread entitled "Golden Gate Recoded Data" Started 27th July)
#' 
#' Script to check the Barcode alleles that Angela extracts for the GoldenGate
#' barcode against those Diego reports in his recoding. 
#' 
#' Notes from email sent by Angela on 29th April 2020 in the thread entitled
#' WGS-to-barcode data for extended Colombia-Ecuador analysis" 
#' 
#' Here are the genotypes for Ecuador. (I have them in a vcf with the Colombia
#' samples so both are actually in this matrix.) I've coded heterozygous calls
#' with 'H' and missing calls as '-1'. I'm also attaching a list of countries,
#' including the 2 we mentioned today that have a probable origin in Venezuela.
#' In the third file, I've put the alleles that we have for each position in the
#' vcf. I thought matching this to the barcode data could be a first quick
#' sanity check. (Happy to do this part myself if you send over the data!)
##################################################################################
rm(list = ls())
library(tidyverse)

# Load barcode alleles provided by Angela
Barcode_alleles <- read.delim(file = '../../OriginalData/Diego_barcode_alleles.txt') 
Barcode_alleles$X.CHROM
Barcode_alleles$POS

# Load file used to map GoldenGate code to 3D7 v3 alleles: 
load(file = '../../RData/GoldenGate_to_3D7_map.RData')

GoldenGate_to_3D7_map$chrom
GoldenGate_to_3D7_map$pos
GoldenGate_to_3D7_map$Ref_nucleotide_3D7

# Check is pos are unqiue: yes
length(GoldenGate_to_3D7_map$pos) == length(unique(GoldenGate_to_3D7_map$pos))
length(Barcode_alleles$POS) == length(unique(Barcode_alleles$POS))

# Check pos match
all(Barcode_alleles$POS %in% GoldenGate_to_3D7_map$pos) # All those from Angela are in GoldenGate_to_3D7_map$pos
sum(!GoldenGate_to_3D7_map$pos %in% Barcode_alleles$POS) # Two missing from Barcode_alleles$POS

# Which positions are missing? 
GoldenGate_to_3D7_map[!GoldenGate_to_3D7_map$pos %in% Barcode_alleles$POS,]

# Add the Alt nucleotide
GoldenGate_to_3D7_map$Alt_nucleotide_3D7 <- sapply(1:nrow(GoldenGate_to_3D7_map), function(i){
  if(GoldenGate_to_3D7_map$Ref_nucleotide_3D7[i] == GoldenGate_to_3D7_map$Nucleotide_3D7_1[i]){
    return(GoldenGate_to_3D7_map$Nucleotide_3D7_2[i])
  } else if (GoldenGate_to_3D7_map$Ref_nucleotide_3D7[i] == GoldenGate_to_3D7_map$Nucleotide_3D7_2[i]) {
    return(GoldenGate_to_3D7_map$Nucleotide_3D7_1[i])
  } else {
    stop("Something wrong")
  }
})

# Rename 
Joined <- Barcode_alleles %>% 
  as.tibble() %>%
  rename(pos = POS) %>%
  full_join(GoldenGate_to_3D7_map, by = "pos")

# All ref match
all(Joined$REF == Joined$Ref_nucleotide_3D7, na.rm = TRUE)

# Not all alt match
all(Joined$ALT == Joined$Alt_nucleotide_3D7, na.rm = TRUE)

# Inspect mismatches (seems to be a problem with Barcode_allleles, not GoldenGate_to_3D7_map) 
Joined %>% 
filter(ALT != Alt_nucleotide_3D7)

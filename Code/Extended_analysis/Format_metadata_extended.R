rm(list = ls())

# Formatted version of agreed upon extended SNP data
load('../../RData/snpdata_extended.RData')

# Meta data used in Taylor et al. 2020 (Plos Genetics)
load('../../RData/SNPData.RData')

# Meta data on Guapi samples: 
metadata <- read.csv("../../recodedgoldengatedata/Guapi_metadata_forAimee.csv")

# Other meta data
metadata_extra <- read.csv("../../recodedgoldengatedata/Colombia-Ecuador_Barcode_Clusters_Missing40_Master.csv")

# Check for overlap
table(metadata$SAMPLE.CODE %in% SNPData$SAMPLE.CODE) # No overlap
table(metadata$SAMPLE.CODE %in% metadata_extra$Sample) # Partial overlap
table(SNPData$SAMPLE.CODE %in% metadata_extra$Sample) # Partial overlap

# Check metadata available for all samples: no
sids <- colnames(snpdata[,-c(1:2)])
sids_missing_from_two <- sids[!sids %in% c(metadata$SAMPLE.CODE, 
                                  as.character(SNPData$SAMPLE.CODE))] # Some missing still 
sids_missing_from_three <- sids_missing_from_two[!sids_missing_from_two %in% metadata_extra$Sample]

# Check: yes 
setequal(sids_missing_from_three, sids[!sids %in% c(metadata$SAMPLE.CODE, 
                  as.character(SNPData$SAMPLE.CODE), 
                  metadata_extra$Sample)])

write.csv(sids_missing_from_two, file = "../../recodedgoldengatedata/sids_missing_v1.csv")
write.csv(sids_missing_from_three, file = "../../recodedgoldengatedata/sids_missing_v2.csv")

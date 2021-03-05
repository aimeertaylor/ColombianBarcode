#################################################################
#' In this script I load the metadata for the samples analysed 
#' in Taylor et al. 2020 and join to additional metadata for samples
#' analysed in upcoming Carrasquilla et al. 2021. 
#' These additional metadata were provided by Manuela and Angela after
#' several rounds of troubleshooting, inc. a problem caused by 
#' a join that introduced some inconsistent city names for five 
#' Ecuadorian samples. This problem is fixed by hand below. As such, 
#' this script needs to be run line-by-line
#################################################################
# Set wrking data to source file location
rm(list = ls())

# Meta data on new samples: 
metadata <- read.csv("../../recodedgoldengatedata/Guapi_metadata_forAimee_v2.csv", 
                     check.names = FALSE)

# Meta data used in Taylor et al. 2020 (Plos Genetics)
load('../../RData/SNPData.RData')
# Convert factors to classes
for(j in 1:ncol(SNPData)){
  if(is.factor(SNPData[, j])){
    SNPData[, j] <- as.character(SNPData[, j])
  }
}

# Check if any samples in Taylor et al. are listed in metadata and vice versa
any(SNPData$SAMPLE.CODE %in% metadata$SAMPLE.CODE)
any(metadata$SAMPLE.CODE %in% SNPData$SAMPLE.CODE)

# Get rid of spaces in nouns
metadata$SAMPLE.CODE <- gsub(" ", "", metadata$SAMPLE.CODE) 
metadata$STATE <- gsub(" ", "", metadata$STATE) 
metadata$City <- gsub(" ", "", metadata$City) 

# Check for spelling inconsistencies and homogenise 
unique(metadata$STATE)
unique(metadata$City) # Problem with "Sucumbios" and "Sucumbíos" 
# metadata$City[metadata$City == "Sucumbos"] <- "Sucumbíos" # This doesn't work for some reason, must be a accent problem
metadata$City[grepl("Sucumb", metadata$City)] <- "Sucumbios"
unique(metadata$City) # Check problem resolved

# Check unique sample per row: no!
nrow(metadata) == length(unique(metadata$SAMPLE.CODE))

# Remove duplicate rows 
nrow(metadata)
metadata <- dplyr::distinct(metadata) # Removes 14 rows. 
nrow(metadata)

# Check unique sample per row: still no
nrow(metadata) == length(unique(metadata$SAMPLE.CODE))

# Pull out non-unique rows to check for differences
rows_per_sid <- table(metadata$SAMPLE.CODE)
inds_not_unique <- metadata$SAMPLE.CODE %in% names(rows_per_sid[rows_per_sid > 1])

# Send back to Angela/Manuela for checking: 
problem_metadata <- dplyr::arrange(metadata[inds_not_unique, ], SAMPLE.CODE)
write.csv(problem_metadata, file = '../../recodedgoldengatedata/problem_meta.csv')

# Manuela checked and there was some error introduced in a merge
# Let's correct by hand: 
metadata <- metadata[!inds_not_unique, ] # Remove the problem data
# Load the corrected problem metadata 
problem_metadata_corrected <- read.csv(file = '../../recodedgoldengatedata/problem_meta_corrected.csv', 
                                       row.names = 1)
problem_metadata_corrected$COLLECTION.DATE <- as.character(problem_metadata_corrected$COLLECTION.DATE)

# Join the corrected problem metadata to the metadata
metadata <- dplyr::full_join(problem_metadata_corrected, metadata)

# Re check metadata: one row per sample :-)
nrow(metadata) == length(unique(metadata$SAMPLE.CODE))

#===========================================================
# Sort COLLECTION.DATE
#===========================================================
# First sort COLLECTION.DATE
class(SNPData$COLLECTION.DATE)
class(metadata$COLLECTION.DATE)

# Separate Year and Collection date for metadata
metadata$Year <- metadata$COLLECTION.DATE
metadata$Year <- sapply(strsplit(metadata$Year, split = "-"), function(x) x[1])

# Convert metadata Collection date into a date
COLLECTION.DATE_ <- as.Date(rep(NA, nrow(metadata)), format = "%Y-%m-%d")
for(i in 1:length(COLLECTION.DATE_)){
  COLLECTION.DATE_[i] <- as.Date(metadata$COLLECTION.DATE[i], format = "%Y-%m-%d")
}

# Check same where match: yes (only non-matching are NAs)
COLLECTION.DATE_[!as.character(COLLECTION.DATE_) %in% as.character(metadata$COLLECTION.DATE)]
unique(as.character(COLLECTION.DATE_) == as.character(metadata$COLLECTION.DATE))

# Replace
metadata$COLLECTION.DATE <- COLLECTION.DATE_

# Join two data sets
metadata_extended <- dplyr::full_join(SNPData[,c(1:5, ncol(SNPData))], metadata)

# Check looks reasonable: yes :-) 
str(metadata_extended) # 526 obs. of  8 variables
lapply(metadata_extended, function(x) unique(x))
head(metadata_extended)

# Add a column to label data origin
metadata_extended$PloSGen2020 <- metadata_extended$SAMPLE.CODE %in% SNPData$SAMPLE.CODE

# Add rownames
rownames(metadata_extended) <- metadata_extended$SAMPLE.CODE

# =============== Save the extended set of metadata  ===============
metadata <- metadata_extended
save(metadata, file = "../../RData/metadata_extended.RData")



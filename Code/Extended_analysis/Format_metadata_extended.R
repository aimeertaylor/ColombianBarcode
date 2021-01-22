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

# Formatted version of agreed upon extended SNP data
load('../../RData/snpdata_extended.RData')

# Meta data used in Taylor et al. 2020 (Plos Genetics)
load('../../RData/SNPData.RData')
# Convert factors to classes
for(j in 1:ncol(SNPData)){
  if(is.factor(SNPData[, j])){
    SNPData[, j] <- as.character(SNPData[, j])
  }
}

# Meta data on Guapi samples: 
metadata <- read.csv("../../recodedgoldengatedata/Guapi_metadata_forAimee_v2.csv", 
                     check.names = FALSE)

# Get rid of spaces in nouns
metadata$SAMPLE.CODE <- gsub(" ", "", metadata$SAMPLE.CODE) 
metadata$STATE <- gsub(" ", "", metadata$STATE) 
metadata$City <- gsub(" ", "", metadata$City) 

# Check for spelling inconsistencies
unique(metadata$STATE)
unique(metadata$City) # Problem with "Sucumbios" and "Sucumbíos" 
# Homogenise 
# metadata$City[metadata$City == "Sucumbos"] <- "Sucumbíos" # This doesn't work for some reason, must be a accent problem
metadata$City[grepl("Sucumb", metadata$City)] <- "Sucumbios"
unique(metadata$City) # Problem with "Sucumbios" and "Sucumbíos" 

# Check unique sample per row: no!
nrow(metadata) == length(unique(metadata$SAMPLE.CODE))

# Remove duplicate rows 
metadata <- dplyr::distinct(metadata) # Removes 40 rows. 

# Check unique sample per row: still no
nrow(metadata) == length(unique(metadata$SAMPLE.CODE))

# Pull out non-unique rows to check for differences
rows_per_sid <- table(metadata$SAMPLE.CODE)
inds_not_unique <- metadata$SAMPLE.CODE %in% names(rows_per_sid[rows_per_sid > 1])

# Send back to Angela/Manuela for checking: 
problem_metadata <- metadata[inds_not_unique, ]
problem_metadata <- dplyr::arrange(problem_metadata, SAMPLE.CODE)
write.csv(problem_metadata, file = '../../recodedgoldengatedata/problem_meta.csv')

# Manuela checked an there was some error introduced in a merge
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

# ==================================================================
#' #' Aside, let's see what the SNP data look like for these samples Summary:
#' #' Although these five samples have poor coverage, they're not among the worst covered.
#' #' Moreover, they don't all belong to the same clone: three samples seem to
#' #' share one genomic sequence across the 250 SNPs (give or take one SNP), while
#' #' the other two samples share another.
#' load("../../RData/snpdata_extended.RData")
#' unique(problem_metadata$SAMPLE.CODE) # Five samples
#' par(mfrow = c(1,1), family = 'serif', mar = c(2,2,1,1))
#' cols <- RColorBrewer::brewer.pal(4,'Dark2')
#' snpdata_problem <- as.matrix(snpdata[,colnames(snpdata) %in% problem_metadata$SAMPLE.CODE])
#' image(snpdata_problem, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols[1:3])
#' 
#' # Are they all the same besides NAs? No
#' polymorphic <- apply(snpdata_problem, 1, function(x){
#'   length(unique(x[!is.na(x)]))>1}) 
#' 
#' # Let's look at the polymorphic SNPs: 
#' # Give or take a SNP, seems to be only two clones:
#' # two samples of one clone, three of the other 
#' image(snpdata_problem[polymorphic, ], 
#'       ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols[1:3])
#' 
#' # Whick ones are they among the masses: 
#' # Change there values so they are differently coloured
#' snpdata[,colnames(snpdata) %in% problem_metadata$SAMPLE.CODE] <-
#'   snpdata[,colnames(snpdata) %in% problem_metadata$SAMPLE.CODE] * 0.5
#' image(as.matrix(snpdata[,-(1:2)]), 
#'       ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols[1:3])
# ==================================================================

# Check for overlap:
table(metadata$SAMPLE.CODE %in% SNPData$SAMPLE.CODE) # No overlap
table(SNPData$SAMPLE.CODE %in% metadata$SAMPLE.CODE) # No overlap

# Check metadata available for all samples: yes
sids_snp <- colnames(snpdata[,-c(1:2)])
sids_meta <- c(metadata$SAMPLE.CODE, as.character(SNPData$SAMPLE.CODE))
all(sids_snp %in% sids_meta)
length(sids_snp); length(sids_meta)
all(sids_snp == unique(sids_snp)) # No duplicates
all(sids_meta == unique(sids_meta)) # No duplicates

#===========================================================
# Stitch together metadata into one file
#===========================================================
# First sort COLLECTION.DATE
class(SNPData$COLLECTION.DATE)
class(metadata$COLLECTION.DATE)

# Separate Year and Collection date for metadata
metadata$Year <- metadata$COLLECTION.DATE
metadata$Year <- sapply(strsplit(metadata$Year, split = "-"), function(x) x[1])

# Convert metadata Collection date into a data
COLLECTION.DATE_ <- as.Date(rep(NA, nrow(metadata)), format = "%Y-%m-%d")
for(i in 1:length(COLLECTION.DATE_)){
  COLLECTION.DATE_[i] <- as.Date(metadata$COLLECTION.DATE[i], format = "%Y-%m-%d")
}

# Check same where match: yes
COLLECTION.DATE_[!as.character(COLLECTION.DATE_) %in% as.character(metadata$COLLECTION.DATE)]

# Replace
metadata$COLLECTION.DATE <- COLLECTION.DATE_

# Join two data sets
metadata_extended <- dplyr::full_join(SNPData[,c(1:5, ncol(SNPData))], metadata)

# Check looks reasonable: yes :-) 
str(metadata_extended) # 526 obs. of  8 variables
lapply(metadata_extended, function(x) unique(x))
head(metadata_extended)

# Add a column to label data origin
metadata_extended$PloSGen2020 <- metadata$SAMPLE.CODE %in% SNPData$SAMPLE.CODE

# =============== Save the extended set of metadata  ===============
metadata <- metadata_extended
save(metadata , file = "../../RData/metadata_extended.RData")



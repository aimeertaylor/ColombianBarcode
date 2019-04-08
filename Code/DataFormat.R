##############################################################################################
# Formatting the Colombian data for relatedness inference and downstream analyses 
# Format for relatedness inference: chrom pos sample
# Time to run: 0.689 sec 
##############################################################################################
rm(list = ls())
library(tictoc)
tic()

# Load data, remove surplus and allocate dimnames
SNPData <- read.csv('../OriginalData/ColombianDataset_from_DiegoEcheverry_with_cities.csv',nrow = 325, row.names = 1)
SNPData <- SNPData[,1:255] # Remove surplus columns
SNPPos <- read.delim('../TxtData/Echeverry2013_3D7v3_posns.txt', skip = 16, sep = '\t', header = FALSE)
nsamp <- nrow(SNPData)
nSNPs <- nrow(SNPPos)
rownames(SNPPos) <- as.character(SNPPos$V1) # Name rows by SNP names for indexing downstream
rownames(SNPData) <- paste('sid', 1:nsamp, sep = '')
colnames(SNPData)[6:255] <- rownames(SNPPos)
         
# Re-order SNPs
pos <- SNPPos$V5
chroms <- as.numeric(do.call(rbind, strsplit(as.character(SNPPos$V2), split = '_'))[,2])
new_order <- rownames(SNPPos)
for(chrom in unique(chroms)){
  ind <- chroms == chrom
  X <- sort(as.numeric(pos[ind]), decreasing = FALSE, ind = TRUE)
  new_order[ind] <- new_order[ind][X$ix]
}

# Reorder and check
SNPData <- SNPData[,c(colnames(SNPData)[1:5], new_order)]
SNPPos <- SNPPos[new_order,]
any(colnames(SNPData)[6:255] != new_order)
any(colnames(SNPData)[6:255] != as.character(SNPPos$V1))

# Record which samples have het calls and how many
X <- as.matrix(SNPData[,6:255])
unique(as.vector(X))
het_count <- rowSums(X == 'AB')
sum(het_count > 0) # 185 with one or mor het calls (all het calls considered mistakes as samples previously deemed single genotype)

# Format data for HMM 
X[X == 'AA'] <- '1' #sample(c('0','1'), sum(X == 'AA'), replace = TRUE)
X[X == 'BB'] <- '0' #sample(c('0','1'), sum(X == 'BB'), replace = TRUE)
X[X == 'AB'] <- NA # Treat mixed signals as missing (since all these sample have previously been deemed single-genotype)
X[X == '--'] <- NA # Treat '--' as missing 
SNPData[,6:255] <- t(apply(X, 1, as.numeric))
SNPData$het_count <- het_count

# Recode as major(1), minor(0) - not strictly necessary
frequencies <- colMeans(SNPData[,6:255], na.rm = TRUE)
MAF_ind <- frequencies < 0.5
To_recode <- SNPData[,6:255][,MAF_ind] # Single out data where 0 currently encodes the minor allele
matrix1s <- array(data = 1, dim = dim(To_recode)) # Create matrix of 1s to add to To_recode
matrix1s[To_recode == 1] <- -1 # Make 1 -> -1 where the goal is to convert 1 -> 0 
SNPData[,6:255][,MAF_ind] <- To_recode + matrix1s # convert 0 + 1 -> 1 and 0 - 1 -> 0
unique(as.vector(as.matrix(SNPData[,6:255]))) # Check only 1, 0 and NA

# Replace Tunaco with Tumaco, which according to googlemaps and the paper is the correct spelling
# and remove accents
SNPData$City <- gsub('Tunaco', 'Tumaco', SNPData$City)
SNPData$City <- gsub('Tad\x97', 'Tado', SNPData$City)
SNPData$City <- gsub('Quibd\x97', 'Quibdo', SNPData$City)
SNPData$STATE <- gsub('Nari\x96o', 'Narino', as.character(SNPData$STATE))
SNPData$STATE <- gsub('Choc\x97', 'Choco', as.character(SNPData$STATE))
SNPData$STATE <- gsub('Nario ', 'Narino', as.character(SNPData$STATE))
SNPData$STATE <- gsub('Narino ', 'Narino', as.character(SNPData$STATE))
SNPData$STATE <- gsub('Valle ', 'Valle', as.character(SNPData$STATE))
unique(SNPData$City)
unique(SNPData$STATE)

# Change date formats
SNPData$COLLECTION.DATE <- as.Date(as.character(SNPData$COLLECTION.DATE), format = "%m/%d/%Y")
SNPData$Year <- do.call(rbind, strsplit(as.character(SNPData$COLLECTION.DATE), split = "-"))[,1]

# Format for HMM 
pos <- SNPPos$V5
chrom <- as.numeric(do.call(rbind, strsplit(as.character(SNPPos$V2), split = '_'))[,2])
HMMData <- cbind(chrom, pos, t(SNPData[,6:255]))

# Save data in format for HMM and RData
write.table(HMMData, file = '/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/TxtData/hmmInput.txt', 
            quote = FALSE, row.names = FALSE, sep = '\t')
save(SNPData, file = '../RData/SNPData.RData')
toc()

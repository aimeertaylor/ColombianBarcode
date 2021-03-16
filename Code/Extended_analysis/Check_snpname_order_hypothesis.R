########################################################################
# Script written to test the hypothesis that the markers were ordered
# alphabetically
########################################################################
require(dplyr)

# Angela's mapping based on sequences across sid matching
fuzzy_match <- read.delim("../../TxtData/GGCall_Version_Marker_Matches_RawCalls.txt")
gsub(":Aimee","",fuzzy_match$Version1)

# Mapping based on numerically reorderin
snp_original_order <- read.csv("../../OriginalData/Colombian snp data and flanking secs.csv")[-1,]
snp_reordered <- data.frame(snp_no = snp_original_order$X.11[1:250], 
                            snp_name_original = snp_original_order$X.12[1:250])
snp_reordered$snp_name_pos_original <- sapply(snp_original_order$X.12[1:250], function(x) as.numeric(strsplit(x, split = "-")[[1]][2]))
snp_reordered$snp_name_chrom_original <- sapply(snp_original_order$X.12[1:250], function(x) as.numeric(gsub("MAL", "", strsplit(x, split = "-")[[1]][1])))
snp_reordered$snp_name_reordered <- snp_reordered %>% 
  arrange(snp_name_chrom_original, snp_name_pos_original) %>% 
  pull(snp_name_original)


pair_sids_recordered <- apply(snp_reordered[, c("snp_name_original", "snp_name_reordered")], 1, function(x) paste(sort(x), collapse = " "))
pair_sids_fuzzymatch <- apply(cbind(gsub(":Aimee","",fuzzy_match$Version1), 
                                    gsub(":Angela","",fuzzy_match$Version2)), 1, function(x) paste(sort(x), collapse = " "))


test_result <- setequal(pair_sids_recordered, pair_sids_fuzzymatch)
if(test_result) {
  writeLines("Is alphabetically miss-ordered hypothesis correct? Yes")
} else {
  writeLines("Is alphabetically miss-ordered hypothesis correct? No")
}

save(snp_reordered, file = "../../RData/snp_reordered.RData")

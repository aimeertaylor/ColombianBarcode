rm(list = ls())
load("../../RData/Clonal_components_extended_all_LCIthrehold_0.75.RData")

# ================ In the "raw" data ================
dataset <- read.delim("../../recodedgoldengatedata/Diego-Vladimir-Fabian_GG3D7_Recode_21Dec2020.txt", 
                      check.names = F, # Stops Xs being added to names beginning with a digit and conversion of "-" to "."
                      stringsAsFactors = F) 

cc_10_and_143_data_unformatted <- dataset[, c(Clonal_components$cc_10, Clonal_components$cc_143)]

# Plot the raw data 
cols <- RColorBrewer::brewer.pal(3,'Dark2')
image(as.matrix(cc_10_and_143_data_unformatted), 
      ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)
legend('topleft', bty = 'n', legend =
         c('0', '1', 'Missing'), fill = c(cols[c(1,3)], adjustcolor('white')))
axis(side = 2, at = seq(0,1,length.out = ncol(cc_10_and_143_data_unformatted)), 
     labels = colnames(cc_10_and_143_data_unformatted), las = 1, cex.axis = 0.5)
axis(side = 1, at = seq(0,1,length.out = nrow(cc_10_and_143_data_unformatted)), 
     labels = 1:nrow(cc_10_and_143_data_unformatted), las = 1, cex.axis = 0.25)
axis(side = 3, at = seq(0,1,length.out = nrow(cc_10_and_143_data_unformatted)), 
     labels = 1:nrow(cc_10_and_143_data_unformatted), las = 1, cex.axis = 0.25)

write.csv(cc_10_and_143_data_unformatted, file = "../../cc_10_and_143_data_unformatted.csv")



# Equality between cc_10 and cc_143 (counts NA as a difference)
NA_snps <- apply(is.na(cc_10_and_143_data_unformatted),1,any)
equal_snps_inc_na <- sapply(1:nrow(dataset), function(i) {
  x <- dataset[i, Clonal_components$cc_10]
  y <- dataset[i, Clonal_components$cc_143]
  setequal(x[!is.na(x)], y[!is.na(y)])
})
diff_snps_inc_na <- sapply(1:nrow(dataset), function(i) {
  x <- dataset[i, Clonal_components$cc_10]
  y <- dataset[i, Clonal_components$cc_143]
  !setequal(x[!is.na(x)], y[!is.na(y)])
})

sum(equal_snps_inc_na)
sum(diff_snps_inc_na)
sum(NA_snps)

# ================ In the formatted data ================
load("../../RData/snpdata_extended.RData")
cc_10_and_143_data_formatted <- snpdata[, c(Clonal_components$cc_10, Clonal_components$cc_143)]
identical(cc_10_and_143_data_formatted, cc_10_and_143_data_unformatted) # Identical




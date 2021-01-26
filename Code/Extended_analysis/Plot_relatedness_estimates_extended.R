##############################################################
# Script to plot relatedness estimates before formatting
# (so that we can plot uninformative etc.)
# This script is a bit redundant: can also find plots in 
#
# Plots also reproduced elsewhere (Format_mles_CIs_extended.R)
##############################################################
rm(list = ls())
PDF <- F # Set to TRUE to plot graphs to pdf
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019
classification_cols = RColorBrewer::brewer.pal(5, "GnBu") # Colours for different thresholds 
names(classification_cols) = c("Unrelated", "Related", "Highly related", "Very highly related", "Clonal")
classification_cols <- c(Other = "gray", classification_cols)
cols <- c("gray", classification_cols["Clonal"]) # Update
if(PDF){pdf("../../Plots/Relatedmess_estimates_extended.pdf")}

# Load relatedness results 
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s.RData', freqs_used)) 

# NAs do not impact dplyr::arrange but they do impact sort
mle_CIs <- mle_CIs[!is.na(mle_CIs$rhat), ] # remove NAs 
# Arrange by rhat, then by CI within rhat
mle_CIs$CI_width <- mle_CIs$r97.5. - mle_CIs$r2.5.
mle_CIs <- dplyr::arrange(mle_CIs, rhat, desc(CI_width)) 

# SNP count below 10 after removing NAs
barplot(table(mle_CIs$snp_count[mle_CIs$snp_count < 10]))

# Load DE meta data file
load('../../RData/metadata_extended.RData') 
DE_sid <- metadata$SAMPLE.CODE[metadata$PloSGen2020]
DE_data_ind <- (mle_CIs$individual1 %in% DE_sid & mle_CIs$individual2 %in% DE_sid)
Between_ind <- !DE_data_ind & (mle_CIs$individual1 %in% DE_sid | mle_CIs$individual2 %in% DE_sid)

# Plot all not NA ====================================
plot(NULL, main = "All estimates", 
     ylim = c(0,1), 
     xlim = c(1,length(mle_CIs$rhat)), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
segments(x0 = 1:nrow(mle_CIs), x1 = 1:nrow(mle_CIs), # Add CIs
         y0 = mle_CIs$r2.5., y1 = mle_CIs$r97.5.,
         lwd = 0.1, col = 'gray')
lines(mle_CIs$rhat, lwd = 1, col = "black") # Add mles
barplot(table(mle_CIs$snp_count[mle_CIs$snp_count < 10]))

# Plot all informative ====================================
mle_CIs <- mle_CIs[!(mle_CIs$r2.5. < eps & mle_CIs$r97.5. > (1-eps)), ]

# SNP count below 10 after removing uninformative
barplot(table(mle_CIs$snp_count[mle_CIs$snp_count < 10]))

plot(NULL, main = "All informative estimates", 
     ylim = c(0,1), 
     xlim = c(1,length(mle_CIs$rhat)), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
segments(x0 = 1:nrow(mle_CIs), x1 = 1:nrow(mle_CIs), # Add CIs
         y0 = mle_CIs$r2.5., y1 = mle_CIs$r97.5.,
         lwd = 0.1, col = 'gray')
lines(mle_CIs$rhat, lwd = 1, col = "black") # Add mles


# Plot Taylor et al. 2020 estimates, coloured by clone =========================
# where clone is defined as having UCI "touching" one.
# Get Diego's data indices from mle_CIs with NAs and uninformative removed
clone_ind <- mle_CIs$r97.5.[DE_data_ind] > (1-eps)

plot(NULL, main = "Re-estimates on data of Taylor et al. 2020", 
     ylim = c(0,1), 
     xlim = c(1,sum(DE_data_ind)), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
segments(x0 = 1:nrow(mle_CIs[DE_data_ind,]), # Add CIs
         x1 = 1:nrow(mle_CIs[DE_data_ind,]),
         y0 = mle_CIs$r2.5.[DE_data_ind], 
         y1 = mle_CIs$r97.5.[DE_data_ind],
         lwd = 0.1, col = cols[clone_ind+1])
lines(mle_CIs$rhat[DE_data_ind], lwd = 1, col = "black") # Add mles
abline(h = 1-0.01, lty = 'dashed', col = cols["Clonal"])


# Plot all new estimates, coloured by clone =======================
# where clone is defined as having LCI
# greater than 0.5 and UCI "touching" one. 
clone_ind <- mle_CIs$r2.5.[!DE_data_ind] > 0.5 & mle_CIs$r97.5.[!DE_data_ind] > (1-eps)

plot(NULL, main = "All estimates on data not in Taylor et al. 2020", 
     ylim = c(0,1), 
     xlim = c(1,sum(!DE_data_ind)), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
segments(x0 = 1:nrow(mle_CIs[!DE_data_ind,]), # Add CIs
         x1 = 1:nrow(mle_CIs[!DE_data_ind,]),
         y0 = mle_CIs$r2.5.[!DE_data_ind], 
         y1 = mle_CIs$r97.5.[!DE_data_ind],
         lwd = 0.1, col = cols[clone_ind+1])
lines(mle_CIs$rhat[!DE_data_ind], lwd = 1, col = "black") # Add mles
abline(h = 1-eps, lty = 'dashed', col = cols["Clonal"])
abline(h = 0.5, lty = 'dashed', col = cols["Clonal"])


# Plot new clonal estimates only ============================
plot(NULL, main = "All estimates on data not in Taylor et al. 2020: Zoom", 
     ylim = c(0,1), 
     xlim = c(73000,sum(!DE_data_ind)), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
segments(x0 = 1:nrow(mle_CIs[!DE_data_ind,]), # Add CIs
         x1 = 1:nrow(mle_CIs[!DE_data_ind,]),
         y0 = mle_CIs$r2.5.[!DE_data_ind], 
         y1 = mle_CIs$r97.5.[!DE_data_ind],
         lwd = 0.5, col = cols[clone_ind+1])
lines(mle_CIs$rhat[!DE_data_ind], lwd = 1, col = "black") # Add mles
abline(h = 1-eps, lty = 'dashed', col = cols["Clonal"])
abline(h = 0.5, lty = 'dashed', col = cols["Clonal"])


# Plot estimates for comparisons beteen DE and new =============
# Redefine clones among Between_ind
clone_ind <- mle_CIs$r2.5.[Between_ind] > 0.5 & mle_CIs$r97.5.[Between_ind] > (1-eps)

plot(NULL, main = "Estimates linking to data in Taylor et al. 2020", 
     ylim = c(0,1), 
     xlim = c(1,sum(Between_ind)), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
segments(x0 = 1:nrow(mle_CIs[Between_ind,]), # Add CIs
         x1 = 1:nrow(mle_CIs[Between_ind,]),
         y0 = mle_CIs$r2.5.[Between_ind], 
         y1 = mle_CIs$r97.5.[Between_ind],
         lwd = 0.1, col = cols[clone_ind+1])
lines(mle_CIs$rhat[Between_ind], lwd = 1, col = "black") # Add mles
abline(h = 1-eps, lty = 'dashed', col = cols["Clonal"]) # Here clone was defined by 
abline(h = 0.5, lty = 'dashed', col = cols["Clonal"]) # Here clone was defined by 



# Zoomed-in plot estimates for comparisons beteen DE and new =============
# NULL plot
plot(NULL, main = "Estimates linking to data in Taylor et al. 2020 (Zoom)", 
     ylim = c(0,1), 
     xlim = c(59800,sum(Between_ind)), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
segments(x0 = 1:nrow(mle_CIs[Between_ind,]), # Add CIs
         x1 = 1:nrow(mle_CIs[Between_ind,]),
         y0 = mle_CIs$r2.5.[Between_ind], 
         y1 = mle_CIs$r97.5.[Between_ind],
         lwd = 1, col = cols[clone_ind+1])
lines(mle_CIs$rhat[Between_ind], lwd = 1, col = "black") # Add mles
abline(h = 1-eps, lty = 'dashed', col = cols["Clonal"]) # Here clone was defined by 
abline(h = 0.5, lty = 'dashed', col = cols["Clonal"]) # Here clone was defined by 



# More classifications are possible? =======
# Differs to Taylor et al. 2020 where we classify additionally
# as highly related with 0.25 and 0.50
unrelated <- mle_CIs$r2.5. < eps & mle_CIs$r97.5. < 0.5
clone <- mle_CIs$r2.5. > 0.5 & mle_CIs$r97.5. > (1-eps)
related <- mle_CIs$r2.5. > eps & mle_CIs$r97.5. < (1-eps)
other <- !(unrelated | related | clone)

Classifications <- rep(NA, nrow(mle_CIs))
Classifications[unrelated] <- "Unrelated"
Classifications[related] <- "Related"
Classifications[clone] <- "Clonal"
Classifications[other] <- "Other"

# All estimates
plot(NULL, 
     main = "All estimates", 
     ylim = c(0,1), xlim = c(1,nrow(mle_CIs)), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
abline(h = 0.5, lty = "dashed")
segments(x0 = 1:nrow(mle_CIs), x1 = 1:nrow(mle_CIs),
         y0 = mle_CIs$r2.5., y1 = mle_CIs$r97.5.,
         lwd = 0.1, col = classification_cols[Classifications])
lines(mle_CIs$rhat, lwd = 1, col = "black") 

# Between estimates
plot(NULL, 
     main = "All estimates linking to data in Taylor et al. 2020", 
     ylim = c(0,1), xlim = c(1,nrow(mle_CIs[Between_ind,])), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
abline(h = 0.5, lty = "dashed")
segments(x0 = 1:nrow(mle_CIs[Between_ind,]), x1 = 1:nrow(mle_CIs[Between_ind,]),
         y0 = mle_CIs$r2.5.[Between_ind], y1 = mle_CIs$r97.5.[Between_ind],
         lwd = 0.1, col = classification_cols[Classifications[Between_ind]])
lines(mle_CIs$rhat[Between_ind], lwd = 1, col = "black") 

# New estimates
plot(NULL, 
     main = "All estimates on data not in Taylor et al. 2020", 
     ylim = c(0,1), xlim = c(1,nrow(mle_CIs[!DE_data_ind,])), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
abline(h = 0.5, lty = "dashed")
segments(x0 = 1:nrow(mle_CIs[!DE_data_ind,]), x1 = 1:nrow(mle_CIs[!DE_data_ind,]),
         y0 = mle_CIs$r2.5.[!DE_data_ind], y1 = mle_CIs$r97.5.[!DE_data_ind],
         lwd = 0.1, col = classification_cols[Classifications[!DE_data_ind]])
lines(mle_CIs$rhat[!DE_data_ind], lwd = 1, col = "black") 

# Original estimates
plot(NULL, 
     main = "Re-estimates on data of Taylor et al. 2020", 
     ylim = c(0,1), xlim = c(1,nrow(mle_CIs[DE_data_ind,])), 
     ylab = 'Relatedness estimate', 
     xlab = "Sample pair index", 
     las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
abline(h = 0.5, lty = "dashed")
segments(x0 = 1:nrow(mle_CIs[DE_data_ind,]), x1 = 1:nrow(mle_CIs[DE_data_ind,]),
         y0 = mle_CIs$r2.5.[DE_data_ind], y1 = mle_CIs$r97.5.[DE_data_ind],
         lwd = 0.1, col = classification_cols[Classifications[DE_data_ind]])
lines(mle_CIs$rhat[DE_data_ind], lwd = 1, col = "black") 

if(PDF)dev.off()


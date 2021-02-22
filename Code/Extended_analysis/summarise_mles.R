# Function to summarise mles to track formatting
summarise_mles <- function(x, metadata_, snp_count_threshold = 10, 
                           zoom = F, eps = 0.01){
  
  x <- dplyr::arrange(x, rhat, CI_width)
  
  # Plot relatedness values
  plot(NULL, main = "", 
       ylim = c(0,1), 
       xlim = c(1,length(x$rhat)), 
       ylab = 'Relatedness estimate', 
       xlab = "Sample pair index", 
       las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
  axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
  segments(x0 = 1:nrow(x), x1 = 1:nrow(x), # Add CIs
           y0 = x$r2.5., y1 = x$r97.5.,
           lwd = 0.1, col = 'gray')
  lines(x$rhat, lwd = 1, col = "black") # Add mles
  abline(h = 1-eps, lty = "dashed", col = 'blue', lwd = 0.5)
  
  if(zoom) {
    # Plot relatedness values
    plot(NULL, main = "", 
         ylim = c(0,1), 
         xlim = c(sum(x$rhat < 0.8, na.rm = T),length(x$rhat)), 
         ylab = 'Relatedness estimate', 
         xlab = "Sample pair index", 
         las = 2, panel.first = grid(), bty = "n", tck = -0.01, xaxt = "n")
    axis(1, cex.axis = 1, las = 2, line = -1.5, tick = FALSE)
    segments(x0 = 1:nrow(x), x1 = 1:nrow(x), # Add CIs
             y0 = x$r2.5., y1 = x$r97.5.,
             lwd = 0.1, col = 'gray')
    lines(x$rhat, lwd = 1, col = "black") # Add mles
    abline(h = 1-eps, lty = "dashed", col = 'blue', lwd = 0.5)
  }
  
  # Summarise data paucity
  snp_counts <- array(0, dim = c(2,snp_count_threshold), 
                      dimnames = list(c("per_sample", "per_sample_pair"), 1:snp_count_threshold))
  snp_counts[1, ] <- table(metadata_[unique(c(x$individual1, x$individual2)), "snp_count"])[colnames(snp_counts)]
  snp_counts[2, ] <- table(x$snp_count)[colnames(snp_counts)]
  
  writeLines(paste(sprintf("No. sample pairs: %s", nrow(x)), 
                   sprintf("No. samples: %s", length(unique(c(x$individual1, x$individual2)))), 
                   "Numbers of samples and sample pairs with sparse data:", sep = "\n"))
  print(snp_counts)
}

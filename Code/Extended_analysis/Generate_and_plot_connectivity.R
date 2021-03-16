###############################################################################
#' Script to explore evidence of connectivity accross sampling locations. Could
#' re-analyse without relatedness estimates whose CIs range from zero to one.
###############################################################################
library(tidymodels)
load("../../RData/mles_CIs_extended_freqsTaylor2020_meta.RData")
load("../../RData/metadata_extended.RData")
nboot <- 50
PDF <- TRUE

#============== Data wrangling =================

# Compute the fraction highly related
high_relatedness_threshold <- 0.25
mle_CIs$highly_related <- as.numeric(mle_CIs$rhat > high_relatedness_threshold)

# Extract the fraction related and number of pairs
fraction_highly_related <- mle_CIs %>%
  group_by(City12) %>%
  summarise(fHR = mean(highly_related), 
            npairs = length(highly_related)) %>%
  arrange(fHR)

# Compute CIs using the bootstrap 
# use group by (instead of argument strata because 
# "strata below 10% of the total are pooled together")
boots <- mle_CIs %>%
  group_by(City12) %>%
  bootstraps(times = nboot) 

# Matrix to collect fractions highly related
boots_fHR_matrix = array(dim = c(nrow(fraction_highly_related), nboot), 
                         dimnames = list(fraction_highly_related$City12, NULL))

# Extract fractions highly related into the matrix  
for(i in 1:length(boots$splits)){ # using forloop because boots$splits not all same length
  x <- analysis(boots$splits[[i]])
  booti <- x %>%
    group_by(City12) %>%
    summarise(fHR = mean(highly_related))  
  boots_fHR_matrix[booti$City12, i] <- booti$fHR
}

# Compute CIs
CIs <- t(apply(boots_fHR_matrix, 1, quantile, probs = c(0.025, 0.975), na.rm = T))

# Extract ind for within same city or not
within_city_ind <- sapply(fraction_highly_related$City12, function(x){
  split_cities <- strsplit(x, split = "_")[[1]]
  if("Iscuande" %in% split_cities) {
    length(unique(split_cities)) == 2
  } else {
    length(unique(split_cities)) == 1
  }
})

# Extract ind for in Taylor et al. 2020 or not
load("../../RData/geo_dist_info_cities.RData")
old_ind <- fraction_highly_related$City12 %in% names(geo_dist_info$pairwise_site_distance_all)

# Check
all(rownames(boots_fHR_matrix) == fraction_highly_related$City12)

#============== Plots =================
if(PDF) pdf("../../Plots/Connectivity_extended.pdf")
for(n_threshold in c(0,100,200)){
  
  plot_inds <- fraction_highly_related$npairs > n_threshold 
  
  # Plot all 
  plot(x = fraction_highly_related$fHR[plot_inds], 
       y = 1:nrow(fraction_highly_related[plot_inds,]),
       pch = c(16,17)[old_ind[plot_inds] + 1], 
       col = c("red","blue")[within_city_ind[plot_inds] + 1],
       xlim = c(0, 1.2), bty = "n", cex.main = 0.75, cex = 0.5, 
       xaxt = "n", yaxt = "n", ylab = "", 
       xlab = "Fraction highly related", 
       main = sprintf("All location pairs with %s or more sample pairs", n_threshold),
       panel.first = abline(v = seq(0,1,0.2), 
                            col = "lightgray", 
                            lty = "dotted"))
  
  text(y = 1:nrow(fraction_highly_related[plot_inds,]), 
       x = fraction_highly_related$fHR[plot_inds], 
       labels = fraction_highly_related$City12[plot_inds], 
       cex = 0.25, pos = 4)
  
  axis(side = 1, at = seq(1,0,-0.2))
  
  # Add boostrapped fractions
  segments(x0 = CIs[plot_inds, "2.5%"], 
           x1 = CIs[plot_inds, "97.5%"],
           y0 = 1:sum(plot_inds),
           y1 = 1:sum(plot_inds),
           col = sapply(c("red","blue")[within_city_ind[plot_inds] + 1], 
                        adjustcolor, alpha.f = 0.5))
  
  # Circle points
  port_ind <- which(fraction_highly_related$City12[plot_inds] == "Esmeraldas_Tumaco" | 
                      fraction_highly_related$City12[plot_inds] == "Buenaventura_Tumaco" |
                      fraction_highly_related$City12[plot_inds] == "Buenaventura_Esmeraldas" |
                      fraction_highly_related$City12[plot_inds] == "Buenaventura_Buenaventura"|
                      fraction_highly_related$City12[plot_inds] == "Tumaco_Tumaco" | 
                      fraction_highly_related$City12[plot_inds] == "Esmeraldas_Esmeraldas") 
  
  points(x = fraction_highly_related$fHR[plot_inds][port_ind], 
         y = (1:nrow(fraction_highly_related[plot_inds,]))[port_ind],
         pch = 1)
  
  # Legends
  legend("right", pch = c(16,17,1), 
         legend = c("Does not feature in Taylor et al. 2020", 
                    "Features in Taylor et al. 2020", 
                    "Comparison across ports"), 
         bty = "n", cex = 0.5, pt.cex = 0.75, title.adj = 0,  
         title = "Point character coding")
  
  legend("bottomright", col = c("blue", "red"), pch = 15,  
         legend = c("Within loction comparison", 
                    "Accross location comparison"), 
         bty = "n", cex = 0.5, pt.cex = 1, title.adj = 0,  
         title = "Colour coding")
}
if(PDF) dev.off()
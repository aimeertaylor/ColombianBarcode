################################################################################
# Generate_components.R and Compare_components.R and Plot_extended_components.R
# Comparing back to original analysis done with CIs 
################################################################################
rm(list = ls())
library(igraph) 
library(RColorBrewer)
source('../igraph_functions.R') # construct_adj_matrix
source("./summarise_mles.R")
load('../../RData/metadata_extended.RData')
load('../../RData/mles_CIs_extended_freqsTaylor2020.RData')
load(file = "../../RData/Clonal_components.RData") # For comparison with CC_extended
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019
a = 0; b = 1; c = 1 # Scaling parameters for plotting edge width and transparancy 
PDF <- TRUE
LCI_threshold <- 0.75
cut_offs <- c(50, 100)
set.seed(1)

# Define city cols for plotting
cities <- unique(metadata$City)
cols_cities <- array(c(rev(brewer.pal(5, 'Spectral')), 
                       brewer.pal(length(cities)-5, 'Dark2')), 
                     dimnames = list(cities))

# Function to progressively remove samples with high NA value counts
remv_function <- function(input){
  A_est_full <- construct_adj_matrix(input, Entry = 'rhat')
  A_est_full[lower.tri(A_est_full)] <- t(A_est_full)[lower.tri(A_est_full)]
  diag(A_est_full) <- 1
  na_count_per_sample <- rowSums(is.na(A_est_full))
  while(max(na_count_per_sample) > 0){
    to_remove <- which.max(na_count_per_sample)
    A_est_full <- A_est_full[-to_remove, -to_remove]
    na_count_per_sample <- rowSums(is.na(A_est_full))
  }
  sids_keep <- unique(colnames(A_est_full))
  return(sids_keep)
}

if(PDF) pdf("../../Plots/CCs_snp_cutoff.pdf")

for(cut_off in cut_offs){
  
  # Filtering by SNP cut offs
  mles <- mle_CIs[mle_CIs$snp_count > cut_off,] # Filter relatedness 
  sids_to_keep <- remv_function(mles)
  inds_to_keep <- mles$individual1 %in% sids_to_keep & mles$individual2 %in% sids_to_keep
  mles <- mles[inds_to_keep,] # Filer samples
  mles$CI_width <- mles$r97.5. - mles$r2.5. # Needed for summarise_mles
  
  summarise_mles(mles, metadata_ = metadata)
  mles <- mles[,-(which(colnames(mles) %in% c("r2.5.", "r97.5.")))] # Remove s.t. cannot be accidentally used
  
  # Using definition of clone where UCI "touches" one and LCI > threshold
  A_est <- construct_adj_matrix(mles, Entry = 'rhat')
  A_clonal <- A_est 
  A_clonal[A_est < 0.9] <- 0
  G_clonal <- graph_from_adjacency_matrix(A_clonal, mode='upper', diag=F, weighted = T) 
  ccompnts <- components(G_clonal)

  # Extract sids per component, including components of size one 
  CC_extended <- lapply(1:length(ccompnts$csize), function(cc){
    names(which(ccompnts$membership == as.numeric(cc)))
  })
  names(CC_extended) <- paste0("cc_", 1:length(CC_extended))
  
  #==================================
  # Comparing to the original
  #==================================
  identical_components <- list()
  extended_components <- list()
  broken_components <- list()
  
  # Cartegorise each of the clonal components reported in Taylor et al. 2020
  for(i in names(Clonal_components)){
    for(j in names(CC_extended)){
      if(setequal(Clonal_components[[i]],CC_extended[[j]])) {
        identical_components[[i]][[j]] <- "Equal"
      } else if (all(Clonal_components[[i]] %in% CC_extended[[j]])) {
        extended_components[[i]][[j]][["original"]] <- Clonal_components[[i]]
        extended_components[[i]][[j]][["additional"]] <- setdiff(CC_extended[[j]],Clonal_components[[i]])
      } else {
        intersecting_sids <- intersect(Clonal_components[[i]],CC_extended[[j]])
        if(length(intersecting_sids) > 0) {
          broken_components[[i]][[j]][["intersect"]] <- intersecting_sids
          broken_components[[i]][[j]][["additional"]] <- setdiff(CC_extended[[j]],Clonal_components[[i]])
        } else {
          next()
        }
      }
    }
  }
  
  sids_extended_components <- unique(unlist(extended_components))
  years_range <- range(sort(unique(as.numeric(metadata[sids_extended_components,"Year"]))))
  
  # Category breakdown
  writeLines(sprintf("Number of broken components: %s", length(broken_components)))
  writeLines(sprintf("Number of identical components: %s", length(identical_components)))
  writeLines(sprintf("Number of extended components: %s", length(extended_components)))
  writeLines(sprintf("Evidence of clonal propagation: %s", paste(years_range, collapse = " ")))
  writeLines(sprintf("Number of samples among extended ccs: %s", length(sids_extended_components)))
  
  # save
  save(CC_extended, extended_components, identical_components, broken_components,
       file = sprintf("../../RData/CCs_snp_cutoff%s", cut_off))

  
  # ============== All components plot ==============
  V(G_clonal)$marker_count <- metadata[V(G_clonal)$name,  "snp_count"]
  V(G_clonal)$city <- metadata[V(G_clonal)$name, "City"]
  V(G_clonal)$color <- cols_cities[V(G_clonal)$city] 
  
  # Rescale edge weights for plotting
  if (identical(round(unique(E(G_clonal)$weight)), 1)){
    weights_rescaled <- E(G_clonal)$weight
  } else {
    weights_rescaled = a + ((b-a) * (E(G_clonal)$weight - min(E(G_clonal)$weight)) / (max(E(G_clonal)$weight) - min(E(G_clonal)$weight)))
  }
  
  # Plot graph 
  plot(G_clonal,
       vertex.label.color = "black", 
       vertex.label.cex = 0.5, 
       vertex.label = NA, 
       vertex.frame.color = NA, 
       vertex.shape = c("square","circle")[metadata[V(G_clonal)$name, "PloSGen2020"] + 1],
       edge.width = weights_rescaled * c, 
       edge.color = sapply(weights_rescaled, function(x)adjustcolor('black', alpha.f = x)),
       vertex.size = 0.5*log(metadata[V(G_clonal)$name, "snp_count"]))
  
  # Title
  title(main = sprintf("Clonal def. UCI >= (1-%s) & LCI > %s", 
                       eps, LCI_threshold), cex.main = 1)
  
  # Legend
  cities_ <- unique(V(G_clonal)$city)
  legend('bottomleft', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, 
         inset = -0.1, col = cols_cities[cities_], legend = cities_)
 
  
  
  # ============== Extended component plot ==============
  cities_inc <- unique(metadata[unlist(lapply(extended_components, unlist)),"City"])
  x <- unique(lapply(extended_components, unlist))
  
  # Check for any converged components
  x_compared <- sapply(x, function(xi) sapply(x, function(xj) setequal(xi,xj)))
  if (any(x_compared[lower.tri(x_compared)])) {
    stop("Some components converged: modify by hand")
  }

  cities_per_year_tab <- array(0, dim = length(cities_inc), dimnames = list(cities_inc))
  y_pos <- seq(-1,1,length.out = length(x)+1)
  
  par(mar = c(7,1,7,1))
  for(i in 1:length(x)){
    
    city_years <- metadata[x[[i]],c("City","Year")]
    years_to_plot <- table(city_years$Year)
    cities_per_year <- lapply(names(years_to_plot), function(y, cities_per_year_tab) {
      z <- table(city_years$City[city_years$Year == y])
      cities_per_year_tab[names(z)] <- z
      return(cities_per_year_tab)}, cities_per_year_tab)
    
    Adjmatrix <- matrix(1, ncol = length(years_to_plot), nrow = length(years_to_plot))
    diag(Adjmatrix) <- 0
    Graph <- igraph::graph_from_adjacency_matrix(Adjmatrix, mode = "undirected")
    
    plot(Graph,
         add = (i!=1), 
         ylim = c(1,length(x)), 
         xlim = years_range, 
         rescale = FALSE, 
         layout = cbind(as.numeric(names(years_to_plot)),rep(i, length(years_to_plot))), 
         vertex.shape = "pie", 
         vertex.pie = cities_per_year,#[V(Graph)$name], 
         vertex.pie.color = list(cols_cities[cities_inc]), 
         vertex.pie.lwd = 0.25, 
         vertex.frame.color = NA, #NB, colours pi outline 
         vertex.size = years_to_plot*7, 
         vertex.label = NA,
         edge.width = 0.5,
         edge.color = "black")
  }
  
  axis(side = 1, at = min(years_range):max(years_range), cex.axis = 1, las = 2)
  legend("top",legend = cities_inc, pch = 16, ncol = 2, 
         col = cols_cities[cities_inc], cex = 0.75, pt.cex = 1, bty = "n")

}

if(PDF) dev.off()




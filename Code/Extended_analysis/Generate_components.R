################################################################################
#' Script to generate and plot clonal components for the extended analysis. To
#' see how FSVC sids (samples that didn't feature in Taylor et al. 2020) cluster
#' with DE sids (samples that did feature in Taylor et al. 2020), generate
#' clonal components from all sids together. To compare 250 SNP based clonal
#' components with WGS-based clusters, generate clonal components from FSVC sids
#' only. Because NA graph edges are ambiguous (setting NAs to one augments
#' connectivity, setting NAs to zero decreases it), we remove samples with one
#' or more NA relatedness estimates; see Generate_sids_remv.R. Because the
#' relatedness estimates for FSVC sample ids are associated with more
#' uncertainty than those for DE sids we consider three different definitions of
#' a clone: statistically indistinguishable from one with UCI >= 1-eps
#' (definition used in Taylor et al. 2020); statistically indistinguishable from
#' one and distinguishable from 0.75 with UCI >= 1-eps and LCI > 0.75
#' (definition equivalent to that used in Taylor et al. 2020 because 
#' 0.75<min(LCI of DE sid clones) - see below); and statistically
#' indistinguishable from and distinguishable from 0.95 (this definition is
#' more stringent than the one used in Taylor et al. 2020 and so breaks up many
#' of the 46 clonal components listed in Taylor et al. 2020 - see
#' Compare_components.R)
################################################################################
rm(list = ls())
library(igraph) 
library(RColorBrewer)
source('../igraph_functions.R') # construct_adj_matrix
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019
a = 0; b = 1; c = 1 # Scaling parameters for plotting edge width and transparancy 
PDF <- TRUE
label_art_pnts <- TRUE
set.seed(1)

# Load relatedness results and metadata
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s_meta.RData', freqs_used)) 
load(sprintf('../../RData/metadata_extended.RData'))
load(file = "../../RData/sids_to_remove_from_graphs.RData")

# Define city cols
cities <- unique(metadata$City)
cols_cities <- array(c(rev(brewer.pal(5, 'Spectral')), 
                       brewer.pal(length(cities)-5, 'Dark2')), 
                     dimnames = list(cities))

# Remove samples with one or more NA estimate
sids_to_keep <- setdiff(metadata$SAMPLE.CODE, sids_remv) 
inds_to_keep <- mle_CIs$individual1 %in% sids_to_keep & mle_CIs$individual2 %in% sids_to_keep
mles_to_keep <- mle_CIs[inds_to_keep,]

# Extract sids ordered by date
DE_sid <- dplyr::arrange(metadata[metadata$PloSGen2020, ], COLLECTION.DATE)$SAMPLE.CODE
FSVC_sid <- dplyr::arrange(metadata[!metadata$PloSGen2020, ], COLLECTION.DATE)$SAMPLE.CODE
DE_pair_ind <- (mles_to_keep$individual1 %in% DE_sid) & (mles_to_keep$individual2 %in% DE_sid)
FSVC_pair_ind <- (mles_to_keep$individual1 %in% FSVC_sid) & (mles_to_keep$individual2 %in% FSVC_sid)

# Find the min(LCI of DE sid clones) and define LCI thresholds accordingly
# Differs between original mles_CIs and new due to randomness of parametric bootstrap
# Original mles_CIs: 
load("../../RData/mles_CIs.RData")
min(mle_CIs$r2.5.[mle_CIs$r97.5. >= 1-eps]) # Min LCI
min(mle_CIs$rhat[mle_CIs$r97.5. >= 1-eps]) # Min relatedness
hist(mle_CIs$r2.5.[mle_CIs$r97.5. >= 1-eps], 
     xlab = "Lower CI limit among clones", main = "")
rm(mle_CIs)

# New DE_sids mles_CIs
min(mles_to_keep$r2.5.[mles_to_keep$r97.5. >= 1-eps & DE_pair_ind])
hist(mles_to_keep$r2.5.[mles_to_keep$r97.5. >= 1-eps & DE_pair_ind], 
     xlab = "Lower CI limit among clones", main = "")

# Define threholds
LCI_thresholds <- c(0,0.75,0.95)

if(PDF) pdf("../../Plots/Clonal_components.pdf")
for(comparisons in c("all", "FSVC")){
  for(LCI_threshold in LCI_thresholds){
    
    if(comparisons == "all") {
      mles <- mles_to_keep
    } else {
      mles <- mles_to_keep[FSVC_pair_ind,]
    }
    
    # Using definition of clone where UCI "touches" one and LCI > threshold
    A_low <- construct_adj_matrix(mles, Entry = 'r2.5.')
    A_high <- construct_adj_matrix(mles, Entry = 'r97.5.')
    A_est <- construct_adj_matrix(mles, Entry = 'rhat')
    
    #' Extract a adjacency matrix for clones In Taylor et al. I used A_clonal <-
    #' A_high. This difference doesn't change the clonal definition and thus
    #' doesn't change the component deinition. However, it does allow
    #' visualisation of the extent of relatedness within clonal components. That
    #' wasn't necessary in Taylor et al. 2020, because the clonal components
    #' were more-or-less all cliques (min rhat = 0.94; see above). However, it
    #' is helpful here because there are some low relatedness point estimates
    #' among those that are statistically indistinguishable from one.
    A_clonal <- A_est 
    A_clonal[A_high < (1-eps) | A_low <= LCI_threshold] <- 0
    writeLines(sprintf("Minimum point estimate among those considered clones: %s", 
                       round(sort(unique(as.vector(A_clonal)))[2],2)))

    # Create graph 
    G_clonal <- graph_from_adjacency_matrix(A_clonal, mode='upper', diag=F, weighted = T) 
    art_pnts <- names(articulation.points(G_clonal))
    ccompnts <- components(G_clonal)
    
    # Add metadata
    V(G_clonal)$marker_count <- metadata[V(G_clonal)$name,  "snp_count"]
    V(G_clonal)$city <- metadata[V(G_clonal)$name, "City"]
    V(G_clonal)$color <- cols_cities[V(G_clonal)$city] 
    
    # Highlight articulation points if any
    if(length(art_pnts) > 0 & label_art_pnts){
      vertex_labels <- V(G_clonal)$name
      vertex_labels[!V(G_clonal)$name %in% art_pnts] <- ""
    }
    
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
    
    # Extract sids per component
    # include components of size one (for Compare_components.R)
    Clonal_components <- lapply(1:length(ccompnts$csize), function(cc){
      names(which(ccompnts$membership == as.numeric(cc)))
    })
    names(Clonal_components) <- paste0("cc_", 1:length(Clonal_components))
    
    # How many non-cliques? 
    # Check to see how many cliques there are per clonal_component
    clique_count_per_component <- sapply(names(Clonal_components), function(cc){
      sids <- Clonal_components[[cc]]
      if(length(sids) > 1){
        subG <- subgraph(G_clonal, sids)
        x <- maximal.cliques.count(subG, min = 2)
        # Plot, remembering this is after setting edges to zero if not clonal
        if(x != 1) plot(subG, 
                        main = cc, 
                        sub = x,
                        vertex.color = cols_cities[metadata[V(subG)$name, "City"]] , 
                        edge.width = 1, 
                        edge.color = "black", 
                        vertex.size = 5, 
                        vertex.label.size = 0.5, 
                        vertex.label = V(subG)$name, #metadata[V(subG)$name, "Year"], 
                        vertex.label.color = "black", 
                        vertex.shape = c("square","circle")[metadata[V(subG)$name, "PloSGen2020"] + 1])
        return(x)
      } else {
        return(0)
      }
    })
    names(clique_count_per_component) <- names(Clonal_components)
  
    # save
    save(Clonal_components, clique_count_per_component, 
         file = sprintf("../../RData/Clonal_components_extended_%s_LCIthrehold_%s.RData", 
                        comparisons, LCI_threshold))
  }
}

if(PDF) dev.off()




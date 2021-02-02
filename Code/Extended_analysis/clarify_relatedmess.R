#=============================================================
#' This script was written to guide structure discovery in relatedmess graphs
#' based on the removal of vertices that belong to many cliques, removing those
#' with sparse data first. It requires a graph with named vertices that also
#' have marker_count, city and colour variables
#'
#' On Feb 2 I discovered that graph_from_adjacency_matrix imputes NAs as "1"
#' when weighted = T (and 0 when weighted = NULL). I thus realised that
#' membership to many cliques is just a correlate of have many NA relatedness
#' comparisons. This raises several points. First, instead of pruning based on
#' clique membership, we should prune based on NA comparisons, rendering this
#' scipt obsolte. Second, is pruning based on NA comparisons a more round-about
#' way of introducing a hard cut off on marker count per sample? Third, could we
#' use the transitivity property of clones etc. to impute NAs? This third point
#' is a project in itself. 
#'
#' Compare any relatedmess approach with thresholding and ground truth
#=============================================================
require(igraph)
source("../Extended_analysis/get_vertex_clique_matrix.R")
PDF <- FALSE
if(PDF) pdf("../../Plots/Clarify_relatedmess.pdf")
par(mfrow = c(1,1))

# Check the graph has name vertices
if (is.null(V(inputG)$name) | 
    is.null(V(inputG)$marker_count) | 
    is.null(V(inputG)$city) |
    is.null(V(inputG)$color)) {
  stop("inputG requires named vertices with marker_counts, city and colour variables")
}

# Extract sample ids and marker counts
sids <- V(inputG)$name
marker_counts <- array(V(inputG)$marker_count, dimnames = list(sids))

# Get the matrix of clique membership per vertex
clique_sample_matrix <- get_vertex_clique_matrix(inputG)

# Extract number of cliques that each sample belongs to and sort
clique_count_per_sample <- sort(colSums(clique_sample_matrix), decreasing = T)

# Get maximums for plotting limits
max_clique_count_per_sample <- max(clique_count_per_sample)
max_marker_count <- max(marker_counts[names(clique_count_per_sample)])

# Plot clique versus marker count for inputG
clique_count <- sum(apply(clique_sample_matrix, 1, any))
plot(y = as.numeric(clique_count_per_sample), 
     x = marker_counts[names(clique_count_per_sample)], 
     ylab = "Per-sample clique membership count", 
     xlab = "Per-sample marker data count",
     bty = "n", pch = 20, cex.main = 1, 
     main = sprintf("Total number of cliques: %s", clique_count)) 

# Initiate pruning
max_acceptable_cliques <- 1
sids_to_maybe_remove <- names(which(clique_count_per_sample > max_acceptable_cliques))
marker_counts_of_sids_to_maybe_remove <- sort(unique(marker_counts[sids_to_maybe_remove]))
threshold <- marker_counts_of_sids_to_maybe_remove[1]
sids_to_keep <- sids
sub_graphs <- list()

# Prune
while(threshold <= max_marker_counts){ 
  
  # Obtain and record sub-graph
  sids_to_remove <- names(which(marker_counts[sids_to_maybe_remove] <= threshold))
  sids_to_keep <- sids_to_keep[!sids_to_keep %in% sids_to_remove]
  subG <- subgraph(inputG, sids_to_keep)
  sub_graphs[[as.character(threshold)]] <- subG 
  
  # Compute the clique counts among the sub-graph 
  # to update vertices for removal
  clique_sample_matrix <- get_vertex_clique_matrix(subG)
  if(!is.null(clique_sample_matrix)) {
    clique_count_per_sample <- sort(colSums(clique_sample_matrix), decreasing = T)
    sids_to_maybe_remove <- names(which(clique_count_per_sample > max_acceptable_cliques))
     if(length(sids_to_maybe_remove) == 0) break()
    marker_counts_of_sids_to_maybe_remove <- sort(unique(marker_counts[sids_to_maybe_remove]))
    threshold <- marker_counts_of_sids_to_maybe_remove[1]
  } else break()
}

# Extract thresholds
thresholds <- names(sub_graphs)

# Plot clique versus marker count for inputG and extract clique counts
clique_counts <- list()
for(threshold in thresholds){
  
  subG <- sub_graphs[[as.character(threshold)]]  
  clique_sample_matrix <- get_vertex_clique_matrix(subG)
  clique_count_per_sample <- sort(colSums(clique_sample_matrix), decreasing = T)
  clique_counts[[as.character(threshold)]] <- sum(apply(clique_sample_matrix, 1, any))
   
  # Plot clique count versus marker count
  plot(y = as.numeric(clique_count_per_sample), 
       x = marker_counts[names(clique_count_per_sample)], 
       ylim = c(0, max_clique_count_per_sample), 
       xlim = c(0, max_marker_count), 
       ylab = "Per-sample clique membership count", 
       xlab = "Per-sample marker data count",
       bty = "n", pch = 20)
}

# Get samples removed on each step
removed_samples_cum <- lapply(sub_graphs, function(x) sids[!sids %in% V(x)$name]) 
removed_samples <- c(removed_samples_cum[1], 
                     sapply(2:length(removed_samples_cum), function(i){
                       setdiff(removed_samples_cum[[i]], removed_samples_cum[[i-1]])
                     }))
names(removed_samples) <- thresholds


# Extract articulation points and component counts
articulation_pts <- list()
compnt_counts <- list()
for(threshold in thresholds){
  subG <- sub_graphs[[as.character(threshold)]]  
  compnt_counts[[as.character(threshold)]] <- sum(components(subG)$csize > 1) # Components greater than one
  articulation_pts[[as.character(threshold)]] <- articulation.points(subG) # Articulation points
}

# Compare overlap between articulation points and samples removed
for(threshold in thresholds){
  print(any(names(articulation_pts[[threshold]]) %in% removed_samples[[threshold]]))
}


# Use plots of clonal component counts and articulation 
# points to select the best threshold for filtering 
par(mfrow = c(1,1))
plot(y = unlist(clique_counts), 
     x = as.numeric(thresholds), type = "l", 
     ylab = "Clique / component count", 
     xlab = "Number of markers with data per-sample")
lines(y = unlist(compnt_counts), 
      x = as.numeric(thresholds), type = "l", lty = "dashed")
abline(v = as.numeric(names(which.max(compnt_counts))))
plot(y = sapply(articulation_pts, length), 
     x = as.numeric(thresholds), type = "l")
plot(y = sapply(removed_samples_cum, length), 
     x = as.numeric(thresholds), type = "l", 
     ylab = "Cumulative number of samples removed")

# Visualise graph and vertex removed 
# To visually discern vertices that are flagged for removal
vertex_removal_ind <- V(inputG)$name %in% removed_samples[[thresholds[1]]]
vertex_articul_ind <- V(inputG)$name %in% names(articulation_pts[[thresholds[1]]])
vertex_labels <- V(inputG)$name
vertex_labels[!vertex_removal_ind] <- NA
vertex_labels[vertex_articul_ind] <- "ART"

plot(inputG, 
     vertex.frame.color = NA,
     vertex.label.color = "black", 
     vertex.label.cex = 0.5, 
     vertex.label = vertex_labels, 
     vertex.size = c(3,6)[vertex_removal_ind+1], 
     vertex.shape = c("circle", "square")[vertex_removal_ind+1])

writeLines(sprintf("%s clonal component/s with two or more samples", 
                   sum(components(inputG)$csize > 1))) 

# Legend
legend('bottom', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, inset = 0.01, 
       col = cols_cities[unique(V(inputG)$city)], 
       legend = unique(V(inputG)$city))

for(i in 1:length(thresholds)){ 
  
  threshold <- thresholds[i]
  subG <- sub_graphs[[threshold]]
  vertex_removal_ind <- V(subG)$name %in% removed_samples[[thresholds[i+1]]]
  vertex_articul_ind <- V(inputG)$name %in% names(articulation_pts[[thresholds[i+1]]])
  vertex_labels <- V(subG)$name
  vertex_labels[!vertex_removal_ind] <- NA 
  vertex_labels[vertex_articul_ind] <- "ART"
  
  plot(subG, main = threshold, 
       vertex.label.color = "black", 
       vertex.label.cex = 0.5, 
       vertex.label = vertex_labels, 
       vertex.frame.color = NA,
       vertex.size = c(3,6)[vertex_removal_ind+1], 
       vertex.shape = c("circle", "square")[vertex_removal_ind+1])
  
  legend('bottom', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, inset = 0.01, 
         col = cols_cities[unique(V(subG)$city)], 
         legend = unique(V(subG)$city))
  
}

if(PDF) dev.off()

# # clique overlap (not used for anything)
# clique_names <- paste0("clique", 1:clique_no)
# names(cliques) <- clique_names
# clique_intersects <- list()
# for(i in 1:(clique_no-1)){
#   for(j in (i+1):clique_no){
#     i_ <- clique_names[i]
#     j_ <- clique_names[j]
#     x <- intersect(names(cliques[[i_]]),
#                    names(cliques[[j_]])) # names returns characters (otherwise returns vertex count)
#     if(length(x) > 0) {
#       clique_intersects[[i_]][[j_]] <- x 
#     } else {
#       next()
#     }
#   }
# }



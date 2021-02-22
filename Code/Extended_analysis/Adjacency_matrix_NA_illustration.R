###########################################################################
#' Obsolte, since graph_from_adjacency_matrix understood
#'
#' A script to demonstrate igraph's treatment of NAs In plots of graphs
#' generated using graph_from_adjacency_matrix with weighted = TRUE vertices
#' whose edge weights are NA are connected. Using graph_from_adjacency_matrix
#' with weighted = NULL (default) vertices whose edge weights are NA are
#' unconnected.
###########################################################################
A <- array(sample(0:2, 100, replace = T), dim = c(10,10))
A[lower.tri(A)] <- NA
A[A == 2] <- NA
par(mfrow = c(2,2))

image(A)

image(as.matrix(as_adjacency_matrix(type = "upper", 
  graph = graph_from_adjacency_matrix(A, 
                                      mode='upper', 
                                      diag=F, 
                                      weighted = T))))

image(as.matrix(as_adjacency_matrix(type = "upper", 
  graph = graph_from_adjacency_matrix(A, 
                                      mode='upper', 
                                      diag=F))))

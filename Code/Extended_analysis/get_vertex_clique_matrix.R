#===============================================================================
#' Function written to get a N by M vertex clique membership matrix where N is
#' the number of maximal cliques of two or more vertices in the input G and
#' M is the number of vertices
#===============================================================================
get_vertex_clique_matrix <- function(G){
  
  # Extract vertex names
  vids <- V(G)$name
  
  # Check the G has named vertices
  if (is.null(vids)) stop("Please provide a G with named vertices")
  
  # Find the maximal cliques of two or more samples
  cliques <- max_cliques(G, min = 2) 
  
  # Number of maximal cliques
  clique_no <- length(cliques)
  
  if (clique_no > 0) {
      
    # Name cliques
    clique_names <- paste0("Clique", 1:clique_no)
    names(cliques) <- clique_names
    
    # Extract clique membership per vertex
    vertex_clique_matrix <- sapply(vids, function(vid){
      sapply(cliques, function(clique){
        as.character(vid) %in% names(clique)
      })
    })
    
    return(vertex_clique_matrix)
    
  } else {
    writeLines("There are no cliques, returning NULL...")
    return(NULL)
  }

}
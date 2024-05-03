if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") 
}
BiocManager::install("tradeSeq")

library(magrittr)
library(RColorBrewer)
library(SingleCellExperiment)
library(igraph)
library(edgeR)
library(dplyr)

# Based on the TradeSeq and Monocle3 workflows provided by the creators of the TradeSeq package

# Input: Any Dynwrap wrapped trajectory with root and global pseudotime (Note: an update to add in the future could be to calculate the pseudotime
# for each cell as the shortest distance from the root to the cell's projection onto the milestone network, this would remove the global pseudotime requirement)
#Outputs: TradeSeq association test object for the trajectory
tradeSeq_Dynwrap_Edge_Method <- function(trajectory, dataset) {
  milestone_graph <- graph_from_data_frame(trajectory$milestone_network, directed = FALSE)
  
  # Initialize an empty list to store the indices of the edge ids
  edges <- list()
  filtered_trajectory <- trajectory$progressions %>%
    group_by(cell_id) %>%
    filter(
      # If there is at least one non-zero percentage value, keep only rows with non-zero percentage
      any(percentage != 0) & percentage != 0 |
        # If all percentage values are zero, keep only one random row
        all(percentage == 0) & row_number() == sample(row_number(), 1)
    ) %>%
    ungroup()
  trajectory$progressions <- filtered_trajectory
  
  ordered_cell_ids <- rownames(dataset$counts)
  
  trajectory$progressions <- trajectory$progressions[match(ordered_cell_ids, trajectory$progressions$cell_id), ]
  
  # Iterate over each cell
  for (i in 1:nrow(trajectory$progressions)) {
    # Get edge 
    edge <- c(trajectory$progressions$from[i], trajectory$progressions$to[i])
    
    # Get edge id
    edge_id <- get.edge.ids(milestone_graph, edge, directed = FALSE)
    
    # Store the index of the edge for the current cell
    edges[[i]] <- edge_id
  }
  
  # Convert the list to a matrix with # cells rows and 1 column
  edges <- matrix(unlist(edges), ncol = 1)
  
  root <- trajectory$root_milestone_id

  # Get the other endpoints
  endpoints <- names(which(igraph::degree(milestone_graph) == 1))
  endpoints <- endpoints[!endpoints %in% root]
  
  # For each endpoint
  cellWeights <- matrix(nrow = nrow(dataset$expression), ncol = length(endpoints))
  rownames(cellWeights) <- rownames(dataset$expression)
  colnames(cellWeights) <- endpoints
  for (endpoint in endpoints) {
    # We find the path between the endpoint and the root
    path <- igraph::shortest_paths(milestone_graph, from=root, to=endpoint)$vpath[[1]]
    path <- as.character(path)
    
    if(length(path) > 2) {
      duplicated_vertices <- rep(path[-c(1, length(path))], each = 2)
      
      path_edges <- as.numeric(c(path[1], duplicated_vertices, tail(path, 1)))
    }
    else {
      path_edges <- as.numeric(path)
    }
    # We find the cells that map along that path
    indices <- which(edges %in% get.edge.ids(milestone_graph, path_edges, directed = FALSE))
    idx_rows <- dataset$expression[indices, ]
    weights <- sapply(rownames(dataset$expression), function(cell_id) {
      if (cell_id %in% idx_rows@Dimnames[[1]]) {
        return(1)
      } else {
        return(0)
      }
    })
    
    cellWeights[, endpoint] <- weights
  }
  
  pseudotime <- matrix(nrow = nrow(dataset$expression), ncol = length(endpoints))
  rownames(pseudotime) <- rownames(dataset$expression)
  colnames(pseudotime) <- endpoints
  
  for (lineage in colnames(cellWeights)) {
    pseudotime[, lineage] <- trajectory$pseudotime
  }
  
  sce <- fitGAM(counts = t(dataset$counts),
                pseudotime = pseudotime,
                cellWeights = cellWeights)
}


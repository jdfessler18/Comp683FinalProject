if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") 
}
BiocManager::install("tradeSeq")

library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(igraph)
library(magrittr)
library(edgeR)

# Based on TradeSeq workflows provided in the TradeSeq package

# Input: Any Dynwrap wrapped trajectory
# Outputs: TradeSeq GAM
tradeSeq_Dynwrap_Vertex_Method <- function(trajectory, dataset) {
  # Get the closest vertice for every cell
  
  # Initialize an empty list to store the indices of the closest milestones
  y_to_cells <- list()
  
  # Iterate over each cell
  for (i in 1:nrow(trajectory$dimred)) {
    # Extract the coordinates of the current cell
    cell_coords <- trajectory$dimred[i, ]
    
    # Calculate the Euclidean distance to all milestones
    distances <- apply(trajectory$dimred_milestones, 1, function(milestone_coords) {
      sqrt(sum((cell_coords - milestone_coords)^2))
    })
    
    # Find the index of the closest milestone
    closest_index <- which.min(distances)
    
    # Store the index of the closest milestone for the current cell
    y_to_cells[[i]] <- closest_index
  }
  
  # Convert the list to a matrix with # cells rows and 1 column
  y_to_cells <- matrix(unlist(y_to_cells), ncol = 1)
  
  root <- trajectory$root_milestone_id
  
  milestone_graph <- graph_from_data_frame(trajectory$milestone_network, directed = FALSE)
  # Get the other endpoints
  endpoints <- names(which(igraph::degree(milestone_graph) == 1))
  endpoints <- endpoints[!endpoints %in% root]
  
  # For each endpoint
  cellWeights <- matrix(nrow = nrow(dataset$expression), ncol = length(endpoints))
  rownames(cellWeights) <- rownames(dataset$expression)
  colnames(cellWeights) <- endpoints
  for (endpoint in endpoints) {
   # We find the path between the endpoint and the root
   path <- igraph::shortest_paths(milestone_graph, root, endpoint)$vpath[[1]]
   path <- as.character(path)
   # We find the cells that map along that path
   indices <- which(y_to_cells %in% path )
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
  # ----eval = FALSE-------------------------------------------------------------
  sce <- fitGAM(counts = t(dataset$counts),
                 pseudotime = pseudotime,
                 cellWeights = cellWeights)
}


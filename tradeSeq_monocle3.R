## ----echo = FALSE-------------------------------------------------------------
library(knitr)

## ----warning=F, message=F-----------------------------------------------------
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(monocle3)

# Derived from TradeSeq package Monocle3 workflow
# Edited to work with Dynwrap

tradeSeq_monocle3 <- function(dataset) {
  counts <- t(dataset$counts)
  cell_metadata = data.frame(celltype = dataset$prior_information$groups_id$group_id, row.names=dataset$prior_information$groups_id$cell_id)
  cell_metadata <- cell_metadata[colnames(counts), , drop=FALSE]
  
  cds <- new_cell_data_set(counts, cell_metadata = cell_metadata,
  gene_metadata = data.frame(gene_short_name = rownames(counts),
                                             row.names = rownames(counts)))
  # Run PCA then UMAP on the data
  cds <- preprocess_cds(cds, method = "PCA")
  cds <- reduce_dimension(cds, preprocess_method = "PCA",
                          reduction_method = "UMAP")

  # Running the clustering method. This is necessary to the construct the graph
  cds <- cluster_cells(cds, reduction_method = "UMAP")

  # Construct the graph
  # Note that, for the rest of the code to run, the graph should be fully connected
  cds <- learn_graph(cds, use_partition = FALSE)

  # We find all the cells that are close to the starting point
  cell_ids <- colnames(cds)[cell_metadata$celltype == dataset$prior_information$start_milestones]
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  closest_vertex <- closest_vertex[cell_ids, ]
  closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
  mst <- principal_graph(cds)$UMAP
  root_pr_nodes <- igraph::V(mst)$name[closest_vertex]

  # We compute the trajectory
  cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
  plot_cells(cds, color_cells_by = "pseudotime")
  
  library(magrittr)
  # Get the closest vertice for every cell
  y_to_cells <-  principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
  y_to_cells$cells <- rownames(y_to_cells)
  y_to_cells$Y <- y_to_cells$V1
  
  # Get the root vertices
  # It is the same node as above
  root <- cds@principal_graph_aux$UMAP$root_pr_nodes
  
  # Get the other endpoints
  endpoints <- names(which(igraph::degree(mst) == 1))
  endpoints <- endpoints[!endpoints %in% root]
  
  # For each endpoint
  cellWeights <- lapply(endpoints, function(endpoint) {
    # We find the path between the endpoint and the root
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    # We find the cells that map along that path
    df <- y_to_cells[y_to_cells$Y %in% path, ]
    df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
    colnames(df) <- endpoint
    return(df)
  }) %>% do.call(what = 'cbind', args = .) %>%  as.matrix()
  rownames(cellWeights) <- colnames(cds)
  pseudotime <- matrix(pseudotime(cds), ncol = ncol(cellWeights),
              nrow = ncol(cds), byrow = FALSE)
  
  sce <- fitGAM(counts = counts,
               pseudotime = pseudotime,
               cellWeights = cellWeights)
}

#tradeSeq_monocle3(data)

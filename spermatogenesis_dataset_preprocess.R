library(tidyverse)
library(dynbenchmark)
library(dplyr)
library(GEOquery)
library(Seurat)

setwd("C:/Users/Jacob/Desktop/Comp683/dynbenchmark")
data <- readRDS("GSE230746_seurat.obj_P15.testis.WT.and.MeiocKO.rds")

# This reference trajectory was derived from https://doi.org/10.1101/2023.09.20.557439
milestone_network <- tribble(
  ~from,    ~to,       ~length, ~directed,
  "Undiff", "A1-4",    1,       TRUE, 
  "A1-4",   "In/B",    1,       TRUE,  
  "In/B",   "B G2/M",       1,       TRUE,  
  "B G2/M",      "pL G1",      1,       TRUE,  
  "pL G1",  "pL eS",    1,       TRUE,  
  "pL eS",   "pL lS",      1,       TRUE,  
  "pL lS",     "L",       1,       TRUE,  
  "L",      "Z",       1,       TRUE,  
  "Z",      "P",       1,       TRUE,  
  "L",      "Mut",       1,       TRUE
)

cell_ids <- names(data@active.ident)
milestone_ids <- unlist(data@active.ident)
cell_info <- tibble(
  cell_id = cell_ids,
  milestone_id = as.character(milestone_ids),
  batch = data@meta.data$Batch
)

# Drop cell types not in the reference trajectory
milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

cell_info <- cell_info[cell_info$milestone_id %in% milestone_ids, ]

cell_ids_to_keep <- cell_info$cell_id
cell_ids_in_data <- colnames(data@assays$RNA$counts)

counts <- data@assays$RNA$counts[, cell_ids_to_keep, drop = FALSE]
rm(data)

# Downsample to 50% due to memory requirements
total_cells <- ncol(counts)
cell_indices <- seq_len(total_cells)  # Indices corresponding to all cells
sampled_indices <- sample(cell_indices, size = round(0.5 * total_cells), replace = FALSE)
counts <- counts[, sampled_indices, drop= FALSE]
counts <- as.matrix(counts)
counts <- t(counts)


cell_ids <- rownames(counts)
cell_info <- cell_info[cell_info$cell_id %in% cell_ids, ]

grouping <- cell_info %>% select(cell_id, milestone_id) %>% deframe()

feature_info = tibble(feature_id = colnames(counts))

# This part is derived from the process_raw_dataset function in the Dynbenchmark package

ensembl_human <- org.Hs.eg.db::org.Hs.egENSEMBL2EG %>% as.list() %>% unlist() %>% tibble(ensembl = names(.), entrez = .)
symbol_human <- org.Hs.eg.db::org.Hs.egSYMBOL %>% as.list() %>% unlist() %>% tibble(entrez = names(.), symbol = .)
human_mapper <- full_join(ensembl_human, symbol_human, by = "entrez")

ensembl_mouse <- org.Mm.eg.db::org.Mm.egENSEMBL2EG %>% as.list() %>% unlist() %>% tibble(ensembl = names(.), entrez = .)
symbol_mouse <- org.Mm.eg.db::org.Mm.egSYMBOL %>% as.list() %>% unlist() %>% tibble(entrez = names(.), symbol = .)
mouse_mapper <- full_join(ensembl_mouse, symbol_mouse, by = "entrez")

id_mapper <- bind_rows(human_mapper, mouse_mapper) %>%
  group_by(entrez) %>%
  filter(row_number() == 1)

colnames(counts) <- tibble(gene_id = colnames(counts)) %>%
  left_join(id_mapper, by = c("gene_id" = "ensembl")) %>%
  mutate(gene_id = ifelse(is.na(symbol), gene_id, symbol)) %>%
  pull(gene_id)
filtered <- colnames(counts) %in% names(which(table(colnames(counts)) == 1))
counts <- counts[, filtered]
conversion_out <- lst(counts, filtered)
counts_prefilter <- conversion_out$counts
feature_info <- feature_info %>% filter(conversion_out$filtered) %>% mutate(feature_id = colnames(counts_prefilter))

# normalise and filter expression
norm_out <- dynnormaliser::normalise_filter_counts(
  counts_prefilter
)

normalisation_info <- norm_out$info
expression <- norm_out$expression
counts <- norm_out$counts

cell_info <- cell_info %>% slice(match(cell_ids, cell_id))
grouping <- grouping[cell_ids] 
feature_info <- feature_info %>% slice(match(colnames(counts), feature_id))

milestone_network$directed <- any(milestone_network$directed)

dataset <-
  dynwrap::wrap_data(
    cell_ids = cell_ids,
    cell_info = cell_info,
    normalisation_info = normalisation_info,
    creation_date = Sys.time()
  )  %>%
  dynwrap::add_cluster_graph(
    milestone_network = milestone_network,
    grouping = grouping
  ) %>%
  dynwrap::add_expression(
    counts = counts,
    expression = expression,
    feature_info = feature_info
  ) %>%
  dynwrap::add_prior_information()
  

dataset <- dynwrap::add_prior_information(dataset) %>%
  dynwrap::add_cell_waypoints()

root_milestone_id <- milestone_network$from[[1]]
dataset <- dynwrap::add_root(dataset, root_milestone_id = root_milestone_id, flip_edges = FALSE)
saveRDS(dataset, file = "spermatogenesis_dynwrap.rds")
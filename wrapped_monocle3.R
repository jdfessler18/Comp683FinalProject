# derived from dynverse/ti_monocle3

# Added batch alignment

run_fun <- function(expression, parameters, priors, verbose, seed) {
  
  requireNamespace("igraph")
  requireNamespace("dplyr")
  requireNamespace("tidyverse")
  
  checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))
  
  # construct metadata
  cell_info <- priors$dataset$cell_info
  if (!is.null(cell_info)) {
    cell_info <- cell_info %>%
      mutate(num_genes_expressed = Matrix::rowSums(expression > 0)) %>%
      as.data.frame() %>%
      tibble::column_to_rownames("cell_id")
  }
  feature_info <- priors$dataset$feature_info
  if (!is.null(feature_info)) {
    if (!"gene_short_name" %in% colnames(feature_info)) {
      if ("feature_name" %in% colnames(feature_info)) {
        feature_info <- feature_info %>% mutate(gene_short_name = feature_name)
      } else {
        feature_info <- feature_info %>% mutate(gene_short_name = feature_id)
      }
    }
    
    feature_info <- feature_info %>%
      as.data.frame() %>%
      tibble::column_to_rownames("feature_id")
  }
  
  cds <- monocle3::new_cell_data_set(
    expression_data = Matrix::t(expression),
    cell_metadata = cell_info,
    gene_metadata = feature_info
  )
  
  # perform data preprocessing
  cds <- monocle3::preprocess_cds(
    cds,
    num_dim = parameters$num_dim,
    norm_method = "size_only"
  )
  
  # perform data preprocessing
  if(!is.null(cell_info$batch)) {
    cds <- monocle3::align_cds(
      cds,
      alignment_group = "batch"
  )
  }
  
  # perform dimensionality reduction
  cds <- monocle3::reduce_dimension(
    cds,
    max_components = parameters$max_components,
    reduction_method = parameters$reduction_method
  )
  
  # perform clustering
  cds <- monocle3::cluster_cells(
    cds,
    k = parameters$k,
    cluster_method = parameters$cluster_method
  )
  
  # calculate trajectory
  cds <- monocle3::learn_graph(cds, use_partition = FALSE)
  
  # process monocle3 output
  dimred <- SingleCellExperiment::reducedDim(cds, parameters$reduction_method) %>%
    magrittr::set_colnames(c("comp_1", "comp_2"))
  dimred_milestones <- t(cds@principal_graph_aux$UMAP$dp_mst) %>%
    magrittr::set_colnames(colnames(dimred))
  milestone_network <-
    igraph::as_data_frame(cds@principal_graph$UMAP) %>%
    transmute(
      from,
      to,
      length = sqrt(rowSums((dimred_milestones[from, ] - dimred_milestones[to, ])^2)),
      directed = FALSE
    )
  
  milestone_percentages <-
    cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex %>%
    magrittr::set_colnames("index") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_id") %>%
    transmute(
      cell_id,
      milestone_id = rownames(dimred_milestones)[index],
      percentage = 1
    )
  
  dimred_segment_progressions <-
    milestone_network %>%
    select(from, to) %>%
    mutate(percentage = map(seq_len(n()), ~ c(0, 1))) %>%
    unnest(percentage)
  
  dsp_names <-
    dimred_segment_progressions %>%
    {ifelse(.$percentage == 0, .$from, .$to)}
  dimred_segment_points <- dimred_milestones[dsp_names, , drop = FALSE]
  
  # TIMING: done with method
  checkpoints$method_aftermethod <- as.numeric(Sys.time())
  
  # return output
  output <-
    dynwrap::wrap_data(
      cell_ids = rownames(expression)
    ) %>%
    dynwrap::add_trajectory(
      milestone_ids = rownames(dimred_milestones),
      milestone_network = milestone_network,
      milestone_percentages = milestone_percentages
    ) %>%
    dynwrap::add_dimred(
      dimred = dimred,
      dimred_milestones = dimred_milestones,
      dimred_segment_progressions = dimred_segment_progressions,
      dimred_segment_points = dimred_segment_points
    ) %>%
    dynwrap::simplify_trajectory() %>%
    dynwrap::add_timings(
      timings = checkpoints
    )
}

#' @importFrom dynwrap convert_definition
definition <- dynwrap::convert_definition(yaml::read_yaml("definition.yml"))


#' Infer a trajectory using Monocle 3
#'
#' @eval dynwrap::generate_parameter_documentation(definition)
#'
#' @examples
#' library(dyntoy)
#' library(dynwrap)
#'
#' priors$dataset <- dyntoy::generate_priors$dataset(
#'   num_cells = 99,
#'   num_features = 101,
#'   model = "tree",
#'   normalise = FALSE
#' )
#' model <- dynwrap::infer_trajectory(priors$dataset, ti_monocle3())
#'
#' @importFrom dynwrap create_ti_method_r
#'
#' @export
ti_monocle3 <- dynwrap::create_ti_method_r(
  definition = definition,
  run_fun = run_fun,
  return_function = TRUE
)
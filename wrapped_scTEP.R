library(scTEP)
library(dynwrap)
library(scDHA)
library(SingleCellExperiment)
library(dplyr)
library(dynparam)


definition <- definition(
  method = def_method(
    id = "sctep",
    # name = "scTEP",
    platform = "R",
    url = "https://github.com/cran/scTEP",
    authors = list(
      def_author(
        given = "Yifang",
        family = "Zhang"
      ),
      def_author(
        given = "Duc",
        family = "Tran"
      ),
      def_author(
        given = "Tin",
        family = "Nguyen"
      ),
      def_author(
        given = "Sergiu",
        family = "Dascalu"
      ),
      def_author(
        given = "Frederick",
        family = "Harris Jr."
      )
    )
  ),
  wrapper = def_wrapper(
    input_required = c("expression","dataset"),
    type = "trajectory",
    trajectory_types = "tree",
    topology_inference = "free"
  ),
  parameters = def_parameters(
    character_parameter(id = "organism", default = "hsa", values = c("hsa", "mmu"))
  )
)

determine_start_milestones <- function(milestone_network, milestone_ids, verbose) {
  # Credit: Function derived from Dynbenchmark package
  
  # convert milestone network to an igraph object
  is_directed <- any(milestone_network$directed)
  gr <- igraph::graph_from_data_frame(
    milestone_network,
    directed = is_directed,
    vertices = milestone_ids
  )
  
  # determine starting milestones
  if (verbose) cat("Computing start milestones\n")
  start_milestones <-
    if (is_directed) {
      deg_in <- igraph::degree(gr, mode = "in")
      deg_out <- igraph::degree(gr, mode = "out")
      names(which(deg_in == 0))
    } else {
      deg <- igraph::degree(gr)
      names(which(deg <= 1))
    }
  
  # if no milestones can be determined as start, pick a random one
  if (length(start_milestones) == 0) {
    start_milestones <- sample(milestone_ids, 1)
  }
  
  # if all milestones are start (ie. cyclic), pick a random one
  if (setequal(start_milestones, milestone_ids)) {
    start_milestones <- sample(start_milestones, 1)
  }
  return(start_milestones)
}


run_fun <- function(expression, parameters, priors, verbose) {
  data('genesets')
  
  #find starting milestone(s) and their indices
  starting_milestones <- determine_start_milestones(priors$dataset$milestone_network, priors$dataset$milestone_ids, verbose=TRUE)
  row_indices <- which(priors$dataset$prior_information$groups_id$group_id %in% starting_milestones)
  
  preprocessed = preprocessing(expression)
  
  data_fa = scTEP.fa(preprocessed, genesets, ncores = 2, data_org = parameters$organism)
  
  allCluster = scTEP::clustering(preprocessed, ncores = 2)
  
  scDHA_res <- scDHA(data_fa, do.clus = T, gen_fil = T, ncores = 2)
  
  out = trajectoryinference(preprocessed, row_indices, scDHA_res, allCluster, ncores = 2)
  
  milestone_ids <- union(out$milestone_network$to, out$milestone_network$from)
  dimred <- scDHA_res$latent
  rownames(dimred) <- rownames(expression)
  
  trajectory <- wrap_data(cell_ids = rownames(expression)) %>%
    dynwrap::add_dimred_projection(milestone_ids, out$milestone_network, dimred, out$data_clus_cent[, 1:15], grouping = out$cluster)

  pseudotime <- setNames(out$pseudotime, rownames(expression))
  trajectory <- dynwrap::add_pseudotime(trajectory, pseudotime)
}

ti_scTEP <- create_ti_method_r(definition, run_fun, package_loaded = "dplyr")




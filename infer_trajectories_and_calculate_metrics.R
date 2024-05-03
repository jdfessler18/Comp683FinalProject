#Assume dataset, ti_scTEP and ti_monocle3 are loaded

# model <- initialise_model(backbone_trifurcating(), distance_metric="euclidean")
# feature_network_default(
#   realnet = "regulatorycircuits_10_lymphocytes",
#   damping = 0.01,
#   target_resampling = Inf,
#   max_in_degree = 5
# )

#dataset <- dyngen::generate_dataset(model, format="list")
#dataset <- as_dyno(dataset$model)
#dataset <- add_prior_information(dataset)
dataset <- readRDS("aging-hsc-young_kowalczyk.rds")

dataset$counts <- Matrix::Matrix(dataset$counts, sparse = TRUE)
dataset$expression <- Matrix::Matrix(dataset$expression, sparse = TRUE)

models <- infer_trajectories(dataset, method=list(ti_monocle3(), dynmethods::ti_slingshot(), ti_scTEP(), dynmethods::ti_paga(), dynmethods::ti_paga_tree()), parameters=list(list(), list(), list(organism="mmu"), list(), list()))

dataset <- add_cell_waypoints(dataset)
models$model <- map(models$model, add_cell_waypoints)

metric_ids <- c("him", "F1_milestones", "F1_branches", "correlation")
metrics <- map_dfr(models$model, dyneval::calculate_metrics, dataset = dataset, metrics = metric_ids)

bind_cols(metrics, models) %>% 
  select(method_id, !!metric_ids) %>% 
  gather("metric_id", "metric_value", -method_id) %>% 
  ggplot(aes(method_id, metric_id, fill = metric_value)) +
  geom_tile()
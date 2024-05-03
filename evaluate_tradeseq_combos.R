library(iCOBRA)
library(scales)
library(slingshot)
library(dynwrap)

data <- dyntoy::generate_dataset(model="bifurcating")
counts <- data$counts

falseGenes <- data$tde_overall$feature_id[data$tde_overall$differentially_expressed]
nullGenes <- data$tde_overall$feature_id[!data$tde_overall$differentially_expressed]

truth <- as.data.frame(matrix(rep(0,nrow(counts)), dimnames=list(rownames(counts),"status")))
truth[falseGenes,"status"] <- 1

models <- infer_trajectories(data, method=list(dynmethods::ti_mfa(), ti_scTEP_without_pathways()), parameters=list(list(), list()))

# Add milestone that was labeled
models$model[[1]] <- add_root(models$model[[1]], root_milestone_id="milestone_begin")

sce <- tradeSeq_Dynwrap_Edge_Method(models$model[[1]], data)

mfa_end <- diffEndTest(sce)$pval
mfa_assoc <- associationTest(sce)$pval
mfa_pattern <- patternTest(sce)$pval

dynplot::plot_dimred(models$model[[2]], label_milestones = TRUE)

#Using plot to assess root
models$model[[2]] <- add_root(models$model[[2]], root_milestone_id="2")

sce_scTEP <- tradeSeq_Dynwrap_Edge_Method(models$model[[2]], data)

scTEP_end <- diffEndTest(sce_scTEP)$pval
scTEP_assoc <- associationTest(sce_scTEP)$pval
scTEP_pattern <- patternTest(sce_scTEP)$pval

sce_slingshot <- tradeSeq_slingshot(data)

# Note that slingshot only produced one lineage for this dataset, resulting in only assoc running correctly
#slingshot_end <- diffEndTest(sce_slingshot)$pval
slingshot_assoc <- associationTest(sce_slingshot)$pval
#slingshot_pattern <- patternTest(sce_slingshot)$pval

sce_monocle3 <- tradeSeq_monocle3(data)

monocle3_end <- diffEndTest(sce_monocle3)$pval
monocle3_assoc <- associationTest(sce_monocle3)$pval
monocle3_pattern <- patternTest(sce_monocle3)$pval

pval <- data.frame( slingshot_assoc=slingshot_assoc, mfa_end=mfa_end,mfa_assoc=mfa_assoc, mfa_pattern=mfa_pattern,
                    scTEP_end=scTEP_end, scTEP_assoc=scTEP_assoc, scTEP_pattern=scTEP_pattern,
                    monocle3_end=monocle3_end, monocle3_assoc=monocle3_assoc, monocle3_pattern=monocle3_pattern,
                    row.names=rownames(t(counts)))
cobra <- COBRAData(pval=pval, truth=truth)

cobra <- calculate_adjp(cobra)
cobra <- calculate_performance(cobra, binary_truth = "status", aspects = c("fdrtpr", "fdrtprcurve", "fdrnbr", "fdrnbrcurve", "tpr", "fpr", "roc", "fpc", "overlap", "corr", "scatter", "deviation", "fsrnbr", "fsrnbrcurve"))

cobraplot <- prepare_data_for_plot(cobra, colorscheme = "Dark2",
                                   facetted = TRUE)
plot_tpr(cobraplot)
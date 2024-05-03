## ----echo = FALSE-------------------------------------------------------------
library(knitr)

## ----eval = FALSE-------------------------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
#  }
#  BiocManager::install("tradeSeq")

## ----warning=F, message=F-----------------------------------------------------
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(umap)

# Derived from TradeSeq package Slingshot workflow
# Edited to work with Dynwrap

tradeSeq_slingshot <- function(dataset) {
  # Get dimred of the counts
  rd <- umap(dataset$counts, n_components=2)$layout
  
  # Get grouping of the cells
  cl <- dataset$prior_information$groups_id$group_id
  crv <- slingshot(rd, cl)
  
  pseudotime <- slingPseudotime(crv, na = FALSE)
  cellWeights <- slingCurveWeights(crv)
  sce <- fitGAM(counts = t(dataset$counts), pseudotime = pseudotime, cellWeights = cellWeights)
}
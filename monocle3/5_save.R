library(tradeSeq)
library(slingshot)

sce <- fitGAM(counts = t(expression_matrix),
              pseudotime = pseudotime,
              cellWeights = cellWeights)

#Save the sce object
saveRDS(sce, "sce_mon2tradeseq.rds")


BiocManager::install("HSMMSingleCell")
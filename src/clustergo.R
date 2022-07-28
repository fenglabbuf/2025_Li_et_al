setwd("~/Desktop/transdiff_GRN")

library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)

gene_ids <- read.csv("data/count_matrix.csv", sep = ",")
gene_ids <- data.frame(gene_ids)
gene_ids <- gene_ids[, c(1, 3)]

mdl <- read.csv("results/temp/community.csv")
entrezmdl <- left_join(mdl, gene_ids, by = c("X" = "Approved.symbol"))
mdls <- unique(mdl$community)

for (m in mdls){
  mdlid <- entrezmdl[entrezmdl$community == m, ]$NCBI.gene.ID
  if (length(mdlid) >= 20) {
    result <- enrichGO(mdlid, OrgDb = "org.Hs.eg.db",
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr",
                       qvalueCutoff = 0.05,
                       minGSSize = 10,
                       maxGSSize = 2000,
                       readable = FALSE)
    resultframe <- data.frame(result)
    if (nrow(resultframe) >= 1) {
        rownames(resultframe) <- NULL
    } else {
      resultframe <- data.frame(0)
    }
  } else {
    resultframe <- data.frame(0)
  }
  write.csv(resultframe, row.names = FALSE,
            paste("results/temp/", as.character(m), ".csv", sep = ""))
}
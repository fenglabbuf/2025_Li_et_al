setwd("~/Desktop/transdiff_GRN")

library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)

gene_ids <- read.csv("data/count_matrix.csv", sep = ",")
gene_ids <- data.frame(gene_ids)
gene_ids <- gene_ids[, c(1, 3)]

mdl <- read.csv("results/rank_community_AMp.txt", sep = "\t")
mdl <- unique(data.frame(mdl))
entrezmdl <- left_join(mdl, gene_ids, by = c("geneid" = "Approved.symbol"))
entrezmdl$thmdls <- paste(entrezmdl$run, entrezmdl$nodes,
                          entrezmdl$community, sep = "_")

mdls <- unique(entrezmdl$thmdls)

for (m in mdls){
  mdlid <- entrezmdl[entrezmdl$thmdls == m, ]$NCBI.gene.ID
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
        write.csv(resultframe, row.names = FALSE,
                  paste("results/GO/AMp_", as.character(m), ".csv", sep = "")
                  )
    }
  }
}

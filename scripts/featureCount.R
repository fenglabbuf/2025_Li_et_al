library(DESeq2)
setwd("~/Desktop/transdiff_GRN")

fcm <- read.table("data/count_matrix.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)
names(fcm)[1] <- "geneid"

row.names(fcm) <- fcm$geneid
fcm[1:3] <- NULL
coldata <- data.frame(colnames(fcm))
dds <- DESeqDataSetFromMatrix(countData = fcm, colData = coldata, design = ~ 1)

# transform with rlog method, output log2 fold change matrix
logm <- rlog(dds, blind=TRUE)
write.csv(assay(logm), "data/rlog_matrix.csv")
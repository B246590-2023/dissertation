library(readxl)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)
library(clusterProfiler)
library(biomaRt)
library(AnnotationDbi)
library(enrichplot)

file_path <- "C:/Users/nem/Desktop/nem1.xlsx"
data <- read_excel(file_path, sheet = "Sheet1", col_names = TRUE)
data <- as.data.frame(data)


rownames(data) <- data[[1]]
data <- data[-1]
data[] <- lapply(data, as.numeric)
data[] <- lapply(data, round)
data[] <- lapply(data, as.integer)
data <- na.omit(data)

head(rownames(data))

str(data)
sample_conditions <- c(rep("high", 8), rep("medium", 11), rep("low", 13))

selected_samples <- colnames(data)[sample_conditions != "low"]
selected_data <- data[, selected_samples]
selected_conditions <- sample_conditions[sample_conditions != "low"]

sample_info <- data.frame(
  row.names = selected_samples,
  condition = factor(selected_conditions, levels = c("high", "medium"))
)

dds <- DESeqDataSetFromMatrix(countData = selected_data, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)

res_high_vs_medium <- results(dds, contrast = c("condition", "high", "medium"))

sig_genes <- res_high_vs_medium[which(res_high_vs_medium$padj < 0.1 | abs(res_high_vs_medium$log2FoldChange) > 0.5), ]

gene_symbols <- rownames(sig_genes)


erich_go_BP <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_CC <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_MF <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db, 
  keyType       = "SYMBOL",
  ont           = "MF",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

head(erich_go_BP)
print(erich_go_BP)
View(erich_go_BP)
write.csv(erich_go_BP, "C:/Users/name/Desktop/erich_go_BP.csv")

write_xlsx(erich_go_BP, path = "C:/Users/name/Desktop/erich_go_BP.xlsx")

library(writexl)


head(erich_go_CC)



dotplot(erich_go_BP, showCategory = 10) + ggtitle("Dotplot for GO BP enrichment analysis(nem/high_vs_medium)")
dotplot(erich_go_CC, showCategory = 10) + ggtitle("Dotplot for GO CC enrichment analysis(nem/high_vs_medium)")
dotplot(erich_go_MF, showCategory = 10) + ggtitle("Dotplot for GO MF enrichment analysis(nem/high_vs_medium)")

barplot(erich_go_BP, showCategory = 10) + ggtitle("Barplot for GO BP enrichment analysis(nem/high_vs_medium)")
barplot(erich_go_CC, showCategory = 10) + ggtitle("Barplot for GO CC enrichment analysis(nem/high_vs_medium)")
barplot(erich_go_MF, showCategory = 10) + ggtitle("Barplot for GO MF enrichment analysis(nem/high_vs_medium)")







selected_samples <- colnames(data)[sample_conditions != "medium"]
selected_data <- data[, selected_samples]
selected_conditions <- sample_conditions[sample_conditions != "medium"]

sample_info <- data.frame(
  row.names = selected_samples,
  condition = factor(selected_conditions, levels = c("high", "low"))
)

dds <- DESeqDataSetFromMatrix(countData = selected_data, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)

res_high_vs_low <- results(dds, contrast = c("condition", "high", "low"))


sig_genes <- res_high_vs_low[which(res_high_vs_low$padj < 0.06 | abs(res_high_vs_low$log2FoldChange) > 0.8), ]
gene_symbols <- rownames(sig_genes)



erich_go_BP <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_CC <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_MF <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "MF",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

head(erich_go_BP)

head(erich_go_CC)



dotplot(erich_go_BP, showCategory = 10) + ggtitle("Dotplot for GO BP enrichment analysis(nem/high_vs_low)")
dotplot(erich_go_CC, showCategory = 10) + ggtitle("Dotplot for GO CC enrichment analysis(nem/high_vs_low)")
dotplot(erich_go_MF, showCategory = 10) + ggtitle("Dotplot for GO MF enrichment analysis(nem/high_vs_low)")

barplot(erich_go_BP, showCategory = 10) + ggtitle("Barplot for GO BP enrichment analysis(nem/high_vs_low)")
barplot(erich_go_CC, showCategory = 10) + ggtitle("Barplot for GO CC enrichment analysis(nem/high_vs_low)")
barplot(erich_go_MF, showCategory = 10) + ggtitle("Barplot for GO MF enrichment analysis(nem/high_vs_low)")







selected_samples <- colnames(data)[sample_conditions != "high"]
selected_data <- data[, selected_samples]
selected_conditions <- sample_conditions[sample_conditions != "high"]

sample_info <- data.frame(
  row.names = selected_samples,
  condition = factor(selected_conditions, levels = c("medium", "low"))
)

dds <- DESeqDataSetFromMatrix(countData = selected_data, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)

res_medium_vs_low <- results(dds, contrast = c("condition", "medium", "low"))

sig_genes <- res_medium_vs_low[which(res_medium_vs_low$padj < 0.06 | abs(res_medium_vs_low$log2FoldChange) > 0.8), ]
gene_symbols <- rownames(sig_genes)



erich_go_BP <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_CC <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_MF <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db, 
  keyType       = "SYMBOL",
  ont           = "MF",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

head(erich_go_BP)

head(erich_go_CC)


dotplot(erich_go_BP, showCategory = 10) + ggtitle("Dotplot for GO BP enrichment analysis(nem/medium_vs_low)")
dotplot(erich_go_CC, showCategory = 10) + ggtitle("Dotplot for GO CC enrichment analysis(nem/medium_vs_low)")
dotplot(erich_go_MF, showCategory = 10) + ggtitle("Dotplot for GO MF enrichment analysis(nem/medium_vs_low)")

barplot(erich_go_BP, showCategory = 10) + ggtitle("Barplot for GO BP enrichment analysis(nem/medium_vs_low)")
barplot(erich_go_CC, showCategory = 10) + ggtitle("Barplot for GO CC enrichment analysis(nem/medium_vs_low)")
barplot(erich_go_MF, showCategory = 10) + ggtitle("Barplot for GO MF enrichment analysis(nem/medium_vs_low)")






file_path <- "C:/Users/name/Desktop/stron1.xlsx"
data <- read_excel(file_path, sheet = "Sheet1", col_names = TRUE)
data <- as.data.frame(data)

rownames(data) <- data[[1]]
data <- data[-1]
data[] <- lapply(data, as.numeric)
data[] <- lapply(data, round)
data[] <- lapply(data, as.integer)
data <- na.omit(data)

head(rownames(data))

str(data)
sample_conditions <- c(rep("high", 9), rep("medium", 14), rep("low", 9))


selected_samples <- colnames(data)[sample_conditions != "low"]
selected_data <- data[, selected_samples]
selected_conditions <- sample_conditions[sample_conditions != "low"]

sample_info <- data.frame(
  row.names = selected_samples,
  condition = factor(selected_conditions, levels = c("high", "medium"))
)
dds <- DESeqDataSetFromMatrix(countData = selected_data, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)

res_high_vs_medium <- results(dds, contrast = c("condition", "high", "medium"))


sig_genes <- res_high_vs_medium[which(res_high_vs_medium$padj < 0.06 | abs(res_high_vs_medium$log2FoldChange) > 0.8), ]
gene_symbols <- rownames(sig_genes)


erich_go_BP <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_CC <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db, 
  keyType       = "SYMBOL",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_MF <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db, 
  keyType       = "SYMBOL",
  ont           = "MF",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

head(erich_go_BP)

head(erich_go_CC)



dotplot(erich_go_BP, showCategory = 10) + ggtitle("Dotplot for GO BP enrichment analysis(stron/high_vs_medium)")
dotplot(erich_go_CC, showCategory = 10) + ggtitle("Dotplot for GO CC enrichment analysis(stron/high_vs_medium)")
dotplot(erich_go_MF, showCategory = 10) + ggtitle("Dotplot for GO MF enrichment analysis(stron/high_vs_medium)")

barplot(erich_go_BP, showCategory = 10) + ggtitle("Barplot for GO BP enrichment analysis(stron/high_vs_medium)")
barplot(erich_go_CC, showCategory = 10) + ggtitle("Barplot for GO CC enrichment analysis(stron/high_vs_medium)")
barplot(erich_go_MF, showCategory = 10) + ggtitle("Barplot for GO MF enrichment analysis(stron/high_vs_medium)")







selected_samples <- colnames(data)[sample_conditions != "medium"]
selected_data <- data[, selected_samples]
selected_conditions <- sample_conditions[sample_conditions != "medium"]

sample_info <- data.frame(
  row.names = selected_samples,
  condition = factor(selected_conditions, levels = c("high", "low"))
)
dds <- DESeqDataSetFromMatrix(countData = selected_data, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)

res_high_vs_low <- results(dds, contrast = c("condition", "high", "low"))


sig_genes <- res_high_vs_low[which(res_high_vs_low$padj < 1 | abs(res_high_vs_low$log2FoldChange) > 0.1), ]
gene_symbols <- rownames(sig_genes)



erich_go_BP <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_CC <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)



erich_go_MF <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db, 
  keyType       = "SYMBOL",
  ont           = "MF",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

head(erich_go_BP)

head(erich_go_CC)


dotplot(erich_go_BP, showCategory = 10) + ggtitle("Dotplot for GO BP enrichment analysis(stron/high_vs_low)")
dotplot(erich_go_CC, showCategory = 10) + ggtitle("Dotplot for GO CC enrichment analysis(stron/high_vs_low)")
dotplot(erich_go_MF, showCategory = 10) + ggtitle("Dotplot for GO MF enrichment analysis(stron/high_vs_low)")

barplot(erich_go_BP, showCategory = 10) + ggtitle("Barplot for GO BP enrichment analysis(stron/high_vs_low)")
barplot(erich_go_CC, showCategory = 10) + ggtitle("Barplot for GO CC enrichment analysis(stron/high_vs_low)")
barplot(erich_go_MF, showCategory = 10) + ggtitle("Barplot for GO MF enrichment analysis(stron/high_vs_low)")








selected_samples <- colnames(data)[sample_conditions != "high"]
selected_data <- data[, selected_samples]
selected_conditions <- sample_conditions[sample_conditions != "high"]

sample_info <- data.frame(
  row.names = selected_samples,
  condition = factor(selected_conditions, levels = c("medium", "low"))
)

dds <- DESeqDataSetFromMatrix(countData = selected_data, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)

res_medium_vs_low <- results(dds, contrast = c("condition", "medium", "low"))

sig_genes <- res_medium_vs_low[which(res_medium_vs_low$padj < 0.5 | abs(res_medium_vs_low$log2FoldChange) > 0.5), ]
gene_symbols <- rownames(sig_genes)



erich_go_BP <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_CC <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db, 
  keyType       = "SYMBOL",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

erich_go_MF <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db, 
  keyType       = "SYMBOL",
  ont           = "MF",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

head(erich_go_BP)

head(erich_go_CC)


dotplot(erich_go_BP, showCategory = 10) + ggtitle("Dotplot for GO BP enrichment analysis(stron/medium_vs_low)")
dotplot(erich_go_CC, showCategory = 10) + ggtitle("Dotplot for GO CC enrichment analysis(stron/medium_vs_low)")
dotplot(erich_go_MF, showCategory = 10) + ggtitle("Dotplot for GO MF enrichment analysis(stron/medium_vs_low)")

barplot(erich_go_BP, showCategory = 10) + ggtitle("Barplot for GO BP enrichment analysis(stron/medium_vs_low)")
barplot(erich_go_CC, showCategory = 10) + ggtitle("Barplot for GO CC enrichment analysis(stron/medium_vs_low)")
barplot(erich_go_MF, showCategory = 10) + ggtitle("Barplot for GO MF enrichment analysis(stron/medium_vs_low)")















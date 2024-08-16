library(readxl)
library(DESeq2)
library(ggplot2)
library(ggrepel)

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

res_high_vs_medium_df <- as.data.frame(res_high_vs_medium)
res_high_vs_medium_df$log10_pvalue <- -log10(res_high_vs_medium_df$pvalue)
res_high_vs_medium_df <- na.omit(res_high_vs_medium_df)
res_high_vs_medium_df <- res_high_vs_medium_df[res_high_vs_medium_df$log10_pvalue >= 0, ]


res_high_vs_medium_df$color <- ifelse(res_high_vs_medium_df$padj < 0.05 & abs(res_high_vs_medium_df$log2FoldChange) > 1, 
                                      "p-value and log2 FC", 
                                      ifelse(res_high_vs_medium_df$padj < 0.05, 
                                             "p-value", 
                                             ifelse(abs(res_high_vs_medium_df$log2FoldChange) > 1, 
                                                    "log2 FC", 
                                                    "NS")))
unique(res_high_vs_medium_df$color)
dev.new()


p <- ggplot(res_high_vs_medium_df, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = color), alpha = 0.5) +
  scale_color_manual(values = c("NS" = "black", "log2 FC" = "blue", "p-value" = "red", "p-value and log2 FC" = "purple"))

print(p)

p <- p + geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed")

print(p)

p <- p + labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot (High vs Medium)") +
  theme_minimal() +
  theme(legend.position = "right")

print(p)

p <- p + geom_text_repel(aes(label = ifelse(color != "NS", rownames(res_high_vs_medium_df), "")), 
                         box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')


p <- p + geom_text_repel(data = res_high_vs_medium_df,
                         aes(label = rownames(res_high_vs_medium_df), color = color),
                         box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50',
                         force = 2, max.overlaps = Inf) 

print(p)

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

res_high_vs_medium_df <- as.data.frame(res_high_vs_medium)
res_high_vs_medium_df$log10_pvalue <- -log10(res_high_vs_medium_df$pvalue)
res_high_vs_medium_df <- na.omit(res_high_vs_medium_df)
res_high_vs_medium_df <- res_high_vs_medium_df[res_high_vs_medium_df$log10_pvalue >= 0, ]


res_high_vs_medium_df$color <- ifelse(res_high_vs_medium_df$padj < 0.05 & abs(res_high_vs_medium_df$log2FoldChange) > 1, 
                                      "p-value and log2 FC", 
                                      ifelse(res_high_vs_medium_df$padj < 0.05, 
                                             "p-value", 
                                             ifelse(abs(res_high_vs_medium_df$log2FoldChange) > 1, 
                                                    "log2 FC", 
                                                    "NS")))
unique(res_high_vs_medium_df$color)
dev.new()

p <- ggplot(res_high_vs_medium_df, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = color), alpha = 0.5) +
  scale_color_manual(values = c("NS" = "black", "log2 FC" = "blue", "p-value" = "red", "p-value and log2 FC" = "purple"))

print(p)

p <- p + geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed")

print(p)

p <- p + labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot (High vs Medium)") +
  theme_minimal() +
  theme(legend.position = "right")

print(p)

p <- p + geom_text_repel(aes(label = ifelse(color == "p-value and log2 FC", rownames(res_high_vs_medium_df), "")), 
                         box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps = 100)

print(p)






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

res_high_vs_low_df <- as.data.frame(res_high_vs_low)
res_high_vs_low_df$log10_pvalue <- -log10(res_high_vs_low_df$pvalue)
res_high_vs_low_df <- na.omit(res_high_vs_low_df)
res_high_vs_low_df <- res_high_vs_low_df[res_high_vs_low_df$log10_pvalue >= 0, ]

res_high_vs_low_df$color <- ifelse(res_high_vs_low_df$padj < 0.05 & abs(res_high_vs_low_df$log2FoldChange) > 1, 
                                   "p-value and log2 FC", 
                                   ifelse(res_high_vs_low_df$padj < 0.05, 
                                          "p-value", 
                                          ifelse(abs(res_high_vs_low_df$log2FoldChange) > 1, 
                                                 "log2 FC", 
                                                 "NS")))
unique(res_high_vs_low_df$color)
dev.new()

p <- ggplot(res_high_vs_low_df, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = color), alpha = 0.5) +
  scale_color_manual(values = c("NS" = "black", "log2 FC" = "blue", "p-value" = "red", "p-value and log2 FC" = "purple"))

print(p)

p <- p + geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed")

print(p)

p <- p + labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot (High vs Low)") +
  theme_minimal() +
  theme(legend.position = "right")

print(p)


p <- p + geom_text_repel(aes(label = ifelse(color != "NS", rownames(res_high_vs_low_df), "")), 
                         box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50')
print(p)











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

res_medium_vs_low_df <- as.data.frame(res_medium_vs_low)
res_medium_vs_low_df$log10_pvalue <- -log10(res_medium_vs_low_df$pvalue)
res_medium_vs_low_df <- na.omit(res_medium_vs_low_df)
res_medium_vs_low_df <- res_medium_vs_low_df[res_medium_vs_low_df$log10_pvalue >= 0, ]

res_medium_vs_low_df$color <- ifelse(res_medium_vs_low_df$padj < 0.05 & abs(res_medium_vs_low_df$log2FoldChange) > 1, 
                                     "p-value and log2 FC", 
                                     ifelse(res_medium_vs_low_df$padj < 0.05, 
                                            "p-value", 
                                            ifelse(abs(res_medium_vs_low_df$log2FoldChange) > 1, 
                                                   "log2 FC", 
                                                   "NS")))
dev.new()

p <- ggplot(res_medium_vs_low_df, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = color), alpha = 0.5) +
  scale_color_manual(values = c("NS" = "black", "log2 FC" = "blue", "p-value" = "red", "p-value and log2 FC" = "purple"))

print(p)

p <- p + geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed")

print(p)

p <- p + labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot (Medium vs Low)") +
  theme_minimal() +
  theme(legend.position = "right")

print(p)

p <- p + geom_text_repel(aes(label = ifelse(color != "NS", rownames(res_medium_vs_low_df), "")), 
                         box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')

p <- p + geom_text_repel(aes(label = ifelse(color != "NS", rownames(res_medium_vs_low_df), "")), 
                         box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50')

print(p)



if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("pheatmap")

install.packages("RColorBrewer")

library(RColorBrewer)

library(DESeq2)
library(pheatmap)
sig_genes <- rownames(res_high_vs_low)[which(res_high_vs_low$padj < 0.05)]

normalized_counts <- counts(dds, normalized = TRUE)
sig_norm_counts <- normalized_counts[sig_genes, ]
annotation_col <- data.frame(
  Condition = selected_conditions
)
rownames(annotation_col) <- selected_samples


pheatmap(sig_norm_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",
         annotation_col = annotation_col,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of Significant DEGsmain of Strongyle(high_vs_low)")

pvalue_threshold <- 0.1
log2fc_threshold <- 0.5



sig_genes <- rownames(res_medium_vs_low)[which(res_medium_vs_low$padj < pvalue_threshold & abs(res_medium_vs_low$log2FoldChange) > log2fc_threshold)]




if (length(sig_genes) == 0) {
  stop("No significant genes found. Check the filtering criteria.")
}

normalized_counts <- counts(dds, normalized = TRUE)
sig_norm_counts <- normalized_counts[sig_genes, ]
print(head(sig_norm_counts))





sig_norm_counts <- sig_norm_counts[apply(sig_norm_counts, 1, function(x) all(is.finite(x))), ]


if (length(sig_genes) == 0) {
  stop("No significant genes found. Check the filtering criteria.")
}

annotation_col <- data.frame(
  Condition = selected_conditions
)
rownames(annotation_col) <- selected_samples

pheatmap(sig_norm_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",
         annotation_col = annotation_col,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of Significant DEGsmain of Strongyle(medium_vs_low)")




sig_genes <- rownames(res_high_vs_medium)[which(res_high_vs_medium$padj < 0.05)]

normalized_counts <- counts(dds, normalized = TRUE)
sig_norm_counts <- normalized_counts[sig_genes, ]

annotation_col <- data.frame(
  Condition = selected_conditions
)
rownames(annotation_col) <- selected_samples

pheatmap(sig_norm_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",
         annotation_col = annotation_col,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of Significant DEGsmain of Strongyle(high_vs_medium)")


library(readxl)
library(DESeq2)
library(ggplot2)
library(ggrepel)

file_path <- "C:/Users/name/Desktop/nem1.xlsx"
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

res_high_vs_medium_df <- as.data.frame(res_high_vs_medium)
res_high_vs_medium_df$log10_pvalue <- -log10(res_high_vs_medium_df$pvalue)
res_high_vs_medium_df <- na.omit(res_high_vs_medium_df)
res_high_vs_medium_df <- res_high_vs_medium_df[res_high_vs_medium_df$log10_pvalue >= 0, ]


res_high_vs_medium_df$color <- ifelse(res_high_vs_medium_df$padj < 0.05 & abs(res_high_vs_medium_df$log2FoldChange) > 1, 
                                      "p-value and log2 FC", 
                                      ifelse(res_high_vs_medium_df$padj < 0.05, 
                                             "p-value", 
                                             ifelse(abs(res_high_vs_medium_df$log2FoldChange) > 1, 
                                                    "log2 FC", 
                                                    "NS")))
unique(res_high_vs_medium_df$color)
dev.new()

p <- ggplot(res_high_vs_medium_df, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = color), alpha = 0.5) +
  scale_color_manual(values = c("NS" = "black", "log2 FC" = "blue", "p-value" = "red", "p-value and log2 FC" = "purple"))

print(p)

p <- p + geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed")

print(p)

p <- p + labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot (High vs Medium)") +
  theme_minimal() +
  theme(legend.position = "right")

print(p)

p <- p + geom_text_repel(aes(label = ifelse(color != "NS", rownames(res_high_vs_medium_df), "")), 
                         box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')

print(p)








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

res_high_vs_low_df <- as.data.frame(res_high_vs_low)
res_high_vs_low_df$log10_pvalue <- -log10(res_high_vs_low_df$pvalue)
res_high_vs_low_df <- na.omit(res_high_vs_low_df)
res_high_vs_low_df <- res_high_vs_low_df[res_high_vs_low_df$log10_pvalue >= 0, ]


res_high_vs_low_df$color <- ifelse(res_high_vs_low_df$padj < 0.05 & abs(res_high_vs_low_df$log2FoldChange) > 1, 
                                   "p-value and log2 FC", 
                                   ifelse(res_high_vs_low_df$padj < 0.05, 
                                          "p-value", 
                                          ifelse(abs(res_high_vs_low_df$log2FoldChange) > 1, 
                                                 "log2 FC", 
                                                 "NS")))
unique(res_high_vs_low_df$color)
dev.new()

p <- ggplot(res_high_vs_low_df, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = color), alpha = 0.5) +
  scale_color_manual(values = c("NS" = "black", "log2 FC" = "blue", "p-value" = "red", "p-value and log2 FC" = "purple"))

print(p)

p <- p + geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed")

print(p)

p <- p + labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot (High vs Low)") +
  theme_minimal() +
  theme(legend.position = "right")

print(p)

p <- p + geom_text_repel(aes(label = ifelse(color != "NS", rownames(res_high_vs_low_df), "")), 
                         box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')

print(p)











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

res_medium_vs_low_df <- as.data.frame(res_medium_vs_low)
res_medium_vs_low_df$log10_pvalue <- -log10(res_medium_vs_low_df$pvalue)
res_medium_vs_low_df <- na.omit(res_medium_vs_low_df)
res_medium_vs_low_df <- res_medium_vs_low_df[res_medium_vs_low_df$log10_pvalue >= 0, ]

res_medium_vs_low_df$color <- ifelse(res_medium_vs_low_df$padj < 0.05 & abs(res_medium_vs_low_df$log2FoldChange) > 1, 
                                     "p-value and log2 FC", 
                                     ifelse(res_medium_vs_low_df$padj < 0.05, 
                                            "p-value", 
                                            ifelse(abs(res_medium_vs_low_df$log2FoldChange) > 1, 
                                                   "log2 FC", 
                                                   "NS")))
dev.new()

p <- ggplot(res_medium_vs_low_df, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = color), alpha = 0.5) +
  scale_color_manual(values = c("NS" = "black", "log2 FC" = "blue", "p-value" = "red", "p-value and log2 FC" = "purple"))

print(p)

p <- p + geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed")

print(p)

p <- p + labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot (Medium vs Low)") +
  theme_minimal() +
  theme(legend.position = "right")

print(p)

p <- p + geom_text_repel(aes(label = ifelse(color != "NS", rownames(res_medium_vs_low_df), "")), 
                         box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')

print(p)



if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("pheatmap")

install.packages("RColorBrewer")


library(RColorBrewer)

library(DESeq2)
library(pheatmap)
sig_genes <- rownames(res_high_vs_low)[which(res_high_vs_low$padj < 0.05)]

normalized_counts <- counts(dds, normalized = TRUE)
sig_norm_counts <- normalized_counts[sig_genes, ]

annotation_col <- data.frame(
  Condition = selected_conditions
)
rownames(annotation_col) <- selected_samples



pheatmap(sig_norm_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",
         annotation_col = annotation_col,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of Significant DEGs of Nematodirs(high_vs_low)")

sig_genes <- rownames(res_medium_vs_low)[which(res_medium_vs_low$padj < 0.05)]

normalized_counts <- counts(dds, normalized = TRUE)
sig_norm_counts <- normalized_counts[sig_genes, ]

sig_genes <- rownames(res_medium_vs_low)[which(res_medium_vs_low$padj < 0.05)]

normalized_counts <- counts(dds, normalized = TRUE)
sig_norm_counts <- normalized_counts[sig_genes, ]

annotation_col <- data.frame(
  Condition = selected_conditions
)
rownames(annotation_col) <- selected_samples



pheatmap(sig_norm_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",
         annotation_col = annotation_col,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of Significant DEGs  of Nematodirs(medium_vs_low)")



sig_genes <- rownames(res_high_vs_medium)[which(res_high_vs_medium$padj < 0.05)]

normalized_counts <- counts(dds, normalized = TRUE)
sig_norm_counts <- normalized_counts[sig_genes, ]



annotation_col <- data.frame(
  Condition = selected_conditions
)
rownames(annotation_col) <- selected_samples



pheatmap(sig_norm_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",
         annotation_col = annotation_col,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of Significant DEGs  of Nematodirs(high_vs_medium)")









if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "biomaRt", "clusterProfiler", "AnnotationDbi", "enrichplot"))

library(DESeq2)
library(biomaRt)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)

res_high_vs_medium <- results(dds, contrast = c("condition", "high", "medium"))

sig_genes <- res_high_vs_medium[which(res_high_vs_medium$padj < 0.05 | abs(res_high_vs_medium$log2FoldChange) > 1), ]
gene_symbols <- rownames(sig_genes)

mart <- useMart(biomart = "ensembl", dataset = "oaries_gene_ensembl")

converted_genes <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = gene_symbols,
  mart = mart
)

print("Converted Ensembl Gene IDs:")
print(head(converted_genes))

converted_genes <- converted_genes[!duplicated(converted_genes$external_gene_name), ]

all_gene_annotation <- getBM(
  attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
  mart = mart
)

print("All gene annotations:")
print(head(all_gene_annotation))

annotated_degs <- all_gene_annotation[all_gene_annotation$ensembl_gene_id %in% converted_genes$ensembl_gene_id, ]

print("Annotated DEGs:")
print(head(annotated_degs))


go_data <- annotated_degs[, c("ensembl_gene_id", "go_id")]
colnames(go_data) <- c("ENSEMBL", "GO")

print("GO data:")
print(head(go_data))

degs_with_go <- converted_genes$ensembl_gene_id[converted_genes$ensembl_gene_id %in% go_data$ENSEMBL]


print(paste("Number of DEGs with GO annotations:", length(degs_with_go)))
print("Example of DEGs with GO annotations:")
print(head(degs_with_go))

print("All Ensembl Gene IDs in GO data:")
print(unique(go_data$ENSEMBL))

ego <- enricher(
  gene = degs_with_go,
  TERM2GENE = go_data,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
ego <- enrichGO(
  gene = degs_with_go,
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID",
  ont = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)


print("GO analysis result:")
print(head(ego))


barplot(ego, showCategory = 20)

dotplot(ego, showCategory = 20)

cnetplot(ego, foldChange = res_high_vs_medium$log2FoldChange)

goplot(ego)












library(DESeq2)
library(biomaRt)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)

res_high_vs_medium <- results(dds, contrast = c("condition", "high", "medium"))

sig_genes <- res_high_vs_medium[which(res_high_vs_medium$padj < 0.05 | abs(res_high_vs_medium$log2FoldChange) > 1), ]
sig_genes <- res_high_vs_medium[which(res_high_vs_medium$padj < 0.05 & abs(res_high_vs_medium$log2FoldChange) > 1), ]
sig_genes <- res_high_vs_medium[which(res_high_vs_medium$padj < 0.1 & abs(res_high_vs_medium$log2FoldChange) > 0.5), ]
sig_genes <- res_high_vs_medium[which(res_high_vs_medium$padj < 0.2 & abs(res_high_vs_medium$log2FoldChange) > 0.3), ]
head(sig_genes)
gene_symbols <- rownames(sig_genes)
head(gene_symbols)



library(DESeq2)

desktop_path <- "C:/Users/name/Desktop/sig_genes3.csv" 



write.csv(sig_genes, file = desktop_path, row.names = TRUE)

if (!requireNamespace("openxlsx", quietly = TRUE))
  install.packages("openxlsx")

library(openxlsx)

excel_path <- "C:/Users/name/Desktop/sig_genes.xlsx" 
write.xlsx(sig_genes, file = excel_path, row.names = TRUE)







mart <- useMart(biomart = "ensembl", dataset = "oaries_gene_ensembl")

converted_genes <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = gene_symbols,
  mart = mart
)

print("Converted Ensembl Gene IDs:")
print(head(converted_genes))

converted_genes <- converted_genes[!duplicated(converted_genes$external_gene_name), ]

ensembl_ids <- converted_genes$ensembl_gene_id

print(paste("Number of DEGs with Ensembl IDs:", length(ensembl_ids)))
print("Example of DEGs with Ensembl IDs:")
print(head(ensembl_ids))

all_gene_annotation <- getBM(
  attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
  mart = mart
)

go_data <- all_gene_annotation[, c("ensembl_gene_id", "go_id")]
colnames(go_data) <- c("ENSEMBL", "GO")

print("GO data:")
print(head(go_data))

degs_with_go <- ensembl_ids[ensembl_ids %in% go_data$ENSEMBL]

print(paste("Number of DEGs with GO annotations:", length(degs_with_go)))
print("Example of DEGs with GO annotations:")
print(head(degs_with_go))

ego <- enrichGO(
  gene = degs_with_go,
  OrgDb = go_data,  
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)


print("GO analysis result:")
print(head(ego))



barplot(ego, showCategory = 20)


dotplot(ego, showCategory = 20)

cnetplot(ego, foldChange = res_high_vs_medium$log2FoldChange)


goplot(ego)










library(DESeq2)
library(biomaRt)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(org.Mm.eg.db)

res_high_vs_medium <- results(dds, contrast = c("condition", "high", "medium"))

sig_genes <- res_high_vs_medium[which(res_high_vs_medium$padj < 0.5 | abs(res_high_vs_medium$log2FoldChange) > 0.5), ]
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

head(erich_go_BP)

head(erich_go_CC)



dotplot(erich_go_BP, showCategory = 10) + ggtitle("Dotplot for GO BP enrichment analysis")
dotplot(erich_go_CC, showCategory = 10) + ggtitle("Dotplot for GO CC enrichment analysis")
barplot(erich_go_BP, showCategory = 10) + ggtitle("Barplot for GO BP enrichment analysis")
barplot(erich_go_CC, showCategory = 10) + ggtitle("Barplot for GO CC enrichment analysis")

































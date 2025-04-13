# Load count matrix (filtered) with genes as rows and samples as columns
dat <- read.table("D:/GitHub/RNA-Seq Analysis of Drug-Resistant PANC1 Cells on Linux/Codes/counts_filtered.txt",
                  header = TRUE, row.names = 1)
head(dat)

# Create metadata with sample IDs and treatment groups
metadata <- data.frame(
  SampleID = c("SRR30861169", "SRR30861170", "SRR30861171", "SRR30861166", "SRR30861167", "SRR30861168"),
  Treatment = c("no treatment", "no treatment", "no treatment", "treatment ", "treatment ", "treatment "))
rownames(metadata) <- metadata$SampleID
metadata$SampleID <- NULL  # Remove the SampleID column, now it’s in rownames

# Load DESeq2 for differential expression analysis
library(DESeq2)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = metadata,
                              design = ~ Treatment)
dds

# Set "no treatment" as the reference level
dds$Treatment <- relevel(dds$Treatment, ref = "no treatment")

# Run differential expression analysis
dds <- DESeq(dds)

# Extract results (log2FoldChange, p-values, adjusted p-values, etc.)
res <- results(dds)

# Remove rows with NA padj values
res <- res[!is.na(res$padj), ]
head(res)

# Summarize number of significantly up/down-regulated genes
summary(res)

# Merge normalized counts with DESeq2 results
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

# Load biomaRt to annotate gene IDs
library(biomaRt)

# Connect to Ensembl and select human gene dataset
human_ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")

# Retrieve gene annotations
gene_ids <- row.names(dat)
annotations <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), 
                     filters = 'ensembl_gene_id', values = gene_ids, mart = human_ensembl)

# Merge annotations with results
resdata_annotated <- merge(resdata, annotations, by.x = "Gene", by.y = "ensembl_gene_id", all.x = TRUE)
head(resdata_annotated)

# Filter significant genes based on adjusted p-value and fold change threshold
sig_genes <- resdata_annotated[resdata_annotated$padj < 0.05 & abs(resdata_annotated$log2FoldChange) > 1, ]
head(sig_genes)
print(nrow(sig_genes))

# PCA analysis for sample visualization
rld <- rlog(dds)  # rlog transformation
plotPCA(rld, intgroup = "Treatment")

# Create sample distance matrix
sampleDists <- dist(t(assay(rld)))  
sampleDistMatrix <- as.matrix(sampleDists)

# Load heatmap visualization packages
library(RColorBrewer)
library(pheatmap)

# Plot sample clustering heatmap
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))

# Plot MA plot of all genes
plotMA(res, ylim = c(-2, 2))

# Shrink log2 fold changes for better visualization
res_shrink <- lfcShrink(dds, coef = 2, type = "apeglm")
plotMA(res_shrink, ylim = c(-2, 2))

# Add significance labels to result table
res$padj[is.na(res$padj)] <- 1
res$significant <- ifelse(res$padj < 0.05, "Significant", "Not significant")
res$largeFoldChange <- ifelse(abs(res$log2FoldChange) > 1, "Large", "Small")

# Basic volcano plot using ggplot2
library(ggplot2)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.8, size = 1.75) + 
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted p-value") +
  theme(plot.title = element_text(hjust = 0.5))

# Enhanced volcano plot with threshold lines
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1.75) +
  scale_color_manual(values = c("Not significant" = "gray", "Significant" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted p-value") +
  theme(plot.title = element_text(hjust = 0.5))

# Load EnhancedVolcano for more advanced volcano plot
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0,
                title = 'Volcano Plot',
                subtitle = 'treatment vs no treatment',
                caption = 'Log2FC and p-value')

# Prepare heatmap for top 20 DEGs
library(tidyverse)
sig_genes_order <- sig_genes[order(sig_genes$padj), ]
top_20_genes <- head(sig_genes_order, 20)

# Select expression columns only
top_genes_HM <- top_20_genes[, c(14, 8, 9, 10, 11, 12, 13)]

# Fix missing gene names
rownames(top_genes_HM) <- NULL
top_genes_HM$external_gene_name[2] <- "UnKnown1"

# Set gene names as row names
top_genes_HM %>% column_to_rownames("external_gene_name") -> heatmap_data

# Basic heatmap
pheatmap(heatmap_data)

# Log2 transform
heatmap_data %>% log2() -> heatmap_data_log2
pheatmap(heatmap_data_log2)

# Subtract row means
heatmap_data_log2 - rowMeans(heatmap_data_log2) -> heatmap_data_meanSubtract
pheatmap(heatmap_data_meanSubtract)

# Calculate Z-scores for each row
heatmap_data_meanSubtract / rowSds(as.matrix(heatmap_data_log2)) -> heatmap_data_zscores
pheatmap(heatmap_data_zscores)

# Alternative way using scale option
pheatmap(heatmap_data_log2, scale = "row")

# Customize heatmap clustering
pheatmap(heatmap_data_meanSubtract, cluster_cols = FALSE)
pheatmap(heatmap_data_meanSubtract, cluster_rows = FALSE)
pheatmap(heatmap_data_meanSubtract, cluster_rows = FALSE, cluster_cols = FALSE)



# Final customized heatmap
pheatmap(heatmap_data_meanSubtract, 
         cluster_cols = FALSE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         fontsize_row = 8,
         fontsize_col = 8,
         treeheight_row = 100,
         legend = FALSE,
         angle_col = 45)


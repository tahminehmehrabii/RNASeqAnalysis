# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "biomaRt"))
library(DESeq2)
library(biomaRt)

# Load count data
dat <- read.table("counts.txt", header = TRUE, row.names = 1)
head(dat)  # Display first few rows

# Create metadata with sample IDs and treatment information
metadata <- data.frame(
  SampleID = c("SRR30861169", "SRR30861170", "SRR30861171", "SRR30861166", "SRR30861167", "SRR30861168"),
  Treatment = factor(c("notreatment", "notreatment", "notreatment", "treatment", "treatment", "treatment"))
)
rownames(metadata) <- metadata$SampleID  # Set row names
metadata$SampleID <- NULL  # Remove redundant column
print(metadata)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = dat, colData = metadata, design = ~Treatment)
dds$Treatment <- relevel(dds$Treatment, ref = "notreatment")  # Set reference level

# Run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)
res <- res[!is.na(res$padj),]  # Remove NA values
head(res)  # Display first few results
summary(res)  # Summarize results

# Filter results with padj < 0.05
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

# Merge results with normalized counts
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), by = "row.names", sort = FALSE)
names(resdata)[1] <- "Gene"  # Rename first column
head(resdata)

# Connect to Ensembl database
human_ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")

# Fetch gene annotations
gene_ids <- resdata$Gene
annotations <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
                     filters = 'ensembl_gene_id', values = gene_ids, mart = human_ensembl)

# Merge annotations with results
resdata_annotated <- merge(resdata, annotations, by.x = "Gene", by.y = "ensembl_gene_id", all.x = TRUE)
head(resdata_annotated)

# Save annotations to file
write.table(annotations, file = "annotations.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

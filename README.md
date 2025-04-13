# RNA-Seq Analysis of Drug-Resistant PANC1 Cells on Linux

Differential gene expression in oxaliplatin-resistant versus untreated PANC1 cells was also analyzed to identify potential biomarkers associated with drug resistance. This RNA-Seq pipeline is designed to process, analyze, and quantify gene expression data efficiently. Most steps were performed using the Linux operating system (as noted in script comments), except for the final stage. The key objectives of this pipeline are as follows:


1: Data Collection:
RNA-Seq analysis begins with the collection of raw data, typically consisting of FASTQ files that contain sequencing reads obtained from RNA samples.

2: Quality Control: 
Use tools like FastQC and MultiQC to assess the quality of raw sequencing data and summarize potential issues such as low-quality reads or adapter contamination.

3: Read Alignment: Map RNA-Seq reads to a reference genome using HISAT2 with high accuracy and summarize alignment statistics using MultiQC to evaluate the efficiency and accuracy of the mapping process.

4: Gene Annotation & Feature Counting: Incorporate genomic annotation data to link sequencing reads to known gene features and use FeatureCounts to quantify gene expression by counting reads mapped to exons and assigning them to corresponding gene IDs.

5: Gene Expression Data Processing: Ensuring high-quality, normalized data for downstream analyses like differential expression and pathway enrichment studies.

# Required software

1. Ubuntu (version 2)
2. java (version 17.0.12)
3. Visual Studio Code (version 1.97.2)
4. FastQC (version 0.12.1)
5. MultiQC (version 1.18)
6. HISAT2 (version 2.2.1)
7. SAMtools (version 1.21)
8. Feature Counts (version 2.1.0)



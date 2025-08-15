if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "pheatmap", "EnhancedVolcano", "apeglm"))

# CRAN packages
install.packages("ggplot2")

#Load packages
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(apeglm)

#Setting working directory
setwd("C:/Users/Robert/Documents/UU/2ndSem/2ndPeriod/Genome Analysis/Project/Project_Lab/Results/03_Transcriptomics/Htseq - Feature Counts/5thCount_flipped_strand=reverse_manymore_options")
#Reads counts of R7 for each sample, without the last 5 lines including the special counters
#R7 - 8538 total features
R7_samp1<-read.delim('R7_SRR24516462_counts.txt', header = FALSE, sep = "")
R7_samp2<-read.delim('R7_SRR24516463_counts.txt', header = FALSE, sep = "")
R7_samp3<-read.delim('R7_SRR24516464_counts.txt', header = FALSE, sep = "")
#HP126 - 7761 total features
HP126_samp1<-read.delim('HP126_SRR24516459_counts.txt', header = FALSE, sep = "")
HP126_samp2<-read.delim('HP126_SRR24516460_counts.txt', header = FALSE, sep = "")
HP126_samp3<-read.delim('HP126_SRR24516461_counts.txt', header = FALSE, sep = "")
#DV3 - 7768 total features
DV3_samp1<-read.delim('DV3_SRR24516456_counts.txt', header = FALSE, sep = "")
DV3_samp2<-read.delim('DV3_SRR24516457_counts.txt', header = FALSE, sep = "")
DV3_samp3<-read.delim('DV3_SRR24516458_counts.txt', header = FALSE, sep = "")

#Special counters 
tail(R7_samp1)
tail(R7_samp2)
tail(R7_samp3)
tail(HP126_samp1)
tail(HP126_samp2)
tail(HP126_samp3)
tail(DV3_samp1)
tail(DV3_samp2)
tail(DV3_samp3)

#Extracting only assigned counts to features
R7_samp1 <- R7_samp1[grep('BOCGMJEH', R7_samp1$V1), ]
R7_samp2 <- R7_samp2[grep('BOCGMJEH', R7_samp2$V1), ]
R7_samp3 <- R7_samp3[grep('BOCGMJEH', R7_samp3$V1), ]
HP126_samp1 <- HP126_samp1[grep('BOCGMJEH', HP126_samp1$V1), ]
HP126_samp2 <- HP126_samp2[grep('BOCGMJEH', HP126_samp2$V1), ]
HP126_samp3 <- HP126_samp3[grep('BOCGMJEH', HP126_samp3$V1), ]
DV3_samp1 <- DV3_samp1[grep('BOCGMJEH', DV3_samp1$V1), ]
DV3_samp2 <- DV3_samp2[grep('BOCGMJEH', DV3_samp2$V1), ]
DV3_samp3 <- DV3_samp3[grep('BOCGMJEH', DV3_samp3$V1), ]

#Total number of assigned read to annotated features (HTSeq) per sample of each strain
print(sum(R7_samp1$V2))
print(sum(R7_samp2$V2))
print(sum(R7_samp3$V2))

print(sum(HP126_samp1$V2))
print(sum(HP126_samp2$V2))
print(sum(HP126_samp3$V2))

print(sum(DV3_samp1$V2))
print(sum(DV3_samp2$V2))
print(sum(DV3_samp3$V2))

#Verifying that both datasets contain the same feature IDs
stopifnot(identical(HP126_samp1$V1, HP126_samp2$V1),
          identical(HP126_samp1$V1, HP126_samp3$V1),
          identical(HP126_samp1$V1, R7_samp1$V1),
          identical(HP126_samp1$V1, R7_samp2$V1),
          identical(HP126_samp1$V1, R7_samp3$V1))

stopifnot(identical(HP126_samp1$V1, HP126_samp2$V1),
          identical(HP126_samp1$V1, HP126_samp3$V1),
          identical(HP126_samp1$V1, DV3_samp1$V1),
          identical(HP126_samp1$V1, DV3_samp2$V1),
          identical(HP126_samp1$V1, DV3_samp3$V1))

###Generating count matrix for the differential gene expression analyses between 
###strains R7 and HP126
every_count_HP126_R7 <- cbind(HP126_samp1$V2, HP126_samp2$V2, HP126_samp3$V2, 
                              R7_samp1$V2, R7_samp2$V2, R7_samp3$V2)

rownames(every_count_HP126_R7) <- R7_samp1$V1
colnames(every_count_HP126_R7) <- c("SRR24516459", "SRR24516460", "SRR24516461",
                                    "SRR24516462", "SRR24516463", "SRR24516464")
every_count_HP126_R7 <- as.matrix(every_count_HP126_R7)
column_data_HP126_R7 <- data.frame(condition = factor(c("HP126", "HP126", "HP126", "R7", "R7", "R7")),
                          row.names=colnames(every_count_HP126_R7))
# Create DESeq2 object
dds_HP126_R7 <- DESeqDataSetFromMatrix(countData = every_count_HP126_R7,
                                       colData = column_data_HP126_R7,
                                       design = ~ condition)

###strains DV3 and HP126
every_count_HP126_DV3 <- cbind(HP126_samp1$V2, HP126_samp2$V2, HP126_samp3$V2,
                               DV3_samp1$V2, DV3_samp2$V2, DV3_samp3$V2)

#Creating the count matrix to compare HP126 with DV3
rownames(every_count_HP126_DV3) <- R7_samp1$V1
colnames(every_count_HP126_DV3) <- c("SRR24516459", "SRR24516460", "SRR24516461",
                                     "SRR24516456", "SRR24516457", "SRR24516458")
every_count_HP126_DV3 <- as.matrix(every_count_HP126_DV3)
column_data_HP126_DV3 <- data.frame(condition = factor(c("HP126", "HP126", "HP126", "DV3", "DV3", "DV3")),
                                   row.names=colnames(every_count_HP126_DV3))

dds_HP126_DV3 <- DESeqDataSetFromMatrix(countData = every_count_HP126_DV3,
                                       colData = column_data_HP126_DV3,
                                       design = ~ condition)

#Verifying matrix dimensions and experimental design set-up
dim((dds_HP126_R7))
dim((dds_HP126_DV3))

head(every_count_HP126_R7,10)
head(every_count_HP126_DV3,10)

levels(dds_HP126_R7$condition)
levels(dds_HP126_DV3$condition)

### Filtering raw counts before differential expression analysis to reduce memory 
### size of the dds objects and increase speed of count modelling with DESeq2. 
### Removing features with counts equal or lower than 10.

dds_HP126_R7 <- dds_HP126_R7[rowSums(counts(dds_HP126_R7)) > 10, ]
dds_HP126_DV3 <- dds_HP126_DV3[rowSums(counts(dds_HP126_DV3)) > 10, ]

### Run DESeq
dds_HP126_R7 <- DESeq(dds_HP126_R7)
dds_HP126_DV3 <- DESeq(dds_HP126_DV3)

# Extract results and shrinking log fold results for better visualization
res_HP126_R7 <- results(dds_HP126_R7)
res_HP126_R7_LFC <- lfcShrink(dds_HP126_R7, coef="condition_R7_vs_HP126", type="apeglm")

res_HP126_DV3 <- results(dds_HP126_DV3)
res_HP126_DV3_LFC <- lfcShrink(dds_HP126_DV3, coef="condition_HP126_vs_DV3", type="apeglm")

# Ordering by adjusted p-value, summarizing the results, and saving to CSV file
resOrdered_HP126_R7 <- res_HP126_R7_LFC[order(res_HP126_R7_LFC$padj), ]
summary(resOrdered_HP126_R7)
write.csv(as.data.frame(resOrdered_HP126_R7), file = "DESeq2_R7_HP126_results.csv")

resOrdered_HP126_DV3 <- res_HP126_DV3_LFC[order(res_HP126_DV3_LFC$padj), ]
summary(resOrdered_HP126_DV3)
write.csv(as.data.frame(resOrdered_HP126_DV3), file = "DESeq2_DV3_HP126_results.csv")

#Metrics of specific genes in the DESeq2 result matrix
resOrdered_HP126_R7["BOCGMJEH_00794", ]

###Filtering the results by most significant features and showing the 10 most
sigGenes_HP126_R7 <- subset(resOrdered_HP126_R7, padj < 0.05)
sigGenes_HP126_DV3 <- subset(resOrdered_HP126_DV3, padj < 0.05)

#Number of genes Downregulated in strain HP126
sigGenes_positive_HP126_R7 <- subset(resOrdered_HP126_R7, padj < 0.05 & log2FoldChange > 0)
#Number of genes Upregulated in strain DV3
sigGenes_positive_HP126_DV3 <- subset(resOrdered_HP126_DV3, padj < 0.05 & log2FoldChange > 0)

#Number of genes Upregulated in strain HP126
sigGenes_negative_HP126_R7 <- subset(resOrdered_HP126_R7, padj < 0.05 & log2FoldChange < 0)
#Number of genes Downregulated in strain DV3
sigGenes_negative_HP126_DV3 <- subset(resOrdered_HP126_DV3, padj < 0.05 & log2FoldChange < 0)

#Number of gene differentially expressed by 2-fold the magnitude
sigGenes_twice_exp_HP126_R7 <- subset(resOrdered_HP126_R7, padj < 0.05 & abs(log2FoldChange) > 1)
sigGenes_twice_exp_HP126_DV3 <- subset(resOrdered_HP126_DV3, padj < 0.05 & abs(log2FoldChange) > 1)

#Top 10 differentially expressed genes from each category
head(sigGenes_HP126_R7, 10)
head(sigGenes_HP126_DV3, 10)
head(sigGenes_twice_exp_HP126_R7, 10)
head(sigGenes_twice_exp_HP126_DV3, 10)

#Total amount of genes driving differential expression, exhibiting statistically
#significant differences in expression levels
nrow(sigGenes_HP126_R7)
nrow(sigGenes_HP126_DV3)

#Proportion of genes driving differential expression, up- and down-regulation.
#A total 8422 CDS were annotated for strain R7
nrow(sigGenes_HP126_R7)/nrow(every_count_HP126_R7)
nrow(sigGenes_HP126_DV3)/nrow(every_count_HP126_DV3)

#Proportion of genes Downregulated in strain HP126
nrow(sigGenes_positive_HP126_R7)/nrow(every_count_HP126_R7)
#Proportion of genes Upregulated in strain DV3
nrow(sigGenes_positive_HP126_DV3)/nrow(every_count_HP126_DV3)
#Proportion of genes Upregulated in strain HP126
nrow(sigGenes_negative_HP126_R7)/nrow(every_count_HP126_R7)
#Proportion of genes Downregulated in strain DV3
nrow(sigGenes_negative_HP126_DV3)/nrow(every_count_HP126_DV3)

# PCA plot
#R7 vs HP126
vsd_HP126_R7 <- vst(dds_HP126_R7, blind = FALSE)
plotPCA(vsd_HP126_R7, intgroup = "condition", pcsToUse = 1:2) + 
  ggplot2::ggtitle("PCA of the log-fold gene expression including three samples per strain") +
  ggplot2::xlim(-40, 40) +
  ggplot2::ylim(-25, 25)

#DV3 vs HP!26
vsd_HP126_DV3 <- vst(dds_HP126_DV3, blind = FALSE)
plotPCA(vsd_HP126_DV3, intgroup = "condition", pcsToUse = 1:2) + 
  ggplot2::ggtitle("PCA of the log-fold gene expression including three samples per strain") +
  ggplot2::xlim(-40, 40) +
  ggplot2::ylim(-25, 25)

#Heatmap
#R7 vs HP126
topVarGenes <- head(order(rowVars(assay(vsd_HP126_R7)), decreasing = TRUE), 50)
# Create heatmap
pheatmap(assay(vsd_HP126_R7)[topVarGenes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
#         labels_row = row.names(res_HP126_R7),
         annotation_col = as.data.frame(colData(dds_HP126_R7))[1],
         main="Top 50 most variable genes")

#DV3 vs HP!26
topVarGenes <- head(order(rowVars(assay(vsd_HP126_DV3)), decreasing = TRUE), 50)
# Create heatmap
pheatmap(assay(vsd_HP126_DV3)[topVarGenes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         #         labels_row = row.names(res_HP126_R7),
         annotation_col = as.data.frame(colData(dds_HP126_DV3))[1],
         main="Top 50 most variable genes")

# Volcano Plot
#R7 vs HP!26
vol_HP126_R7 <- EnhancedVolcano::EnhancedVolcano(res_HP126_R7_LFC,
                                 lab = rownames(res_HP126_R7_LFC),
                                 x = 'log2FoldChange',
                                 y = 'pvalue',
                                 title = 'Differential Expression between R7 and HP126',
                                 pCutoff = 0.05,
                                 FCcutoff = 1)
vol_HP126_R7

#DV3 vs HP!26
vol_HP126_DV3 <- EnhancedVolcano::EnhancedVolcano(res_HP126_DV3_LFC,
                                        lab = rownames(res_HP126_DV3_LFC),
                                        x = 'log2FoldChange',
                                        y = 'pvalue',
                                        title = 'Differential Expression between DV3 and HP126',
                                        pCutoff = 0.05,
                                        FCcutoff = 1)
vol_HP126_DV3

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "pheatmap", "EnhancedVolcano", "apeglm"))

# CRAN packages
install.packages("ggplot2")

library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(apeglm)


setwd("C:/Users/Robert/Documents/UU/2ndSem/2ndPeriod/Genome Analysis/Project+Labs/
      Project_Lab/Results/DESeq2 - Differential Expression")
#Reads counts of R7 for each sample
#8538 CDS
R7_samp1<-read.delim('R7_SRR24516462_counts.txt', header = FALSE, sep = "")
R7_samp2<-read.delim('R7_SRR24516463_counts.txt', header = FALSE, sep = "")
R7_samp3<-read.delim('R7_SRR24516464_counts.txt', header = FALSE, sep = "")

#Extracting only assigned counts to features
R7_samp1 <- R7_samp1[grep('BOCGMJEH', R7_samp1$V1), ]
R7_samp2 <- R7_samp2[grep('BOCGMJEH', R7_samp2$V1), ]
R7_samp3 <- R7_samp3[grep('BOCGMJEH', R7_samp3$V1), ]

#Read counts of HP126
#7649 CDS
HP126_samp1<-read.delim('HP126_SRR24516459_counts.txt', header = FALSE, sep = "")
HP126_samp2<-read.delim('HP126_SRR24516460_counts.txt', header = FALSE, sep = "")
HP126_samp3<-read.delim('HP126_SRR24516461_counts.txt', header = FALSE, sep = "")

#Extracting only assigned counts to features
HP126_samp1 <- HP126_samp1[grep('BOCGMJEH', HP126_samp1$V1), ]
HP126_samp2 <- HP126_samp2[grep('BOCGMJEH', HP126_samp2$V1), ]
HP126_samp3 <- HP126_samp3[grep('BOCGMJEH', HP126_samp3$V1), ]

#Read counts of DV3
#7656 CDS
DV3_samp1<-read.delim('DV3_SRR24516456_counts.txt', header = FALSE, sep = "")
DV3_samp2<-read.delim('DV3_SRR24516457_counts.txt', header = FALSE, sep = "")
DV3_samp3<-read.delim('DV3_SRR24516458_counts.txt', header = FALSE, sep = "")

#Extracting only assigned counts to features
DV3_samp1 <- DV3_samp1[grep('BOCGMJEH', DV3_samp1$V1), ]
DV3_samp2 <- DV3_samp2[grep('BOCGMJEH', DV3_samp2$V1), ]
DV3_samp3 <- DV3_samp3[grep('BOCGMJEH', DV3_samp3$V1), ]

#Total number of assigned counts per sample of each strain
print(count(R7_samp1$V2))
print(count(R7_samp2$V2))
print(count(R7_samp3$V2))

print(count(HP126_samp1$V2))
print(count(HP126_samp2$V2))
print(count(HP126_samp3$V2))

print(count(DV3_samp1$V2))
print(count(DV3_samp2$V2))
print(count(DV3_samp3$V2))

#Verifying that both datasets contain the same number of feature IDs
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

#Count matrix for gene expression comparisons between R7 and HP126
every_count_HP126_R7 <- cbind(HP126_samp1$V2, HP126_samp2$V2, HP126_samp3$V2, 
                              R7_samp1$V2, R7_samp2$V2, R7_samp3$V2)

#Creating the count matrix for comparison HP126 and R7
rownames(every_count_HP126_R7) <- R7_samp1$V1
colnames(every_count_HP126_R7) <- c("SRR24516459", "SRR24516460", "SRR24516461",
                                    "SRR24516462", "SRR24516463", "SRR24516464")
every_count_HP126_R7 <- as.matrix(every_count_HP126_R7)
column_data_HP126_R7 <- data.frame(condition = factor(c("HP126", "HP126", "HP126", "R7", "R7", "R7")),
                          row.names=colnames(every_count_HP126_R7))

#Count matrix for gene expression comparisons between DV3 and HP126
every_count_HP126_DV3 <- cbind(DV3_samp1$V2, DV3_samp2$V2, DV3_samp3$V2,
                               HP126_samp1$V2, HP126_samp2$V2, HP126_samp3$V2)

#Creating the count matrix for comparison HP126 and DV3
rownames(every_count_HP126_DV3) <- R7_samp1$V1
colnames(every_count_HP126_DV3) <- c("SRR24516456", "SRR24516457", "SRR24516458",
                                     "SRR24516459", "SRR24516460", "SRR24516461")
every_count_HP126_DV3 <- as.matrix(every_count_HP126_DV3)
column_data_HP126_DV3 <- data.frame(condition = factor(c("DV3", "DV3", "DV3","HP126", "HP126", "HP126")),
                                   row.names=colnames(every_count_HP126_DV3))
# Create DESeq2 object
dds_HP126_R7 <- DESeqDataSetFromMatrix(countData = every_count_HP126_R7,
                              colData = column_data_HP126_R7,
                              design = ~ condition)
dds_HP126_DV3 <- DESeqDataSetFromMatrix(countData = every_count_HP126_DV3,
                                       colData = column_data_HP126_DV3,
                                       design = ~ condition)
dim((dds_HP126_DV3))
# Filter low count genes to reduce memory size of the dds data object and increase 
#the speed of count modelling within DESeq2
dds_HP126_R7 <- dds_HP126_R7[rowSums(counts(dds_HP126_R7)) > 10, ]
dds_HP126_DV3 <- dds_HP126_DV3[rowSums(counts(dds_HP126_DV3)) > 10, ]

# Run DESeq
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

# PCA plot
#R7 vs HP126
vsd_HP126_R7 <- vst(dds_HP126_R7, blind = FALSE)
plotPCA(vsd_HP126_R7, intgroup = "condition", pcsToUse = 1:2) + 
  ggplot2::ggtitle("PCA of the log-fold gene expression including three samples per strain")

#DV3 vs HP!26
vsd_HP126_DV3 <- vst(dds_HP126_DV3, blind = FALSE)
plotPCA(vsd_HP126_DV3, intgroup = "condition", pcsToUse = 1:2) + 
  ggplot2::ggtitle("PCA of the log-fold gene expression including three samples per strain")
?plotPCA
#Heatmap
#R7 vs HP126
topVarGenes <- head(order(rowVars(assay(vsd_HP126_R7)), decreasing = TRUE), 50)
# Create heatmap
pheatmap(assay(vsd_HP126_R7)[topVarGenes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = as.data.frame(colData(dds_HP126_R7)))

#DV3 vs HP!26
topVarGenes <- head(order(rowVars(assay(vsd_HP126_DV3)), decreasing = TRUE), 50)
# Create heatmap
pheatmap(assay(vsd_HP126_DV3)[topVarGenes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = as.data.frame(colData(dds_HP126_DV3)))

# Volcano Plot
#R7 vs HP!26
vol_HP126_R7 <- EnhancedVolcano::EnhancedVolcano(res_HP126_R7_LFC,
                                 lab = rownames(res_HP126_R7_LFC),
                                 x = 'log2FoldChange',
                                 y = 'pvalue',
                                 title = 'Differential Expression between R7 and HP126',
                                 pCutoff = 0.05,
                                 FCcutoff = 1)
#DV3 vs HP!26
vol_HP126_DV3 <- EnhancedVolcano::EnhancedVolcano(res_HP126_DV3_LFC,
                                        lab = rownames(res_HP126_DV3_LFC),
                                        x = 'log2FoldChange',
                                        y = 'pvalue',
                                        title = 'Differential Expression between DV3 and HP126',
                                        pCutoff = 0.05,
                                        FCcutoff = 1)


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

print(count(R7_samp1$V2))
print(count(R7_samp2$V2))
print(count(R7_samp3$V2))

print(count(HP126_samp1$V2))
print(count(HP126_samp2$V2))
print(count(HP126_samp3$V2))

print(count(DV3_samp1$V2))
print(count(DV3_samp2$V2))
print(count(DV3_samp3$V2))


stopifnot(identical(HP126_samp1$V1, HP126_samp2$V1),
          identical(HP126_samp1$V1, HP126_samp3$V1),
          identical(HP126_samp1$V1, R7_samp1$V1),
          identical(HP126_samp1$V1, R7_samp2$V1),
          identical(HP126_samp1$V1, R7_samp3$V1))

#Count matrix for gene expression compairsons between R7 and HP126
every_count <- cbind(
  HP126_samp1$V2, HP126_samp2$V2, HP126_samp3$V2,
  R7_samp1$V2, R7_samp2$V2, R7_samp3$V2
)

rownames(every_count) <- HP126_samp1$V1
colnames(every_count) <- c("SRR24516459", "SRR24516460", "SRR24516461",
                         "SRR24516462", "SRR24516463", "SRR24516464")
every_count <- as.matrix(every_count)

column_data <- data.frame(condition = factor(c("HP126", "HP126", "HP126", "R7", "R7", "R7")),
                          row.names=colnames(every_count))
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = every_count,
                              colData = column_data,
                              design = ~ condition)


# Filter low count genes to reduce memory size of the dds data object and increase 
#the speed of count modelling within DESeq2
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq
dds <- DESeq(dds)

# Extract results
res <- results(dds)
#Shrinking the log fold results for better visualization
resLFC <- lfcShrink(dds, coef="condition_R7_vs_HP126", type="apeglm")

# Order by adjusted p-value
resOrdered <- resLFC[order(resLFC$padj), ]

# Getting a summary of the results
summary(resLFC)
# Saving the results to a CSV file
write.csv(as.data.frame(resOrdered), file = "DESeq2_R7_HP126_results.csv")

# PCA plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition") + 
  ggplot2::ggtitle("PCA of the log-fold gene expression including three samples per strain")

dev.off()
#Heatmap
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
# Create heatmap
pheatmap(assay(vsd)[topVarGenes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = as.data.frame(colData(dds)))

# Volcano Plot
vol <- EnhancedVolcano::EnhancedVolcano(resLFC,
                                 lab = rownames(resLFC),
                                 x = 'log2FoldChange',
                                 y = 'pvalue',
                                 title = 'Differential Expression between R7 and HP126',
                                 pCutoff = 0.05,
                                 FCcutoff = 1)
png("DESeq2_R7_HP126_results.png", width = 800, height = 600)
dev.off()

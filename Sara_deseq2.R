
## Distance Matrix Heatmap

library(pheatmap)
library(RColorBrewer)

data <- read.table("raw_counts.txt", comment.char = "#", header = TRUE, row.names = 1, sep= "\t")
data <- data[-c(1:5)]
names(data) <- gsub(x=names(data), pattern=("X.home.hugorodrigues.data.sara.corkoak_reads.aligned_reads."), replacement = "")
names(data) <- gsub(x=names(data), pattern=("Aligned.out.bam"), replacement = "")
names(data) <- c("ERR12072747_WD-I","ERR12072748_WD-I", "ERR12072750_WD-P","ERR12072751_WD-P","ERR12072752_WD-P","ERR12072753_WD-X","ERR12072756_WW-I","ERR12072757_WW-I","ERR12072759_WW-P","ERR12072760_WW-P","ERR12072762_WW-X","ERR12072763_WW-X","ERR12072764_WW-X")

sample_distance <- dist(t(data))
dist_mat <- as.matrix(sample_distance)
pheatmap(dist_mat,clustering_distance_rows = sample_distance, clustering_distance_cols = sample_distance)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

pheatmap(dist_mat,clustering_distance_rows = sample_distance, clustering_distance_cols = sample_distance, col=colors)

## Differential Gene Expression

library(DESeq2)
library(tidyverse)

# Pedro stuff

keep <- rowMeans(data) > 10
countdata <- data[keep, ]

coldata <- data.frame(row.names = colnames(countdata), Treatment = c(rep("WD",6),rep("WW",7)), Tissue = c("I","I","P","P","P","X","I","I","P","P","X","X","X"))
coldata$Treatment <-as.factor(coldata$Treatment)
coldata$Tissue <-as.factor(coldata$Tissue)

#coldata$Tissue <- factor(coldata$Tissue, levels = c("Phellem", "InnerBark", "Xylem"))
#coldata$Treatment <- relevel(coldata$Treatment, ref="WW")

ann_colors = list(Tissue = c(Phellem="#E6AB02", InnerBark="#009E73", Xylem ="#1F78B4"),Treatment = c(WW="#BDC3C7", WD="#2C3E50"))

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~Tissue + Tissue:Treatment)
dds

# Principal Component Analysis (PCA) with ggplot2

vsd <- vst(dds, blind = TRUE) 

data <- plotPCA(vsd, intgroup = c( "Treatment", "Tissue"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pcaDATA <- plotPCA(vsd,ntop = 1000,intgroup = c("Treatment","Tissue"),
                   returnData = TRUE)

percentvar <- round(100*attr(pcaDATA,"percentVar"))

pca_plot <- ggplot(pcaDATA,aes(PC1,PC2,colour = Tissue, 
                               shape = Treatment)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ",percentvar[1],"% variance")) +
  ylab(paste0("PC2: ",percentvar[2],"% variance")) +
  scale_color_manual(values=c( "#E6AB02","#009E73", "#1F78B4")) +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        legend.text = element_text(size = 14))

pca_plot

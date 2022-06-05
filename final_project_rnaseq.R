library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(stringr)

# Data process

samples_Mus <- list.files(path = "./quant_salmon", full.names = T, pattern="EV|N2NL")
files_Mus <- file.path(samples_Mus, "quant.sf")
names(files_Mus) <- str_replace(samples_Mus, "./quant_salmon/", "") %>% str_replace("_quant", "")

df_Mus <- data.frame(names(files_Mus), files_Mus)
colnames(df_Mus)
names(df_Mus)[names(df_Mus) == "names.files_Mus."] <- "names"
names(df_Mus)[names(df_Mus) == "files_Mus"] <- "files"

samples_Homo <- list.files(path = "./quant_salmon", full.names = T, pattern="NOTCH2NL")
files_Homo <- file.path(samples_Homo, "quant.sf")
names(files_Homo) <- str_replace(samples_Homo, "./quant_salmon/", "") %>% str_replace("_quant", "")

df_Homo <- data.frame(names(files_Homo), files_Homo)
colnames(df_Homo)
names(df_Homo)[names(df_Homo) == "names.files_Homo."] <- "names"
names(df_Homo)[names(df_Homo) == "files_Homo"] <- "files"

# Create annotation

library(GenomicFeatures)
library(AnnotationDbi)

txdb_Homo = makeTxDbFromGFF('Homo_sapiens.GRCh38.106.gtf.gz')
k_Homo <- keys(txdb_Homo, keytype = "TXNAME")
tx2gene_Homo <- select(txdb_Homo, k_Homo, "GENEID", "TXNAME")

txdb_Mus = makeTxDbFromGFF('Mus_musculus.GRCm39.106.gtf.gz')
k_Mus <- keys(txdb_Mus, keytype = "TXNAME")
tx2gene_Mus <- select(txdb_Mus, k_Mus, "GENEID", "TXNAME")

# Tximport

library(tximport)

txi_Mus <- tximport(files_Mus, type = "salmon", tx2gene = tx2gene_Mus, countsFromAbundance="lengthScaledTPM", ignoreTxVersion=TRUE)
colnames(txi_Mus$counts)[1] <- "EV_ctrl_1"
colnames(txi_Mus$counts)[2] <- "EV_ctrl_2"
colnames(txi_Mus$counts)[3] <- "EV_ctrl_3"
colnames(txi_Mus$counts)[4] <- "N2NL-Sh_1"
colnames(txi_Mus$counts)[5] <- "N2NL-Sh_2"
colnames(txi_Mus$counts)[6] <- "N2NL-Sh_3"

txi_Homo <- tximport(files_Homo, type = "salmon", tx2gene = tx2gene_Homo, countsFromAbundance="lengthScaledTPM", ignoreTxVersion=TRUE)
colnames(txi_Homo$counts)[1] <- "N2NL_KO_1"
colnames(txi_Homo$counts)[2] <- "N2NL_KO_2"
colnames(txi_Homo$counts)[3] <- "N2NL_KO_3"
colnames(txi_Homo$counts)[4] <- "N2NL_WT_1"
colnames(txi_Homo$counts)[5] <- "N2NL_WT_2"
colnames(txi_Homo$counts)[6] <- "N2NL_WT_3"

# Mouse samples with human reference

samples_M_H <- list.files(path = "./quant_m_to_h_salmon", full.names = T, pattern="EV|N2NL")
files_M_H <- file.path(samples_M_H, "quant.sf")
names(files_M_H) <- str_replace(samples_M_H, "./quant_m_to_h_salmon/", "") %>% str_replace("_quant", "")

df_M_H <- data.frame(names(files_M_H), files_M_H)
colnames(df_M_H)
names(df_M_H)[names(df_M_H) == "names.files_M_H."] <- "names"
names(df_M_H)[names(df_M_H) == "files_M_H"] <- "files"

txi_M_H <- tximport(files_M_H, type = "salmon", tx2gene = tx2gene_Homo, countsFromAbundance="lengthScaledTPM", ignoreTxVersion=TRUE)
colnames(txi_M_H$counts)[1] <- "EV_ctrl_1"
colnames(txi_M_H$counts)[2] <- "EV_ctrl_2"
colnames(txi_M_H$counts)[3] <- "EV_ctrl_3"
colnames(txi_M_H$counts)[4] <- "N2NL-Sh_1"
colnames(txi_M_H$counts)[5] <- "N2NL-Sh_2"
colnames(txi_M_H$counts)[6] <- "N2NL-Sh_3"

# DESeq2

coldata_Mus <- data.frame(condition=c(rep("Control",3),rep("N2NL_Sh",3)), 
                          row.names= c("EV_ctrl_1", 
                                       "EV_ctrl_2", 
                                       "EV_ctrl_3", 
                                       "N2NL-Sh_1", 
                                       "N2NL-Sh_2", 
                                       "N2NL-Sh_3"))

coldata_Homo <- data.frame(condition=c(rep("KO",3),rep("Control",3)), 
                           row.names= c("N2NL_KO_1", 
                                        "N2NL_KO_2", 
                                        "N2NL_KO_3", 
                                        "N2NL_WT_1", 
                                        "N2NL_WT_2", 
                                        "N2NL_WT_3"))

coldata_M_H <- data.frame(condition=c(rep("Control",3),rep("N2NL_Sh",3)), 
                          row.names= c("EV_ctrl_1", 
                                       "EV_ctrl_2", 
                                       "EV_ctrl_3", 
                                       "N2NL-Sh_1", 
                                       "N2NL-Sh_2", 
                                       "N2NL-Sh_3"))

dds_Homo <- DESeqDataSetFromTximport(txi_Homo, coldata_Homo, ~condition)
dim(dds_Homo)
# pre-filtering
dim(dds_Homo[rowSums(counts(dds_Homo)) > 5, ])
dds_Homo <- dds_Homo[rowSums(counts(dds_Homo)) > 5, ]
dds_Homo$condition <- relevel(dds_Homo$condition, ref = "Control")
dds_Homo <- estimateSizeFactors(dds_Homo)
sizeFactors(dds_Homo)
# normalize read counts
normalized_counts_Homo <- counts(dds_Homo, normalized=TRUE)
# quality control for Tximport results
rld_Homo <- rlog(dds_Homo, blind=TRUE)
plotPCA(rld_Homo, intgroup="condition")
rld_mat_Homo <- assay(rld_Homo) 
rld_cor_Homo <- cor(rld_mat_Homo)
head(rld_cor_Homo)
pheatmap(rld_cor_Homo, annotation = coldata_Homo)
dds_Homo$condition
# find differential expression
dds_Homo <- DESeq(dds_Homo)
res0.05_Homo <- results(dds_Homo, alpha = 0.05)

View(data.frame(res0.05_Homo))
summary(res0.05_Homo)

dds_Mus <- DESeqDataSetFromTximport(txi_Mus, coldata_Mus, ~condition)
dim(dds_Mus)
# pre-filtering
dim(dds_Mus[rowSums(counts(dds_Mus)) > 5, ])
dds_Mus <- dds_Mus[rowSums(counts(dds_Mus)) > 5, ]
dds_Mus$condition <- relevel(dds_Mus$condition, ref = "Control")
dds_Mus <- estimateSizeFactors(dds_Mus)
sizeFactors(dds_Mus)
# normalize read counts
normalized_counts_Mus <- counts(dds_Mus, normalized=TRUE)
# quality control for Tximport results
rld_Mus <- rlog(dds_Mus, blind=TRUE)
plotPCA(rld_Mus, intgroup="condition")
rld_mat_Mus <- assay(rld_Mus) 
rld_cor_Mus <- cor(rld_mat_Mus)
head(rld_cor_Mus)
pheatmap(rld_cor_Mus, annotation = coldata_Mus)
dds_Mus$condition
# find differential expression
dds_Mus <- DESeq(dds_Mus)
res0.05_Mus <- results(dds_Mus, alpha = 0.05)

View(data.frame(res0.05_Mus))
summary(res0.05_Mus)

dds_M_H <- DESeqDataSetFromTximport(txi_M_H, coldata_M_H, ~condition)
dim(dds_M_H)
# pre-filtering
dim(dds_M_H[rowSums(counts(dds_M_H)) > 5, ])
dds_M_H <- dds_M_H[rowSums(counts(dds_M_H)) > 5, ] 
dds_M_H$condition <- relevel(dds_M_H$condition, ref = "Control")
dds_M_H <- estimateSizeFactors(dds_M_H)
sizeFactors(dds_M_H)
# normalize read counts
normalized_counts_M_H <- counts(dds_M_H, normalized=TRUE)
# quality control for Tximport results
rld_M_H <- rlog(dds_M_H, blind=TRUE)
plotPCA(rld_M_H, intgroup="condition")
rld_mat_M_H <- assay(rld_M_H) 
rld_cor_M_H <- cor(rld_mat_M_H)
head(rld_cor_M_H)
pheatmap(rld_cor_M_H, annotation = coldata_M_H)
dds_M_H$condition
# find differential expression
dds_M_H <- DESeq(dds_M_H)
res0.05_M_H <- results(dds_M_H, alpha = 0.05)

View(data.frame(res0.05_M_H))
summary(res0.05_M_H)

# Plot single gene count
Homo_A <- plotCounts(dds_Homo, gene="ENSG00000264343", 
                     intgroup="condition", returnData=TRUE)
Homo_A %>% View()
ggplot(Homo_A, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(Homo_A))) + 
  theme_bw() +
  ggtitle("NOTCH2NLA") +
  theme(plot.title = element_text(hjust = 0.5))

Homo_B <- plotCounts(dds_Homo, gene="ENSG00000286019", 
                     intgroup="condition", returnData=TRUE)
Homo_B %>% View()
ggplot(Homo_B, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(Homo_B))) + 
  theme_bw() +
  ggtitle("NOTCH2NLB") +
  theme(plot.title = element_text(hjust = 0.5))

Homo_C <- plotCounts(dds_Homo, gene="ENSG00000286019", 
                     intgroup="condition", returnData=TRUE)
Homo_C %>% View()
ggplot(Homo_C, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(Homo_C))) + 
  theme_bw() +
  ggtitle("NOTCH2NLC") +
  theme(plot.title = element_text(hjust = 0.5))


M_H_A <- plotCounts(dds_M_H, gene="ENSG00000264343", 
                    intgroup="condition", returnData=TRUE)
M_H_A %>% View()
ggplot(M_H_A, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(M_H_A))) + 
  theme_bw() +
  ggtitle("NOTCH2NLA") +
  theme(plot.title = element_text(hjust = 0.5))

M_H_B <- plotCounts(dds_M_H, gene="ENSG00000286019", 
                    intgroup="condition", returnData=TRUE)
M_H_B %>% View()
ggplot(M_H_B, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(M_H_B))) + 
  theme_bw() +
  ggtitle("NOTCH2NLB") +
  theme(plot.title = element_text(hjust = 0.5))

M_H_C <- plotCounts(dds_M_H, gene="ENSG00000286019", 
                    intgroup="condition", returnData=TRUE)
M_H_C %>% View()
ggplot(M_H_C, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(M_H_C))) + 
  theme_bw() +
  ggtitle("NOTCH2NLC") +
  theme(plot.title = element_text(hjust = 0.5))

# MA pot
plotMA(res0.05_Homo)
plotMA(res0.05_Mus)
plotMA(res0.05_M_H)

# Quality control using DEGreport
# Size factor QC
counts_Homo <- counts(dds_Homo, normalized = TRUE)
design_Homo <- as.data.frame(colData(dds_Homo))
degCheckFactors(counts_Homo[, 1:6])

counts_Mus <- counts(dds_Mus, normalized = TRUE)
design_Mus <- as.data.frame(colData(dds_Mus))
degCheckFactors(counts_Mus[, 1:6])

counts_M_H <- counts(dds_M_H, normalized = TRUE)
design_M_H <- as.data.frame(colData(dds_M_H))
degCheckFactors(counts_M_H[, 1:6])

# Mean-Variance QC plots
degQC(counts_Homo, colData(dds_Homo)[["condition"]], pvalue = res0.05_Homo[["pvalue"]])
degQC(counts_Mus, colData(dds_Mus)[["condition"]], pvalue = res0.05_Mus[["pvalue"]])
degQC(counts_M_H, colData(dds_M_H)[["condition"]], pvalue = res0.05_M_H[["pvalue"]])

# Covariates correlation with metrics
cor_Homo <- degCorCov(colData(dds_Homo))
cor_Mus <- degCorCov(colData(dds_Mus))
cor_M_H <- degCorCov(colData(dds_M_H))

# Extract significants

padj.cutoff <- 0.05

res_sig_Homo <- res0.05_Homo %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sig_Homo <- res_sig_Homo %>%
  filter(padj < padj.cutoff)
sig_Homo
sig_Homo_500 <- sig_Homo[order(sig_Homo$padj), ]
sig_Homo_500[1:500, ]
write.table(sig_Homo[sig_Homo[, 3] > 0, 1], file="sig_Homo_up.txt", 
            sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
write.table(sig_Homo[sig_Homo[, 3] < 0, 1], file="sig_Homo_down.txt", 
            sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
write.table(sig_Homo[, 1], file="sig_Homo_all.txt", 
            sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
write.table(sig_Homo_500[1:500, 1], file="sig_Homo_500.txt", 
            sep="\t", quote=F, row.names = FALSE, col.names = FALSE)


res_sig_Mus <- res0.05_Mus %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sig_Mus <- res_sig_Mus %>%
  filter(padj < padj.cutoff)
sig_Mus
write.table(sig_Mus[sig_Mus[, 3] > 0, 1], file="sig_Mus_up.txt", 
            sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
write.table(sig_Mus[sig_Mus[, 3] < 0, 1], file="sig_Mus_down.txt", 
            sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
write.table(sig_Mus[, 1], file="sig_Mus_all.txt", 
            sep="\t", quote=F, row.names = FALSE, col.names = FALSE)

# Volcano plots
degVolcano(res0.05_Homo[,c("log2FoldChange", "padj")])
degVolcano(res0.05_Mus[,c("log2FoldChange", "padj")])

par(mfrow=c(1,1))
with(res0.05_Homo, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-10,10)))
with(subset(res0.05_Homo, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))

par(mfrow=c(1,1))
with(res0.05_Mus, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot"))
with(subset(res0.05_Mus, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))

# DGE heatmap
library(reshape2)
library(viridis)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)

Homo_VST <- vst(dds_Homo)
Homo_VST <- assay(Homo_VST)
Homo_VST <- as.data.frame(Homo_VST)
Homo_VST$Gene <- rownames(Homo_VST)
head(Homo_VST)
res0.05_Homo_DF <- as.data.frame(res0.05_Homo)
Homo_sigGenes <- rownames(res0.05_Homo_DF[res0.05_Homo_DF$padj < .05 & abs(res0.05_Homo_DF$log2FoldChange) > 1,])
Homo_VST <- Homo_VST[Homo_VST$Gene %in% Homo_sigGenes,]
Homo_VST <- melt(Homo_VST, id.vars=c("Gene"))
Homo_heatmap <- ggplot(Homo_VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
Homo_heatmap

Homo_VSTMatrix <- dcast(Homo_VST, Gene ~ variable)
rownames(Homo_VSTMatrix) <- Homo_VSTMatrix$Gene
Homo_VSTMatrix$Gene <- NULL
Homo_distanceGene <- dist(Homo_VSTMatrix)
Homo_distanceSample <- dist(t(Homo_VSTMatrix))
Homo_clusterGene <- hclust(Homo_distanceGene, method="average")
Homo_clusterSample <- hclust(Homo_distanceSample, method="average")
Homo_sampleModel <- as.dendrogram(Homo_clusterSample)
Homo_sampleDendrogramData <- segment(dendro_data(Homo_sampleModel, type = "rectangle"))
Homo_sampleDendrogram <- ggplot(Homo_sampleDendrogramData) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  theme_dendro()
Homo_VST$variable <- factor(Homo_VST$variable, levels=Homo_clusterSample$labels[Homo_clusterSample$order])
Homo_heatmap <- ggplot(Homo_VST, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
Homo_heatmap
grid.arrange(Homo_sampleDendrogram, Homo_heatmap, ncol=1, heights=c(1,5))

Homo_sampleDendrogram_1 <- Homo_sampleDendrogram + 
  scale_x_continuous(expand=c(.0085, .0085)) + 
  scale_y_continuous(expand=c(0, 0))
Homo_heatmap_1 <- Homo_heatmap + scale_x_discrete(expand=c(0, 0)) + 
  scale_y_discrete(expand=c(0, 0))
Homo_sampleDendrogramGrob <- ggplotGrob(Homo_sampleDendrogram_1)
Homo_heatmapGrob <- ggplotGrob(Homo_heatmap_1)
Homo_sampleDendrogramGrob$widths
Homo_heatmapGrob$widths
Homo_sampleDendrogramGrob <- gtable_add_cols(Homo_sampleDendrogramGrob, 
                                             Homo_heatmapGrob$widths[7], 6)
Homo_sampleDendrogramGrob <- gtable_add_cols(Homo_sampleDendrogramGrob, 
                                             Homo_heatmapGrob$widths[8], 7)
Homo_maxWidth <- unit.pmax(Homo_sampleDendrogramGrob$widths, Homo_heatmapGrob$widths)
Homo_sampleDendrogramGrob$widths <- as.list(Homo_maxWidth)
Homo_heatmapGrob$widths <- as.list(Homo_maxWidth)
Homo_finalGrob <- arrangeGrob(Homo_sampleDendrogramGrob, Homo_heatmapGrob, 
                              ncol=1, heights=c(2,5))
grid.draw(Homo_finalGrob)


Mus_VST <- vst(dds_Mus)
Mus_VST <- assay(Mus_VST)
Mus_VST <- as.data.frame(Mus_VST)
Mus_VST$Gene <- rownames(Mus_VST)
head(Mus_VST)
res0.05_Mus_DF <- as.data.frame(res0.05_Mus)
Mus_sigGenes <- rownames(res0.05_Mus_DF[res0.05_Mus_DF$padj <= .05 & abs(res0.05_Mus_DF$log2FoldChange) > 1,])
Mus_VST <- Mus_VST[Mus_VST$Gene %in% Mus_sigGenes,]
Mus_VST <- melt(Mus_VST, id.vars=c("Gene"))
Mus_heatmap <- ggplot(Mus_VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
Mus_heatmap

Mus_VSTMatrix <- dcast(Mus_VST, Gene ~ variable)
rownames(Mus_VSTMatrix) <- Mus_VSTMatrix$Gene
Mus_VSTMatrix$Gene <- NULL
Mus_distanceGene <- dist(Mus_VSTMatrix)
Mus_distanceSample <- dist(t(Mus_VSTMatrix))
Mus_clusterGene <- hclust(Mus_distanceGene, method="average")
Mus_clusterSample <- hclust(Mus_distanceSample, method="average")
Mus_sampleModel <- as.dendrogram(Mus_clusterSample)
Mus_sampleDendrogramData <- segment(dendro_data(Mus_sampleModel, type = "rectangle"))
Mus_sampleDendrogram <- ggplot(Mus_sampleDendrogramData) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  theme_dendro()
Mus_VST$variable <- factor(Mus_VST$variable, levels=Mus_clusterSample$labels[Mus_clusterSample$order])
Mus_heatmap <- ggplot(Mus_VST, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
Mus_heatmap
grid.arrange(Mus_sampleDendrogram, Mus_heatmap, ncol=1, heights=c(1,5))

Mus_sampleDendrogram_1 <- Mus_sampleDendrogram + 
  scale_x_continuous(expand=c(.0085, .0085)) + 
  scale_y_continuous(expand=c(0, 0))
Mus_heatmap_1 <- Mus_heatmap + scale_x_discrete(expand=c(0, 0)) + 
  scale_y_discrete(expand=c(0, 0))
Mus_sampleDendrogramGrob <- ggplotGrob(Mus_sampleDendrogram_1)
Mus_heatmapGrob <- ggplotGrob(Mus_heatmap_1)
Mus_sampleDendrogramGrob$widths
Mus_heatmapGrob$widths
Mus_sampleDendrogramGrob <- gtable_add_cols(Mus_sampleDendrogramGrob, 
                                            Mus_heatmapGrob$widths[7], 6)
Mus_sampleDendrogramGrob <- gtable_add_cols(Mus_sampleDendrogramGrob, 
                                            Mus_heatmapGrob$widths[8], 7)
Mus_maxWidth <- unit.pmax(Mus_sampleDendrogramGrob$widths, Mus_heatmapGrob$widths)
Mus_sampleDendrogramGrob$widths <- as.list(Mus_maxWidth)
Mus_heatmapGrob$widths <- as.list(Mus_maxWidth)
Mus_finalGrob <- arrangeGrob(Mus_sampleDendrogramGrob, Mus_heatmapGrob, 
                             ncol=1, heights=c(2,5))
grid.draw(Mus_finalGrob)

# GO results

Homo_go_up <- rbind("ENSG00000108231", 
                    "ENSG00000148704", 
                    "ENSG00000162374", 
                    "ENSG00000149571", 
                    "ENSG00000080824", 
                    "ENSG00000162761", 
                    "ENSG00000101558", 
                    "ENSG00000109062", 
                    "ENSG00000187821", 
                    "ENSG00000196277", 
                    "ENSG00000116106", 
                    "ENSG00000127152", 
                    "ENSG00000156687", 
                    "ENSG00000144355", 
                    "ENSG00000016082", 
                    "ENSG00000167081", 
                    "ENSG00000134198", 
                    "ENSG00000185070")

Homo_go_down <- rbind("ENSG00000055163", 
                      "ENSG00000204580", 
                      "ENSG00000121966", 
                      "ENSG00000083290", 
                      "ENSG00000160888", 
                      "ENSG00000105245", 
                      "ENSG00000168505", 
                      "ENSG00000135541", 
                      "ENSG00000092421", 
                      "ENSG00000104884", 
                      "ENSG00000131094", 
                      "ENSG00000128482", 
                      "ENSG00000074181", 
                      "ENSG00000187122", 
                      "ENSG00000240771", 
                      "ENSG00000183317", 
                      "ENSG00000128342", 
                      "ENSG00000188064", 
                      "ENSG00000168610", 
                      "ENSG00000205336")

Mus_go_up <- rbind("ENSMUSG00000022382", 
                   "ENSMUSG00000042372", 
                   "ENSMUSG00000021743", 
                   "ENSMUSG00000021994", 
                   "ENSMUSG00000020950", 
                   "ENSMUSG00000068748", 
                   "ENSMUSG00000067786")

Homo_go_down_vec <- c("CYFIP2", 
                      "DDR1", 
                      "CXCR4", 
                      "ULK2", 
                      "IER2", 
                      "NUMBL", 
                      "GBX2", 
                      "AHI1", 
                      "SEMA6A", 
                      "ERCC2", 
                      "C1QL1", 
                      "RNF112", 
                      "NOTCH3", 
                      "SLIT1", 
                      "ARHGEF25", 
                      "EPHA10", 
                      "LIF", 
                      "WNT7B", 
                      "STAT3", 
                      "ADGRG1")

Homo_go_up_vec <- c("LGI1", 
                    "VAX1", 
                    "ELAVL4", 
                    "KIRREL3", 
                    "HSP90AA1", 
                    "LMX1A", 
                    "VAPA", 
                    "SLC9A3R1", 
                    "HELT", 
                    "GRM7", 
                    "EPHA4", 
                    "BCL11B", 
                    "UNC5D", 
                    "DLX1", 
                    "ISL1", 
                    "PBX3", 
                    "TSPAN2", 
                    "FLRT2")

Mus_go_up_vec <- c("Wnt7b", 
                   "Dmrt3", 
                   "Fezf2", 
                   "Wnt5a", 
                   "Foxg1", 
                   "Ptprz1", 
                   "Nnat")

# Z-score normalization
Homo_go_VST <- vst(dds_Homo)
Homo_go_VST <- assay(Homo_go_VST)
Homo_go_VST <- t(scale(t(Homo_go_VST)))
Homo_go_VST <- as.data.frame(Homo_go_VST)
Homo_go_VST$Gene <- rownames(Homo_go_VST)
head(Homo_go_VST)
Homo_go_VST <- Homo_go_VST[Homo_go_VST$Gene %in% Homo_go_up,]
Homo_go_VST <- melt(Homo_go_VST, id.vars=c("Gene"))
Homo_go_heatmap <- ggplot(Homo_go_VST, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + scale_fill_viridis() + 
  scale_y_discrete(labels=Homo_go_up_vec)
theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
Homo_go_heatmap

# Z-score normalization
Homo_go_VST2 <- vst(dds_Homo)
Homo_go_VST2 <- assay(Homo_go_VST2)
Homo_go_VST2 <- t(scale(t(Homo_go_VST2)))
Homo_go_VST2 <- as.data.frame(Homo_go_VST2)
Homo_go_VST2$Gene <- rownames(Homo_go_VST2)
head(Homo_go_VST2)
Homo_go_VST2 <- Homo_go_VST2[Homo_go_VST2$Gene %in% Homo_go_down,]
Homo_go_VST2 <- melt(Homo_go_VST2, id.vars=c("Gene"))
Homo_go_heatmap2 <- ggplot(Homo_go_VST2, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + scale_fill_viridis() + 
  scale_y_discrete(labels=Homo_go_down_vec)
theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
Homo_go_heatmap2

# Z-score normalization
Mus_go_VST <- vst(dds_Mus)
Mus_go_VST <- assay(Mus_go_VST)
Mus_go_VST <- t(scale(t(Mus_go_VST)))
Mus_go_VST <- as.data.frame(Mus_go_VST)
Mus_go_VST$Gene <- rownames(Mus_go_VST)
head(Mus_go_VST)
Mus_go_VST <- Mus_go_VST[Mus_go_VST$Gene %in% Mus_go_up,]
Mus_go_VST <- melt(Mus_go_VST, id.vars=c("Gene"))
Mus_go_heatmap <- ggplot(Mus_go_VST, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + scale_fill_viridis() + 
  scale_y_discrete(labels=Mus_go_up_vec) 
theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
Mus_go_heatmap
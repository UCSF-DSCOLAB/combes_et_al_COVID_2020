## Author: Arja Ray

library(assertthat)  # For sanity checks
library(ggplot2)     # For pretty plots
library(cowplot)     # For pretty plots
library(dittoSeq)    # For pretty, colorblind friendly plots
library(dplyr)       # For inline modification of matrices
library(grid)        # For plotting multiple plots in one frame
library(gridExtra)   # For plotting multiple plots in one frame
library(reshape2)    # For "melting" dataframesb
library(scales)      # To access break formatting functions
library(reshape)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(Seurat)


# Post loading the original data with the subsetted "TNK" object

# Looking at data at the original resolution
DimPlot(TNK, reduction = 'umap',group.by = "subset_louvain_res0.4", label = T)

tnk_Markers_0.4 <- read.table(file = "merged_colossal_paper_final_res0.6_TNK_subset_louvain_res0.4_markers.tsv", sep = '\t', header = TRUE,
                              stringsAsFactors = FALSE)


tnk_Markers_0.4_padj0.1 <- tnk_Markers_0.4[which(tnk_Markers_0.4$p_val_adj<0.1),]
tnk_Markers_0.4_padj0.1 <- tnk_Markers_0.4_padj0.1[order(tnk_Markers_0.4_padj0.1$avg_logFC,decreasing = TRUE),]
tnk_Markers_0.4_padj0.1 <- tnk_Markers_0.4_padj0.1[order(tnk_Markers_0.4_padj0.1$cluster,decreasing = FALSE),]

tnk_Markers_0.4_padj0.1_Top5 <- tnk_Markers_0.4_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tnk_Markers_0.4_padj0.1_Top5_unique <- unique(tnk_Markers_0.4_padj0.1_Top5$gene)

pdf("DotPlot_TNK_res0.4.pdf", width = 10, height = 12)
DotPlot(TNK, features= tnk_Markers_0.4_padj0.1_Top5_unique, cols= 'RdYlBu', group.by = "subset_louvain_res0.4",
        assay = "RNA", dot.scale = 5) +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10),
                                               axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
dev.off()

# cleaning up Neut and RBC contaminants
tnk_clean <- subset(TNK, subset = subset_louvain_res0.4 %in% c(0:8, 10))

tnk_clean <- RunUMAP(tnk_clean, dims = 1:30, reduction='harmony',
                     n.neighbors = 30, min.dist = 0.3, spread = 1, verbose = FALSE)

tnk_clean <- FindNeighbors(tnk_clean,  dims = 1:30, k.param = 20,
                           verbose = FALSE, reduction = 'harmony')

tnk_clean_0.8 <- FindClusters(tnk_clean, verbose = TRUE,
                              algorithm = 1, resolution = 0.8, random.seed = 21212)

pdf("UMAP_tnk_clean_res0.8.pdf", width = 10, height = 10)
DimPlot(tnk_clean_0.8, reduction = 'umap', label = T)
dev.off()

tnk_clean_Markers_0.8 <- FindAllMarkers(tnk_clean_0.8, only.pos = TRUE, 
                                        min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson",
                                        assay = "RNA", latent.vars = "LIBRARY")
tnk_clean_Markers_0.8_padj0.1 <- tnk_clean_Markers_0.8[which(tnk_clean_Markers_0.8$p_val_adj<0.1),]
tnk_clean_Markers_0.8_padj0.1 <- tnk_clean_Markers_0.8_padj0.1[order(tnk_clean_Markers_0.8_padj0.1$avg_logFC,decreasing = TRUE),]
tnk_clean_Markers_0.8_padj0.1 <- tnk_clean_Markers_0.8_padj0.1[order(tnk_clean_Markers_0.8_padj0.1$cluster,decreasing = FALSE),]

tnk_clean_Markers_0.8_padj0.1_Top5 <- tnk_clean_Markers_0.8_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tnk_clean_Markers_0.8_padj0.1_Top5_unique <- unique(tnk_clean_Markers_0.8_padj0.1_Top5$gene)

pdf("DotPlot_tnk_clean_res0.6.pdf", width = 10, height = 12)
DotPlot(tnk_clean_0.6, features= tnk_clean_Markers_0.6_padj0.1_Top5_unique, cols= 'RdYlBu',
        assay = "RNA", dot.scale = 5) +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10),
                                               axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
dev.off()


# Adding the original identities at this resolution to the metadata
tnk_clean@meta.data$clean_louvain_res0.8 <- Idents(tnk_clean_0.8)
tnk_clean@meta.data$clean_louvain_res0.6 <- Idents(tnk_clean_0.6)

tnk_clean2 <- subset(tnk_clean, subset = clean_louvain_res0.8 %in% c(0:16))

tnk_clean2 <- RunUMAP(tnk_clean2, dims = 1:30, reduction='harmony',
                      n.neighbors = 30, min.dist = 0.3, spread = 1, verbose = FALSE)

tnk_clean2 <- FindNeighbors(tnk_clean2,  dims = 1:30, k.param = 20,
                            verbose = FALSE, reduction = 'harmony')

tnk_clean2_0.8 <- FindClusters(tnk_clean2, verbose = TRUE,
                               algorithm = 1, resolution = 0.8, random.seed = 21212)

pdf("UMAP_tnk_clean2_res0.8.pdf", width = 10, height = 10)
DimPlot(tnk_clean2_0.8, reduction = 'umap', label = T)
dev.off()

tnk_clean2_Markers_0.8 <- FindAllMarkers(tnk_clean2_0.8, only.pos = TRUE, 
                                         min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson",
                                         assay = "RNA", latent.vars = "LIBRARY")
tnk_clean2_Markers_0.8_padj0.1 <- tnk_clean2_Markers_0.8[which(tnk_clean2_Markers_0.8$p_val_adj<0.1),]
tnk_clean2_Markers_0.8_padj0.1 <- tnk_clean2_Markers_0.8_padj0.1[order(tnk_clean2_Markers_0.8_padj0.1$avg_logFC,decreasing = TRUE),]
tnk_clean2_Markers_0.8_padj0.1 <- tnk_clean2_Markers_0.8_padj0.1[order(tnk_clean2_Markers_0.8_padj0.1$cluster,decreasing = FALSE),]

tnk_clean2_Markers_0.8_padj0.1_Top5 <- tnk_clean2_Markers_0.8_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tnk_clean2_Markers_0.8_padj0.1_Top5_unique <- unique(tnk_clean2_Markers_0.8_padj0.1_Top5$gene)

pdf("DotPlot_tnk_clean2_res0.8.pdf", width = 10, height = 12)
DotPlot(tnk_clean2_0.8, features= tnk_clean_Markers_0.6_padj0.1_Top5_unique, cols= 'RdYlBu',
        assay = "RNA", dot.scale = 5) +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10),
                                               axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
dev.off()


# annotating object for removal and re-Harmony
tnk_clean2@meta.data$clean_louvain_res0.8 <- Idents(tnk_clean2_0.8)
tnk_clean3 <- subset(tnk_clean2, subset = clean_louvain_res0.8 %in% c(0:4, 6:16))

write.table(tnk_clean3@meta.data, file = "TNK_clean_metadata.tsv", sep = '\t', qoute = FALSE)


# Post-Harmony reclean object

load("merged_colossal_paper_final_res0.6_TNK.RData")

tnk_reclean <- RunUMAP(TNK, dims = 1:30, reduction='harmony',
                       n.neighbors = 30, min.dist = 0.3, spread = 1, verbose = FALSE)

tnk_reclean <- FindNeighbors(tnk_reclean,  dims = 1:30, k.param = 20,
                             verbose = FALSE, reduction = 'harmony')

tnk_reclean_0.8 <- FindClusters(tnk_reclean, verbose = TRUE,
                                algorithm = 1, resolution = 0.8, random.seed = 21212)

pdf("UMAP_tnk_reclean_res0.5.pdf", width = 10, height = 10)
DimPlot(tnk_reclean_0.8, reduction = 'umap', label = T)
dev.off()

tnk_reclean_Markers_0.8 <- FindAllMarkers(tnk_reclean_0.8, only.pos = TRUE, 
                                          min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson",
                                          assay = "RNA", latent.vars = "LIBRARY")
tnk_reclean_Markers_0.8_padj0.1 <- tnk_reclean_Markers_0.8[which(tnk_reclean_Markers_0.8$p_val_adj<0.1),]
tnk_reclean_Markers_0.8_padj0.1 <- tnk_reclean_Markers_0.8_padj0.1[order(tnk_reclean_Markers_0.8_padj0.1$avg_logFC,decreasing = TRUE),]
tnk_reclean_Markers_0.8_padj0.1 <- tnk_reclean_Markers_0.8_padj0.1[order(tnk_reclean_Markers_0.8_padj0.1$cluster,decreasing = FALSE),]

tnk_reclean_Markers_0.8_padj0.1_Top5 <- tnk_reclean_Markers_0.8_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tnk_reclean_Markers_0.8_padj0.1_Top5_unique <- unique(tnk_reclean_Markers_0.8_padj0.1_Top5$gene)

pdf("DotPlot_tnk_reclean_res0.8.pdf", width = 10, height = 12)
DotPlot(tnk_reclean_0.8, features= tnk_reclean_Markers_0.8_padj0.1_Top5_unique, cols= 'RdYlBu',
        assay = "RNA", dot.scale = 5) +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10),
                                               axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
dev.off()

# Final Cleanup-removing one platelet cluster

tnk_reclean_0.8@meta.data$clean_louvain_res0.8 <- Idents(tnk_reclean_0.8)
tnk_reclean2 <- subset(tnk_reclean_0.8, subset = clean_louvain_res0.8 %in% c(0:12, 14:16))  # removing one platelet cluster

pdf("UMAP_tnk_reclean2_res0.8.pdf", width = 10, height = 10)
DimPlot(tnk_reclean2, reduction = 'umap', label = T)
dev.off()

tnk_reclean2 <- RunUMAP(tnk_reclean2, dims = 1:30, reduction='harmony',
                        n.neighbors = 30, min.dist = 0.5, spread = 1, verbose = FALSE)

tnk_reclean2@meta.data$reclean_louvain_0.8 <- Idents(tnk_reclean2)

annotations_all <- c("Naive CD4", "CD56hi NK", "Cytotoxic CD8", "Memory CD4", "Memory CD8", "Memory CD4", "CD56lo NK", "Naive CD8", "Treg", 
                     "MAIT", "GammaDelta", "Cycling CD8", "Cycling CD4", "Cytotoxic CD8", "ISG+", "Naive CD4")
names(annotations_all) <- levels(tnk_reclean2@meta.data$reclean_louvain_0.8)
tnk_reclean2 <- RenameIdents(tnk_reclean2, annotations_all)
tnk_reclean2@meta.data$annotated_clusters <- Idents(tnk_reclean2)

DimPlot(tnk_reclean2, reduction = 'umap', group.by = "annotated_clusters", label = T)

pdf("ViolinPlot_tnk_reclean2_Ident_Markers_Original.pdf", width = 10, height = 6)
VlnPlot(tnk_reclean2, features = c("IL7R", "LTB", "CD3E","CD3D", "CD8A","CD8B" ,"PRF1", "GNLY"), group.by = "reclean_louvain_0.8", ncol = 4, pt.size = 0)
dev.off()


write.table(tnk_reclean2@meta.data, "TNK_Reclean2_metadata.tsv", sep = '\t', quote = FALSE)
save (tnk_reclean2, tnk_reclean_Markers_0.8, file = "tnk_cleaned.RData")

mylevels_reclean2 <- c("Naive CD4", "Memory CD4", "Treg", "Cycling CD4", "ISG+", "Naive CD8", 
                       "Memory CD8", "Cytotoxic CD8", "Cycling CD8", "CD56lo NK", "CD56hi NK", "MAIT", "GammaDelta")

mygenelist_reclean2 <- c("LTB", "IL7R", "VIM", "CCR7", "LEF1", "SELL", "TCF7", "JUNB", "FYB",  
                         "AQP3", "PLP2", "S100A4", "CD52", "NOSIP",
                         "FOXP3", "RTKN2", "CD27", "IL32", "IL2RA",
                         "ACTG1", "COTL1", "HIST1H4C", "STMN1", "HMGB2", "TUBB", "MKI67",
                         "ISG15", "MX1", "IFI6", "IFIT3", "MT2A",
                         "CD8A", "CD8B", 
                         "CCL5", "GZMK","DUSP2",
                         "GZMH","FGFBP2","GNLY",
                         "NKG7","TYROBP","SPON2",
                         "NCAM1", "PTGDS", "FCER1G",
                         "KLRB1", "TRAV1-2","DUSP1",
                         "TRDV2", "TRGV9", "TRDC")

tnk_reclean2@meta.data$annotated_clusters <- factor(tnk_reclean2@meta.data$annotated_clusters, levels = mylevels_reclean2)
Idents(tnk_reclean2) <- factor (Idents(tnk_reclean2), levels = mylevels_reclean2)

pdf("DotPlot_tnk_reclean2_reorder_mygenelist.pdf", width = 9, height = 12)
DotPlot(tnk_reclean2, features= mygenelist_reclean2, cols= 'RdYlBu', group.by = "annotated_clusters",
        assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10),
                                axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
dev.off()

COMET_10X_CLINICAL_SCORES = read.csv("COMET_10X_CLINICAL_SCORES_PAPER.csv")

tnk_reclean2@meta.data$Age <- COMET_10X_CLINICAL_SCORES$Age[match(tnk_reclean2@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
tnk_reclean2@meta.data$Day_after_onset <- COMET_10X_CLINICAL_SCORES$Day_after_onset[match(tnk_reclean2@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
tnk_reclean2@meta.data$ICU_vs_FLOOR <- COMET_10X_CLINICAL_SCORES$ICU_vs_FLOOR[match(tnk_reclean2@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
tnk_reclean2@meta.data$Other_infection_type <- COMET_10X_CLINICAL_SCORES$Other_infection_type[match(tnk_reclean2@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
tnk_reclean2@meta.data$Death <- COMET_10X_CLINICAL_SCORES$Death[match(tnk_reclean2@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
tnk_reclean2@meta.data$NIH_score <- COMET_10X_CLINICAL_SCORES$NIH_score[match(tnk_reclean2@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
tnk_reclean2@meta.data$Sampling_score <- COMET_10X_CLINICAL_SCORES$Sampling_score[match(tnk_reclean2@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
tnk_reclean2@meta.data$Overall_score <- COMET_10X_CLINICAL_SCORES$Overall_score[match(tnk_reclean2@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
tnk_reclean2@meta.data$Qualitative_score <- COMET_10X_CLINICAL_SCORES$Qualitative_score[match(tnk_reclean2@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]

# Adding IF-responsive gene signature scores and plotting

ISG_genes = c('MT2A', 'ISG15', 'LY6E', 'IFIT1', 'IFIT2', 'IFIT3', 'IFITM1', 'IFITM3', 'IFI44L', "IFI6", "MX1", "IFI27") #final ISG signature

# Get the cell values from the Seurat object
gene_df <- get_genes_s3(ISG_genes, tnk_reclean2)
# Define the score vector. In this case, we call the score for a cell the mean value of all genes in the vector
score <- colMeans(gene_df)
tnk_reclean2@meta.data$ISG_score <- score[match(colnames(tnk_reclean2), names(score))]

pdf("ViolinPlot_ISG_score_Covid_All.pdf", width = 8, height = 10)
VlnPlot(tnk_reclean2, features = "ISG_score", group.by = "covid_status", pt.size = 0, cols = mycolorscovid)
dev.off()

pdf("ViolinPlot_ISG_score_Severity_All.pdf", width = 8, height = 10)
VlnPlot(tnk_reclean2, features = "ISG_score", group.by = "Qualitative_score", pt.size = 0, cols = mycolorsseverity)
dev.off()

pdf("ViolinPlot_ISG_score_Severity_splitBycovid.pdf", width = 8, height = 10)
VlnPlot(tnk_reclean2, features = "ISG_score", group.by = "covid_status", split.by = "Qualitative_score", pt.size = 0, cols = mycolorsseverity)
dev.off()

pdf("UMAP_tnk_reclean2_Severity.pdf", width = 10, height = 10)
DimPlot(tnk_reclean2, reduction = 'umap', group.by = "Qualitative_score", label = F, cols = mycolorsseverity)
#dittoDimPlot(tnk_reclean, var = 'reclean_louvain_previous', reduction.use = 'umap', do.label = F)
dev.off()

pdf("UMAP_tnk_reclean2_Severity.pdf", width = 10, height = 10)
dittoDimPlot(tnk_reclean2, var = 'Qualitative_score', reduction.use = 'umap', do.label = F, colors = mycolorsseverity)
dev.off()

# annotating for Correlation/phEMD analysis

annotations_final <- c("TNK_Naive_CD4", "TNK_Memory_CD4", "TNK_Treg", "TNK_Cycling_CD4",
                       "TNK_ISG+", "TNK_Naive_CD8", "TNK_Memory_CD8", "TNK_Cytotoxic_CD8",  "TNK_Cycling_CD8",
                       "TNK_CD56lo_NK", "TNK_CD56hi_NK", "TNK_MAIT", "TNK_GammaDelta")
names(annotations_final) <- levels(tnk_reclean2@meta.data$annotated_clusters)
tnk_reclean2@meta.data$harmony_final_clusters <- annotations_final[match(tnk_reclean2@meta.data$annotated_clusters, names(annotations_final))]

write.table(tnk_reclean2@meta.data, "TNK_Final_metadata.tsv", sep = '\t', quote = FALSE)
save (tnk_reclean2, file = "TNK_Final.RData")

## Calculating frequencies across subtypes with Covid status and Severity

# By Covid Status

y <- table(tnk_reclean2@meta.data[,c('SAMPLE.by.SNPs', 'annotated_clusters')])
#melt table into one column by id='SAMPLE.by.SNPs; call 'x'
ClustersByCovidStatus <- melt(y, id='SAMPLE.by.SNPs')
#make data.frame from tnk that contains SAMPLE.by.SNPs event #s; call 'temp'
temp <- table(tnk_reclean2@meta.data$SAMPLE.by.SNPs)
temp <- data.frame(temp)
#add a column to x called 'total' that contains total # of cells per SAMPLE.by.SNPs
ClustersByCovidStatus$total <- temp$Freq[match(ClustersByCovidStatus$SAMPLE.by.SNPs, temp$Var1)]
#add a column to x called 'freq' that gives the fraction of cells of each seurat_clusters per total # of cells in sample
ClustersByCovidStatus$freq <- ClustersByCovidStatus$value/ClustersByCovidStatus$total

#generate temp2 dataframe that has covid_status matched to SAMPLE.by.SNPs
temp2 <- tnk_reclean2@meta.data[,c('SAMPLE.by.SNPs', 'covid_status')]
temp2 <- data.frame(temp2)
#add column to x called 'covid_status' that adds pos/neg/ctrl to x matched to SAMPLE.by.SNPs
ClustersByCovidStatus$covid_status <- temp2$covid_status[match(ClustersByCovidStatus$SAMPLE.by.SNPs, temp2$SAMPLE.by.SNPs)]
#make 'seurat_clusters' numbers factors, not values; or else ggplot gets confused
ClustersByCovidStatus$annotated_clusters <- factor(ClustersByCovidStatus$annotated_clusters)

write.csv(ClustersByCovidStatus, file = "Freq_Covid_Status.csv")

# By Severity

z <- table(tnk_reclean2@meta.data[,c('SAMPLE.by.SNPs', 'annotated_clusters')])
#melt table into one column by id='SAMPLE.by.SNPs; call 'x'
ClustersBySeverity <- melt(z, id='SAMPLE.by.SNPs')
#make data.frame from tnk that contains SAMPLE.by.SNPs event #s; call 'temp'
temp3 <- table(tnk_reclean2@meta.data$SAMPLE.by.SNPs)
temp3 <- data.frame(temp3)
#add a column to x called 'total' that contains total # of cells per SAMPLE.by.SNPs
ClustersBySeverity$total <- temp3$Freq[match(ClustersBySeverity$SAMPLE.by.SNPs, temp3$Var1)]
#add a column to x called 'freq' that gives the fraction of cells of each seurat_clusters per total # of cells in sample
ClustersBySeverity$freq <- ClustersBySeverity$value/ClustersBySeverity$total

#ggplot(ClustersByICU, aes(x=broad_clusters, y=freq))+geom_point(stat = 'identity')
temp4 <- tnk_reclean2@meta.data[,c('SAMPLE.by.SNPs', 'Qualitative_score')]
temp4 <- data.frame(temp4)

ClustersBySeverity$Qualitative_score <- temp4$Qualitative_score[match(ClustersBySeverity$SAMPLE.by.SNPs, temp4$SAMPLE.by.SNPs)]

ClustersBySeverity$annotated_clusters <- factor(ClustersBySeverity$annotated_clusters)

write.csv(ClustersBySeverity, file = "Freq_Severity.csv")


## Calculating ISG scores for different groups

ISG_score_ctrl <- tnk_reclean2@meta.data$ISG_score[tnk_reclean2@meta.data$covid_status == 'ctrl']
ISG_score_neg <- tnk_reclean2@meta.data$ISG_score[tnk_reclean2@meta.data$covid_status == 'neg']
ISG_score_pos <- tnk_reclean2@meta.data$ISG_score[tnk_reclean2@meta.data$covid_status == 'pos']
ISG_score_Covid <- cbind(ISG_score_ctrl, ISG_score_neg, ISG_score_pos)

write.csv(ISG_score_Covid, file = "ISG_score_Covid.csv")

ISG_score_neg_mild <- tnk_reclean2@meta.data$ISG_score[tnk_reclean2@meta.data$covid_status == 'neg' & tnk_reclean2@meta.data$Qualitative_score == 'MILD']
ISG_score_neg_severe <- tnk_reclean2@meta.data$ISG_score[tnk_reclean2@meta.data$covid_status == 'neg' & tnk_reclean2@meta.data$Qualitative_score == 'SEVERE']
ISG_score_pos_mild <- tnk_reclean2@meta.data$ISG_score[tnk_reclean2@meta.data$covid_status == 'pos' & tnk_reclean2@meta.data$Qualitative_score == 'MILD']
ISG_score_pos_severe <- tnk_reclean2@meta.data$ISG_score[tnk_reclean2@meta.data$covid_status == 'pos' & tnk_reclean2@meta.data$Qualitative_score == 'SEVERE']
ISG_score_Severity <- cbind(ISG_score_neg_mild, ISG_score_neg_severe, ISG_score_pos_mild, ISG_score_pos_severe)

write.csv(ISG_score_Severity, file = "ISG_score_Severity.csv")

##--END




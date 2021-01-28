## Author: Tristan Courau

# Call libraries
library(assertthat)
library(ggplot2)
library(cowplot)
library(dittoSeq)
library(dplyr)
library(grid)
library(gridExtra)
library(reshape2)
library(scales)
library(reshape)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(Seurat)
library(beepr)
library(ggpubr)
library(EnhancedVolcano)

# Set working directory
setwd("~/Tristan/Labo MK/Manips/COVID/200428 patient_data_SO/201001 SEURAT WD/")

##################################################################################################
####################################### CLINICAL DATA ############################################
##################################################################################################

# Call clinical scores csv file
COMET_10X_CLINICAL_SCORES <- read.csv("~/Tristan/Labo MK/Manips/COVID/COMET_10X_CLINICAL_SCORES_PAPER.csv", sep=',', header=T)

# Apply values to patients in the neutrophil dataset
merged_colossal_full_harmony_Neuts@meta.data$covid_status <- COMET_10X_CLINICAL_SCORES$covid_status[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$Age <- COMET_10X_CLINICAL_SCORES$Age[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$Day_onset_category <- COMET_10X_CLINICAL_SCORES$Day_onset_category[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$ICU_vs_FLOOR <- COMET_10X_CLINICAL_SCORES$ICU_vs_FLOOR[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$Other_infection_type <- COMET_10X_CLINICAL_SCORES$Other_infection_type[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$Death <- COMET_10X_CLINICAL_SCORES$Death[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$NIH_score <- COMET_10X_CLINICAL_SCORES$NIH_score[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$Sampling_score <- COMET_10X_CLINICAL_SCORES$Sampling_score[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$Overall_score <- COMET_10X_CLINICAL_SCORES$Overall_score[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$Qualitative_score <- COMET_10X_CLINICAL_SCORES$Qualitative_score[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$Spike_score <- COMET_10X_CLINICAL_SCORES$Spike_score[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$Nucleocapsid_score <- COMET_10X_CLINICAL_SCORES$Nucleocapsid_score[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$ORF3a_score <- COMET_10X_CLINICAL_SCORES$ORF3a_score[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
merged_colossal_full_harmony_Neuts@meta.data$RBD_score <- COMET_10X_CLINICAL_SCORES$RBD_score[match(merged_colossal_full_harmony_Neuts@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]

# Create color scheme for COVID-19 status and disease severity
mycolorsseverity <- setNames(c("grey40", "orange", "orangered2"), c('CTRL', 'MILD', 'SEVERE'))
mycolorsstatus <- setNames(c("grey40", "dodgerblue3", "firebrick3"), c('CTRL', 'NEG', 'POS'))

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
################################################## FINAL NEUTS ########################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

#Load the data
load("~/Tristan/Labo MK/Manips/COVID/200428 patient_data_SO/merged/merged_colossal_paper_final/annotated_res_0.6/neuts/merged_colossal_paper_final_res0.6_neuts.RData")
merged_colossal_full_harmony_Neuts <- neuts
rm(neuts)

#Reset clustering and try resolution 1
merged_colossal_full_harmony_Neuts_1 <- FindClusters(merged_colossal_full_harmony_Neuts, verbose = TRUE, algorithm = 1, resolution = 1, method = "igraph", random.seed = 21212)  
DimPlot(merged_colossal_full_harmony_Neuts_1 , reduction='umap', label = T, label.size = 6, repel = T) + labs(color = "Resolution 1")

#Explore resolution 1 (calculate DEG between clusters, rank them, create a top10 list of them and plot it in a dotplot)
merged_colossal_full_harmony_Neuts_1_Markers <- FindAllMarkers(merged_colossal_full_harmony_Neuts_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson", latent.vars = "LIBRARY", assay = "RNA")
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers[which(merged_colossal_full_harmony_Neuts_1_Markers$p_val_adj<0.1),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$cluster,decreasing = FALSE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DotPlot(merged_colossal_full_harmony_Neuts_1, features = unique(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10$gene), cols = "RdBu", assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()

############################################################## 1st CLEANUP ###########################################################

# Annotate all cells as Neuts
merged_colossal_full_harmony_Neuts_1@meta.data$All_Neuts_annotations <- ifelse(merged_colossal_full_harmony_Neuts_1@meta.data$seurat_clusters %in% c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21'), 'NEUTS', '?')

# According to previous dotplot, annotate cells that we KEEP and those that we REMOVE (= RBC/T/B/Platelets contaminants)
merged_colossal_full_harmony_Neuts_1@meta.data$Subsetting_annotations <- ifelse(merged_colossal_full_harmony_Neuts_1@meta.data$seurat_clusters %in% c('0', '1', '2', '3', '4', '5', '6', '10', '11', '12', '13', '15', '16', '17', '18', '19', '21'), 'KEEP',
                                                                           ifelse(merged_colossal_full_harmony_Neuts_1@meta.data$seurat_clusters %in% c('7', '8', '9', '14', '20'), 'REMOVE', '?'))

# Observe annotation to confirm its accuracy
DimPlot(merged_colossal_full_harmony_Neuts_1, reduction='umap', label = F, group.by = "Subsetting_annotations")
ggplot(merged_colossal_full_harmony_Neuts_1@meta.data) + geom_bar(aes(x=All_Neuts_annotations, fill=Subsetting_annotations) , stat="count" , position="fill")

# Subset the cells that we KEEP
Idents(merged_colossal_full_harmony_Neuts_1) <- "Subsetting_annotations"
merged_colossal_full_harmony_Neuts <- subset(merged_colossal_full_harmony_Neuts_1, idents = c('KEEP'), invert = FALSE)
DimPlot(merged_colossal_full_harmony_Neuts , reduction='umap')

# Explore 20 first PCs to see if they are all relevant to use
png("NEUTS_DimHeatmap_20_PCs.png" , width = 14 , height = 10, units = "in", res = 200)
DimHeatmap(merged_colossal_full_harmony_Neuts, dims = 1:20, nfeatures = 5)
dev.off()

# Use the 20 first PCs to recalculate the UMAP space and the clusters neighboring
merged_colossal_full_harmony_Neuts <- RunUMAP(merged_colossal_full_harmony_Neuts, dims = 1:20, n.neighbors = 30, min.dist = 0.3, spread = 1, verbose = FALSE, seed.use = 21212, reduction = 'harmony')
merged_colossal_full_harmony_Neuts<- FindNeighbors(merged_colossal_full_harmony_Neuts, dims = 1:20, k.param = 20,  verbose = FALSE, reduction = 'harmony')

#Reset clustering, try and explore resolution 1
merged_colossal_full_harmony_Neuts_1 <- FindClusters(merged_colossal_full_harmony_Neuts, verbose = TRUE, algorithm = 1, resolution = 1, method = "igraph", random.seed = 21212)  
DimPlot(merged_colossal_full_harmony_Neuts_1 , reduction='umap', label = T, label.size = 6, repel = T) + labs(color = "Resolution 1")

merged_colossal_full_harmony_Neuts_1_Markers <- FindAllMarkers(merged_colossal_full_harmony_Neuts_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson", latent.vars = "LIBRARY", assay = "RNA")
beep(2)
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers[which(merged_colossal_full_harmony_Neuts_1_Markers$p_val_adj<0.1),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$cluster,decreasing = FALSE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DotPlot(merged_colossal_full_harmony_Neuts_1, features = unique(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10$gene), cols = "RdBu", assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()

############################################################## 2nd CLEANUP ###########################################################

merged_colossal_full_harmony_Neuts_1@meta.data$Subsetting_annotations <- ifelse(merged_colossal_full_harmony_Neuts_1@meta.data$seurat_clusters %in% c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '11', '13', '14', '15', '16', '17', '20', '21'), 'KEEP',
                                                                         ifelse(merged_colossal_full_harmony_Neuts_1@meta.data$seurat_clusters %in% c('10', '12', '18', '19'), 'REMOVE', '?'))

DimPlot(merged_colossal_full_harmony_Neuts_1, reduction='umap', label = F, group.by = "Subsetting_annotations")
ggplot(merged_colossal_full_harmony_Neuts_1@meta.data) + geom_bar(aes(x=All_Neuts_annotations, fill=Subsetting_annotations) , stat="count" , position="fill")

Idents(merged_colossal_full_harmony_Neuts_1) <- "Subsetting_annotations"
merged_colossal_full_harmony_Neuts <- subset(merged_colossal_full_harmony_Neuts_1, idents = c('KEEP'), invert = FALSE)
DimPlot(merged_colossal_full_harmony_Neuts , reduction='umap')

merged_colossal_full_harmony_Neuts <- RunUMAP(merged_colossal_full_harmony_Neuts, dims = 1:20, n.neighbors = 30, min.dist = 0.3, spread = 1, verbose = FALSE, seed.use = 21212, reduction = 'harmony')
merged_colossal_full_harmony_Neuts <- FindNeighbors(merged_colossal_full_harmony_Neuts, dims = 1:20, k.param = 20,  verbose = FALSE, reduction = 'harmony')
merged_colossal_full_harmony_Neuts_1 <- FindClusters(merged_colossal_full_harmony_Neuts, verbose = TRUE, algorithm = 1, resolution = 1, method = "igraph", random.seed = 21212)  
DimPlot(merged_colossal_full_harmony_Neuts_1 , reduction='umap', label = T, label.size = 6, repel = T) + labs(color = "Resolution 1")

merged_colossal_full_harmony_Neuts_1_Markers <- FindAllMarkers(merged_colossal_full_harmony_Neuts_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson", latent.vars = "LIBRARY", assay = "RNA")
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers[which(merged_colossal_full_harmony_Neuts_1_Markers$p_val_adj<0.1),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$cluster,decreasing = FALSE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DotPlot(merged_colossal_full_harmony_Neuts_1, features = unique(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10$gene), cols = "RdBu", assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()

############################################################## 3rd CLEANUP ###########################################################

merged_colossal_full_harmony_Neuts_1@meta.data$Subsetting_annotations <- ifelse(merged_colossal_full_harmony_Neuts_1@meta.data$seurat_clusters %in% c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '12', '14'), 'KEEP',
                                                                                ifelse(merged_colossal_full_harmony_Neuts_1@meta.data$seurat_clusters %in% c('11', '13', '15'), 'REMOVE', '?'))

DimPlot(merged_colossal_full_harmony_Neuts_1, reduction='umap', label = F, group.by = "Subsetting_annotations")
ggplot(merged_colossal_full_harmony_Neuts_1@meta.data) + geom_bar(aes(x=All_Neuts_annotations, fill=Subsetting_annotations) , stat="count" , position="fill")

Idents(merged_colossal_full_harmony_Neuts_1) <- "Subsetting_annotations"
merged_colossal_full_harmony_Neuts <- subset(merged_colossal_full_harmony_Neuts_1, idents = c('KEEP'), invert = FALSE)
DimPlot(merged_colossal_full_harmony_Neuts , reduction='umap')

merged_colossal_full_harmony_Neuts <- RunUMAP(merged_colossal_full_harmony_Neuts, dims = 1:20, n.neighbors = 30, min.dist = 0.3, spread = 1, verbose = FALSE, seed.use = 21212, reduction = 'harmony')
merged_colossal_full_harmony_Neuts <- FindNeighbors(merged_colossal_full_harmony_Neuts, dims = 1:20, k.param = 20,  verbose = FALSE, reduction = 'harmony')
merged_colossal_full_harmony_Neuts_1 <- FindClusters(merged_colossal_full_harmony_Neuts, verbose = TRUE, algorithm = 1, resolution = 1, method = "igraph", random.seed = 21212)  
DimPlot(merged_colossal_full_harmony_Neuts_1 , reduction='umap', label = T, label.size = 6, repel = T) + labs(color = "Resolution 1")

merged_colossal_full_harmony_Neuts_1_Markers <- FindAllMarkers(merged_colossal_full_harmony_Neuts_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson", latent.vars = "LIBRARY", assay = "RNA")
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers[which(merged_colossal_full_harmony_Neuts_1_Markers$p_val_adj<0.1),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$cluster,decreasing = FALSE),]
write.table(format(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1, digits=2), file="merged_colossal_full_harmony_Neuts_1_Markers_padj0.1.tsv", row.names=T, col.names=T, quote=F, sep="\t")
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DotPlot(merged_colossal_full_harmony_Neuts_1, features = unique(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10$gene), cols = "RdBu", assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()

# Write metadata file and save clean object, to be re-harmonized before final clustering and annotations
write.table(format(merged_colossal_full_harmony_Neuts_1@meta.data, digits=2), file="201001_merged_colossal_full_harmony_Neuts_1_metadata.tsv", row.names=T, col.names=T, quote=F, sep="\t")
save(merged_colossal_full_harmony_Neuts_1, file="201001_merged_colossal_full_harmony_Neuts_1_1.Robj")

############################################################## 4th CLEANUP ###########################################################

load("~/Tristan/Labo MK/Manips/COVID/200428 patient_data_SO/merged/merged_colossal_paper_final/annotated_res_0.6/neuts/CLEANED/merged_colossal_paper_final_res0.6_neuts.RData")
merged_colossal_full_harmony_Neuts <- neuts
rm(neuts)

merged_colossal_full_harmony_Neuts_1 <- FindClusters(merged_colossal_full_harmony_Neuts, verbose = TRUE, algorithm = 1, resolution = 1, method = "igraph", random.seed = 21212)  
DimPlot(merged_colossal_full_harmony_Neuts_1 , reduction='umap', label = T, label.size = 6, repel = T) + labs(color = "Resolution 1")

merged_colossal_full_harmony_Neuts_1_Markers <- FindAllMarkers(merged_colossal_full_harmony_Neuts_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson", latent.vars = "LIBRARY", assay = "RNA")
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers[which(merged_colossal_full_harmony_Neuts_1_Markers$p_val_adj<0.1),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$cluster,decreasing = FALSE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DotPlot(merged_colossal_full_harmony_Neuts_1, features = unique(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10$gene), cols = "RdBu", assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()

merged_colossal_full_harmony_Neuts_1@meta.data$All_Neuts_annotations <- ifelse(merged_colossal_full_harmony_Neuts_1@meta.data$seurat_clusters %in% c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19'), 'NEUTS', '?')

merged_colossal_full_harmony_Neuts_1@meta.data$Subsetting_annotations <- ifelse(merged_colossal_full_harmony_Neuts_1@meta.data$seurat_clusters %in% c('0', '1', '2', '3', '4', '5', '6', '10', '11', '12', '14', '15'), 'KEEP',
                                                                                ifelse(merged_colossal_full_harmony_Neuts_1@meta.data$seurat_clusters %in% c('7', '8', '9', '13', '16', '17', '18', '19'), 'REMOVE', '?'))

DimPlot(merged_colossal_full_harmony_Neuts_1, reduction='umap', label = F, group.by = "Subsetting_annotations")
ggplot(merged_colossal_full_harmony_Neuts_1@meta.data) + geom_bar(aes(x=All_Neuts_annotations, fill=Subsetting_annotations) , stat="count" , position="fill")

Idents(merged_colossal_full_harmony_Neuts_1) <- "Subsetting_annotations"
merged_colossal_full_harmony_Neuts <- subset(merged_colossal_full_harmony_Neuts_1, idents = c('KEEP'), invert = FALSE)
DimPlot(merged_colossal_full_harmony_Neuts , reduction='umap')

merged_colossal_full_harmony_Neuts <- RunUMAP(merged_colossal_full_harmony_Neuts, dims = 1:20, n.neighbors = 30, min.dist = 0.3, spread = 1, verbose = FALSE, seed.use = 21212, reduction = 'harmony')
merged_colossal_full_harmony_Neuts <- FindNeighbors(merged_colossal_full_harmony_Neuts, dims = 1:20, k.param = 20,  verbose = FALSE, reduction = 'harmony')
merged_colossal_full_harmony_Neuts_1 <- FindClusters(merged_colossal_full_harmony_Neuts, verbose = TRUE, algorithm = 1, resolution = 1, method = "igraph", random.seed = 21212)  
DimPlot(merged_colossal_full_harmony_Neuts_1 , reduction='umap', label = T, label.size = 6, repel = T) + labs(color = "Resolution 1")
merged_colossal_full_harmony_Neuts_2 <- FindClusters(merged_colossal_full_harmony_Neuts, verbose = TRUE, algorithm = 1, resolution = 2, method = "igraph", random.seed = 21212)  
DimPlot(merged_colossal_full_harmony_Neuts_2 , reduction='umap', label = T, label.size = 6, repel = T) + labs(color = "Resolution 2")

merged_colossal_full_harmony_Neuts_1_Markers <- FindAllMarkers(merged_colossal_full_harmony_Neuts_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson", latent.vars = "LIBRARY", assay = "RNA")
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers[which(merged_colossal_full_harmony_Neuts_1_Markers$p_val_adj<0.1),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1$cluster,decreasing = FALSE),]
merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10 <- merged_colossal_full_harmony_Neuts_1_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DotPlot(merged_colossal_full_harmony_Neuts_1, features = unique(merged_colossal_full_harmony_Neuts_1_Markers_padj0.1_Top10$gene), cols = "RdBu", assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()

merged_colossal_full_harmony_Neuts_2_Markers <- FindAllMarkers(merged_colossal_full_harmony_Neuts_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson", latent.vars = "LIBRARY", assay = "RNA")
merged_colossal_full_harmony_Neuts_2_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_2_Markers[which(merged_colossal_full_harmony_Neuts_2_Markers$p_val_adj<0.1),]
merged_colossal_full_harmony_Neuts_2_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_2_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_2_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
merged_colossal_full_harmony_Neuts_2_Markers_padj0.1 <- merged_colossal_full_harmony_Neuts_2_Markers_padj0.1[order(merged_colossal_full_harmony_Neuts_2_Markers_padj0.1$cluster,decreasing = FALSE),]
merged_colossal_full_harmony_Neuts_2_Markers_padj0.1_Top10 <- merged_colossal_full_harmony_Neuts_2_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DotPlot(merged_colossal_full_harmony_Neuts_2, features = unique(merged_colossal_full_harmony_Neuts_2_Markers_padj0.1_Top10$gene), cols = "RdBu", assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()

######################################################### FINAL ANNOTATIONS RES 2 ###########################################################

#Annotate the objects with different levels of resolutions, going from ALL cells to Coarse subtypes to Fine subtypes
merged_colossal_full_harmony_Neuts_2@meta.data$ALL_Neuts <- ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'), 'NEUTS', '?')
DimPlot(merged_colossal_full_harmony_Neuts_2, reduction='umap', label = F, group.by = "ALL_Neuts")

merged_colossal_full_harmony_Neuts_2@meta.data$Fine_Subs_annotations <- ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('0'), 'S100A12_1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('1'), 'RIBO_1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('2'), 'G0S2_1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('3'), 'ISG_1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('4'), 'ISG_2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('5'), 'RIBO_2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('6'), 'NEAT1_1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('7'), 'SLPI_1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('8'), 'S100A12_2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('9'), 'NEAT1_2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('10'), 'S100A12_3',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('11'), 'S100A12_4',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('12'), 'S100A12_5',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('13'), 'NEAT1_3',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('14'), 'G0S2_2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('15'), 'S100A12_6',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('16'), 'SLPI_2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('17'), 'ISG_3',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('18'), 'NEAT1_4',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('19'), 'LCN2_1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('20'), 'ISG_4',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('21'), 'S100A12_7',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('22'), 'S100A12_8',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('23'), 'ISG_5',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('24'), 'LCN2_2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('25'), 'S100A12_10',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('26'), 'S100A12_11', '?')))))))))))))))))))))))))))

DimPlot(merged_colossal_full_harmony_Neuts_2 , reduction='umap', label = T, repel = T, group.by = "Fine_Subs_annotations")

merged_colossal_full_harmony_Neuts_2@meta.data$Coarse_Subs_annotations <- ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('0'), 'S100A12',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('1'), 'RIBO.',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('2'), 'G0S2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('3'), 'ISG',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('4'), 'ISG',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('5'), 'RIBO.',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('6'), 'NEAT1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('7'), 'SLPI',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('8'), 'S100A12',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('9'), 'NEAT1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('10'), 'S100A12',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('11'), 'S100A12',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('12'), 'S100A12',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('13'), 'NEAT1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('14'), 'G0S2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('15'), 'S100A12',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('16'), 'SLPI',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('17'), 'ISG',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('18'), 'NEAT1',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('19'), 'LCN2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('20'), 'ISG',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('21'), 'S100A12',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('22'), 'S100A12',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('23'), 'ISG',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('24'), 'LCN2',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('25'), 'S100A12',
                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('26'), 'S100A12', '?')))))))))))))))))))))))))))

DimPlot(merged_colossal_full_harmony_Neuts_2 , reduction='umap', label = T, repel = T, group.by = "Coarse_Subs_annotations")

#Organize the annotations for plotting purposes
Idents(merged_colossal_full_harmony_Neuts_2) <- "seurat_clusters"
My_level_seurat<- c('0', '8', '10', '11', '12', '15', '21', '22', '25', '26', '1', '5', '2', '14', '3', '4', '17', '20', '23', '6', '9', '13', '18', '7', '16', '19', '24')
levels(merged_colossal_full_harmony_Neuts_2) <- My_level_seurat
DotPlot(merged_colossal_full_harmony_Neuts_2, features = unique(merged_colossal_full_harmony_Neuts_2_Markers_padj0.1_Top10$gene), cols = "RdBu", assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
DotPlot(merged_colossal_full_harmony_Neuts_2, features = c('RPS9', 'PI3', 'SLPI', 'NFKBIA', 'CXCL8', 'G0S2', 'FTH1', 'MALAT1', 'NEAT1', 'FCGR3B', 'IFITM3', 'IFITM1', 'IFIT3', 'IFIT2', 'IFIT1', 'LY6E', 'ISG15', 'MT2A', 'S100A11', 'GCA', 'CST7', 'ACTB', 'S100A4', 'S100A6', 'S100A9', 'MYL6', 'TSPO', 'S100A8', 'S100A12', 'PGLYRP1', 'MMP9', 'CAMP', 'RETN', 'LTF', 'LCN2'), cols = "RdBu") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
Idents(merged_colossal_full_harmony_Neuts_2) <- "Coarse_Subs_annotations"
My_level_coarse<- c('LCN2', 'S100A12', 'ISG', 'NEAT1', 'G0S2', 'SLPI', 'RIBO.')
levels(merged_colossal_full_harmony_Neuts_2) <- My_level_coarse
pdf('1C Neuts DotPlot.pdf', width=5, height=6, useDingbats=F)
DotPlot(merged_colossal_full_harmony_Neuts_2, features = c('RPS9', 'PI3', 'SLPI', 'NFKBIA', 'CXCL8', 'G0S2', 'FTH1', 'MALAT1', 'NEAT1', 'FCGR3B', 'IFITM3', 'IFITM1', 'IFIT3', 'IFIT2', 'IFIT1', 'LY6E', 'ISG15', 'MT2A', 'S100A11', 'GCA', 'CST7', 'ACTB', 'S100A4', 'S100A6', 'S100A9', 'MYL6', 'TSPO', 'S100A8', 'S100A12', 'PGLYRP1', 'MMP9', 'CAMP', 'RETN', 'LTF', 'LCN2'), cols = "RdBu") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
dev.off()

#Export counts used to create prism graphs
merged_colossal_full_harmony_Neuts_2_Counts <- table(merged_colossal_full_harmony_Neuts_2@meta.data[,c('SAMPLE.by.SNPs', 'Coarse_Subs_annotations')])
write.table(merged_colossal_full_harmony_Neuts_2_Counts, file="merged_colossal_full_harmony_Neuts_2_Counts.tsv", row.names=T, col.names=T, sep="\t")

# Write metadata file and save clean object
write.table(format(merged_colossal_full_harmony_Neuts_2@meta.data, digits=2), file="201001_merged_colossal_full_harmony_Neuts_2_metadata_FINAL.tsv", row.names=T, col.names=T, quote=F, sep="\t")
save(merged_colossal_full_harmony_Neuts_2, file="201001_merged_colossal_full_harmony_Neuts_2_FINAL.Robj")

# Produce UMAPs of Figure 1
pdf('1D Neuts UMAP.pdf', width=6, height=5, useDingbats=F)
DimPlot(merged_colossal_full_harmony_Neuts_2 , reduction='umap', label = F, group.by = "Coarse_Subs_annotations", cols = "Paired", order = c('RIBO.', 'SLPI', 'G0S2', 'NEAT1', 'ISG', 'S100A12', 'LCN2'))
dev.off()

pdf('1E Neuts UMAP status.pdf', width=6, height=5, useDingbats=F)
DimPlot(merged_colossal_full_harmony_Neuts_2 , reduction='umap', label = F, group.by = "covid_status", cols = mycolorsstatus)
dev.off()

pdf('1F Neuts UMAP severity.pdf', width=6, height=5, useDingbats=F)
DimPlot(merged_colossal_full_harmony_Neuts_2 , reduction='umap', label = F, group.by = "Qualitative_score", cols = mycolorsseverity)
dev.off()

############################################## ANNOTATIONS FOR PHEMD GRAPHS ########################################################

merged_colossal_full_harmony_Neuts_2@meta.data$harmony_cluster_final <- ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('0'), 'NEUTS_S100A12',
                                                                                 ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('1'), 'NEUTS_RIBO',
                                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('2'), 'NEUTS_G0S2',
                                                                                               ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('3'), 'NEUTS_ISG',
                                                                                                      ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('4'), 'NEUTS_ISG',
                                                                                                             ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('5'), 'NEUTS_RIBO',
                                                                                                                    ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('6'), 'NEUTS_NEAT1',
                                                                                                                           ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('7'), 'NEUTS_SLPI',
                                                                                                                                  ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('8'), 'NEUTS_S100A12',
                                                                                                                                         ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('9'), 'NEUTS_NEAT1',
                                                                                                                                                ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('10'), 'NEUTS_S100A12',
                                                                                                                                                       ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('11'), 'NEUTS_S100A12',
                                                                                                                                                              ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('12'), 'NEUTS_S100A12',
                                                                                                                                                                     ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('13'), 'NEUTS_NEAT1',
                                                                                                                                                                            ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('14'), 'NEUTS_G0S2',
                                                                                                                                                                                   ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('15'), 'NEUTS_S100A12',
                                                                                                                                                                                          ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('16'), 'NEUTS_SLPI',
                                                                                                                                                                                                 ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('17'), 'NEUTS_ISG',
                                                                                                                                                                                                        ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('18'), 'NEUTS_NEAT1',
                                                                                                                                                                                                               ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('19'), 'NEUTS_LCN2',
                                                                                                                                                                                                                      ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('20'), 'NEUTS_ISG',
                                                                                                                                                                                                                             ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('21'), 'NEUTS_S100A12',
                                                                                                                                                                                                                                    ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('22'), 'NEUTS_S100A12',
                                                                                                                                                                                                                                           ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('23'), 'NEUTS_ISG',
                                                                                                                                                                                                                                                  ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('24'), 'NEUTS_LCN2',
                                                                                                                                                                                                                                                         ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('25'), 'NEUTS_S100A12',
                                                                                                                                                                                                                                                                ifelse(merged_colossal_full_harmony_Neuts_2@meta.data$seurat_clusters %in% c('26'), 'NEUTS_S100A12', '?')))))))))))))))))))))))))))
DimPlot(merged_colossal_full_harmony_Neuts_2, group.by = "harmony_cluster_final")

write.table(merged_colossal_full_harmony_Neuts_2@meta.data, file="201001_NEUTS_METADATA_FOR_WILL.tsv", row.names=T, col.names=T, sep="\t")

##################################################################################################
##################################### NEUTS SIGNATURES ###########################################
##################################################################################################

# Call neutrophil degranulation genes and use them to produce a degranulation score
Neuts_Degranulation_genes <- scan("~/Tristan/Labo MK/Manips/COVID/NEUTROPHIL_DEGRANULATION.csv", what = "character", sep = NULL)
gene_df <- get_genes_s3(Neuts_Degranulation_genes, merged_colossal_full_harmony_Neuts_2, drop = T)
Score_Degranulation <- colMeans(gene_df)
merged_colossal_full_harmony_Neuts_2@meta.data$DegranuScore_0 <- saturate(vec=Score_Degranulation, sat=0, binary=FALSE)

# Call ISG genes and use them to produce an ISG score
Sig_IFNNeuts = c('MT2A', 'ISG15', 'LY6E', 'IFIT1', 'IFIT2', 'IFIT3', 'IFITM1', 'IFITM3', 'IFI44L', 'IFI6', 'MX1', 'IFI27')
gene_df <- get_genes_s3(Sig_IFNNeuts, merged_colossal_full_harmony_Neuts_2, drop = T)
Score_ISG <- colMeans(gene_df)
merged_colossal_full_harmony_Neuts_2@meta.data$ISGScore_0 <- saturate(vec=Score_ISG, sat=0, binary=FALSE)

# Subset neutrophils according to status and/or severity to produce graphs and stats
Idents(merged_colossal_full_harmony_Neuts_2) <- merged_colossal_full_harmony_Neuts_2@meta.data$covid_status
merged_colossal_full_harmony_Neuts_2_CTRL <- subset(merged_colossal_full_harmony_Neuts_2, idents = c('CTRL'), invert = FALSE)
merged_colossal_full_harmony_Neuts_2_NEG <- subset(merged_colossal_full_harmony_Neuts_2, idents = c('NEG'), invert = FALSE)
merged_colossal_full_harmony_Neuts_2_POS <- subset(merged_colossal_full_harmony_Neuts_2, idents = c('POS'), invert = FALSE)
Idents(merged_colossal_full_harmony_Neuts_2) <- merged_colossal_full_harmony_Neuts_2@meta.data$Qualitative_score
merged_colossal_full_harmony_Neuts_2_MILD <- subset(merged_colossal_full_harmony_Neuts_2, idents = c('MILD'), invert = FALSE)
merged_colossal_full_harmony_Neuts_2_SEVERE <- subset(merged_colossal_full_harmony_Neuts_2, idents = c('SEVERE'), invert = FALSE)
Idents(merged_colossal_full_harmony_Neuts_2_NEG) <- merged_colossal_full_harmony_Neuts_2_NEG@meta.data$Qualitative_score
Idents(merged_colossal_full_harmony_Neuts_2_POS) <- merged_colossal_full_harmony_Neuts_2_POS@meta.data$Qualitative_score
merged_colossal_full_harmony_Neuts_2_NEG_MILD <- subset(merged_colossal_full_harmony_Neuts_2_NEG, idents = c('MILD'), invert = FALSE)
merged_colossal_full_harmony_Neuts_2_POS_MILD <- subset(merged_colossal_full_harmony_Neuts_2_POS, idents = c('MILD'), invert = FALSE)
merged_colossal_full_harmony_Neuts_2_NEG_SEVERE <- subset(merged_colossal_full_harmony_Neuts_2_NEG, idents = c('SEVERE'), invert = FALSE)
merged_colossal_full_harmony_Neuts_2_POS_SEVERE <- subset(merged_colossal_full_harmony_Neuts_2_POS, idents = c('SEVERE'), invert = FALSE)

# Produce violin plots of Figure 1 and S1
pdf('S1M ISG_VIOLIN.pdf', width=4.5, height=5, useDingbats=F)
VlnPlot(merged_colossal_full_harmony_Neuts_2, features = "ISGScore_0", group.by = "covid_status", pt.size = 0.0, split.by = "Qualitative_score", cols = mycolorsseverity)
dev.off()

pdf('S1N ISG VIOLIN ALL SUBS ALL.pdf', width=6.5, height=5, useDingbats=F)
VlnPlot(merged_colossal_full_harmony_Neuts_2, features = "ISGScore_0", group.by = "Coarse_Subs_annotations", pt.size = 0.0, split.by = "Qualitative_score", cols = mycolorsseverity)
dev.off()

pdf('S1O ISG VIOLIN ALL SUBS NEG.pdf', width=6.5, height=5, useDingbats=F)
VlnPlot(merged_colossal_full_harmony_Neuts_2_NEG, features = "ISGScore_0", group.by = "Coarse_Subs_annotations", pt.size = 0.0, split.by = "Qualitative_score", cols = mycolorsseverity)
dev.off()

pdf('1K ISG VIOLIN ALL SUBS POS.pdf', width=6.5, height=5, useDingbats=F)
VlnPlot(merged_colossal_full_harmony_Neuts_2_POS, features = "ISGScore_0", group.by = "Coarse_Subs_annotations", pt.size = 0.0, split.by = "Qualitative_score", cols = mycolorsseverity)
dev.off()

pdf('S1P DEGRANU_VIOLIN_ALL.pdf', width=3, height=5, useDingbats=F)
VlnPlot(merged_colossal_full_harmony_Neuts_2, features = "DegranuScore_0", group.by = "Qualitative_score", pt.size = 0.0, cols = mycolorsseverity)
dev.off()

pdf('S1Q DEGRANU_VIOLIN_SPLIT.pdf', width=4.5, height=5, useDingbats=F)
VlnPlot(merged_colossal_full_harmony_Neuts_2, features = "DegranuScore_0", group.by = "covid_status", pt.size = 0.0, split.by = "Qualitative_score", cols = mycolorsseverity)
dev.off()

# Wilcoxon tests
CTRL_ISG <- merged_colossal_full_harmony_Neuts_2_CTRL@meta.data$ISGScore_0
NEG_ISG <- merged_colossal_full_harmony_Neuts_2_NEG@meta.data$ISGScore_0
POS_ISG <- merged_colossal_full_harmony_Neuts_2_POS@meta.data$ISGScore_0
MILD_ISG <- merged_colossal_full_harmony_Neuts_2_MILD@meta.data$ISGScore_0
SEV_ISG <- merged_colossal_full_harmony_Neuts_2_SEVERE@meta.data$ISGScore_0
PM_ISG <- merged_colossal_full_harmony_Neuts_2_POS_MILD@meta.data$ISGScore_0
PS_ISG <- merged_colossal_full_harmony_Neuts_2_POS_SEVERE@meta.data$ISGScore_0
NM_ISG <- merged_colossal_full_harmony_Neuts_2_NEG_MILD@meta.data$ISGScore_0
NS_ISG <- merged_colossal_full_harmony_Neuts_2_NEG_SEVERE@meta.data$ISGScore_0

my_data <- data.frame(group = c(rep("NEG_MILD", length(NM_ISG)), rep("NEG_SEVERE", length(NS_ISG))),
                      ISG_score = c(NM_ISG, NS_ISG))

ggboxplot(my_data, x = "group", y = "ISG_score", 
          color = "group", palette = c("#00AFBB", "#E7B800"),
          ylab = "ISG_score", xlab = "Groups")

wilcox.test(ISG_score ~ group, data = my_data, exact = FALSE)

CTRL_DEGR <- merged_colossal_full_harmony_Neuts_2_CTRL@meta.data$DegranuScore_0
NEG_DEGR <- merged_colossal_full_harmony_Neuts_2_NEG@meta.data$DegranuScore_0
POS_DEGR <- merged_colossal_full_harmony_Neuts_2_POS@meta.data$DegranuScore_0
MILD_DEGR <- merged_colossal_full_harmony_Neuts_2_MILD@meta.data$DegranuScore_0
SEV_DEGR <- merged_colossal_full_harmony_Neuts_2_SEVERE@meta.data$DegranuScore_0
PM_DEGR <- merged_colossal_full_harmony_Neuts_2_POS_MILD@meta.data$DegranuScore_0
PS_DEGR <- merged_colossal_full_harmony_Neuts_2_POS_SEVERE@meta.data$DegranuScore_0
NM_DEGR <- merged_colossal_full_harmony_Neuts_2_NEG_MILD@meta.data$DegranuScore_0
NS_DEGR <- merged_colossal_full_harmony_Neuts_2_NEG_SEVERE@meta.data$DegranuScore_0

my_data <- data.frame(group = c(rep("NEG_MILD", length(NM_DEGR)), rep("POS_MILD", length(PM_DEGR))),
                      DEGR_score = c(NM_DEGR, PM_DEGR))

ggboxplot(my_data, x = "group", y = "DEGR_score", 
          color = "group", palette = c("#00AFBB", "#E7B800"),
          ylab = "DEGR_score", xlab = "Groups")

wilcox.test(DEGR_score ~ group, data = my_data, exact = FALSE)

Idents(merged_colossal_full_harmony_Neuts_2) <- "Coarse_Subs_annotations"
LCN2_Neuts <- subset(merged_colossal_full_harmony_Neuts_2, idents = c('LCN2'), invert = FALSE)
DimPlot(LCN2_Neuts , reduction='umap')

Idents(LCN2_Neuts) <- LCN2_Neuts@meta.data$covid_status
LCN2_Neuts_POS <- subset(LCN2_Neuts, idents = c('POS'), invert = FALSE)

Idents(LCN2_Neuts_POS) <- "Qualitative_score"
LCN2_Neuts_POS_MILD <- subset(LCN2_Neuts_POS, idents = c('MILD'), invert = FALSE)
LCN2_Neuts_POS_SEVERE <- subset(LCN2_Neuts_POS, idents = c('SEVERE'), invert = FALSE)

PM_ISG <- LCN2_Neuts_POS_MILD@meta.data$ISGScore_0
PS_ISG <- LCN2_Neuts_POS_SEVERE@meta.data$ISGScore_0

my_data <- data.frame(group = c(rep("POS_MILD", length(PM_ISG)), rep("POS_SEVERE", length(PS_ISG))),
                      ISG_score = c(PM_ISG, PS_ISG))

ggboxplot(my_data, x = "group", y = "ISG_score", 
          color = "group", palette = c("#00AFBB", "#E7B800"),
          ylab = "ISG_score", xlab = "Groups")

wilcox.test(ISG_score ~ group, data = my_data, exact = FALSE)

##################################################################################################
####################################### DEG + VOLCANO ############################################
##################################################################################################

# VOLCANO PLOT ALL NEUTS POS vs NEG
Idents(merged_colossal_full_harmony_Neuts_2) <- merged_colossal_full_harmony_Neuts_2@meta.data$covid_status
merged_colossal_full_harmony_Neuts_2_DEG_STATUS_Markers <- FindMarkers(merged_colossal_full_harmony_Neuts_2, ident.1="POS", ident.2="NEG", min.pct = 0.25, test.use="poisson", latent.vars = "LIBRARY", assay = "RNA")
write.table(format(merged_colossal_full_harmony_Neuts_2_DEG_STATUS_Markers, digits=2), file="merged_colossal_full_harmony_Neuts_2_DEG_STATUS_Markers_prevolcano.tsv", row.names=T, col.names=T, quote=F, sep="\t")
#merged_colossal_full_harmony_Neuts_2_DEG_STATUS_Markers = read.table("~/Tristan/Labo MK/Manips/COVID/200428 patient_data_SO/200822 SEURAT WD/merged_colossal_full_harmony_Neuts_2_DEG_STATUS_Markers_prevolcano.tsv", sep='\t', header=T)

# VOLCANO PLOT ALL NEUTS MILD vs SEVERE
Idents(merged_colossal_full_harmony_Neuts_2) <- merged_colossal_full_harmony_Neuts_2@meta.data$Qualitative_score
merged_colossal_full_harmony_Neuts_2_DEG_SEVERITY_Markers <- FindMarkers(merged_colossal_full_harmony_Neuts_2, ident.1="MILD", ident.2="SEVERE", min.pct = 0.25, test.use="poisson", latent.vars = "LIBRARY", assay = "RNA")
write.table(format(merged_colossal_full_harmony_Neuts_2_DEG_SEVERITY_Markers, digits=2), file="merged_colossal_full_harmony_Neuts_2_DEG_SEVERITY_Markers_prevolcano.tsv", row.names=T, col.names=T, quote=F, sep="\t")
#merged_colossal_full_harmony_Neuts_2_DEG_SEVERITY_Markers = read.table("~/Tristan/Labo MK/Manips/COVID/200428 patient_data_SO/200822 SEURAT WD/merged_colossal_full_harmony_Neuts_2_DEG_SEVERITY_Markers_prevolcano.tsv", sep='\t', header=T)

pdf('S1I Volcano STATUS.pdf', width=6, height=6, useDingbats=F)
EnhancedVolcano(merged_colossal_full_harmony_Neuts_2_DEG_STATUS_Markers,
                lab = rownames(merged_colossal_full_harmony_Neuts_2_DEG_STATUS_Markers),
                x = 'avg_logFC',
                y = 'p_val',
                xlim = c(-1,0.75),
                title = 'ALL NEG VS POS',
                pCutoff = 10e-20,
                FCcutoff = 0.2,
                pointSize = 2.0,
                labSize = 2,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1,
                drawConnectors = T,
                widthConnectors = 0.25)
dev.off()

pdf('S1K Volcano SEVERITY.pdf', width=6, height=6, useDingbats=F)
EnhancedVolcano(merged_colossal_full_harmony_Neuts_2_DEG_SEVERITY_Markers,
                lab = rownames(merged_colossal_full_harmony_Neuts_2_DEG_SEVERITY_Markers),
                x = 'avg_logFC',
                y = 'p_val',
                xlim = c(-1.5,1),
                title = 'ALL MILD VS SEVERE',
                pCutoff = 10e-20,
                FCcutoff = 0.2,
                pointSize = 2.0,
                labSize = 2,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1,
                drawConnectors = T,
                widthConnectors = 0.25)
dev.off()

# Enjoy
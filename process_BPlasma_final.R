## Author: Kenneth Hu

#assigning severity and covid status to metadata
BPlasma@meta.data$severity = character(nrow(BPlasma@meta.data))
blah<-COMET_INTEGRATED_FINAL_CLEAN$Qual_score[match(BPlasma@meta.data$SAMPLE.by.SNPs,COMET_INTEGRATED_FINAL_CLEAN$Sample)]
BPlasma@meta.data$severity[blah==0]<-"Healthy"
BPlasma@meta.data$severity[blah==1]<-"Mild"
BPlasma@meta.data$severity[blah==2]<-"Severe"
BPlasma@meta.data$covid_status = character(nrow(BPlasma@meta.data))
BPlasma$covid_status<-COMET_INTEGRATED_FINAL_CLEAN$covid_status[match(BPlasma@meta.data$SAMPLE.by.SNPs,COMET_INTEGRATED_FINAL_CLEAN$Sample)]
#start with b_plasma_harmony, object subsetted from the whole object following harmony
BPlasma <- FindNeighbors(b_plasma_harmony, dims = 1:30, k.param = 20,  verbose = FALSE) #Arjun: "k depends but sqrt (total cells) should be good to start"
BPlasma <- FindClusters(BPlasma, verbose = TRUE, algorithm = 1, resolution = 0.2, random.seed = 21212)
BPlasma <- RunUMAP(BPlasma, dims = 1:15, n.neighbors = 30, min.dist = 0.1,  verbose = FALSE,assay = "RNA")
#visualize
DimPlot(BPlasma,label=TRUE,reduction = 'umap')
#find markers for clusters
markers0.1<- FindAllMarkers(BPlasma,
                          test.use='poisson',
                          only.pos=TRUE,
                          min.pct=0.25,
                          logfc.threshold=0.4,
                          assay='RNA',
                          latent.vars = 'LIBRARY'
)
b_plasma_harmony_0.2_Markers_padj0.1 <- markers0.1[which(markers0.1$p_val_adj<0.1),]
b_plasma_harmony_0.2_Markers_padj0.1 <- b_plasma_harmony_0.2_Markers_padj0.1[order(b_plasma_harmony_0.2_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
b_plasma_harmony_0.2_Markers_padj0.1 <- b_plasma_harmony_0.2_Markers_padj0.1[order(b_plasma_harmony_0.2_Markers_padj0.1$cluster,decreasing = FALSE),]
b_plasma_harmony_0.2_Markers_padj0.1_Top10 <- b_plasma_harmony_0.2_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
b_plasma_harmony_0.2_Markers_padj0.1_Top5 <- b_plasma_harmony_0.2_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#Dotplot to visualize marker genes
DotPlot(BPlasma,
        features = unique(b_plasma_harmony_0.2_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
#cleaning time, remove clusters 5/6/10 as contaminating clusters based on RBC/platelet/neutrophil markers
BPlasma_clean = SubsetData(BPlasma,ident.remove = c(5,6,10))
#redo cluster+umap
BPlasma_clean <- FindNeighbors(BPlasma_clean, dims = 1:15, k.param = 20,  verbose = FALSE,reduction = "harmony") #Arjun: "k depends but sqrt (total cells) should be good to start"
BPlasma_clean <- FindClusters(BPlasma_clean, verbose = TRUE, algorithm = 1, resolution = 0.6, random.seed = 21212,force.recalc=TRUE)
BPlasma_clean <- RunUMAP(BPlasma_clean, dims = 1:20, n.neighbors = 30, min.dist = 0.05,  verbose = FALSE,assay = "RNA",spread = 0.5)
DimPlot(BPlasma_clean,label=TRUE,reduction = 'umap')
#find markers here
markers0.4<- FindAllMarkers(BPlasma_clean,
                            test.use='poisson',
                            only.pos=TRUE,
                            min.pct=0.25,
                            logfc.threshold=0.4,
                            assay='RNA',
                            latent.vars = 'LIBRARY'
)

b_plasma_harmony_0.2_Markers_padj0.1 <- markers0.4[which(markers0.4$p_val_adj<0.1),]
b_plasma_harmony_0.2_Markers_padj0.1 <- b_plasma_harmony_0.2_Markers_padj0.1[order(b_plasma_harmony_0.2_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
b_plasma_harmony_0.2_Markers_padj0.1 <- b_plasma_harmony_0.2_Markers_padj0.1[order(b_plasma_harmony_0.2_Markers_padj0.1$cluster,decreasing = FALSE),]
b_plasma_harmony_0.2_Markers_padj0.1_Top10 <- b_plasma_harmony_0.2_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
b_plasma_harmony_0.2_Markers_padj0.1_Top5 <- b_plasma_harmony_0.2_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DotPlot(BPlasma_clean,
        features = unique(b_plasma_harmony_0.2_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
#cleaning again, removing clusters 8 and 9 based on RBC and platelet contaminating markers
BPlasma_clean2 = SubsetData(BPlasma_clean,ident.remove = c(8,9))
#Redo clustering and umap
BPlasma_clean2 <- FindNeighbors(BPlasma_clean2, dims = 1:15, k.param = 30,  verbose = FALSE,reduction = "harmony") #Arjun: "k depends but sqrt (total cells) should be good to start"
BPlasma_clean2 <- FindClusters(BPlasma_clean2, verbose = TRUE, algorithm = 1, resolution = 0.8, random.seed = 21212)
BPlasma_clean2 <- RunUMAP(BPlasma_clean2, dims = 1:20, n.neighbors = 30, min.dist = 0.05,  verbose = FALSE,assay = "RNA")
DimPlot(BPlasma_clean2,label=TRUE)
markers0.8<- FindAllMarkers(BPlasma_clean2,
                            test.use='poisson',
                            only.pos=TRUE,
                            min.pct=0.25,
                            logfc.threshold=0.4,
                            assay='RNA',
                            latent.vars = 'LIBRARY'
)
b_plasma_harmony_0.2_Markers_padj0.1 <- markers0.8[which(markers0.8$p_val_adj<0.1),]
b_plasma_harmony_0.2_Markers_padj0.1 <- b_plasma_harmony_0.2_Markers_padj0.1[order(b_plasma_harmony_0.2_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
b_plasma_harmony_0.2_Markers_padj0.1 <- b_plasma_harmony_0.2_Markers_padj0.1[order(b_plasma_harmony_0.2_Markers_padj0.1$cluster,decreasing = FALSE),]
b_plasma_harmony_0.2_Markers_padj0.1_Top10 <- b_plasma_harmony_0.2_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
b_plasma_harmony_0.2_Markers_padj0.1_Top5 <- b_plasma_harmony_0.2_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DotPlot(BPlasma_clean2,
        features = unique(b_plasma_harmony_0.2_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
# at this point, the subsetted object was re-harmonized
# Reclustered and UMAP
BPlasma <- FindNeighbors(BPlasma, dims = 1:20, k.param = 30,  verbose = FALSE) #Arjun: "k depends but sqrt (total cells) should be good to start"
BPlasma <- FindClusters(BPlasma, verbose = TRUE, algorithm = 1, resolution = 0.65, random.seed = 21212)
BPlasma <- RunUMAP(BPlasma, dims = 1:20, n.neighbors = 30, min.dist = 0.1,  verbose = FALSE,assay = "RNA")
DimPlot(BPlasma,label = TRUE)
markers0.8<- FindAllMarkers(BPlasma,
                            test.use='poisson',
                            only.pos=TRUE,
                            min.pct=0.2,
                            logfc.threshold=0.25,
                            assay='RNA',
                            latent.vars = 'LIBRARY'
)
markers <- markers0.8[which(markers0.8$p_val_adj<0.1),]
markers <- markers[order(markers$avg_logFC,decreasing = TRUE),]
markers <- markers[order(markers$cluster,decreasing = FALSE),]
markersTop10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
markersTop5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DotPlot(BPlasma,
        features = unique(c(markersTop5$gene,"CD74")),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
##remove cluster 8 which shows signs of cytotoxic cd8 t cell contamination
BPlasma <- SubsetData(BPlasma,ident.remove = c(8))
BPlasma <- FindNeighbors(BPlasma, dims = 1:15, k.param = 40,  verbose = FALSE) #Arjun: "k depends but sqrt (total cells) should be good to start"
BPlasma <- FindClusters(BPlasma, verbose = TRUE, algorithm = 1, resolution = 0.63, random.seed = 21212)
BPlasma <- RunUMAP(BPlasma, dims = 1:20, n.neighbors = 30, min.dist = 0.2,  verbose = FALSE,assay = "RNA")
DimPlot(BPlasma,label=TRUE)
markers0.8<- FindAllMarkers(BPlasma,
                            test.use='poisson',
                            only.pos=TRUE,
                            min.pct=0.2,
                            logfc.threshold=0.2,
                            assay='RNA',
                            latent.vars = 'LIBRARY'
)
markers <- markers0.8[which(markers0.8$p_val_adj<0.1),]
markers <- markers[order(markers$avg_logFC,decreasing = TRUE),]
markers <- markers[order(markers$cluster,decreasing = FALSE),]
markersTop10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
markersTop5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
###Final Plotting####
pdf("~/Krummel lab/COVID/redo/b_plasma_harmony_0.63_DotPlot.pdf", width = 7, height = 8,useDingbats = FALSE)
DotPlot(BPlasma,
        features = unique(c(markersTop5$gene)),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
dev.off()
#annotation of the clusters
BPlasma$subtype_final_annotation <- BPlasma$seurat_clusters
annotations <- c(
  "B_naive_1",
  "B_mature",
  "B_plasma_1",
  "B_naive_2",
  "B_memory",
  "B_ISG",
  "B_plasma_2" 
  )
names(annotations) <- levels(BPlasma$seurat_clusters)
BPlasma <- RenameIdents(BPlasma, annotations)
BPlasma@meta.data$seurat_cluster_final <- Idents(BPlasma)
#get counts of the clusters
table(BPlasma@meta.data$seurat_cluster_final)
write.table(b_plasma_harmony_0.4@meta.data, file="b_plasma_harmony_0.4_metadata.tsv", row.names=T, col.names=T, quote=F, sep="\t")
Idents(BPlasma) <- 'seurat_cluster_final'
reorder_clusters <- c('B_naive_1','B_naive_2','B_ISG','B_memory','B_mature','B_plasma_1','B_plasma_2')
levels(BPlasma) <- reorder_clusters
BPlasma@meta.data$harmony_cluster_final<-BPlasma$seurat_cluster_final
#Violin plot for selected marker genes by clusters
pdf("~/Krummel lab/COVID/redo/vln_final.pdf" , width = 6 , height = 16)
VlnPlot(BPlasma, features = c('MS4A1', 'CD74', 'CD79A', 'TCL1A', 'JCHAIN', 'CD38', 'IGHD', 'CD27',"IFIT3","IFITM3"), pt.size = 0, assay = 'RNA', ncol = 2,cols =  brewer.pal(n=7,"Paired")) + NoLegend()
dev.off()
#Final UMAP with annotated clusters
pdf("~/Krummel lab/COVID/redo/b_plasma_harmony_UMAP.pdf" , width = 5 , height = 5,useDingbats = FALSE)
DimPlot(BPlasma, reduction = "umap")+scale_color_manual(values = brewer.pal(n=7,"Paired")) + NoLegend()
dev.off()
##################################################################################
#final dotplot 
pdf("b_plasma_harmony_0.4_DotPlot.pdf", width = 7, height = 7)
Idents(b_plasma_harmony_0.2) <- 'seurat_cluster_final'
reorder_clusters <- c('naive B cell_1', 'naive B cell_2', 'naive B cell_3', 'memory B cell', 'unidentified', 'Plasma_1', 'Plasma_2', 'Plasma_3')
levels(b_plasma_harmony_0.4) <- reorder_clusters
DotPlot(b_plasma_harmony_0.4,
        features = unique(b_plasma_harmony_0.4_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
dev.off()

## Author: Alexis Combes

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
setwd("~/Alexis/Manips/COVID/200428 patient_data_SO/201001 SEURAT WD/")

## Call clinical scores csv file and add clinical value as metadata in the object
COMET_10X_CLINICAL_SCORES <- read.csv("~/Desktop/COMET_10X_CLINICAL_SCORES_PAPER.csv", sep=',', header=T)
monodc_clean_1.5NWUMAP@meta.data$covid_status <- COMET_10X_CLINICAL_SCORES$covid_status[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$Age <- COMET_10X_CLINICAL_SCORES$Age[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$Day_onset_category <- COMET_10X_CLINICAL_SCORES$Day_onset_category[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$ICU_vs_FLOOR <- COMET_10X_CLINICAL_SCORES$ICU_vs_FLOOR[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$Other_infection_type <- COMET_10X_CLINICAL_SCORES$Other_infection_type[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$Death <- COMET_10X_CLINICAL_SCORES$Death[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$NIH_score <- COMET_10X_CLINICAL_SCORES$NIH_score[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$Sampling_score <- COMET_10X_CLINICAL_SCORES$Sampling_score[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$Overall_score <- COMET_10X_CLINICAL_SCORES$Overall_score[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$Qualitative_score <- COMET_10X_CLINICAL_SCORES$Qualitative_score[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$Spike_score <- COMET_10X_CLINICAL_SCORES$Spike_score[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$Nucleocapsid_score <- COMET_10X_CLINICAL_SCORES$Nucleocapsid_score[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$ORF3a_score <- COMET_10X_CLINICAL_SCORES$ORF3a_score[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]
monodc_clean_1.5NWUMAP@meta.data$RBD_score <- COMET_10X_CLINICAL_SCORES$RBD_score[match(monodc_clean_1.5NWUMAP@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES$SAMPLE.by.SNPs)]

# Create color scheme for COVID-19 status and disease severity
mycolorsseverity <- setNames(c("grey40", "orange", "orangered2"), c('CTRL', 'MILD', 'SEVERE'))
mycolorsstatus <- setNames(c("grey40", "dodgerblue3", "firebrick3"), c('CTRL', 'NEG', 'POS'))

#load data
load("~/Alexis/Manips/COVID/200428 patient_data_SO/merged/merged_colossal_paper_final/annotated_res_0.6/monos/merged_colossal_paper_final_res0.6_monodc.RData")

sobj <- RunUMAP(sobj, 
                dims = pcs_to_use,
                reduction='harmony',
                n.neighbors = 20,
                min.dist = 0.1, 
                spread = 2,
                verbose = FALSE)
sobj <- FindNeighbors(sobj,
                      dims = pcs_to_use,  # Num PCs to use
                      reduction='harmony',
                      k.param = 20,  # k for the knn algorithm
                      verbose = FALSE)

pcs_to_use = c()
for (pc in 1:30) {
  top_10 <- names(sort(monodc_harmony_subset@reductions$harmony@feature.loadings[,pc]))[1:10]
  bottom_10 <- names(sort(monodc_harmony_subset@reductions$harmony@feature.loadings[,pc], decreasing=T))[1:10]
  if ((sum(top_10 %in% mito_genes) + sum(top_10 %in% ribo_genes)) >= 6){
    next
  } else if ((sum(bottom_10 %in% mito_genes) + sum(bottom_10 %in% ribo_genes))>= 6){
    next
  } else {
    pcs_to_use <- c(pcs_to_use, pc)
  }
}
sobj <- RunUMAP(sobj, 
                dims = pcs_to_use,
                reduction='harmony',
                n.neighbors = 30,
                min.dist = 0.3, 
                spread = 1,
                verbose = FALSE)


Marker_monoDC_refresh_0.2 <- FindAllMarkers(monodc_select, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.40, test.use="poisson", latent.vars = "LIBRARY", assay = "RNA")
Marker_monoDC_refresh_0.2_padj0.001 <- Marker_monoDC_refresh_0.2[which(Marker_monoDC_refresh_0.2$p_val_adj<0.001),]
Marker_monoDC_refresh_0.2_padj0.001<- Marker_monoDC_refresh_0.2_padj0.001[order(Marker_monoDC_refresh_0.2_padj0.001$avg_logFC,decreasing = TRUE),]
Marker_monoDC_refresh_0.2_padj0.001<- Marker_monoDC_refresh_0.2_padj0.001[order(Marker_monoDC_refresh_0.2_padj0.001$cluster,decreasing = FALSE),]
Idents(monodc_select) <- monodc_select@meta.data$integrated_snn_res.0.4
Marker_monoDC_refresh_0.4<- FindAllMarkers(monodc_select,logfc.threshold = 0.4,min.pct = 0.25,only.pos = TRUE,test.use = "poisson",latent.vars = "LIBRARY",assay = "RNA")
Marker_monoDC_refresh_0.4_padj0.001 <- Marker_monoDC_refresh_0.4[which(Marker_monoDC_refresh_0.4$p_val_adj<0.001),]
Marker_monoDC_refresh_0.4_padj0.001<- Marker_monoDC_refresh_0.4_padj0.001[order(Marker_monoDC_refresh_0.4_padj0.001$avg_logFC,decreasing = TRUE),]
Marker_monoDC_refresh_0.4_padj0.001<- Marker_monoDC_refresh_0.4_padj0.001[order(Marker_monoDC_refresh_0.4_padj0.001$cluster,decreasing = FALSE),]
Idents(monodc_select) <- monodc_select@meta.data$integrated_snn_res.0.6
Marker_monoDC_refresh_0.6<- FindAllMarkers(monodc_select,logfc.threshold = 0.4,min.pct = 0.25,only.pos = TRUE,test.use = "poisson",latent.vars = "LIBRARY",assay = "RNA")
Marker_monoDC_refresh_0.6_padj0.001 <- Marker_monoDC_refresh_0.6[which(Marker_monoDC_refresh_0.6$p_val_adj<0.001),]
Marker_monoDC_refresh_0.6_padj0.001<- Marker_monoDC_refresh_0.6_padj0.001[order(Marker_monoDC_refresh_0.6_padj0.001$avg_logFC,decreasing = TRUE),]
Marker_monoDC_refresh_0.6_padj0.001<- Marker_monoDC_refresh_0.6_padj0.001[order(Marker_monoDC_refresh_0.6_padj0.001$cluster,decreasing = FALSE),]
Idents(monodc_select) <- monodc_select@meta.data$integrated_snn_res.0.8
Marker_monoDC_refresh_0.8<- FindAllMarkers(monodc_select,logfc.threshold = 0.4,min.pct = 0.25,only.pos = TRUE,test.use = "poisson",latent.vars = "LIBRARY",assay = "RNA")
Marker_monoDC_refresh_0.8_padj0.001 <- Marker_monoDC_refresh_0.8[which(Marker_monoDC_refresh_0.8$p_val_adj<0.001),]
Marker_monoDC_refresh_0.8_padj0.001<- Marker_monoDC_refresh_0.8_padj0.001[order(Marker_monoDC_refresh_0.8_padj0.001$avg_logFC,decreasing = TRUE),]
Marker_monoDC_refresh_0.8_padj0.001<- Marker_monoDC_refresh_0.8_padj0.001[order(Marker_monoDC_refresh_0.8_padj0.001$cluster,decreasing = FALSE),]

mono_subset_final_0.6<- FindAllMarkers(mono_subset_final_0.6,logfc.threshold = 0.25,min.pct = 0.25,assay = "RNA",latent.vars = "LIBRARY",test.use = "poisson",only.pos = TRUE)
mono_subset_final_0.6_marker_0.001 <- mono_subset_final_0.6_marker[which(mono_subset_final_0.6_marker$p_val_adj<0.001),]
mono_subset_final_0.6_marker_0.001<- mono_subset_final_0.6_marker_0.001[order(mono_subset_final_0.6_marker_0.001$avg_logFC,decreasing = TRUE),]
mono_subset_final_0.6_marker_0.001<- mono_subset_final_0.6_marker_0.001[order(mono_subset_final_0.6_marker_0.001$cluster,decreasing = FALSE),]

monodc_clean_1.5 <- RunUMAP(monodc_clean_1.5, dims = 1:30, n.neighbors = 30, min.dist = 0.3, spread = 1, verbose = FALSE, seed.use = 21212, reduction = 'harmony')
monodc_clean<- FindNeighbors(monodc_clean, dims = 1:30, k.param = 20,  verbose = FALSE, reduction = 'harmony')
DimPlot(monodc_clean_1.5)

#Change Spread and and min dist if you want to play with UMAP if increase min dist cluster will get closer same for spread
monodc_clean_1.5NWUMAP<- RunUMAP(monodc_clean_1.5,dims = 1:30, n.neighbors = 30, min.dist = 3, spread = 3, verbose = FALSE, seed.use = 21212, reduction = 'harmony')
my_level_2<- c("MONO_Classical monocytes","MONO_Classical monocytes ISG","MONO_non Classical monocytes","MONO_non Classical monocytes C1Q+","MONO_DCs","MONO_pDCs","MONO_Classical monocytes S100A12")

DimPlot(monodc_clean_1.5NWUMAP, reduction='umap', label = F, cols = "Paired", order = c('MONO_Classical monocytes S100A12', 'MONO_pDCs', 'MONO_DCs', 'MONO_non classical monocytes C1Q+', 'MONO_Classical monocytes ISG', 'MONO_non classical monocytes', 'MONO_Classical monocytes'))

#LOADING GENE SET OF INTEREST
OXPHOS_genesig<- scan(file = "~/Box/Alexis Data/Gene Set used/OXPHOS.txt",what = "character")
Glycolysis_genesig<-scan(file = "~/Box/Alexis Data/Gene Set used/Glycolysis_Zenith.txt",what = "character") 
Mono_clasical_genesig<-scan(file = "~/Box/Alexis Data/Gene Set used/Mono1_Blood_vilanni.txt",what = "character") 
Mono_non_clasical_genesig<-scan(file = "~/Box/Alexis Data/Gene Set used/Mono2_Blood_vilanni.txt",what = "character")
DCs_genesig<-scan(file = "~/Box/Alexis Data/Gene Set used/DCs_Blood_vilanni.txt",what = "character")
pDCs_genesig<-scan(file = "~/Box/Alexis Data/Gene Set used/pDCs_Blood_vilanni.txt",what = "character")
#CALCULATING MONO CLASSICAL FROM VILlani et al 2017
gene_df <- get_genes_s3(Mono_clasical_genesig, monodc_clean_1.5NWUMAP,drop = TRUE)
# Define the score vector. In this case, we call the score for a cell the mean value of all genes in the vector
score_monoclassical <- colMeans(gene_df)
# load in the metadata
monodc_clean_1.5NWUMAP@meta.data$scoremonoClass <- saturate(vec=score_monoclassical, sat=0.70, binary=FALSE)
# Generate a feature plot with the score
FeaturePlot(monodc_clean_1.5NWUMAP,features = "scoremonoClass",cols = c("grey","red"))

#Export counts used to create prism graphs
merged_colossal_full_harmony_MP_Counts <- table(monodc_clean_1.5NWUMAP@meta.data[,c('SAMPLE.by.SNPs', 'harmony_cluster_final','covid_status','Qualitative_score')])
write.table(merged_colossal_full_harmony_MP_Counts, file="~/Desktop/merged_colossal_full_harmony_MP.tsv", row.names=T, col.names=T, sep="\t")
write.table(ClustersByCovidStatus,file = "~/Desktop/merged_colossal_full_harmony_MP.tsv",row.names=T, col.names=T, sep="\t")


#CALCULATING MONO Non CLASSICAL FROM VILlani et al 2017
gene_df <- get_genes_s3(Mono_non_clasical_genesig, monodc_clean_1.5NWUMAP,drop = TRUE)
# Define the score vector. In this case, we call the score for a cell the mean value of all genes in the vector
score_monononclassical <- colMeans(gene_df)
# load in the metadata
monodc_clean_1.5NWUMAP@meta.data$scoremonoNonClass <- saturate(vec=score_monononclassical, sat=0.70, binary=FALSE)
# Generate a feature plot with the score
FeaturePlot(monodc_clean_1.5NWUMAP,features = "scoremonoNonClass",cols = c("grey","red"))

#CALCULATING DCs score FROM VILlani et al 2017
gene_df <- get_genes_s3(DCs_genesig, monodc_clean_1.5NWUMAP,drop = TRUE)
# Define the score vector. In this case, we call the score for a cell the mean value of all genes in the vector
score_DCs <- colMeans(gene_df)
# load in the metadata
monodc_clean_1.5NWUMAP@meta.data$scoreDCs <- saturate(vec=score_DCs, sat=0.70, binary=FALSE)
# Generate a feature plot with the score
FeaturePlot(monodc_clean_1.5NWUMAP,features = "scoreDCs",cols = c("grey","red"))

#CALCULATING pDCs score FROM VILlani et al 2017
gene_df <- get_genes_s3(pDCs_genesig, monodc_clean_1.5NWUMAP,drop = TRUE)
# Define the score vector. In this case, we call the score for a cell the mean value of all genes in the vector
score_pDCs <- colMeans(gene_df)
# load in the metadata
monodc_clean_1.5NWUMAP@meta.data$scorepDCs <- saturate(vec=score_pDCs, sat=0.70, binary=FALSE)
# Generate a feature plot with the score
FeaturePlot(monodc_clean_1.5NWUMAP,features = "scorepDCs",cols = c("grey","red"))

#CALCULATING Glycolysis score FROM Arguello et al 2020
gene_df <- get_genes_s3(Glycolysis_genesig, monodc_clean_1.5NWUMAP,drop = TRUE)
# Define the score vector. In this case, we call the score for a cell the mean value of all genes in the vector
score_Glyco <- colMeans(gene_df)
# load in the metadata
monodc_clean_1.5NWUMAP@meta.data$scoreGlyco <- saturate(vec=score_Glyco, sat=0.70, binary=FALSE)
# Generate a feature plot with the score
FeaturePlot(monodc_clean_1.5NWUMAP,features = "scoreGlyco",cols = c("grey","red"))

###   BiocManager::install('EnhancedVolcano')
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
#Volcano_plot
EnhancedVolcano(gene_list,
                lab = rownames(gene_list),
                x = 'avg_logFC',
                y = 'p_val',
                xlim = c(-1.5,1.5),
                title = 'APC MILD VS SEVERE',
                pCutoff = 10e-20,
                FCcutoff = 0.2,
                pointSize = 2.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                selectLab = c('IFITM1','IFITM3','ISG15',
                              'S100A12','LYZ','LY6E','S100A9','IFI27','OAS1','IRF7','STAT1','TNFSF10','S100A8','CD52'),
                colAlpha = 1)



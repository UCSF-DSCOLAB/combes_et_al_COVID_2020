## Author: Nicholas Kuhn

# Call libraries
library(assertthat)  # For sanity checks
library(ggplot2)     # For pretty plots
library(cowplot)     # For pretty plots
library(dittoSeq)    # For pretty, colorblind friendly plots
library(dplyr)       # For inline modification of matrices
library(grid)        # For plotting multiple plots in one frame
library(gridExtra)   # For plotting multiple plots in one frame
library(reshape2)    # For "melting" dataframesb
library(scales)      # To access break formatting functions
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(Seurat)

# Load the data
load("~/merged_colossal_paper_final_res0.6_platelets_reHarmony.RData")

platelet_reHarmony <- platelets

platelet_reHarmony <- RunUMAP(platelet_reHarmony, dims = 1:30, n.neighbors = 30, min.dist = 0.3, spread = 1, verbose = FALSE, seed.use = 21212, reduction = 'harmony')
platelet_reHarmony <- FindNeighbors(platelet_reHarmony, dims = 1:30, k.param = 20,  verbose = FALSE, reduction = 'harmony')
platelet_reHarmony_0.4 <- FindClusters(platelet_reHarmony, verbose = TRUE, algorithm = 1, resolution = 0.4, random.seed = 21212)  
##### Figure 3B
DimPlot(platelet_reHarmony_0.4 , reduction='umap', label = T) + labs(color = "0.4")

platelet_reHarmony_0.4_Markers <- FindAllMarkers(platelet_reHarmony_0.4,
                                               test.use='poisson',
                                               only.pos=TRUE,
                                               min.pct=0.25,
                                               logfc.threshold=0.4,
                                               assay='RNA',
                                               latent.vars = 'LIBRARY'
)
platelet_reHarmony_0.4_Markers_padj0.1 <- platelet_reHarmony_0.4_Markers[which(platelet_reHarmony_0.4_Markers$p_val_adj<0.1),]
platelet_reHarmony_0.4_Markers_padj0.1 <- platelet_reHarmony_0.4_Markers_padj0.1[order(platelet_reHarmony_0.4_Markers_padj0.1$avg_logFC,decreasing = TRUE),]
platelet_reHarmony_0.4_Markers_padj0.1 <- platelet_reHarmony_0.4_Markers_padj0.1[order(platelet_reHarmony_0.4_Markers_padj0.1$cluster,decreasing = FALSE),]
platelet_reHarmony_0.4_Markers_padj0.1_Top10 <- platelet_reHarmony_0.4_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

platelet_reHarmony_0.4_Markers_padj0.1_Top5 <- platelet_reHarmony_0.4_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
Idents(platelet_reHarmony_0.4) <- 'seurat_cluster_final'
DotPlot(platelet_reHarmony_0.4,
        features = unique(platelet_reHarmony_0.4_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 8), text = element_text(size = 14)) + coord_flip()

### add column harmony_cluster_final with annotation
Idents(platelet_reHarmony_0.4) <- 'seurat_clusters'
platelet_reHarmony_0.4$harmony_cluster_final <- platelet_reHarmony_0.4$seurat_clusters
### rename clusters in column "seurat_cluster_final"
annotations <- c(
  "Platelet_ACTB", # 0
  "Platelet_HIST1H2AC", # 1
  "Platelet_PPBP", # 2
  "Platelet_H3F3B", # 3
  "Platelet_TMSB4X", #4
  "Platelet_RGS18" #5
)
names(annotations) <- levels(platelet_reHarmony_0.4)
platelet_reHarmony_0.4 <- RenameIdents(platelet_reHarmony_0.4, annotations)
platelet_reHarmony_0.4@meta.data$harmony_cluster_final <- Idents(platelet_reHarmony_0.4)
table(platelet_reHarmony_0.4@meta.data$harmony_cluster_final)

Idents(platelet_reHarmony_0.4) <- 'harmony_cluster_final'
reorder_clusters <- c('Platelet_H3F3B', 'Platelet_HIST1H2AC', 'Platelet_ACTB', 'Platelet_PPBP', 'Platelet_TMSB4X', 'Platelet_RGS18')
levels(platelet_reHarmony_0.4) <- reorder_clusters
##### Figure 3A
DotPlot(platelet_reHarmony_0.4,
        features = unique(platelet_reHarmony_0.4_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 8), text = element_text(size = 14)) + coord_flip()

##### Figures 3D and S5A
VlnPlot(platelet_reHarmony_0.4, features = c('BCL2L1'), pt.size = 0, assay = 'RNA', ncol = 1) + NoLegend()
VlnPlot(platelet_reHarmony_0.4, features = c('TBXA2R'), pt.size = 0, assay = 'RNA', ncol = 1) + NoLegend()
VlnPlot(platelet_reHarmony_0.4, features = c('P2RX1'), pt.size = 0, assay = 'RNA', ncol = 1) + NoLegend()
VlnPlot(platelet_reHarmony_0.4, features = c('PTGIR'), pt.size = 0, assay = 'RNA', ncol = 1) + NoLegend()

FeaturePlot(platelet_reHarmony_0.4, features = c('BCL2L1'), cols = c("grey", "red"), pt.size = 0.3, min.cutoff = 0.2)

########################################################
##################plot frequencies of clusters by Qualitative score
########################################################
######## add annotation based on clinical scores
COMET_10X_CLINICAL_SCORES_PAPER <- read.csv("~/COMET_10X_CLINICAL_SCORES_PAPER.csv", sep=',', header=T)
platelet_reHarmony_0.4@meta.data$Qualitative_score <- COMET_10X_CLINICAL_SCORES_PAPER$Qualitative_score[match(platelet_reHarmony_0.4@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES_PAPER$SAMPLE.by.SNPs)]

#create file 'y' with table containing SAMPLE.by.SNPs as rows and seurat_clusters as columns with corresponding event #s
y <- table(platelet_reHarmony_0.4@meta.data[,c('SAMPLE.by.SNPs', 'harmony_cluster_final')])

#load reshape2
library(reshape2)

#melt table into one column by id='SAMPLE.by.SNPs; call 'ClustersByCovidStatus'
ClustersByQualitativeScore <- melt(y, id='SAMPLE.by.SNPs')

#make data.frame from sobj that contains SAMPLE.by.SNPs event #s; call 'temp'
temp <- table(platelet_reHarmony_0.4@meta.data$SAMPLE.by.SNPs)
temp <- data.frame(temp)

#add a column to ClustersByQualitativeScore called 'total' that contains total # of cells per SAMPLE.by.SNPs pulled from 'temp' data.frame in column 'Freq'
ClustersByQualitativeScore$total <- temp$Freq[match(ClustersByQualitativeScore$SAMPLE.by.SNPs, temp$Var1)]
#add a column to ClustersByQualitativeScore called 'freq' that gives the fraction of cells of each harmony_cluster_final per total # of cells in sample
ClustersByQualitativeScore$freq <- ClustersByQualitativeScore$value/ClustersByQualitativeScore$total
ggplot(ClustersByQualitativeScore, aes(x=harmony_cluster_final, y=freq))+geom_point(stat = 'identity')

#generate temp2 dataframe that has Qualitative_score matched to SAMPLE.by.SNPs
temp2 <- platelet_reHarmony_0.4@meta.data[,c('SAMPLE.by.SNPs', 'Qualitative_score')]
temp2 <- data.frame(temp2)

#add column to ClustersByQualitativeScore called 'Qualitative_score' that adds Ctrl/moderate/severe/critical to ClustersByQualitativeScore matched to SAMPLE.by.SNPs pulled from 'temp2' data.frame in column 'Qualitative_score'
ClustersByQualitativeScore$Qualitative_score <- temp2$Qualitative_score[match(ClustersByQualitativeScore$SAMPLE.by.SNPs, temp2$SAMPLE.by.SNPs)]
#make 'harmony_cluster_final' numbers factors, not values; or else ggplot gets confused
#and reorder clusters with "levels" command
ClustersByQualitativeScore$harmony_cluster_final <- factor(ClustersByQualitativeScore$harmony_cluster_final, levels=c('Platelet_ACTB', 'Platelet_HIST1H2AC', 'Platelet_PPBP', 'Platelet_H3F3B', 'Platelet_TMSB4X', 'Platelet_RGS18'))
write.table(ClustersByQualitativeScore, file="ClustersByQualitativeScore.tsv", row.names=T, col.names=T, quote=F, sep="\t")

#plot
mycolorsseverity <- setNames(c("grey40", "orange", "orangered2"), c('CTRL', 'MILD', 'SEVERE'))
e <- ggplot(ClustersByQualitativeScore, aes(x=harmony_cluster_final, y=freq))

##### Figure 3C
e + geom_boxplot(
  aes(fill = Qualitative_score),
  outlier.size = 0.5,
  position = position_dodge(0.8)) +
  labs(title = NULL,x="harmony_cluster_final", y = "frequency") +
  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14), panel.background = element_rect(fill = "white",size = 2, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid'), axis.line = element_line(colour = "black")) +
  theme_classic() +
  scale_fill_manual(values=mycolorsseverity)

########################################################
##################plot frequencies of clusters by status
########################################################
#create file 'y' with table containing SAMPLE.by.SNPs as rows and seurat_clusters as columns with corresponding event #s
y <- table(plateletCLEAN_harmony_0.4@meta.data[,c('SAMPLE.by.SNPs', 'seurat_cluster_final')])

#load reshape2
library(reshape2)

#melt table into one column by id='SAMPLE.by.SNPs; call 'ClustersByCovidStatus'
ClustersByCovidStatus <- melt(y, id='SAMPLE.by.SNPs')

#make data.frame from sobj that contains SAMPLE.by.SNPs event #s; call 'temp'
temp <- table(plateletCLEAN_harmony_0.4@meta.data$SAMPLE.by.SNPs)
temp <- data.frame(temp)

#add a column to ClustersByCovidStatus called 'total' that contains total # of cells per SAMPLE.by.SNPs pulled from 'temp' data.frame in column 'Freq'
ClustersByCovidStatus$total <- temp$Freq[match(ClustersByCovidStatus$SAMPLE.by.SNPs, temp$Var1)]
#add a column to ClustersByCovidStatus called 'freq' that gives the fraction of cells of each seurat_cluster_final per total # of cells in sample
ClustersByCovidStatus$freq <- ClustersByCovidStatus$value/ClustersByCovidStatus$total
ggplot(ClustersByCovidStatus, aes(x=seurat_cluster_final, y=freq))+geom_point(stat = 'identity')

#generate temp2 dataframe that has covid_status matched to SAMPLE.by.SNPs
temp2 <- plateletCLEAN_harmony_0.4@meta.data[,c('SAMPLE.by.SNPs', 'covid_status')]
temp2 <- data.frame(temp2)

#add column to ClustersByCovidStatus called 'covid_status' that adds neg/neg/ctrl to ClustersByCovidStatus matched to SAMPLE.by.SNPs pulled from 'temp2' data.frame in column 'covid_status'
ClustersByCovidStatus$covid_status <- temp2$covid_status[match(ClustersByCovidStatus$SAMPLE.by.SNPs, temp2$SAMPLE.by.SNPs)]
#make 'seurat_cluster_final' numbers factors, not values; or else ggplot gets confused
#and reorder clusters with "levels" command
ClustersByCovidStatus$seurat_cluster_final <- factor(ClustersByCovidStatus$seurat_cluster_final, levels=c('H3F3B', 'PPBP', 'ACTB', 'TMSB4X', 'HIST1H2AC', 'RGS18'))

#plot
mycolorsstatus <- setNames(c("grey40", "dodgerblue4", "firebrick3"), c('ctrl', 'neg', 'pos'))
e <- ggplot(ClustersByCovidStatus, aes(x=seurat_cluster_final, y=freq))

e + geom_boxplot(
  aes(fill = covid_status),
  outlier.size = 0.5,
  position = position_dodge(0.8)) +
  labs(title = NULL,x="seurat_cluster_final", y = "frequency") +
  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14), panel.background = element_rect(fill = "white",size = 2, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid'), axis.line = element_line(colour = "black")) +
  theme_classic() +
  scale_fill_manual(values=mycolorsstatus)
dev.off()


######## ISG Score
### load functions 'get_genes_s3' and 'saturate'

get_genes_s3 <- function(genes, sobj, drop=FALSE) {
  ###
  # DESCRIPTION
  # Get a dataframe of gene values per cell from the input seurat version3 object
  #
  # INPUTS
  # genes: A character vector containing gene names
  # sobj: The Seurat object
  # OUTPUT
  # A dataframe of input genes as rows and all cells as columns
  ###
  gene_idx <- sapply(genes, function(x) { match(x, rownames(sobj@assays[["RNA"]]@data)) })
  if (sum(is.na(gene_idx)) > 0) {
    print("The following genes are not in the gene list.")
    print(names(gene_idx[is.na(gene_idx)]))
    if (drop){
      drop_names <- names(gene_idx[is.na(gene_idx)])
      genes <- genes[!genes%in%drop_names]
    } else {
      return(1)  
    }
  }
  genes <- as.data.frame(as.matrix(sobj@assays[[sobj@active.assay]]@data[genes,]))
}

saturate <- function(vec, sat=0, binary=FALSE){
  ###
  # DESCRIPTION
  # A Function to convert a vector of scores into a saturated vectore of scores. A saturated vector is one where all values below the 
  # provided "saturation" (percentile of data) are set to 0. If the binary flag is specified, all values greater than or equal to the
  # saturation will be set to 1.
  #
  # INPUTS
  # vec: A numeric vector of scores
  # sat: A value to saturate the vector at (float (0.0-1.0) or percent (1.0-100.0))
  # binary: A flag to indicate if we should make the output vector a binary one.
  #
  # OUTPUT
  # A vector of saturated scores
  ###
  sat = if (sat > 1.0) sat/100 else sat
  z <- quantile(vec, sat)
  for (i in 1:length(vec)){
    if (vec[i] < z) {
      vec[i] = 0
    } else if(binary) {
      vec[i] = 1
    }
  }
  vec
}

platelet_reHarmony_0.4@meta.data$Qualitative_score <- COMET_10X_CLINICAL_SCORES_PAPER$Qualitative_score[match(platelet_reHarmony_0.4@meta.data$SAMPLE.by.SNPs,COMET_10X_CLINICAL_SCORES_PAPER$SAMPLE.by.SNPs)]

# Adding IF-responsive gene signature scores and plotting
# A list of genes in the signature
ISG_genes = c('MT2A', 'ISG15', 'LY6E', 'IFIT1', 'IFIT2', 'IFIT3', 'IFITM1', 'IFITM3', 'IFI44L', 'IFI6', 'MX1', 'IFI27')
# Get the cell values from the Seurat object
gene_df <- get_genes_s3(ISG_genes, platelet_reHarmony_0.4, drop = T)
# Define the score vector. In this case, we call the score for a cell the mean value of all genes in the vector
score <- colMeans(gene_df)
# Add score columns to the metadata column of the Seurat object
platelet_reHarmony_0.4@meta.data$ISG_Score_0 <- saturate(vec=score, sat=0, binary=FALSE) # saturate at 0%
# Make a featureplot
##### Figure S5C
FeaturePlot(platelet_reHarmony_0.4,  features = "ISG_Score_0",
            cols =c("light grey", "red"), pt.size = 0.1, min.cutoff = 0.1)
##### Figure 3G
VlnPlot(platelet_reHarmony_0.4, features = "ISG_Score_0", group.by = "covid_status", split.by = 'Qualitative_score', pt.size = 0, col = mycolorsseverity) + NoLegend()

###Wilcoxon rank sum test on ISG score
#subset all groups
Idents(platelet_reHarmony_0.4) <- 'covid_status'
platelet_reHarmony_0.4_pos <- subset(platelet_reHarmony_0.4, idents = 'pos')
platelet_reHarmony_0.4_neg <- subset(platelet_reHarmony_0.4, idents = 'neg')
Idents(platelet_reHarmony_0.4_pos) <- 'Qualitative_score'
platelet_reHarmony_0.4_pos_mild <- subset(platelet_reHarmony_0.4_pos, idents = 'MILD')
platelet_reHarmony_0.4_pos_severe <- subset(platelet_reHarmony_0.4_pos, idents = 'SEVERE')
Idents(platelet_reHarmony_0.4_neg) <- 'Qualitative_score'
platelet_reHarmony_0.4_neg_mild <- subset(platelet_reHarmony_0.4_neg, idents = 'MILD')
platelet_reHarmony_0.4_neg_severe <- subset(platelet_reHarmony_0.4_neg, idents = 'SEVERE')
#stats
NM_ISG <- platelet_reHarmony_0.4_neg_mild@meta.data$ISG_Score_0
PM_ISG <- platelet_reHarmony_0.4_pos_mild@meta.data$ISG_Score_0
my_data <- data.frame(group = c(rep("NEG_MILD",length(NM_ISG)),
                                rep("POS_MILD",length(PM_ISG))),
                      ISG_score = c(NM_ISG,  PM_ISG))
group_by(my_data, group) %>%
  summarise(
    count = n(),
    median = median(ISG_score, na.rm = TRUE),
    IQR = IQR(ISG_score, na.rm = TRUE))
p <- ggplot(my_data, aes(x=group, y=ISG_score)) +
  geom_violin()
p
wilcox.test(ISG_score ~ group, data = my_data,
            exact = FALSE)



## Author: Will Chen

savedest <- 'C:/Users/wilyu/OneDrive/Desktop/Research/ImmunoX/COVID/Results/2020-10-01/'
library(Seurat)

# Load single cell cluster assignments and sample of origin data
bplasma.metadata <- read.delim(file='C:/Users/wilyu/OneDrive/Desktop/Research/Datasets/Krummel_COVID_2020/celltype_labels_subset/2020-10-01/BPlasma_metadata.tsv', header=T, sep='\t', stringsAsFactors=F)
tnk.metadata <- read.delim(file='C:/Users/wilyu/OneDrive/Desktop/Research/Datasets/Krummel_COVID_2020/celltype_labels_subset/2020-10-01/TNK_Final_metadata.tsv', header=T, sep='\t', stringsAsFactors=F)
monomac.metadata <- read.delim(file='C:/Users/wilyu/OneDrive/Desktop/Research/Datasets/Krummel_COVID_2020/celltype_labels_subset/2020-10-01/merged_monodc_metadata.tsv', header=T, sep='\t', stringsAsFactors=F)
neut.metadata <- read.delim(file='C:/Users/wilyu/OneDrive/Desktop/Research/Datasets/Krummel_COVID_2020/celltype_labels_subset/2020-10-01/201001_NEUTS_METADATA_FOR_KEN_WILL.tsv', header=T, sep='\t', stringsAsFactors=F)
plt.metadata <- read.delim(file='C:/Users/wilyu/OneDrive/Desktop/Research/Datasets/Krummel_COVID_2020/celltype_labels_subset/2020-10-01/platelet_reHarmony_0.4_metadata.tsv', header=T, sep='\t', stringsAsFactors=F)

rownames(bplasma.metadata) <- bplasma.metadata$X
bplasma.metadata <- bplasma.metadata[,-1]
tnk.metadata$harmony_cluster_final <- tnk.metadata$harmony_final_clusters

all.metadata <- read.delim(file='C:/Users/wilyu/OneDrive/Desktop/Research/Datasets/Krummel_COVID_2020/celltype_labels_subset/2020-10-01/merged_colossal_paper_final_merged_annotated_metadata.tsv', header=T, sep='\t', stringsAsFactors=F)

head(all.metadata)

# Load UMAP coordinate data
load('C:/Users/wilyu/OneDrive/Desktop/Research/Datasets/Krummel_COVID_2020/celltype_labels_subset/2020-10-01/merged_colossal_paper_final_merged_annotated.RData')
umap_coords <- Embeddings(merged_data[['umap']])
harmony_coords <- Embeddings(merged_data[['harmony']])

# reorder all.metadata rows to match that of the RData object
all.metadata <- all.metadata[match(rownames(umap_coords), rownames(all.metadata)),]
stopifnot(identical(match(rownames(all.metadata), rownames(umap_coords)), 1:nrow(umap_coords)))

all.metadata <- cbind(all.metadata, harmony_coords)

# Map from split datasets to master dataset
neut2full_map <- match(rownames(neut.metadata), rownames(all.metadata))
tnk2full_map <- match(rownames(tnk.metadata), rownames(all.metadata))
monomac2full_map <-match(rownames(monomac.metadata), rownames(all.metadata))
bplasma2full_map <- match(rownames(bplasma.metadata), rownames(all.metadata))
plt2full_map <- match(rownames(plt.metadata), rownames(all.metadata))

# Assert that all cells were successfully mapped
stopifnot(sum(is.na(neut2full_map)) == 0)
stopifnot(sum(is.na(tnk2full_map)) == 0)
stopifnot(sum(is.na(monomac2full_map)) == 0)
stopifnot(sum(is.na(bplasma2full_map)) == 0)
stopifnot(sum(is.na(plt2full_map)) == 0)


allcells.mapped <- c(neut2full_map, tnk2full_map, monomac2full_map, bplasma2full_map, plt2full_map)
nrow(all.metadata) #total cells before further filtering
length(allcells.mapped) # number of final cells post-filtering


# final set of all cells post-filtering
sc.metadata.final <- all.metadata[allcells.mapped,]
sc.metadata.final$cluster_final <- c(neut.metadata$harmony_cluster_final, tnk.metadata$harmony_cluster_final, monomac.metadata$harmony_cluster_final, bplasma.metadata$harmony_cluster_final, plt.metadata$harmony_cluster_final)

####################
## 1. Plot cell state UMAP for all cell types (final clusters, post-filtering)
##############
# 2D UMAP from Seurat
umap.coords <- sc.metadata.final[,c('UMAP_1', 'UMAP_2')]
ndim.umap <- 3 # use intrinsic dimension of 3

# immune cells only i.e. excluding platelets
sc.metadata.final <- sc.metadata.final[-which(sc.metadata.final$cluster_final %in% c("Platelet_ACTB","Platelet_PPBP", "Platelet_H3F3B" , "Platelet_TMSB4X" ,"Platelet_HIST1H2AC","Platelet_RGS18" )),]

# Run 3D umap on filtered data
merged_colossal_full_harmony_final <- subset(merged_data[,rownames(sc.metadata.final)])
merged_colossal_full_harmony_final <- RunUMAP(merged_colossal_full_harmony_final, reduction='harmony', dims = 1:30, min.dist = 0.4, n.components=ndim.umap)
save(merged_colossal_full_harmony_final, file=paste(savedest, 'UMAP3d.RData', sep=''))

#################
# # 3D UMAP
#############
load(paste(savedest, 'UMAP3d.RData', sep=''))
library(scatterplot3d)
umap_coords_3d <- Embeddings(merged_colossal_full_harmony_final[['umap']])
colnames(umap_coords_3d) <- paste('3D', colnames(umap_coords_3d), sep='')

# Add 3D umap coordinates to master table
sc.metadata.final <- cbind(sc.metadata.final, umap_coords_3d)

# Get per-patient cell count matrix 
sample.ids <- unique(sc.metadata.final$SAMPLE.by.SNPs)
master.features <- names(table(sc.metadata.final$cluster_final))
pt.cellcount <- matrix(NA,nrow=length(sample.ids), ncol=length(master.features))
for(curidx in seq_len(length(sample.ids))) {
  cur.id <- sample.ids[curidx]
  sc.metadata.curpt <- sc.metadata.final[which(sc.metadata.final$SAMPLE.by.SNPs == cur.id),]
  curpt.counts <- table(sc.metadata.curpt$cluster_final)
  curpt2master.order <- match(master.features, names(curpt.counts))
  curpt.counts <- as.numeric(curpt.counts)[curpt2master.order]
  curpt.counts[is.na(curpt.counts)] <- 0
  pt.cellcount[curidx,] <- curpt.counts
}
colnames(pt.cellcount) <- master.features
rownames(pt.cellcount) <- sample.ids


# Convert counts to cell fractions
pt.cellfrac <- pt.cellcount
pt.cellfrac <- apply(pt.cellfrac, 1, function(x) {
  return(x / sum(x))
})
pt.cellfrac <- t(pt.cellfrac)
celltypes.final <- colnames(pt.cellfrac)

# Assert that all cell fractions sum to 100% for each patient
stopifnot(sum(rowSums(pt.cellfrac) != 1) == 0)

# Get cluster centroids
umap.coords <- umap_coords_3d # use 3D cell state umap
centroids <- matrix(data=NA, nrow=length(celltypes.final), ncol=ndim.umap)
celltypes.final <- colnames(pt.cellfrac)
for(curidx in seq_len(ncol(pt.cellfrac))) {
  cur.celltype <- celltypes.final[curidx]
  cur_centroid <- colMeans(umap.coords[which(sc.metadata.final$cluster_final == cur.celltype),])
  #cur_centroid <- apply(umap.coords[which(sc.metadata.final$cluster_final == cur.celltype),], 2, median)
  centroids[curidx,] <- cur_centroid
}
rownames(centroids) <- celltypes.final

# Generate GDM for EMD
gdm <- as.matrix(dist(centroids))

pt_cellfrac_mat <- pt.cellfrac
# Compute pairwise EMD
library(transport)
library(pracma)
Y <- rep(0, (nrow(pt_cellfrac_mat)-1)*nrow(pt_cellfrac_mat)/2)
counter <- 1
for(i in seq_len(nrow(pt_cellfrac_mat))) {
  cur_f1_weights <- pt_cellfrac_mat[i,]
  for(j in (i+1):nrow(pt_cellfrac_mat)) {
    if(j > nrow(pt_cellfrac_mat)) break
    cur_f2_weights <- pt_cellfrac_mat[j,]
    flow <- transport(cur_f1_weights, cur_f2_weights, gdm[1:ncol(pt_cellfrac_mat),1:ncol(pt_cellfrac_mat)], method='primaldual')
    curdist <- 0
    for(k in seq_len(nrow(flow))) {
      cur_penalty <- gdm[flow[k,1], flow[k,2]]
      curdist <- curdist + cur_penalty*flow[k,3]
    }
    Y[counter] <- curdist
    counter <- counter + 1
  }
}
emd_distmat_final <- squareform(Y)
rownames(emd_distmat_final) <- rownames(pt_cellfrac_mat)

##################
# Cluster samples
ncluster <- 8 #define 8 clusters
hclust.out <- hclust(as.dist(emd_distmat_final), method='ward.D2')
sample.clusters <- cutree(hclust.out, k=ncluster) 

# Renumber clusters to be ordered in a more interpretable way
cluster_assignments_factored <- factor(sample.clusters)
levels(cluster_assignments_factored) <- c(3,5,7,2,6,8,4,1)
cluster_assignments_factored <- factor(cluster_assignments_factored,levels(cluster_assignments_factored)[c(8,4,1,7,2,5,3,6)])
cluster_assignments_final <- as.numeric(cluster_assignments_factored)
sample.clusters <- cluster_assignments_final

#####################
# Generate color maps
#####################
library(RColorBrewer)
getPalette <- colorRampPalette(brewer.pal(11, 'Spectral'))
cmap <- getPalette(11)
cmap <- cmap[-6]
cmap <- sample(cmap)
sample.clusters.letter <- sapply(cluster_assignments_final, function (x) {
  intToUtf8(64 + x)
})

##############
# Clinical data
###############
myclinicaldata<-  read.delim(file='C:/Users/wilyu/OneDrive/Desktop/Research/Datasets/Krummel_COVID_2020/celltype_labels_subset/2020-10-01/COMET_10X_CLINICAL_SCORES_PAPER.csv', sep=',', header=T, stringsAsFactors=F)
head(myclinicaldata)

seurat2clin_map <- match(rownames(pt_cellfrac_mat),myclinicaldata$SAMPLE.by.SNPs)
myclinicaldata.reordered <- myclinicaldata[seurat2clin_map,]

##############
# Compute PHATE embedding
##########
# Install phateR
#reticulate::py_install("phate", pip=TRUE)
#devtools::install_github("KrishnaswamyLab/phateR")
library(phateR)
phateobj <- phate(emd_distmat_final, ndim=3, knn.dist.method='precomputed_distance', gamma=0.5, knn=3, t=12)

covidstatus.cmap <- c('grey40', 'dodgerblue3', 'firebrick3') #ctrl, neg, pos
disease.severity.cmap <- c("grey40", "orange", "orangered2") # 'CTRL', 'MILD', 'SEVERE'


# Rescale point size to be proportional to y-value to improve depth perception, per reviewer suggestion
pt.szs <- phateobj$embedding[,2] - min(phateobj$embedding[,2])
pt.szs <- pt.szs / max(pt.szs)
pt.szs.final <- pt.szs / 0.25 + 2

# Save PhEMD diffusion map colored by de novo clusters pdf  
library(scatterplot3d)
pdf(file=paste(savedest, 'dm.46pts.clusters.unlabeled.pdf', sep=''), width=7, height=6)
s3d <- scatterplot3d(phateobj$embedding[,1], phateobj$embedding[,2], phateobj$embedding[,3], color=cmap[sample.clusters], pch=20, box=F, tick.marks = F, xlab = 'PHATE1', ylab='PHATE2', zlab='PHATE3', angle=60, scale.y = .8,cex.symbols=pt.szs.final)
legend(s3d$xyz.convert(0.26, -0.08, 0.02), col= cmap[1:ncluster], bg="white", lty=NA, lwd=0, yjust=0, pch=19,legend = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J')[1:ncluster], cex = 0.6)
dev.off()

################
# Save PhEMD PHATE plot colored by covid / disease status without sample labels
library(scales)
# Save plot colored by covid status pdf
legend.labs.covid <- c('Control', 'Neg', 'Pos')
pdf(file=paste(savedest, 'dm.46pts.covidstatus.unlabeled.pdf', sep=''), width=7, height=6)
s3d <- scatterplot3d(phateobj$embedding[,1], phateobj$embedding[,2], phateobj$embedding[,3], color=covidstatus.cmap[as.numeric(factor(myclinicaldata.reordered$covid_status))], pch=20, box=F, tick.marks = F, xlab = 'PHATE1', ylab='PHATE2', zlab='PHATE3', angle=60, scale.y = .8, cex.symbols=pt.szs.final)
legend(s3d$xyz.convert(-0.18, -0.08, 0.07), col= covidstatus.cmap, bg="white", lty=NA, lwd=0, yjust=0, pch=19,legend = legend.labs.covid, cex = 0.8)
dev.off()

# Save PhEMD PHATE plot colored by disease severity pdf
legend.labs.severity <- c('Control', 'Mild', 'Severe')
pdf(file=paste(savedest, 'dm.46pts.diseaseseverity.unlabeled.pdf', sep=''), width=7, height=6)
s3d <- scatterplot3d(phateobj$embedding[,1], phateobj$embedding[,2], phateobj$embedding[,3], color=disease.severity.cmap[as.numeric(factor(myclinicaldata.reordered$Qualitative_score))], pch=20, box=F, tick.marks = F, xlab = 'PHATE1', ylab='PHATE2', zlab='PHATE3', angle=60, scale.y = .8, cex.symbols=pt.szs.final)
legend(s3d$xyz.convert(-0.18, -0.08, 0.07), col= disease.severity.cmap, bg="white", lty=NA, lwd=0, yjust=0, legend = legend.labs.severity, pch = 19, cex = 1)
dev.off()  


########
# Plot celltype fractions (granular binning) as barplot
b_idx <- 1:5
plasma_idx <- 6:7
monomac_idx <- 8:14
neut_idx <- 15:21
nk_idx <- 22:23
t_idx <- 24:34

cluster_assignments <- sample.clusters
cluster_weights <- pt_cellfrac_mat
colnames(cluster_weights) <- c('B_ISG', 'B_mature', 'B_mem', 'B_naive_1', 'B_naive_2', 'Plasma_1', 'Plasma_2', 'Mono_classical_1', 'Mono_classical_ISG', 'Mono_classical_S100A12', 'Mono_DCs', 'Mono_non-classical', 'Mono_non-classical_C1Q', 'Mono_pDCs', 'Neut_G0S2', 'Neut_ISG', 'Neut_LCN2', 'Neut_NEAT1', 'Neut_Ribo', 'Neut_S100A12', 'Neut_SLPI', 'NK_CD56hi', 'NK_CD56lo', 'T_CD4_cyc', 'T_CD8_cyc', 'T_CD8_cytotoxic', 'T_GammaDelta', 'T_ISG', 'T_MAIT', 'T_CD4_mem', 'T_CD8_mem', 'T_CD4_naive', 'T_CD8_naive',   'T_Treg')

cmap.celltype <- c(hue_pal(h = c(0, 45))(length(b_idx)),  hue_pal(h = c(215, 230))(length(plasma_idx)), hue_pal(h = c(45, 90))(length(monomac_idx)), hue_pal(h = c(135, 180))(length(neut_idx)), hue_pal(h = c(250,280))(length(nk_idx)), hue_pal(h = c(290, 360))(length(t_idx)))

proto_inhibs <- matrix(0, max(cluster_assignments), ncol(cluster_weights))
for (i in seq_len(max(cluster_assignments))) {
  if (sum(cluster_assignments == i) == 1) {
    proto_inhibs[i, ] <- cluster_weights[which(cluster_assignments == i), ]
  }
  else {
    proto_inhibs[i, ] <- colMeans(cluster_weights[which(cluster_assignments == i), ])
  }
}
colnames(proto_inhibs) <- colnames(cluster_weights)

# Color legend
cols.celltype <- c(cmap.celltype[1], cmap.celltype[6], cmap.celltype[11], cmap.celltype[15], cmap.celltype[22], cmap.celltype[3])
cmap.celltype <- c(rep(cols.celltype[1], length(b_idx)), rep(cols.celltype[2],length(plasma_idx)), rep(cols.celltype[3], length(monomac_idx)), rep(cols.celltype[4], length(neut_idx)), rep(cols.celltype[5], length(nk_idx)), rep(cols.celltype[6], length(t_idx)))

# Reorder x-axis (cell type ordering) for plotting
proto_inhibs <- proto_inhibs[,c(1:7, 12:13, 8:11, 14:33)]


library(scales)
for (i in seq_len(max(cluster_assignments))) {
  pdf(file=paste(savedest, sprintf('summary_hists/GranularGroup%s_legend.pdf',intToUtf8(64 + i)), sep=''), width=9, height=6)
  par(mar=c(7.1, 5.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  if(max(proto_inhibs[i, ]) > 0.4) {
    ymax <- max(proto_inhibs[i, ]) + 0.1
  } else {
    ymax <- 0.4
  }
  barplot(proto_inhibs[i, ], col = cmap.celltype, main = "", 
          xlab = "", ylab = "", ylim = c(0,  ymax), cex.axis = 1, cex.names = 0.6, cex.lab = 2, 
          names.arg = colnames(proto_inhibs), las=2)
  title(xlab = "Cell subtype", line = 5.5, cex.lab = 2)
  title(ylab = "Frequency (%)", line = 3.5, cex.lab = 2)
  title(main = sprintf("Group %s", intToUtf8(64 + i)), line = 0, cex.main = 3)
  legend(x=32,y=0.25, col= cols.celltype, bg="white", lty=NA, lwd=0, yjust=0, pch=15,legend = c('B cell', 'Plasma cell','Monocyte/Macrophage', 'Neutrophil',  'NK cell', 'T cell'), cex = 0.7)
  dev.off()
}



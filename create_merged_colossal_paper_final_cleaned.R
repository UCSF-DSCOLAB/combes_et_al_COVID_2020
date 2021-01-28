library(cowplot)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(Seurat)
library(harmony)

source('/krummellab/data1/arrao/aarao_scripts/R/BootstrappedAddModuleScore.R')
source('/krummellab/data1/arrao/aarao_scripts/R/signature_dotplot_from_metadata.R')

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

file_name = 'merged_colossal_paper_final'
prefix = file.path('/krummellab/data1/arrao/projects/MVIR1/patient_data_SO/merged/', file_name, file_name)
rdata_prefix = file.path('/krummellab/data1/arrao/projects/MVIR1/patient_data_SO/merged/', file_name, 'rdatas', file_name)
setwd(dirname(prefix))

load(paste0(rdata_prefix, '_merged_temp.RData'))

cleaned_objects_dir = '/krummellab//data1/arrao/projects/MVIR1/patient_data_SO/merged/merged_colossal_paper_final/annotated_res_0.6/monocles/inputs'


rds = Sys.glob(paste0(cleaned_objects_dir, '/*rds'))
cells_to_keep = colnames(readRDS(rds))
sobj = NULL
cls = c(ls(), 'cls')

for (rdata in c(Sys.glob(paste0(cleaned_objects_dir, '/*Robj')), 
                Sys.glob(paste0(cleaned_objects_dir, '/*RData')))) {
  load(rdata)
  sobj = setdiff(ls(), cls)
  stopifnot(length(sobj) == 1)
  cells_to_keep = c(cells_to_keep, colnames(get(sobj)))
  rm(list=sobj)
}


merged_data <- subset(merged_data, cells=cells_to_keep)
prefix = paste0(prefix, '_cleaned')

pdf(paste0(prefix, '_harmony_convergence.pdf'))
merged_data <- RunHarmony(merged_data, "LIBRARY", plot_convergence = TRUE, max.iter.harmony=30, max.iter.cluster=30)
dev.off()

source('/krummellab/data1/arrao/projects/MVIR1/auxiliary_files/five_gene_sigs_final.R')

for (sig in names(five_gene_sigs_final)) {
    nb = 30
    while (TRUE) {
        try(merged_data <- BootstrappedAddModuleScore(merged_data,
                                                      assay='RNA',
                                                      features=list(five_gene_sigs_final[[sig]]),
                                                      name=sig,
                                                      nbin = nb,
                                                      num_iters=100))
        if (paste0(sig, '1') %in% colnames(merged_data@meta.data)) {
            merged_data@meta.data[sig] <- merged_data@meta.data[paste0(sig, '1')]
            merged_data@meta.data[paste0(sig, '1')] <- NULL
            merged_data@meta.data[paste0(sig, '_75sat')] <- saturate(merged_data@meta.data[[sig]], sat=0.75, binary=T)
            merged_data@meta.data[paste0(sig, '_75nbsat')] <- saturate(merged_data@meta.data[[sig]], sat=0.75, binary=F)
            break
        } else if (nb == 5) {
          stop()
        } else {
          nb <- nb - 1
        }
    }
}

source('/krummellab/data1/arrao/projects/MVIR1/auxiliary_files/covid_status_all_consent.R')
if (!all(names(covid_status)[!names(covid_status) %in% merged_data@meta.data$SAMPLE.by.SNPs])){
  stop('Not all sample in the covid status file')
}
names(covid_status)[!names(covid_status) %in% merged_data@meta.data$SAMPLE.by.SNPs]

merged_data@meta.data$covid_status = covid_status[merged_data@meta.data$SAMPLE.by.SNPs]

save(merged_data, file=paste0(rdata_prefix, '_merged_temp.RData'))

ribo_genes <- read.table("/krummellab/data1/ipi/data/refs/10x/genesets/GRCh38/ribo_genes.tsv",
                         sep = "\t", 
                         header=TRUE,
                         stringsAsFactors = FALSE)
ribo_genes <- ribo_genes[ribo_genes[["HUGO"]] %in% rownames(merged_data), 'HUGO']
mito_genes <- rownames(merged_data)[grep('^MT-', rownames(merged_data))]

pcs_to_use = c()
for (pc in 1:30) {
  top_10 <- names(sort(merged_data@reductions$harmony@feature.loadings[,pc]))[1:10]
  bottom_10 <- names(sort(merged_data@reductions$harmony@feature.loadings[,pc], decreasing=T))[1:10]
  if ((sum(top_10 %in% mito_genes) + sum(top_10 %in% ribo_genes)) >= 6){
    next
  } else if ((sum(bottom_10 %in% mito_genes) + sum(bottom_10 %in% ribo_genes))>= 6){
    next
  } else {
    pcs_to_use <- c(pcs_to_use, pc)
  }
}

# Use 20, 0.2, 1.5
merged_data <- RunUMAP(merged_data,
                       dims = pcs_to_use,
                       reduction = 'harmony',
                       n.neighbors = 20,
                       min.dist = 0.2,
                       spread = 1.5,
                       a = NULL,
                       b = NULL,
                       verbose = FALSE)

# Calculate the neighborhood graph
merged_data <- FindNeighbors(merged_data,
                             dims = pcs_to_use,  # Num PCs to use
                             reduction = 'harmony',
                             k.param = 20,  # k for the knn algorithm
                             verbose = FALSE
                             )

md <- read.table('rdatas/merged_colossal_paper_final_merged_annotated_metadata.tsv', 
                 sep='\t', 
                 header=T, row.names=1)


png(filename=paste0(prefix, '_sample_umap.png'), width = 10, height = 10, units = "in", res = 300)
print(DimPlot(merged_data, group.by='LIBRARY'))
dev.off()


png(filename=paste0(prefix, '_split_sample_umap.png'), width = 21, height = 18, units = "in", res = 300)
print(DimPlot(merged_data, split.by='LIBRARY', ncol=7) + theme(legend.position="none", axis.title=element_blank(), axis.text=element_blank()))
dev.off()


png(filename=paste0(prefix, '_covid_status_umap.png'), width = 10, height = 10, units = "in", res = 300)
print(DimPlot(merged_data, group.by='covid_status'))
dev.off()


png(filename=paste0(prefix, '_split_covid_status_umap.png'), width = 15, height = 5, units = "in", res = 300)
print(DimPlot(merged_data, split.by='covid_status', ncol=7) + theme(legend.position="none", axis.title=element_blank(), axis.text=element_blank()))
dev.off()


signatures = paste(names(five_gene_sigs_final), '_75sat', sep='')
signatures2 = paste(names(five_gene_sigs_final), '_75nbsat', sep='')
genes_to_consider = c()
for (gs in names(five_gene_sigs_final)){
  for (gene in five_gene_sigs_final[[gs]]){
    if (!gene %in% genes_to_consider){
      genes_to_consider <- c(genes_to_consider, gene)
    }
  }
}

save(merged_data, file=paste0(rdata_prefix, '_merged_unprocessed.RData'))

for (res in c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2)){
  prefix2 <- paste0(dirname(prefix), '/res_', res, '/', basename(prefix))
  merged_data <- FindClusters(merged_data, verbose = TRUE,
                            algorithm = 1,
                            resolution = res)
  merged_data@meta.data[[paste0('louvain_res', res)]] <- merged_data@meta.data$seurat_clusters

  png(filename=paste0(prefix2, '_clusters_louvain_res', res,'_umap.png'), width = 10, height = 10, units = "in", res = 300)
  print(DimPlot(merged_data, group.by=paste0('louvain_res', res)))
  dev.off()

  png(filename=paste0(prefix2, '_clusters_louvain_res', res,'_labeled_umap.png'), width = 10, height = 10, units = "in", res = 300)
  print(DimPlot(merged_data, group.by=paste0('louvain_res', res), label=T))
  dev.off()

  png(paste(prefix2, 'res', res, 'SigDotPlot.png', sep='_'), height=100 * length(signatures), 
      width=100 * length(unique(as.vector(merged_data@meta.data[[paste0('louvain_res', res)]]))), units='px', res=150)
  print(signature_dotplot_from_metadata(merged_data@meta.data, signatures=signatures, group.by=paste0('louvain_res', res)))
  dev.off()
  
  png(paste(prefix2, 'res', res, 'NBSigDotPlot.png', sep='_'), height=100 * length(signatures2), 
      width=100 * length(unique(as.vector(merged_data@meta.data[[paste0('louvain_res', res)]]))), units='px', res=150)
  print(signature_dotplot_from_metadata(merged_data@meta.data, signatures=signatures2, group.by=paste0('louvain_res', res)))
  dev.off()

  png(paste(prefix2, 'res', res, 'FiveGeneDotPlot.png', sep='_'), height=100 * length(genes_to_consider), 
      width=100 * length(unique(as.vector(merged_data@meta.data[[paste0('louvain_res', res)]]))), units='px', res=150)
  print(DotPlot(merged_data, 
                assay='RNA', 
                group.by=paste0('louvain_res', res), 
                features = genes_to_consider,
                cols='RdYlBu') + coord_flip())
  dev.off()

  Idents(merged_data) <- merged_data@meta.data[[paste0('louvain_res', res)]]
  plots1 <- list()
  plots2 <- list()
  plots3 <- list()
  plots1[['clusters']] <- print(DimPlot(object=merged_data)+ NoLegend())
  plots2[['clusters']] <- print(DimPlot(object=merged_data)+ NoLegend())
  plots3[['clusters']] <- print(DimPlot(object=merged_data)+ NoLegend())

  for (sig in names(five_gene_sigs_final)) {
      plots1[[sig]] <- print(FeaturePlot(merged_data, features=sig))
      plots2[[sig]] <- print(FeaturePlot(merged_data, features=paste0(sig, '_75sat')))
      plots3[[sig]] <- print(FeaturePlot(merged_data, features=paste0(sig, '_75nbsat')))
  }

  png(filename=paste(prefix2, 'res', res, 'fivegene_umap.png', sep='_'), width = 15, height = 15, units = "in", res = 300)
  print(plot_grid(plotlist = plots1, ncol=4))
  dev.off()

  png(filename=paste(prefix2, 'res', res, 'fivegene_sat_umap.png', sep='_'), width = 15, height = 15, units = "in", res = 300)
  print(plot_grid(plotlist = plots2, ncol=4))
  dev.off()

  png(filename=paste(prefix2, 'res', res, 'fivegene_nbsat_umap.png', sep='_'), width = 15, height = 15, units = "in", res = 300)
  print(plot_grid(plotlist = plots3, ncol=4))
  dev.off()


  png(filename=paste(prefix2, 'res', res, 'fivegene_violins.png', sep='_'), width = 15, height = 15, units = "in", res = 300)
  print(VlnPlot(merged_data, features=names(five_gene_sigs_final), pt.size=0))
  dev.off()

  x <- reshape2::melt(merged_data@meta.data[, c(paste0('louvain_res', res), paste(names(five_gene_sigs_final), '_75sat', sep=''))])
  y <- x %>% 
          mutate(cluster_temp= eval(parse(text=paste0('louvain_res', res)))) %>% 
          group_by(cluster_temp, variable) %>% 
          mutate(mean=mean(value), cluster=as.numeric(as.vector(cluster_temp))) %>% 
          ungroup() %>%
          distinct(cluster, variable, mean)


  png(paste(prefix2, 'res', res, 'mean_score_per_cluster.png', sep='_'), width = 15, height = 15, units = "in", res = 300)
  print(ggplot(y) +
            geom_point(aes(x=cluster, y=mean, col=variable)) +
            geom_line(aes(x=cluster, y=mean, col=variable)) +
            scale_x_continuous(breaks=sort(unique(y$cluster)), labels=as.character(sort(unique(y$cluster)))) +
            theme_bw())
  dev.off()
}

save(merged_data, file=paste0(rdata_prefix, '_merged_processed.RData'))

metadata <- merged_data@meta.data
for (redn in c('pca', 'umap')){
  cell_embeddings <- as.data.frame(merged_data@reductions[[redn]]@cell.embeddings[, c(1: min(5, dim(merged_data@reductions[[redn]]@cell.embeddings)[2]))])
  metadata <- merge(metadata, 
                    cell_embeddings,
                    by=0)
  rownames(metadata) <- metadata$Row.names
  metadata$Row.names <- NULL
}

write.table(metadata,
            file=paste(rdata_prefix, 'metadata.tsv', sep='_'),
            sep='\t',
            row.names=T,
            col.names=T,
            quote=F)

stop('manual from here')
plots <- list()
for (l in levels(metadata$louvain_res0.6)){
  plots[[l]] <- ggplot() + 
                  geom_point(data=metadata, aes(x=UMAP_1, y=UMAP_2), col='grey', size=0.15) + 
                  geom_point(data=metadata[metadata$louvain_res0.6 == l, ], aes(x=UMAP_1, y=UMAP_2), col='red', size=0.15) + 
                  annotate("text", label = l, x = -Inf, y = Inf, size = 10, colour = "red", hjust = 0, vjust = 1)+
                  theme_bw() +
                  theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
}

png(paste0(prefix, '_res_0.6_per_cluster_dimplots.png'), width=15, height=(3 * ceiling(length(plots)/5)), units='in', res=150)
plot_grid(plotlist=plots, ncol=5)
dev.off()

plots <- list()
metadata$LIBRARY <- as.factor(metadata$LIBRARY)
libraries <- setNames(1:length(levels(metadata$LIBRARY)), levels(metadata$LIBRARY))
for (l in names(libraries)){
  plots[[l]] <- ggplot() + 
                  geom_point(data=metadata, aes(x=UMAP_1, y=UMAP_2), col='grey', size=0.15) + 
                  geom_point(data=metadata[metadata$LIBRARY == l, ], aes(x=UMAP_1, y=UMAP_2), col='red', size=0.15) + 
                  annotate("text", label = l, x = -Inf, y = Inf, size = 5, colour = "red", hjust = 0, vjust = 1)+
                  theme_bw() +
                  theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
}

png(paste0(prefix, '_per_LIBRARY_dimplots.png'), width=15, height=(3 * ceiling(length(plots)/5)), units='in', res=150)
plot_grid(plotlist=plots, ncol=5)
dev.off()

Idents(merged_data) <- merged_data@meta.data$louvain_res0.6
markers <- FindAllMarkers(merged_data, 
                          test.use='poisson', 
                          latent.vars='LIBRARY', 
                          assay='RNA', 
                          logfc.threshold=0.4, 
                          min.pct=0.2, 
                          only.pos=TRUE)
write.table(markers, file=paste0(prefix, '_res_0.6_markers.tsv'), sep='\t', col.names=T, row.names=F, quote=F)


top_markers = (markers %>% group_by(cluster) %>% top_n(5, wt=avg_logFC))$gene
genes_to_consider = c()
for (gene in top_markers){
  if (!gene %in% genes_to_consider){
    genes_to_consider <- c(genes_to_consider, gene)
  }
}

png(paste(prefix, 'res0.6_TopMarkerDotPlot.png', sep='_'), 
    height=50 * length(genes_to_consider), 
    width=100 * length(unique(as.vector(merged_data@meta.data$louvain_res0.6))), 
    units='px', 
    res=150)
print(DotPlot(merged_data, 
              features=genes_to_consider, 
              group.by='louvain_res0.6',
              assay='RNA',
              cols='RdYlBu') + coord_flip())
dev.off()
## Author: Arjun Arkal Rao

library(cowplot)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(Seurat)
library(harmony)

# These files are in auxiliary_files/
source('/krummellab/data1/arrao/aarao_scripts/R/DoMultiBarHeatmap.R')
source('/krummellab/data1/arrao/aarao_scripts/R/plot_cluster_diffs_by_sample.R')

# Input file
source('/krummellab/data1/arrao/projects/MVIR1/auxiliary_files/merged_colossal_full_paper_final_res0.6_annotations.R')


if (!all(c('file_name', 'res', 'annotations') %in% ls())){
  stop('Cannot run without file_name, res, and annotations.')
}

prefix = paste0('/krummellab/data1/arrao/projects/MVIR1/patient_data_SO/merged/', file_name, '/annotated_res_', res, '/', file_name)
rdata_prefix = paste0('/krummellab/data1/arrao/projects/MVIR1/patient_data_SO/merged/', file_name, '/annotated_res_', res, '/rdatas/', file_name)

dir.create(dirname(rdata_prefix), recursive=TRUE, showWarnings=TRUE)
if (!dir.exists(dirname(rdata_prefix))){
  stop(paste0('Could not create output dir recursively: ', dirname(rdata_prefix)))
}
setwd(dirname(prefix))

load(file=paste0(dirname(dirname(prefix)), '/rdatas/', file_name, '_merged_processed.RData'))

source('/krummellab/data1/arrao/projects/MVIR1/auxiliary_files/five_gene_sigs_final.R')

Idents(merged_data) <- merged_data@meta.data[[paste0('louvain_res', res)]]
names(annotations) <- c(0: (length(annotations)-1))
annotations <- annotations[levels(Idents(merged_data))]

merged_data <- RenameIdents(merged_data, annotations)
merged_data@meta.data$fine_annotations <- as.vector(Idents(merged_data))

merged_data@meta.data$coarse_annotations <- gsub("_+[0-9]+$", "" , merged_data@meta.data$fine_annotations)

for (annot_type in c('fine_annotations', 'coarse_annotations')){
  Idents(merged_data) <- merged_data@meta.data[[annot_type]]
  prefix2 <- paste0(dirname(prefix), '/', annot_type, '/', basename(prefix))
  dir.create(dirname(prefix2), recursive=FALSE, showWarnings=TRUE)
  png(filename=paste0(prefix2, '_res', res, '_', annot_type, '_umap.png'), width = 10, height = 10, units = "in", res = 300)
  DimPlot(merged_data, group.by=annot_type, label=T)
  dev.off()

  for (id_type in c('SAMPLE.by.SNPs', 'LIBRARY')){
    outs <- plot_cluster_diffs_by_sample(merged_data, 
                                         ann.col=id_type, 
                                         ident.col=annot_type)
    ggsave(filename = paste0(prefix2, '_res', res, '_', annot_type, '_', id_type, '_fractions_per_CellType_hist.png'),
           plot=outs[['plots']][['ann_fractions_per_ident_hist']],
           height=5, 
           width=10, 
           units='in',
           dpi=150)
    ggsave(filename = paste0(prefix2, '_res', res, '_', annot_type, '_CellType_fractions_per_', id_type, '_hist.png'),
           plot=outs[['plots']][['ident_fractions_per_ann_hist']],
           height=5, 
           width=10, 
           units='in',
           dpi=150)
    ggsave(filename = paste0(prefix2, '_res', res, '_', annot_type, '_CellType_fractions_per_', id_type, '_line.png'),
           plot=outs[['plots']][['ident_fractions_per_ann_line']],
           height=5, 
           width=10, 
           units='in',
           dpi=150)
    write.table(format(as.data.frame(outs[['fraction_table']]),  digits=2),
                file=paste0(prefix2, '_res', res, '_', annot_type, '_', id_type, '_counts.tsv'),
                sep='\t',
                col.names=T,
                row.names=F,
                quote=F)
  }
}

save(merged_data, file=paste0(dirname(dirname(prefix)), '/rdatas/', file_name, '_merged_annotated.RData'))
metadata <- merged_data@meta.data
for (redn in c('harmony', 'pca', 'umap')){
  cell_embeddings <- as.data.frame(merged_data@reductions[[redn]]@cell.embeddings[, c(1, 2)])
  metadata <- merge(metadata, 
                    cell_embeddings, 
                    by=0)
  rownames(metadata) <- metadata$Row.names
  metadata$Row.names <- NULL
}
write.table(metadata,
            file=paste0(dirname(dirname(prefix)), '/rdatas/', file_name, '_merged_annotated_metadata.tsv'),
            sep='\t',
            row.names=T,
            col.names=T,
            quote=F)

# Do these separately since the first half is considerably faster than the second
for (annot_type in c('fine_annotations', 'coarse_annotations')){
  prefix2 <- paste0(dirname(prefix), '/', annot_type, '/', basename(prefix))
  Idents(merged_data) <- merged_data@meta.data[[annot_type]]
  markers <- FindAllMarkers(merged_data, 
                            only.pos = TRUE, 
                            min.pct = 0.25, 
                            logfc.threshold = 0.40, 
                            test.use="poisson",
                            latent.vars = 'LIBRARY',
                            assay='RNA')
  write.table(markers, 
              file=paste0(prefix2, '_res', res, '_', annot_type, '_markers.tsv'),
              sep='\t',
              col.names=T,
              row.names=T,
              quote=F)
  
  gene_lists = list(
    top5=(markers %>% 
              group_by(cluster) %>% 
                top_n(5, wt=avg_logFC))$gene,
    five_gene=unname(unlist(five_gene_sigs_final))
  )

  for (gl in names(gene_lists)){
    # Can't be png see https://github.com/satijalab/seurat/issues/1277
    plot <- DoHeatmap(merged_data, 
                      group.by=annot_type, 
                      features=gene_lists[[gl]])
    ggsave(filename = paste0(prefix2, '_res', res, '_', annot_type, '_', gl, '_heatmap.pdf'), 
                             plot = plot, 
                             width=20, 
                             height=10, 
                             units='in', 
                             dpi=300)

    plot <- DoMultiBarHeatmap(merged_data,
                              group.by=annot_type, 
                              additional.group.by=c('LIBRARY', 'SAMPLE.by.SNPs'),
                              additional.group.sort.by=c('LIBRARY', 'SAMPLE.by.SNPs'),
                              features=gene_lists[[gl]])
    ggsave(filename = paste0(prefix2, '_res', res, '_', annot_type, '_', gl, '_multibar_heatmap.pdf'), 
                             plot = plot, 
                             width=20, 
                             height=10, 
                             units='in', 
                             dpi=300)
  }
}

ribo_genes <- read.table("/krummellab/data1/ipi/data/refs/10x/genesets/GRCh38/ribo_genes.tsv",
                         sep = "\t", 
                         header=TRUE,
                         stringsAsFactors = FALSE)
ribo_genes <- ribo_genes[ribo_genes[["HUGO"]] %in% rownames(merged_data), 'HUGO']
mito_genes <- rownames(merged_data)[grep('^MT-', rownames(merged_data))]

Idents(merged_data) <- merged_data@meta.data$coarse_annotations
for (cluster in levels(Idents(merged_data))) {
  if (cluster %in% c('rbc', 'junk')){
    next
  }

  prefix3 <- paste0(dirname(prefix), '/', cluster, '/', basename(prefix))
  dir.create(dirname(prefix3), recursive=FALSE, showWarnings=TRUE)

  sobj <- subset(merged_data, idents=cluster)
  
  # This will happen on a repeat after QC
  if (file.exists(paste0(dirname(prefix3), '/cleaned_', cluster, '_cellnames.list'))){
    cleaned_cells <- readLines(paste0(dirname(prefix3), '/cleaned_', cluster, '_cellnames.list'))
    sobj <- subset(merged_data, cells=cleaned_cells)
  }

  sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
  sobj <- ScaleData(object = sobj, 
                    vars.to.regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"))
  
  sobj <- RunPCA(object = sobj)

  pdf(paste0(prefix3, '_', cluster, '_harmony_convergence.pdf'))
  sobj <- RunHarmony(sobj, "LIBRARY", plot_convergence = TRUE, max.iter.harmony=30, max.iter.cluster=30)
  dev.off()

  pcs_to_use = c()
  for (pc in 1:30) {
    top_10 <- names(sort(sobj@reductions$harmony@feature.loadings[,pc]))[1:10]
    bottom_10 <- names(sort(sobj@reductions$harmony@feature.loadings[,pc], decreasing=T))[1:10]
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
                  n.neighbors = 20,
                  min.dist = 0.1, 
                  spread = 2,
                  verbose = FALSE)

  sobj <- FindNeighbors(sobj,
                        dims = pcs_to_use,  # Num PCs to use
                        reduction='harmony',
                        k.param = 20,  # k for the knn algorithm
                        verbose = FALSE)

  png(filename=paste0(prefix3, '_', cluster, '_sample_umap.png'), 
      width = 10, height = 10, units = "in", res = 300)
  print(DimPlot(sobj, group.by='SAMPLE.by.SNPs') + NoLegend())
  dev.off()

  png(filename=paste0(prefix3, '_', cluster, '_covid_status_umap.png'), 
      width = 10, height = 10, units = "in", res = 300)
  print(DimPlot(sobj, group.by='covid_status'))
  dev.off()

  for (res2 in c(0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0)){
    sobj <- FindClusters(sobj, 
                         verbose = TRUE, 
                         algorithm = 1,
                         resolution = res2)
    sobj@meta.data[[paste0('subset_louvain_res', res2)]] <- sobj@meta.data$seurat_clusters
    png(filename=paste0(prefix3, '_res', res, '_', cluster, '_subset_louvain_res', res2,'_umap.png'), 
        width = 10, height = 10, units = "in", res = 300)
    print(DimPlot(sobj, group.by=paste0('subset_louvain_res', res2), label=T))
    dev.off()
  }

  assign(cluster, sobj)
  save(list=cluster, file=paste0(prefix3, '_res', res, '_', cluster, '.RData'))
  rm(list=cluster)
}

Idents(merged_data) <- merged_data@meta.data$coarse_annotations
for (cluster in levels(Idents(merged_data))) {
  if (cluster %in% c('rbc', 'junk')){
    next
  }
  # Just use the lowest resolution
  Idents(sobj) <- sobj@meta.data[['subset_louvain_res0.4']]
  markers <- FindAllMarkers(sobj, 
                            test.use='poisson', 
                            latent.vars='LIBRARY', 
                            assay='RNA', 
                            logfc.threshold=0.4, 
                            min.pct=0.2, 
                            only.pos=TRUE)
  write.table(markers, 
              file=paste0(prefix3, '_res', res, '_', cluster, '_subset_louvain_res0.4_markers.tsv'), 
              sep='\t', 
              col.names=T, 
              row.names=F, 
              quote=F)
  
  top_markers = (markers %>% group_by(cluster) %>% top_n(5, wt=avg_logFC))$gene
  genes_to_consider = c()
  for (gene in top_markers){
    if (!gene %in% genes_to_consider){
      genes_to_consider <- c(genes_to_consider, gene)
    }
  }

  png(paste0(prefix3, '_res', res, '_', cluster, '_subset_louvain_res0.4_TopMarkerDotPlot.png'), 
      height=50 * length(genes_to_consider), 
      width=100 * length(unique(as.vector(sobj@meta.data$subset_louvain_res0.4))), 
      units='px', 
      res=150)
  print(DotPlot(sobj, 
                features=genes_to_consider, 
                group.by='subset_louvain_res0.4',
                assay='RNA',
                cols='RdYlBu') + coord_flip())
  dev.off()

  Idents(sobj) <- sobj@meta.data[['covid_status']]
  markers <- FindAllMarkers(sobj, 
                            test.use='poisson', 
                            latent.vars='LIBRARY', 
                            assay='RNA', 
                            logfc.threshold=0.4, 
                            min.pct=0.2, 
                            only.pos=TRUE)
  write.table(markers, 
              file=paste0(prefix3, '_res', res, '_', cluster, '_subset_louvain_res0.4_covid_status_markers.tsv'), 
              sep='\t', 
              col.names=T, 
              row.names=F, 
              quote=F)
  
  top_markers = (markers %>% group_by(cluster) %>% top_n(10, wt=avg_logFC))$gene
  genes_to_consider = c()
  for (gene in top_markers){
    if (!gene %in% genes_to_consider){
      genes_to_consider <- c(genes_to_consider, gene)
    }
  }

  png(paste0(prefix3, '_res', res, '_', cluster, '_subset_louvain_res0.4_CovidStatusTopMarkerDotPlot.png'), 
      height=100 * length(genes_to_consider), 
      width=300 * length(unique(as.vector(sobj@meta.data$covid_status))), 
      units='px', 
      res=150)
  print(DotPlot(sobj, 
                features=genes_to_consider, 
                group.by='covid_status',
                assay='RNA',
                cols='RdYlBu') + coord_flip())
  dev.off()
}

for (cluster in levels(Idents(merged_data))) { 
  if (prefix == 'neuts') {
    md <- read.table('merged_colossal_full_seurat_Neuts_0.4_metadata.tsv',
                     sep='\t', 
                     header=T,
                     row.names=1,
                     stringsAsFactors=F)
  } else if (prefix == 'monodc') {
    md <- read.table('merged_colossal_full_res1.2_monodc_selected_noNonConsent_metadata.tsv',
                     sep='\t', 
                     header=T,
                     row.names=1,
                     stringsAsFactors=F)
  } else if (prefix == 'tnk') {
    md <- read.table('TNK_metadata_updated.tsv',
                     sep='\t', 
                     header=T,
                     row.names=1,
                     stringsAsFactors=F)
    md$Coarse_Subs_annotations <- md$broad_clusters
  } else {
    stop('No metadata received')
  }

  md <- md[rownames(md) %in% colnames(sobj), 'Coarse_Subs_annotations', drop=F]
  missing_cells <- colnames(sobj)[!colnames(sobj) %in% rownames(md)]
  md2 <- data.frame(Coarse_Subs_annotations=rep('MISSING', length(missing_cells)), row.names=missing_cells)
  md_final <- rbind(md, md2)
  sobj <- AddMetaData(sobj, md_final)

  md <- cbind(sobj@meta.data, sobj@reductions$umap@cell.embeddings)
  png(filename=paste0(prefix, '_', cluster, '_seurat_clusters_umap.png'), 
      width = 10, height = 10, units = "in", res = 300)
  print(DimPlot(sobj, group.by='Coarse_Subs_annotations'))
  dev.off()

  plots <- list()
  for (l in unique(sobj$Coarse_Subs_annotations)) {
  plots[[l]] <- ggplot() + 
                  geom_point(data=md, aes(x=UMAP_1, y=UMAP_2), col='grey', size=0.15) + 
                  geom_point(data=md[md$Coarse_Subs_annotations == l, ], aes(x=UMAP_1, y=UMAP_2), col='red', size=0.15) + 
                  annotate("text", label = l, x = -Inf, y = Inf, size = 5, colour = "red", hjust = 0, vjust = 1)+
                  theme_bw() +
                  theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
  }
  png(paste0(prefix, '_', cluster, '_per_seurat_clusters_dimplots.png'), width=15, height=(3 * ceiling(length(plots)/5)), units='in', res=150)
  plot_grid(plotlist=plots, ncol=5)
  dev.off()

} 

stop("auxiliary stuff")
# subset Neuts
md <- read.table('merged_colossal_full_harmony2_Neuts_1_metadata.tsv',
                 sep='\t',
                 header=T,
                 row.names=1,
                 stringsAsFactors=F)
stopifnot(all(rownames(md) == colnames(sobj)))
sobj <- subset(sobj, cells=colnames(sobj)[md$Subsetting_annotations == 'KEEP' & sobj@reductions$umap@cell.embeddings[,1] < 10])

stop('neuts end')
## subset Monos
md <- read.table('merged_colossal_full_harmony_res0.6_monodc_metadata.txt',
                 sep='\t',
                 header=T,
                 row.names=1,
                 stringsAsFactors=F)
stopifnot(all(rownames(md) == colnames(sobj)))
sobj <- subset(sobj, cells=colnames(sobj)[md$Subsetting_annotations == 'KEEP' & sobj@reductions$umap@cell.embeddings[,1] < 10])



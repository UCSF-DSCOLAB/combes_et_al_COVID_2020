## Author: Arjun Arkal Rao

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


plot_cluster_diffs_by_sample <- function(sobj, sample.col='SAMPLE', cluster.col='seurat_clusters', palette='Dark2', prefix='S') {
  samples <- levels(as.factor(sobj@meta.data[[sample.col]]))
  suppressWarnings({
    if (any(!is.na(as.numeric(samples)))){
      cat(paste0('Not all samples are strictly non-numeric. Prefixing with `', prefix, '`.'))
      sobj@meta.data[[sample.col]] <- gsub('^', prefix, as.vector(sobj@meta.data[[sample.col]]))
      samples <- gsub('^', prefix, as.vector(samples))
    }
  })

  colors <- palette
  #try(colors <- brewer.pal(name = palette, n=length(samples)), silent = TRUE)

  #if (is.null(colors)){
  #  stop("Palette `", palette, "` supports only (", brewer.pal.info[palette, 'maxcolors'], ") colors but you have (", length(samples), ") annotations. please choose a different palette.")
  #}

  #names(colors) <- samples
  outs <- list()
  outs[['plots']] <- list()
  outs[['plots']][['dimplot']] <- DimPlot(sobj, split.by = sample.col) +
                                      guides(colour = guide_legend(override.aes = list(colour='white'))) +
                                      theme(legend.title = element_blank(),
                                            legend.key=element_blank(),
                                            legend.text = element_blank())

  outs[['sample_fractions']] <- sobj@meta.data %>%
      count((!!sym(cluster.col)), (!!sym(sample.col))) %>%
      group_by((!!sym(cluster.col))) %>%
      mutate(cluster_size=sum(n)) %>%
      mutate(fraction_of_cluster = n/cluster_size) %>%
      group_by((!!sym(sample.col))) %>%
      mutate(sample_size=sum(n)) %>%
      mutate(fraction_of_sample = n/sample_size) %>%
      arrange(-cluster_size, .by_group = TRUE) %>%
      mutate(idx=row_number())
  outs[['plots']][['sample_fractions_per_cluster_hist']] <- ggplot(outs[['sample_fractions']],
                                                                   aes_string(y="fraction_of_cluster",
                                                                              x=cluster.col,
                                                                              fill=sample.col)) +
                                                               geom_col(position = "dodge") +
                                                               theme(legend.key.size = unit(0.3, "cm"))

  # convenience
  temp <- outs[['sample_fractions']]
  outs[['plots']][['cluster_fractions_per_sample_hist']] <- ggplot(temp, aes_string(y="fraction_of_sample", x=sample.col, fill=cluster.col)) +
                                                               geom_col(position = "dodge") +
                                                               theme(legend.key.size = unit(0.3, "cm"))

  outs[['plots']][['cluster_fractions_per_sample_line']] <- ggplot(temp, aes_string(x="idx", y="fraction_of_sample", col=sample.col)) +
                                                               geom_line() +
                                                               geom_point(shape=1, size=3) +
                                                               #scale_color_brewer(palette=palette)
                                                               xlab(cluster.col) +
                                                               ylab('percent') +
                                                               scale_x_continuous(breaks=temp[temp[[sample.col]]==temp[[sample.col]][1], 'idx', drop=T],
                                                                                  labels=as.vector(temp[temp[[sample.col]]==temp[[sample.col]][1], cluster.col, drop=T]))
  outs
}


file_name = 'merged_colossal_paper_final_platelets_first'
prefix = file.path('/krummellab/data1/arrao/projects/MVIR1/patient_data_SO/merged/', file_name, file_name)
rdata_prefix = file.path('/krummellab/data1/arrao/projects/MVIR1/patient_data_SO/merged/', file_name, 'rdatas', file_name)
rbc_cluster <- read.csv('/krummellab/data1/arrao/projects/MVIR1/auxiliary_files/merged_colossal_rbc_clusters_only_consent.csv', row.names=1)
colnames(rbc_cluster) <- gsub('[.]', '-', colnames(rbc_cluster))

source('/krummellab/data1/arrao/projects/MVIR1/auxiliary_files/freemuxlet_assignments.R')
dir.create(dirname(prefix), showWarnings=F)
dir.create(dirname(rdata_prefix), showWarnings=F)
setwd(dirname(prefix))

sobjs <- list()
for (sample in colnames(rbc_cluster)){
  load(paste0('/krummellab/data1/arrao/projects/MVIR1/patient_data_SO/', sample, '/', sample, '_scTransformed_processed.RData'))
  sobjs[[sample]] <- get(sample)
  rm(list=sample)

  sobjs[[sample]] <- NormalizeData(object = sobjs[[sample]], assay="RNA")
  sobjs[[sample]]@active.assay = 'RNA'
  sobjs[[sample]][["SCT"]] <- NULL

  # Remove non-consented sample from the pooled object.
  if (sample == 'MVIR1-POOL-SCG4'){
      sobjs[[sample]] <- subset(sobjs[[sample]], cells=colnames(sobjs[[sample]])[sobjs[[sample]]$SAMPLE.by.SNPs == 1])
  }

  md <- read.table(paste0('/krummellab/data1/arrao/projects/MVIR1/auxiliary_files/doublets/',
                          sample,
                          '/',
                          sample,
                          '_doubletfinder_metadata.tsv'),
                   sep='\t',
                   header=TRUE,
                   row.names=1,
                   stringsAsFactors=FALSE)

  md <- md[, 2, drop=F]
  colnames(md) <- ('doublet_status')

  sobjs[[sample]] <- AddMetaData(sobjs[[sample]], md)

  rbc_cluster_numbers <- rownames(rbc_cluster[rbc_cluster[[sample]]=="R" & !is.na(rbc_cluster[[sample]]),sample, drop=F])
  if (length(rbc_cluster_numbers) > 0){
    non_rbc_indices = ! as.character(as.vector(sobjs[[sample]]$seurat_clusters)) %in% rbc_cluster_numbers
  } else {
    non_rbc_indices = rep(TRUE, dim(sobjs[[sample]])[2])
  }

  sobjs[[sample]] <- subset(sobjs[[sample]], cells=colnames(sobjs[[sample]])[non_rbc_indices])

  sobjs[[sample]]@meta.data$LIBRARY <- sample
  sobjs[[sample]]@meta.data$temp <- gsub('MVIR', 'M',
                                         gsub('XHLT', 'X',
                                              gsub('POOL', 'P',
                                                   gsub('HS', 'H',
                                                        gsub('BLD', 'B',
                                                             gsub('D0BLD', 'D0B',
                                                                  gsub('SCG', 'S', 
                                                                       gsub('-', '.', sample))))))))
  if ('SAMPLE.by.SNPs' %in% colnames(sobjs[[sample]]@meta.data)) {
    sobjs[[sample]]@meta.data$SAMPLE.by.SNPs <- freemuxlet_assignments[paste(sobjs[[sample]]@meta.data$temp, sobjs[[sample]]@meta.data$SAMPLE.by.SNPs, sep='.C')]
    stopifnot(any(is.na(sobjs[[sample]]@meta.data$SAMPLE.by.SNPs)) == FALSE)
  } else {
    sobjs[[sample]]@meta.data$SAMPLE.by.SNPs <- sobjs[[sample]]@meta.data$temp
  }  
  sobjs[[sample]]@meta.data$temp <- NULL
}

source('/krummellab/data1/arrao/projects/MVIR1/auxiliary_files/five_gene_sigs_final.R')
for (s in names(sobjs)){
  nb = 30
  while (TRUE) {
      try(sobjs[[s]] <- AddModuleScore(sobjs[[s]],
                                       assay='RNA',
                                       features=list(five_gene_sigs_final[['NewPlatelet']]),
                                       name='NewPlatelet',
                                       nbin = nb))
      if ('NewPlatelet1' %in% colnames(sobjs[[s]]@meta.data)) {
          sobjs[[s]]@meta.data['NewPlatelet_per_obj'] <- sobjs[[s]]@meta.data['NewPlatelet1']
          sobjs[[s]]@meta.data['NewPlatelet1'] <- NULL
          sobjs[[s]]@meta.data['NewPlatelet1_per_obj_75sat'] <- saturate(sobjs[[s]]@meta.data[['NewPlatelet_per_obj']], sat=0.75, binary=T)
          sobjs[[s]]@meta.data['NewPlatelet1_per_obj_75nbsat'] <- saturate(sobjs[[s]]@meta.data[['NewPlatelet_per_obj']], sat=0.75, binary=F)
          break
      } else if (nb == 5) {
        stop()
      } else {
        nb <- nb - 1
      }
  }
}

x <- do.call(rbind, 
             lapply(names(sobjs), 
                    function(s){sobjs[[s]]@meta.data[, c('NewPlatelet_per_obj', 'LIBRARY')]}))
x$LIBRARY <- as.factor(x$LIBRARY)
pdf(paste0(prefix, '_per_object_platelet_thresholds.pdf'), width=20, height=20)
print(ggplot(x) + 
  geom_density(aes(x=NewPlatelet_per_obj)) +
  theme(axis.title = element_blank()) +
  facet_wrap(~LIBRARY, ncol=5, scales='free'))
dev.off()

sobjs.orig <- sobjs
sobjs <- sapply(names(sobjs), function(x){ subset(sobjs[[x]], subset= PF4 > 1 | PPBP > 1)})

gene_intersection <- lapply(sobjs, row.names) %>% Reduce(intersect, .)

first_sobj <- sobjs[[names(sobjs)[1]]]
sobjs[[names(sobjs)[1]]] <- NULL
merged_data <- merge(x = first_sobj, y = unname(sobjs))


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

merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 3000)
merged_data <- ScaleData(merged_data, vars.to.regress=c('percent.mt', 'percent.ribo', 'S.Score', 'G2M.Score'), verbose = FALSE)
merged_data <- RunPCA(merged_data, npcs = 30, verbose = FALSE)

save(merged_data, file=paste0(prefix, '_merged_temp.RData'))

pdf(paste0(prefix, '_harmony_convergence.pdf'))
merged_data <- RunHarmony(merged_data, "LIBRARY", plot_convergence = TRUE, max.iter.harmony=30, max.iter.cluster=30)
dev.off()

source('/krummellab/data1/arrao/projects/MVIR1/auxiliary_files/covid_status_all_consent.R')
if (!all(names(covid_status)[!names(covid_status) %in% merged_data@meta.data$SAMPLE.by.SNPs])){
  stop('Not all sample in the covid status file')
}
names(covid_status)[!names(covid_status) %in% merged_data@meta.data$SAMPLE.by.SNPs]

merged_data@meta.data$covid_status = covid_status[merged_data@meta.data$SAMPLE.by.SNPs]

save(merged_data, file=paste0(prefix, '_merged_temp.RData'))

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
                       dims = pcs_to_use,  # Num PCs to use
                       reduction = 'harmony',
                       n.neighbors = 20,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                       min.dist = 0.2,   # Default. Controls the size of the clusters. Should be smaller than spread
                       spread = 3,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                       a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                       b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                       verbose = FALSE)

# Calculate the neighborhood graph
merged_data <- FindNeighbors(merged_data,
                             dims = pcs_to_use,  # Num PCs to use
                             reduction = 'harmony',
                             k.param = 20,  # k for the knn algorithm
                             verbose = FALSE
                             )

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

save(merged_data, file=paste0(prefix, '_merged_unprocessed.RData'))

for (res in c(0.2, 0.4, 0.6)){
  merged_data <- FindClusters(merged_data, verbose = TRUE,
                            algorithm = 1,
                            resolution = res)
  merged_data@meta.data[[paste0('louvain_res', res)]] <- merged_data@meta.data$seurat_clusters

  png(filename=paste0(prefix, '_clusters_louvain_res', res,'_umap.png'), width = 10, height = 10, units = "in", res = 300)
  print(DimPlot(merged_data, group.by=paste0('louvain_res', res)))
  dev.off()
  
  png(filename=paste0(prefix, '_labeled_clusters_louvain_res', res,'_umap.png'), width = 10, height = 10, units = "in", res = 300)
  print(DimPlot(merged_data, group.by=paste0('louvain_res', res), label=T))
  dev.off()

  png(paste(prefix, 'res', res, 'FiveGeneDotPlot.png', sep='_'), height=100 * length(genes_to_consider), 
      width=100 * length(unique(as.vector(merged_data@meta.data[[paste0('louvain_res', res)]]))), units='px', res=150)
  print(DotPlot(merged_data, 
                assay='RNA', 
                #group.by=paste0('louvain_res', res), 
                features = genes_to_consider,
                cols='RdYlBu') + coord_flip())
  dev.off()
}

for (res in c(0.2, 0.4, 0.6)){
  png(paste(prefix, 'res', res, 'SigDotPlot.png', sep='_'), height=100 * length(signatures), 
      width=100 * length(unique(as.vector(merged_data@meta.data[[paste0('louvain_res', res)]]))), units='px', res=150)
  print(signature_dotplot_from_metadata(merged_data@meta.data, signatures=signatures, group.by=paste0('louvain_res', res)))
  dev.off()
  
  png(paste(prefix, 'res', res, 'NBSigDotPlot.png', sep='_'), height=100 * length(signatures2), 
      width=100 * length(unique(as.vector(merged_data@meta.data[[paste0('louvain_res', res)]]))), units='px', res=150)
  print(signature_dotplot_from_metadata(merged_data@meta.data, signatures=signatures2, group.by=paste0('louvain_res', res)))
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

  png(filename=paste(prefix, 'res', res, 'fivegene_umap.png', sep='_'), width = 15, height = 15, units = "in", res = 300)
  print(plot_grid(plotlist = plots1, ncol=4))
  dev.off()

  png(filename=paste(prefix, 'res', res, 'fivegene_sat_umap.png', sep='_'), width = 15, height = 15, units = "in", res = 300)
  print(plot_grid(plotlist = plots2, ncol=4))
  dev.off()

  png(filename=paste(prefix, 'res', res, 'fivegene_nbsat_umap.png', sep='_'), width = 15, height = 15, units = "in", res = 300)
  print(plot_grid(plotlist = plots3, ncol=4))
  dev.off()


  png(filename=paste(prefix, 'res', res, 'fivegene_violins.png', sep='_'), width = 15, height = 15, units = "in", res = 300)
  print(VlnPlot(merged_data, features=names(five_gene_sigs_final), pt.size=0))
  dev.off()

  x <- reshape2::melt(merged_data@meta.data[, c(paste0('louvain_res', res), paste(names(five_gene_sigs_final), '_75sat', sep=''))])
  y <- x %>% 
          mutate(cluster_temp= eval(parse(text=paste0('louvain_res', res)))) %>% 
          group_by(cluster_temp, variable) %>% 
          mutate(mean=mean(value), cluster=as.numeric(as.vector(cluster_temp))) %>% 
          ungroup() %>%
          distinct(cluster, variable, mean)


  png(paste(prefix, 'res', res, 'mean_score_per_cluster.png', sep='_'), width = 15, height = 15, units = "in", res = 300)
  print(ggplot(y) +
            geom_point(aes(x=cluster, y=mean, col=variable)) +
            geom_line(aes(x=cluster, y=mean, col=variable)) +
            scale_x_continuous(breaks=sort(unique(y$cluster)), labels=as.character(sort(unique(y$cluster)))) +
            theme_bw())
  dev.off()
}

# Look at plots and confirm with main plot too
md <- read.table('../merged_colossal_paper_final/rdatas/merged_colossal_paper_final_merged_annotated_metadata.tsv', sep='\t', header=T, row.names=1, stringsAsFactors=F)
md <- md[, c('SAMPLE.by.SNPs', 'coarse_annotations')]
md$rn <- gsub('[^A-Z]', '', rownames(md))
md$rn <- paste(md$rn, md$SAMPLE.by.SNPs, sep='_')
rownames(md) <- md$rn
md$rn <- NULL


clint <- merged_data@meta.data[, c('louvain_res0.2', 'SAMPLE.by.SNPs'), drop=F]

clint$rn <- gsub('[^A-Z]', '', rownames(clint))
clint$rn <- paste(clint$rn, clint$SAMPLE.by.SNPs, sep='_')
rownames(clint) <- clint$rn
clint$rn <- NULL

pooled <- merge(md, clint, by=0)

table(pooled[, c('louvain_res0.2', 'coarse_annotations')])
#              coarse_annotations
#louvain_res0.2 BPlasma eosinophils monodc neuts platelets   RBC   TNK
#            0        8           4      8   174     11132    60   249
#            1       48          31    108  4505       772   138   107
#            2        0           1      0     3      5182     0    18
#            3       37           0    159  1636       887   112    93
#            4        4           5    681    82        96    50    26
#            5        2           1      0    16        85     2   905
#            6        0           1      0    14       132     2   966
#            7      328           1     11     5        11     1     2
#            8        0           0      0     2       164    79     9
#            9       19           0      5    82        60    97    14
#            10       0           0      5     0       291     2    22

annot_table <- setNames(
  c('platelets',
    'neuts',
    'platelets',
    'neuts',
    'myeloid',
    'tnk',
    'tnk',
    'b_pl',
    'rbc',
    'neuts',
    'platelets'),
  c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10')
  )

merged_data$coarse_annots <- factor(annot_table[as.vector(merged_data$louvain_res0.2)],
                                    levels=c('myeloid',
                                             'b_pl',
                                             'neuts',
                                             'platelets',
                                             'rbc',
                                             'tnk'))
save(merged_data, file=paste0(prefix, '_merged_processed.RData'))

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
            file=paste(prefix, 'metadata.tsv', sep='_'),
            sep='\t',
            row.names=T,
            col.names=T,
            quote=F)

Idents(merged_data) <- merged_data$coarse_annots

genes_to_consider = c()
for (gs in names(five_gene_sigs_final)){
  for (gene in five_gene_sigs_final[[gs]]){
    if (!gene %in% genes_to_consider){
      genes_to_consider <- c(genes_to_consider, gene)
    }
  }
}


pdf('F3H_annotated.pdf', width = 10, height = 10)
print(DimPlot(merged_data, group.by='coarse_annots', label=T) + NoLegend())
dev.off()


pdf('F3H.pdf', width = 10, height = 10)
print(DimPlot(merged_data, group.by='coarse_annots') + NoLegend())
dev.off()

pdf('SF5D.pdf', 
    width=10, 
    height=7,
    useDingbats=F)
print(DotPlot(merged_data, 
            assay='RNA',
            features = genes_to_consider,
            cols='RdYlBu') + theme(axis.text.x=element_text(angle=90, hjust=1)))
dev.off()


pdf('SF5Ea.pdf', 
    width=10,
    height=10,
    useDingbats=F)
print(VlnPlot(merged_data, 
              assay='RNA',
              features = 'percent.mt',
              pt.size=0) + theme(axis.text.x=element_text(angle=45, hjust=1)))
dev.off()

pdf('SF5Eb.pdf', 
    width=10,
    height=10,
    useDingbats=F)
print(VlnPlot(merged_data, 
              assay='RNA',
              features = 'percent.ribo',
              pt.size=0) + theme(axis.text.x=element_text(angle=45, hjust=1)))
dev.off()


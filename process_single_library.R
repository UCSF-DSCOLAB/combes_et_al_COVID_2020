## Author: Arjun Arkal Rao
library(Seurat)

data <- Read10X(data.dir=GEX_datadir)
sobj <- CreateSeuratObject(counts = data,
                           project = s,  # The name of the library
                           min.cells = 3,
                           min.features = 100)

# If this was a pooled library
freemuxlet_metadata <- read.table(freemuxlet_file, header=True, row.names=True,
                                  sep="\t", stringsAsFactors=False)
sobj <- AddMetaData(sobj, freemuxlet_metadata)

# Get profiles for Mito, ribo and cell cycling
sobj <- PercentageFeatureSet(sobj,
                             pattern = "^MT-",
                             col.name = "percent.mt")
# These files are in single_cell/auxiliary_files
ribo_genes <- read.table("ribo_genes.tsv", sep = "\t", header=TRUE,
                         stringsAsFactors = FALSE)
ribo_genes <- ribo_genes[ribo_genes[["HUGO"]] %in% rownames(sobj), ]
sobj <- PercentageFeatureSet(sobj,
                             features = ribo_genes[["HUGO"]],
                             col.name = "percent.ribo")

#Store cell cycle state in the object metadata
cc_genes <- read.table("cell_cycle_genes.tsv", sep = "\t", header=TRUE,
                       stringsAsFactors = FALSE)
sobj <- CellCycleScoring(sobj,
                         s.features = cc_genes[cc_genes$stage=="G1-S", "HUGO"],
                         g2m.features = cc_genes[cc_genes$stage=="G2-M", "HUGO"],
                         nbin = 12)

# We would stop at this point, look at some plots for the metrics and save info to
# cutoffs.yml. Those values would be used to filter each library
# The format of cutoffs.yml is
# ---
# percent.mt.high: 20
# percent.ribo.high: 50
# nFeature_RNA.low: 100
# DROPLET.TYPE.keep: SNG

# Every library was filtered with 
# percent.mt.high=[15, 20]
# percent.ribo.high=50
# nFeature_RNA.low=100
# DROPLET.TYPE.keep: SNG

# For the platelet first analyses, we redid this without DROPLET.TYPE.keep=SNG 
# to allow for doublets 

cutoffs <- read_yaml("cutoffs.yml")
# Filter cells with high mito content (dying/dead cells)
keep_rownames = rep(TRUE, dim(sobj@meta.data)[1])
total_cells = c(dim(sobj@meta.data)[1], 100)
additional_text = ""

for (filter_key in names(cutoffs)){
  filter_feature = strsplit(filter_key, split = ".",
                            fixed = TRUE)[[1]]
  suppmsg <- assert_that(length(filter_feature) >= 2,
                         msg=paste0("Names in cutoffs.yml must ",
                                    "be of the form `feature.high`, ",
                                    "`feature.low`, or `feature.keep`. ",
                                    "Got ", filter_key))
  filter_cat <- filter_feature[length(filter_feature)]
  filter_feature <- paste(filter_feature[1:(length(filter_feature)-1)],
                          collapse=".")
  suppmsg <- assert_that(filter_cat %in% c("high", "low", "keep"),
                         msg=paste0("Names in cutoffs.yml must ",
                                    "be of the form `feature.high`, ",
                                    "`feature.low`, or `feature.keep`. ",
                                    "Got ", filter_key))
  if (!filter_feature %in% colnames(sobj@meta.data)){
    print(paste0("Could not find ", filter_feature, " in the metadata ",
                 "for ", s))
    next
  }
  if(!is.na(cutoffs[[filter_key]])){
    if (filter_cat == "high") {
      keep_rownames <- keep_rownames & sobj@meta.data[[filter_feature]] <= cutoffs[[filter_key]]
      filter_cat_text <- " > "
    } else if (filter_cat == "low") {
      keep_rownames <- keep_rownames & sobj@meta.data[[filter_feature]] >= cutoffs[[filter_key]]
      filter_cat_text <- " < "
    } else {
      keep_rownames <- keep_rownames &
        sobj@meta.data[[filter_feature]] %in% strsplit(cutoffs[[filter_key]], ',')[[1]]
      filter_cat_text <- " not in "
    }

    print(paste0("Dropping ", additional_text,
                 total_cells[1] - sum(keep_rownames),
                 " cells (",
                 round((total_cells[1] - sum(keep_rownames))/length(keep_rownames) * 100, 2),
                 "%) for having `", filter_feature,
                 filter_cat_text,
                 cutoffs[[filter_key]], "`."))
    additional_text = "an additional "
    total_cells = c(sum(keep_rownames),
                    round(sum(keep_rownames)/length(keep_rownames) * 100, 2))
  }
}
sobj <- subset(sobj, cells = rownames(sobj@meta.data[keep_rownames, ]))
keep_features <- c()
for (assay in names(sobj@assays)){
  if (assay == "RNA"){
    next
  }
  keep_features <- c(keep_features, rownames(sobj@assays[[assay]]@counts))
}
keep_genes <- rownames(sobj@assays$RNA@counts)[Matrix::rowSums(sobj@assays$RNA@counts>0)>3]
print(paste0("Dropping ",
             (length(rownames(sobj@assays$RNA@counts)) -
              length(keep_genes)),
             " genes for being in fewer than 3 cells"))
sobj <- subset(sobj, features = c(keep_genes, keep_features))

SCTransform(sobj, vars.to.regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"),
            return.only.var.genes = FALSE,
            verbose = FALSE)
RunPCA(sobj,
       verbose = FALSE)
RunPCA(sobj,
       verbose = FALSE)

sobj <- RunUMAP(sobj,
                dims = 1:30,
                n.neighbors = 30,
                min.dist = 0.3,
                spread = 1,
                a = NULL,
                b = NULL,
                verbose = FALSE)

sobj <- FindNeighbors(sobj,
                      dims = 1:30,
                      k.param = 20,
                      verbose = FALSE)

sobj <- FindClusters(sobj, 
                     verbose = FALSE,
                     algorithm = 4)
save(sobj, file=paste0(sample, '_scTransformed_processed.RData'))
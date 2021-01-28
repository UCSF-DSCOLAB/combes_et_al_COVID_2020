## Author: Arjun Arkal Rao

#Load the data
library(Seurat)
setwd('/krummellab/data1/arrao/projects/MVIR1/patient_data_SO/merged/merged_colossal_paper_final')

load('rdatas/merged_colossal_paper_final_merged_annotated.RData')
bpl <- readRDS('annotated_res_0.6/monocles/inputs/b_plasma_cleaned.rds')

table(merged_data@meta.data[colnames(bpl)[bpl$subtype_final_annotation %in%c('plasma cell 1', 'plasma cell 2')],'fine_annotations'])
#BPlasma_1 BPlasma_2 BPlasma_3 
#        5       663         4


translation = c(
  'eosinophils'='Eosinophils',
  'BPlasma'='B cells',
  'TNK'='T/NK cells',
  'monodc'='MPC',
  'platelets'='Platelets',
  'neuts'='Neutrophils')


Idents(merged_data) <- merged_data$coarse_annotations
merged_data <- subset(merged_data, idents = c('RBC'), invert = TRUE)

merged_data$final_annotations <- translation[as.vector(merged_data$coarse_annotations)]
merged_data$final_annotations[merged_data$fine_annotations == 'BPlasma_2'] = 'Plasma cells'
merged_data$final_annotations <- factor(merged_data$final_annotations,
                                        levels=c('Eosinophils','Plasma cells','B cells','T/NK cells',
                                                 'MPC','Platelets','Neutrophils'))

Idents(merged_data) <- merged_data$final_annotations
pdf('Fig_1Ba.pdf', width=4, height=6, useDingbats=F)
DimPlot(merged_data, reduction = "umap") + 
  scale_colour_manual(values = rev(brewer.pal(7, "Paired"))) +
  NoLegend()
dev.off()

pdf('Fig_1Ba_labeled.pdf', width=4, height=6, useDingbats=F)
DimPlot(merged_data, reduction = "umap", label=T) + 
  scale_colour_manual(values = rev(brewer.pal(7, "Paired"))) +
  NoLegend()
dev.off()


merged_data$patient_name <- gsub('[.][BD].*', '',
                                 gsub('X1.H', 'Donor ', 
                                      gsub('M1.H', 'Patient ', merged_data$SAMPLE.by.SNPs)))

patient_levels <-  unique(merged_data$patient_name)
patient_levels <- patient_levels[order(gsub('[0-9]*', '', patient_levels), -as.numeric(gsub('[^0-9]*', '', patient_levels)), decreasing=T)]

merged_data$patient_name <- factor(merged_data$patient_name, levels=patient_levels)

pdf('SFig_1Ba_c.pdf', width=10, height=6, useDingbats=F)
DimPlot(merged_data, reduction = "umap" , group.by='patient_name')
dev.off()


pdf('SFig_1C.pdf', width=8, height=3, useDingbats=F)
DotPlot(merged_data, features = c('G0S2', 'FCGR3B', 'CSF3R', 'S100A9', 'S100A8',
                                  'CLEC1B', 'RGS18', 'RGS10', 'PF4', 'PPBP',
                                  'MS4A7', 'VCAN', 'CD14', 'LYZ',
                                  'GNLY', 'PRF1', 'IL7R', 'CD3D', 'CD3E',
                                  'HLA-DQA1', 'CD79A', 'MS4A1', 'IGHA1', 'JCHAIN',
                                  'CLC', 'FCER1A', 'GATA2', 'HDC'),
cols = "RdBu", assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + NoLegend()
dev.off()

pdf('SFig_1C_labeled.pdf', width=15, height=5, useDingbats=F)
DotPlot(merged_data, features = c('G0S2', 'FCGR3B', 'CSF3R', 'S100A9', 'S100A8',
                                  'CLEC1B', 'RGS18', 'RGS10', 'PF4', 'PPBP',
                                  'MS4A7', 'VCAN', 'CD14', 'LYZ',
                                  'GNLY', 'PRF1', 'IL7R', 'CD3D', 'CD3E',
                                  'HLA-DQA1', 'CD79A', 'MS4A1', 'IGHA1', 'JCHAIN',
                                  'CLC', 'FCER1A', 'GATA2', 'HDC'),
cols = "RdBu", assay = "RNA") +  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14))
dev.off()
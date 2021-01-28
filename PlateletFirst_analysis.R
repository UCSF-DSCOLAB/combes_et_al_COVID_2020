##### Figure 5I
#### Create correlation plots between PlateletFirst and merged_colossal only of Covid pos cells

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

setwd("~/PlateletFirst")
library(readr)
plateletFirst <- read_csv("merged_colossal_paper_final_platelets_first_metadata.csv")

########################################################
##################Create correlation plots between PlateletFirst and merged_colossal
########################################################
colnames(plateletFirst)

#create file 'y' with table containing SAMPLE.by.SNPs as rows and seurat_clusters as columns with corresponding event #s
y <- table(plateletFirst[,c('SAMPLE.by.SNPs', 'coarse_annots')])

#load reshape2
library(reshape2)

#melt table into one column by id='SAMPLE.by.SNPs
PlateletFirstFreq <- melt(y, id='SAMPLE.by.SNPs')

#make data.frame from sobj that contains SAMPLE.by.SNPs event #s; call 'temp'
temp <- table(plateletFirst$SAMPLE.by.SNPs)
temp <- data.frame(temp)

#add a column to PlateletFirstFreq called 'total' that contains total # of cells per SAMPLE.by.SNPs pulled from 'temp' data.frame in column 'Freq'
PlateletFirstFreq$total <- temp$Freq[match(PlateletFirstFreq$SAMPLE.by.SNPs, temp$Var1)]
#add a column to PlateletFirstFreq called 'freq' that gives the fraction of cells of each coarse_annots per total # of cells in sample
PlateletFirstFreq$freqPlat <- PlateletFirstFreq$value/PlateletFirstFreq$total

#create cell type-specific subsets (PF = plateletFirst)
PF_b_pl <- PlateletFirstFreq[PlateletFirstFreq$coarse_annots == 'b_pl', ]
PF_myeloid <- PlateletFirstFreq[PlateletFirstFreq$coarse_annots == 'myeloid', ]
PF_neuts <- PlateletFirstFreq[PlateletFirstFreq$coarse_annots == 'neuts', ]
PF_platelets <- PlateletFirstFreq[PlateletFirstFreq$coarse_annots == 'platelets', ]
PF_rbc <- PlateletFirstFreq[PlateletFirstFreq$coarse_annots == 'rbc', ]
PF_tnk <- PlateletFirstFreq[PlateletFirstFreq$coarse_annots == 'tnk', ]

##### use GLOBAL_COUNT_mod from Tristan
library(readxl)
MC <- read_excel("GLOBAL_COUNTS_mod.xlsx")
#add PF values to MC (= merged_colossal)
MC$PF_neuts <- PF_neuts$freqPlat[match(MC$SAMPLE.by.SNPs, PF_neuts$SAMPLE.by.SNPs)]
MC$PF_myeloid <- PF_myeloid$freqPlat[match(MC$SAMPLE.by.SNPs, PF_myeloid$SAMPLE.by.SNPs)]
MC$PF_tnk <- PF_tnk$freqPlat[match(MC$SAMPLE.by.SNPs, PF_tnk$SAMPLE.by.SNPs)]
MC$PF_b_pl <- PF_b_pl$freqPlat[match(MC$SAMPLE.by.SNPs, PF_b_pl$SAMPLE.by.SNPs)]

mycolorsseverity <- setNames(c("grey40", "orange", "orangered2"), c('CTRL', 'MILD', 'SEVERE'))

##################
###### plot only covid_pos patients
#################
MC_pos <- MC[MC$covid_status == 'POS', ]
MC_ctrl <- MC[MC$covid_status == 'CTRL', ]
MC_ctrl_pos <- rbind(MC_ctrl, MC_pos)

#### make ratios between plateletFirst (PF) to merge_colossal (MC) by cell type
MC_ctrl_pos$PFtoMC_neuts <- MC_ctrl_pos$PF_neuts/MC_ctrl_pos$neuts...14
MC_ctrl_pos$PFtoMC_myeloid <- MC_ctrl_pos$PF_myeloid/MC_ctrl_pos$monodc...17
MC_ctrl_pos$PFtoMC_bplasma <- MC_ctrl_pos$PF_b_pl/MC_ctrl_pos$BPlasma...18
MC_ctrl_pos$PFtoMC_tnk <- MC_ctrl_pos$PF_tnk/MC_ctrl_pos$TNK...16

##################
###### plot only covid_pos patients; without platelet in denominator
#################

#### make PlateletFirst fractions without platelet in denominator; 'PFnP' -> PlateletFirst no plateles
temp <- table(plateletFirst$SAMPLE.by.SNPs)
temp <- data.frame(temp)
MC_ctrl_pos$totalPF <- temp$Freq[match(MC_ctrl_pos$SAMPLE.by.SNPs, temp$Var1)]
MC_ctrl_pos$PF_platelets <- PF_platelets$value[match(MC_ctrl_pos$SAMPLE.by.SNPs, PF_platelets$SAMPLE.by.SNPs)]
MC_ctrl_pos$totalPFnP <- MC_ctrl_pos$totalPF - MC_ctrl_pos$PF_platelets
MC_ctrl_pos$PF_neuts_total <- PF_neuts$value[match(MC_ctrl_pos$SAMPLE.by.SNPs, PF_neuts$SAMPLE.by.SNPs)]
MC_ctrl_pos$PF_myeloid_total <- PF_myeloid$value[match(MC_ctrl_pos$SAMPLE.by.SNPs, PF_myeloid$SAMPLE.by.SNPs)]
MC_ctrl_pos$PF_tnk_total <- PF_tnk$value[match(MC_ctrl_pos$SAMPLE.by.SNPs, PF_tnk$SAMPLE.by.SNPs)]
MC_ctrl_pos$PF_b_pl_total <- PF_b_pl$value[match(MC_ctrl_pos$SAMPLE.by.SNPs, PF_b_pl$SAMPLE.by.SNPs)]
MC_ctrl_pos$PF_rbc_total <- PF_rbc$value[match(MC_ctrl_pos$SAMPLE.by.SNPs, PF_rbc$SAMPLE.by.SNPs)]

MC_ctrl_pos$PFnP_neuts <- MC_ctrl_pos$PF_neuts_total/MC_ctrl_pos$totalPFnP
MC_ctrl_pos$PFnP_myeloid <- MC_ctrl_pos$PF_myeloid_total/MC_ctrl_pos$totalPFnP
MC_ctrl_pos$PFnP_tnk <- MC_ctrl_pos$PF_tnk_total/MC_ctrl_pos$totalPFnP
MC_ctrl_pos$PFnP_b_pl <- MC_ctrl_pos$PF_b_pl_total/MC_ctrl_pos$totalPFnP

#### make MC fractions without platelet in denominator; MCnP -> merged colossal no platelets
MC_ctrl_pos$totalMCnP <- MC_ctrl_pos$TOTAL...12 - MC_ctrl_pos$platelets...5

MC_ctrl_pos$MCnP_neuts <- MC_ctrl_pos$neuts...4/MC_ctrl_pos$totalMCnP
MC_ctrl_pos$MCnP_myeloid <- MC_ctrl_pos$monodc...7/MC_ctrl_pos$totalMCnP
MC_ctrl_pos$MCnP_tnk <- MC_ctrl_pos$TNK...6/MC_ctrl_pos$totalMCnP
MC_ctrl_pos$MCnP_b_pl <- MC_ctrl_pos$BPlasma...8/MC_ctrl_pos$totalMCnP

ggplot(MC_ctrl_pos, aes(x=MCnP_neuts, y=PFnP_neuts, col = Qualitative_score)) +
  geom_point(size = 4) + 
  xlim(0, 1) + ylim(0, 1) +
  geom_abline(slope=1, intercept =0) +
  #geom_smooth(method=lm, aes(fill=Qual_score), se = F) +
  theme_classic() +
  scale_color_manual(values=mycolorsseverity) +
  scale_fill_manual(values=mycolorsseverity) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(legend.position = "none")
ggplot(MC_ctrl_pos, aes(x=MCnP_myeloid, y=PFnP_myeloid, col = Qualitative_score)) +
  geom_point(size = 4) + 
  xlim(0, 0.4) + ylim(0, 0.4) +
  geom_abline(slope=1, intercept =0) +
  #geom_smooth(method=lm, aes(fill=Qual_score), se = F) +
  theme_classic() +
  scale_color_manual(values=mycolorsseverity) +
  scale_fill_manual(values=mycolorsseverity) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(legend.position = "none")
ggplot(MC_ctrl_pos, aes(x=MCnP_tnk, y=PFnP_tnk, col = Qualitative_score)) +
  geom_point(size = 4) + 
  xlim(0, 0.75) + ylim(0, 0.75) +
  geom_abline(slope=1, intercept =0) +
  #geom_smooth(method=lm, aes(fill=Qual_score), se = F) +
  theme_classic() +
  scale_color_manual(values=mycolorsseverity) +
  scale_fill_manual(values=mycolorsseverity) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(legend.position = "none")
ggplot(MC_ctrl_pos, aes(x=MCnP_b_pl, y=PFnP_b_pl, col = Qualitative_score)) +
  geom_point(size = 4) + 
  xlim(0, 0.2) + ylim(0, 0.2) +
  geom_abline(slope=1, intercept =0) +
  #geom_smooth(method=lm, aes(fill=Qual_score), se = F) +
  theme_classic() +
  scale_color_manual(values=mycolorsseverity) +
  scale_fill_manual(values=mycolorsseverity) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(legend.position = "none")

#### make ratios between plateletFirst (PF) to merge_colossal (MC) by cell type
MC_ctrl_pos$PFnPtoMCnP_neuts <- MC_ctrl_pos$PFnP_neuts/MC_ctrl_pos$MCnP_neuts
MC_ctrl_pos$PFnPtoMCnP_myeloid <- MC_ctrl_pos$PFnP_myeloid/MC_ctrl_pos$MCnP_myeloid
MC_ctrl_pos$PFnPtoMCnP_bplasma <- MC_ctrl_pos$PFnP_b_pl/MC_ctrl_pos$MCnP_b_pl
MC_ctrl_pos$PFnPtoMCnP_tnk <- MC_ctrl_pos$PFnP_tnk/MC_ctrl_pos$MCnP_tnk

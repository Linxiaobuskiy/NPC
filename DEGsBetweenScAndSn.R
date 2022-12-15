rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)

library(scater)
library(patchwork)


load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse



## DEG
NPC.tmp <- subset(NPC.combined, subset = Annotation.Coarse %in% c("B"))
DefaultAssay(NPC.tmp) <- "RNA"
Idents(NPC.tmp) <- 'SeqTech'
markers <- FindAllMarkers(object = NPC.tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


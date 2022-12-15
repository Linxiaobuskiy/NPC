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

#############################################################
#########################Figure S1a##########################
#############################################################

p1 <- DimPlot(NPC.combined, reduction = "umap", group.by = 'orig.ident', raster=FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
p2 <- DimPlot(NPC.combined, reduction = "umap", group.by = 'orig.ident',label=TRUE)

# save p1 as 8inch*8inch pdf
# save p2 as 8inch*8inch pdf
# Transfer p1 to semi-png figure
# Adjust p1 to 50mm*50mm
# Add labels and legend from p2 to p1
# Adjust word to 5pt, label gap 2mm, square-word left-left 2.5mm
# Add Axes

###############################################################
#########################Figure S1b,c##########################
###############################################################

rm(list = ls(all = TRUE))
library(imcRtools)
library(cytomapper)
library(dittoSeq)
library(scater)
library(patchwork)
library(batchelor)
library(cowplot)
library(dittoSeq)
library(viridis)
library(Rphenograph)
library(igraph)
library(bluster)
library(BiocParallel)
library(ggplot2)
library(scran)
library(scuttle)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

spe <- readRDS("G:/IMC/NPC/Rdata/step2.spe.k30.annotated.rds")

##UMAP Figure b

dittoDimPlot(spe, var = "celltype.fine", 
             reduction.use = "UMAP", size = 0.01,
             do.label = T) + NoLegend() + NoAxes() + labs(title = NULL)

dittoColors(1)[1:23]
# save p as 8inch*8inch pdf
# Transfer p1 to semi-png figure
# Adjust p to 50mm*50mm

##UMAP Figure c
dittoDimPlot(spe, var = "patient_id", 
             reduction.use = "UMAP", size = 0.01,
             do.label = F) + NoLegend() + NoAxes() + labs(title = NULL)

dittoColors(1)[1:9]
# save p as 8inch*8inch pdf
# Transfer p1 to semi-png figure
# Adjust p to 50mm*50mm
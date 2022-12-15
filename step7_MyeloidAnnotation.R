rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)
library(scDblFinder)

Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.1")
library(infercnv)

load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse
NPC.combined@meta.data$Annotation.Immune <- NPC.combined@meta.data$Annotation.Coarse
NPC.combined@meta.data$number <- 1:length(NPC.combined@meta.data$orig.ident)

`%notin%` <- Negate(`%in%`)
# table(NPC.combined@meta.data$Annotation.Coarse,NPC.combined@meta.data$orig.ident)

###############################################################
##########################   scRNA   ##########################
###############################################################
# First order annotation
NPC.immune <- subset(NPC.combined, subset = SeqTech %in% c("Single-cell"))
#NPC.immune <- subset(NPC.immune, subset = Annotation.Coarse %notin% c("Epithelial","Fibroblast"))
table(NPC.immune@meta.data$Annotation.Coarse)

###############################################################
## Myeloid
NPC.tmp <- subset(NPC.immune, subset = Annotation.Coarse %in% c("Myeloid"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.15)#resoulution越小分的越粗
p1 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE)
p2 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE,group.by = "orig.ident")
p1+p2

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("Mon/Mac_FN1+","Mon/Mac_C1Q+","Mon/Mac_C1Q+","DC_CD52+","DC_CCL5+","DC_CCR7+","DC_CLEC9A+","Neutrophil")

names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Immune <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Immune")
NPC.combined@meta.data$Annotation.Immune[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Immune)


###############################################################
##########################   snRNA   ##########################
###############################################################
# First order annotation
NPC.immune <- subset(NPC.combined, subset = SeqTech %in% c("Single-nucleus"))
#NPC.immune <- subset(NPC.immune, subset = Annotation.Coarse %notin% c("Epithelial","Fibroblast"))
table(NPC.immune@meta.data$Annotation.Coarse)
###############################################################
## Myeloid
NPC.tmp <- subset(NPC.immune, subset = Annotation.Coarse %in% c("Myeloid"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.1)#resoulution越小分的越粗
p1 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE)
p2 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE,group.by = "orig.ident")
p1+p2





Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("Mon/Mac_CTSB+","Mon/Mac_IL18+","Mon/Mac_FCGBP+","DC_CD83+","Mon/Mac_CLNK+","Mon/Mac_PTPRS+")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Immune <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Immune")
NPC.combined@meta.data$Annotation.Immune[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Immune)


#############################################################
##########################   End   ##########################
#############################################################
Annotation.Immune <- NPC.combined@meta.data$Annotation.Immune
save(Annotation.Immune, file = "G:/scRNA data/DS0/step7.AnnotatedImmune.RData")
# DimPlot(NPC.combined, reduction = "umap", group.by = 'Annotation.Coarse',label = TRUE)
# DimPlot(NPC.combined, reduction = "umap", group.by = 'orig.ident',raster = FALSE)
# new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10")



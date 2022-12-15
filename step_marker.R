rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(rhdf5)
library(scater)

load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse
#DimPlot(NPC.combined, reduction = "umap", group.by = 'Annotation.Coarse',label = TRUE)


# tmp <- NPC.combined
# DefaultAssay(tmp) <- "RNA"
# Idents(tmp) <- 'Annotation.Coarse'
# markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# save(markers, file = "G:/scRNA data/DS0/step.test.marker.RData")
load("G:/scRNA data/DS0/step.test.marker.RData")


t <- subset(markers,subset = avg_log2FC>1.5)
t <- t$gene
t <- unique(t)
t <- t[!grepl("^RP[SL]", t, ignore.case = F)]
NPC.combined@assays$RNA@var.features <- t
#NPC.combined <- FindVariableFeatures(NPC.combined, selection.method = "vst", nfeatures = 2000)

#NPC.combined <- ScaleData(NPC.combined,  vars.to.regress=c("S.Score", "G2M.Score"))
NPC.combined <- ScaleData(NPC.combined)

NPC.combined <- RunPCA(NPC.combined)
ElbowPlot(NPC.combined)
dims_parameter <- 15
NPC.combined <- RunUMAP(NPC.combined, dims = 1:dims_parameter)
#NPC.combined <- FindNeighbors(NPC.combined, dims = 1:dims_parameter)
#NPC.combined <- FindClusters(NPC.combined, resolution = 1)#resoulution越小分的越粗
#DimPlot(NPC.combined, reduction = "umap",label = TRUE,raster=FALSE)

NPC.combined@meta.data$tmp <- paste(NPC.combined@meta.data$orig.ident, NPC.combined@meta.data$Annotation.Coarse,sep=': ')

DimPlot(NPC.combined, reduction = "umap", group.by = 'tmp',label = TRUE,raster=FALSE)
DimPlot(NPC.combined, reduction = "umap", group.by = 'SeqTech',label = TRUE,raster=FALSE)
p1 <- DimPlot(NPC.combined, reduction = "umap", group.by = 'Annotation.Coarse',label = TRUE,raster=FALSE)
p2 <- DimPlot(NPC.combined, reduction = "umap", group.by = 'orig.ident',label = TRUE,raster=FALSE)
p1+p2



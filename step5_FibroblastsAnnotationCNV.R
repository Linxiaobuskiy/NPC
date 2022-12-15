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


NPC.Fib <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Fibroblast"))
NPC.Fib@meta.data$number <- 1:length(NPC.Fib@meta.data$orig.ident)
NPC.Fib@meta.data$Annotation.Fib <- "Fibroblast"
DimPlot(NPC.Fib, reduction = "umap", group.by = "Annotation.Fib",label = TRUE)

#############################################################
##########################   N03   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.Fib, subset = orig.ident %in% c("N03"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 10
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.4)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("0","1","2","3","4")
#new.cluster.ids <- c("Matrix Fibroblast","Matrix Fibroblast","Inflammatory Fibroblast","Matrix Fibroblast","Vascular Fibroblast")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Fib <- NPC.tmp@active.ident
# NPC.Fib@meta.data$Annotation.Fib[NPC.tmp@meta.data$number] <- paste("Normal_",NPC.tmp@meta.data$Annotation.Fib,sep="")
NPC.Fib@meta.data$Annotation.Fib[NPC.tmp@meta.data$number] <- "non-CNV"
DimPlot(NPC.Fib, reduction = "umap", group.by = "Annotation.Fib")

#############################################################
##########################   P04   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.combined, subset = orig.ident %in% c("P04"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 10
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.2)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

# Fibroblasts
NPC.tmp <- subset(NPC.tmp, subset = seurat_clusters %in% c(2,5))


# Annotation
Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("non-CNV","CNV")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Fib <- NPC.tmp@active.ident

# a <- colnames(NPC.tmp)
# b <- colnames(NPC.Fib)
# NPC.Fib@meta.data$Annotation.Fib[match(a,b)]
# NPC.Fib@meta.data$Annotation.Fib[match(colnames(NPC.tmp),colnames(NPC.Fib))] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Fib,sep="")
NPC.Fib@meta.data$Annotation.Fib[match(colnames(NPC.tmp),colnames(NPC.Fib))] <- as.character(NPC.tmp@meta.data$Annotation.Fib)
DimPlot(NPC.Fib, reduction = "umap", group.by = "Annotation.Fib")

#############################################################
##########################   P05   ##########################
#############################################################

# First order annotation
NPC.tmp <- subset(NPC.combined, subset = orig.ident %in% c("P05"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 1)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

# Fibroblasts
NPC.tmp <- subset(NPC.tmp, subset = seurat_clusters %in% c(7,12))

# Annotation
Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("CNV","non-CNV")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Fib <- NPC.tmp@active.ident

# a <- colnames(NPC.tmp)
# b <- colnames(NPC.Fib)
# NPC.Fib@meta.data$Annotation.Fib[match(a,b)]
# NPC.Fib@meta.data$Annotation.Fib[match(colnames(NPC.tmp),colnames(NPC.Fib))] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Fib,sep="")
NPC.Fib@meta.data$Annotation.Fib[match(colnames(NPC.tmp),colnames(NPC.Fib))] <- as.character(NPC.tmp@meta.data$Annotation.Fib)
DimPlot(NPC.Fib, reduction = "umap", group.by = "Annotation.Fib")

#############################################################
##########################   P06   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.combined, subset = orig.ident %in% c("P06"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 10
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.2)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

# Fibroblasts
NPC.tmp <- subset(NPC.tmp, subset = seurat_clusters %in% c(5))
NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 10
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.3)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("non-CNV","CNV")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Fib <- NPC.tmp@active.ident
# a <- colnames(NPC.tmp)
# b <- colnames(NPC.Fib)
# NPC.Fib@meta.data$Annotation.Fib[match(a,b)]
# NPC.Fib@meta.data$Annotation.Fib[match(colnames(NPC.tmp),colnames(NPC.Fib))] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Fib,sep="")
NPC.Fib@meta.data$Annotation.Fib[match(colnames(NPC.tmp),colnames(NPC.Fib))] <- as.character(NPC.tmp@meta.data$Annotation.Fib)
DimPlot(NPC.Fib, reduction = "umap", group.by = "Annotation.Fib")

#############################################################
##########################   P07   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.combined, subset = orig.ident %in% c("P07"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.2)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

# Fibroblasts
NPC.tmp <- subset(NPC.tmp, subset = seurat_clusters %in% c(5))
NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 30
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.8)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE,group.by = "seurat_clusters")

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("CNV","non-CNV")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Fib <- NPC.tmp@active.ident
# a <- colnames(NPC.tmp)
# b <- colnames(NPC.Fib)
# NPC.Fib@meta.data$Annotation.Fib[match(a,b)]
# NPC.Fib@meta.data$Annotation.Fib[match(colnames(NPC.tmp),colnames(NPC.Fib))] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Fib,sep="")
NPC.Fib@meta.data$Annotation.Fib[match(colnames(NPC.tmp),colnames(NPC.Fib))] <- as.character(NPC.tmp@meta.data$Annotation.Fib)
DimPlot(NPC.Fib, reduction = "umap", group.by = "Annotation.Fib")

#############################################################
##########################   End   ##########################
#############################################################
Annotation.Fib <- NPC.Fib@meta.data$Annotation.Fib
save(Annotation.Fib, file = "G:/scRNA data/DS0/step5.FibroblastAnnotationCNV.RData")




rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(rhdf5)
library(scater)

load("G:/scRNA data/DS0/step1.NPC.Merge.RData")

NPC.combined@meta.data$number <- 1:length(NPC.combined@meta.data$orig.ident)
# Coarse
NPC.combined@meta.data$Annotation.Coarse <- "Immune Cell"

load("G:/scRNA data/DS0/step.test.marker.RData")
Ref = markers
topRef <- Ref %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#############################################################
########################## N01 N02 ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.combined, subset = orig.ident %in% c("N01","N02"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 10
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.3)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("B","NK/T","NK/T","NK/T","B","B","B","Plasma B","Myeloid","Epithelial","Epithelial","NK/T")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)

# Second order annotation
NPC.tmp <- subset(NPC.tmp, subset = Annotation.Coarse %in% c("Epithelial"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.3)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("Epithelial","Epithelial","Epithelial","Epithelial","pDC","Mast")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)


#############################################################
##########################   N03   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.combined, subset = orig.ident %in% c("N03"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.3)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("B","NK/T","B","Plasma B","NK/T","NK/T","Myeloid","7","Fibroblast","Epithelial","Epithelial")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)

# Second order annotation
NPC.tmp <- subset(NPC.tmp, subset = Annotation.Coarse %in% c("7"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.3)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("B","B","NK/T")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)


#############################################################
##########################   P01   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.combined, subset = orig.ident %in% c("P01"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.4)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("B","NK/T","NK/T","NK/T","NK/T","NK/T","Myeloid","tmp","B","tmp","Myeloid","pDC","Mast")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)

# Second order annotation
NPC.tmp <- subset(NPC.tmp, subset = Annotation.Coarse %in% c("tmp"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 2)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("B","Plasma B","B","NK/T","B","NK/T","NK/T","B","Epithelial","B","Myeloid")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)


#############################################################
##########################   P02   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.combined, subset = orig.ident %in% c("P02"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.2)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("NK/T","NK/T","B","NK/T","Myeloid","Plasma B","tmp","Myeloid","Epithelial","Myeloid","pDC")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)

# Second order annotation
NPC.tmp <- subset(NPC.tmp, subset = Annotation.Coarse %in% c("tmp"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.3)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("NK/T","NK/T","NK/T","B")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)


#############################################################
##########################   P03   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.combined, subset = orig.ident %in% c("P03"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.1)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("NK/T","NK/T","Myeloid","NK/T","B","tmp","Epithelial","NK/T","Myeloid","pDC")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)

# Second order annotation
NPC.tmp <- subset(NPC.tmp, subset = Annotation.Coarse %in% c("tmp"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.3)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("Myeloid","Epithelial","Myeloid","Plasma B") # cluster 0 is epithelial-myeloid dul-like cells
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)

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

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("Epithelial","Epithelial","Fibroblast","Epithelial","tmp","Fibroblast","Epithelial","Epithelial","Epithelial")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)

# Second order annotation
NPC.tmp <- subset(NPC.tmp, subset = Annotation.Coarse %in% c("tmp"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.3)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("Myeloid","Epithelial","B","Epithelial","B","Myeloid")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)


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

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Fibroblast","Epithelial","Epithelial","Epithelial","Epithelial","Fibroblast","Epithelial","Epithelial")
#new.cluster.ids <- c("Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Fibroblast","Epithelial","B","Epithelial","Epithelial","Fibroblast","Epithelial","Fibroblast")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)


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

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("NK/T","Epithelial","Myeloid","Plasma B","Epithelial","Fibroblast","Endothelial")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)


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

Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("Epithelial","NK/T","Plasma B","Myeloid","Epithelial","Fibroblast","Myeloid","Endothelial")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Coarse <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Coarse")
NPC.combined@meta.data$Annotation.Coarse[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Coarse)


#############################################################
##########################   End   ##########################
#############################################################
Annotation.Coarse <- NPC.combined$Annotation.Coarse
save(Annotation.Coarse, file = "G:/scRNA data/DS0/step2.Annotation.RData")
# load("G:/scRNA data/DS0/step2.Annotation.RData")
# NPC.combined$Annotation.Coarse <- Annotation.Coarse
# DimPlot(NPC.combined, reduction = "umap", group.by = 'Annotation.Coarse',label = TRUE)
# DimPlot(NPC.combined, reduction = "umap", group.by = 'orig.ident',raster = FALSE)
# new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10")

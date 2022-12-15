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


NPC.Epi <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Epithelial"))
NPC.Epi@meta.data$number <- 1:length(NPC.Epi@meta.data$orig.ident)
NPC.Epi@meta.data$Annotation.Epi <- "Epithelial"
DimPlot(NPC.Epi, reduction = "umap", group.by = "Annotation.Epi",label = TRUE)

#############################################################
######################## N01 and N02 ########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.Epi, subset = orig.ident %in% c("N01","N02"))
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
new.cluster.ids <- c("Ciliated","Unciliated","Unciliated")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Epi <- NPC.tmp@active.ident
NPC.Epi@meta.data$Annotation.Epi[NPC.tmp@meta.data$number] <- paste("scNormal_",NPC.tmp@meta.data$Annotation.Epi,sep="")
DimPlot(NPC.Epi, reduction = "umap", group.by = "Annotation.Epi")

#############################################################
##########################   N03   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.Epi, subset = orig.ident %in% c("N03"))
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
new.cluster.ids <- c("Ciliated","Ciliated","Unciliated","Unciliated")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Epi <- NPC.tmp@active.ident
NPC.Epi@meta.data$Annotation.Epi[NPC.tmp@meta.data$number] <- paste("snNormal_",NPC.tmp@meta.data$Annotation.Epi,sep="")
DimPlot(NPC.Epi, reduction = "umap", group.by = "Annotation.Epi")

#############################################################
######################## P01 and P02 ########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.Epi, subset = orig.ident %in% c("P01","P02"))

NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 1)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

NPC.tmp@meta.data$Annotation.Epi <- "Unciliated"
NPC.Epi@meta.data$Annotation.Epi[NPC.tmp@meta.data$number] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Epi,sep="")
DimPlot(NPC.Epi, reduction = "umap", group.by = "Annotation.Epi")

#############################################################
##########################   P03   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.Epi, subset = orig.ident %in% c("P03"))
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
new.cluster.ids <- c("Unciliated_MT+","Unciliated_MT-")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Epi <- NPC.tmp@active.ident
NPC.Epi@meta.data$Annotation.Epi[NPC.tmp@meta.data$number] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Epi,sep="")
DimPlot(NPC.Epi, reduction = "umap", group.by = "Annotation.Epi")

#############################################################
##########################   P04   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.Epi, subset = orig.ident %in% c("P04"))
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
new.cluster.ids <- c("Unciliated_IG+","Unciliated_IG+","Unciliated_IG+BPIFA1+","Unciliated_IG-BPIFA1+","Ciliated_TSPAN1+","Unciliated_ADAM23+","Ciliated_AGBL4+")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Epi <- NPC.tmp@active.ident
NPC.Epi@meta.data$Annotation.Epi[NPC.tmp@meta.data$number] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Epi,sep="")
DimPlot(NPC.Epi, reduction = "umap", group.by = "Annotation.Epi")

#############################################################
##########################   P05   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.Epi, subset = orig.ident %in% c("P05"))
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
# new.cluster.ids <- c("Ciliated_IG+","Ciliated_IG+","Ciliated_IG+","Ciliated_AGBL4+","Ciliated_BPIFA1+","Unciliated","Ciliated_XKR4+","Ciliated_IG+")
# new.cluster.ids <- c("Ciliated_IG+","Ciliated_IG+","Ciliated_IG+","Ciliated_AGBL4+","Ciliated_BPIFA1+","Unciliated_BPIFA1+","Ciliated_XKR4+","Ciliated_IG+","Ciliated_VWF+")
new.cluster.ids <- c("Ciliated_IG+","Ciliated_IG+","Ciliated_IG+","Ciliated_BPIFA1+",
                     "Ciliated_AGBL4+","Unciliated_BPIFA1+","Ciliated_XKR4+","Ciliated_VWF+")


names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Epi <- NPC.tmp@active.ident
NPC.Epi@meta.data$Annotation.Epi[NPC.tmp@meta.data$number] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Epi,sep="")
DimPlot(NPC.Epi, reduction = "umap", group.by = "Annotation.Epi")


#############################################################
##########################   P06   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.Epi, subset = orig.ident %in% c("P06"))
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
new.cluster.ids <- c("Unciliated_IG+","Unciliated_ABCA13+","Unciliated_MT+","Unciliated","Unciliated_ZMAT4+")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Epi <- NPC.tmp@active.ident
NPC.Epi@meta.data$Annotation.Epi[NPC.tmp@meta.data$number] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Epi,sep="")
DimPlot(NPC.Epi, reduction = "umap", group.by = "Annotation.Epi")

#############################################################
##########################   P07   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.Epi, subset = orig.ident %in% c("P07"))
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
new.cluster.ids <- c("Unciliated","Unciliated","Unciliated_MT+")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Epi <- NPC.tmp@active.ident
NPC.Epi@meta.data$Annotation.Epi[NPC.tmp@meta.data$number] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Epi,sep="")
DimPlot(NPC.Epi, reduction = "umap", group.by = "Annotation.Epi")


#############################################################
##########################   End   ##########################
#############################################################
Annotation.Epi <- NPC.Epi@meta.data$Annotation.Epi
save(Annotation.Epi, file = "G:/scRNA data/DS0/step4.AnnotatedEpithelial.RData")
# DimPlot(NPC.combined, reduction = "umap", group.by = 'Annotation.Coarse',label = TRUE)
# DimPlot(NPC.combined, reduction = "umap", group.by = 'orig.ident',raster = FALSE)
# new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10")






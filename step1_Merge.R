rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)
library(scDblFinder)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

####################################################################
#########################10x combined data##########################
####################################################################
load.list <- c("N01","N02","N03","P01","P02","P03","P04","P05","P06","P07")

for (i in load.list){
  path_name <- paste("G:/scRNA data/DS0/",i,".10x.RData", sep="")
  #path_name <- paste("/media/linxiaobuskiy/NPC Project/scRNA data/DS0/",i,".10x.RData", sep="")
  load(path_name)
  NPC.10x <- RenameCells(NPC.10x, add.cell.id = i)
  assign(i,NPC.10x)
  # NPC.10x <- FindVariableFeatures(NPC.10x, selection.method = "vst", nfeatures = 2000)
  rm(NPC.10x)
}

NPC = list(N01,N02,N03,P01,P02,P03,P04,P05,P06,P07)
names(NPC) <- load.list

# normalize and identify variable features for each dataset independently
NPC <- lapply(X = NPC, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  #x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 6000)
})


NPC.combined <- merge(N01, 
                      y = c(N02,N03,P01, P02, P03, P04, P05, P06, P07),
                      project = "NPC",
                      merge.data = TRUE)
#save(NPC.combined, file = "D:/LQY/NPC script/PhD/NPC.combined.Human.RData")
#load("D:/LQY/NPC script/new.keepCC/NPC.combined.DS0.RData")


#NPC.combined <- SCTransform(NPC.combined, variable.features.n = 2000)
NPC.combined <- FindVariableFeatures(NPC.combined, selection.method = "vst", nfeatures = 2000)
#NPC.combined@assays$RNA@var.features <- NPC.combined@assays$RNA@var.features[!grepl("^RP[SL]", NPC.combined@assays$RNA@var.features, ignore.case = F)]

NPC.combined <- ScaleData(NPC.combined,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.combined <- RunPCA(NPC.combined)
ElbowPlot(NPC.combined)
dims_parameter <- 10
NPC.combined <- RunUMAP(NPC.combined, dims = 1:dims_parameter)
NPC.combined <- FindNeighbors(NPC.combined, dims = 1:dims_parameter)
NPC.combined <- FindClusters(NPC.combined, resolution = 1)#resoulution越小分的越粗
DimPlot(NPC.combined, reduction = "umap",raster = FALSE,label = TRUE)

DimPlot(NPC.combined, reduction = "umap",raster = FALSE,label = TRUE,group.by = 'seurat_clusters')
DimPlot(NPC.combined, reduction = "umap",group.by = "orig.ident",raster = FALSE)


## Other information
#NPC.combined@meta.data$orig.ident <-  substring(NPC.combined@meta.data$orig.ident,4,6)

# Sequencing Technology information
NPC.combined@meta.data$SeqTech <- NPC.combined@meta.data$orig.ident
Idents(NPC.combined) <- 'SeqTech'
SeqTech <- c("Single-cell", "Single-cell", "Single-nucleus",
             "Single-cell", "Single-cell", "Single-cell",
             "Single-nucleus","Single-nucleus","Single-nucleus","Single-nucleus")
Patient <- c("Normal","Normal","Normal","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor")
names(SeqTech) <- levels(NPC.combined)
NPC.combined <- RenameIdents(NPC.combined, SeqTech)
NPC.combined@meta.data$SeqTech <- NPC.combined@active.ident

## save
# save(NPC.combined, file = "/media/linxiaobuskiy/NPC Project/scRNA data/DS0/step1.NPC.Merge.RData")
save(NPC.combined, file = "G:/scRNA data/DS0/step1.NPC.Merge.RData")

rm(NPC.combined)

####################################################################
#load("G:/scRNA data/DS0/NPC.RNA.10x.RData")
#DimPlot(NPC.combined, reduction = "umap", group.by = "orig.ident",raster=FALSE)

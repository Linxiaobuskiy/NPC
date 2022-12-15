rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)
library(scDblFinder)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

datalist <- c("N01","N02","N03","P01","P02","P03","P04","P05","P06","P07")


for (i in datalist){
  
  print(paste("processing data ",i))
  
  ## load NPC.10X
  tmp.path <- paste("/media/linxiaobuskiy/Storage_4T/DS0/DS0",i,"/outs/filtered_feature_bc_matrix/",sep = "")
  # tmp.path <- paste("F:/",tmp1,"/",tmp1,tmp2,"/outs/filtered_feature_bc_matrix/",sep = "")
  NPC.10x <- Read10X(data.dir = tmp.path)
  NPC.10x <- CreateSeuratObject(NPC.10x, project = i, min.cells = 10, min.features = 200)
  
  NPC.10x$log10GenesPerUMI <- log10(NPC.10x$nFeature_RNA) / log10(NPC.10x$nCount_RNA)
  
  ## Normalization
  print("Normalization")
  QC <- NPC.10x
  QC <- NormalizeData(QC)
  QC <- CellCycleScoring(QC, 
                         g2m.features = g2m.genes, 
                         s.features = s.genes)
  NPC.10x@meta.data <- QC@meta.data
  rm(QC)
  
  ## Remove MT over expressed cells
  print("Remove MT over expressed cells")
  NPC.10x <- PercentageFeatureSet(NPC.10x, pattern = "^MT-", col.name = "percent.mt")
  
  # V2
  #NPC.10x <- subset(NPC.10x, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
  # V3
  # NPC.10x <- subset(NPC.10x, subset = nFeature_RNA > 250 & percent.mt < 20 & nCount_RNA > 500 & log10GenesPerUMI > 0.8)
  NPC.10x <- subset(NPC.10x, subset = nFeature_RNA > 250 & percent.mt < 20 & nCount_RNA > 500 & log10GenesPerUMI > 0.8)
  
  #NPC.10x <- SCTransform(NPC.10x, variable.features.n = 6000, vars.to.regress=c("S.Score", "G2M.Score","percent.mt"))
  NPC.10x <- NormalizeData(NPC.10x, normalization.method = "LogNormalize", scale.factor = 10000)
  # NPC.10x <- FindVariableFeatures(NPC.10x, selection.method = "vst", nfeatures = 6000)
  # all.genes <- rownames(NPC.10x)
  # NPC.10x <- ScaleData(NPC.10x, features = all.genes, vars.to.regress=c("S.Score", "G2M.Score","percent.mt"))
  
  
  ## Doublet identification
  print("Doublet identification")
  tmp <- as.SingleCellExperiment(NPC.10x)
  tmp <- scDblFinder(tmp)
  tmp <- as.Seurat(tmp)
  tmp_doublet <- tmp@meta.data
  rm(tmp)
  
  NPC.10x@meta.data <- tmp_doublet
  #table(P01@meta.data$scDblFinder.class)
  NPC.10x <- subset(NPC.10x, subset = scDblFinder.class == "singlet")
  VlnPlot(NPC.10x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  rm(tmp_doublet)
  
  
  ## save
  path <- paste("/media/linxiaobuskiy/NPC Project/scRNA data/DS0/",i,".10x.RData",sep = "")
  #path <- paste("G:/scRNA data/DS0/",i,".10x.RData",sep = "")
  save(NPC.10x, file = path)
  
  assign(i,NPC.10x)
  rm(NPC.10x)
  
}


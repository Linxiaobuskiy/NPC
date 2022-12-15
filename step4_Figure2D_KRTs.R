rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)
library(scDblFinder)
library(scater)
library(CytoTRACE)
#############################################################
# loading Epithelial cells subtypes annotation

load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse

NPC.Epi <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Epithelial"))
load("G:/scRNA data/DS0/step4.AnnotatedEpithelial.RData")
NPC.Epi@meta.data$Annotation.Epi <- Annotation.Epi
###############################################################
# UMAP
`%notin%` <- Negate(`%in%`)
colors1 <- c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87", "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658", "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398", "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B", "#712820", "#DCC1DD", "#CCE0F5",  "#CCC9E6", "#625D9E", "#68A180", "#3A6963", "#968175")#颜色设置

NPC.tmp <- NPC.Epi

NPC.tmp@meta.data$Annotation.Epi <- factor(NPC.tmp@meta.data$Annotation.Epi, levels = c(
  "scNormal_Ciliated","scNormal_Unciliated",
  "snNormal_Ciliated","snNormal_Unciliated",
  "P01_Unciliated","P02_Unciliated",
  "P03_Unciliated_MT-","P03_Unciliated_MT+",
  "P04_Ciliated_AGBL4+","P04_Ciliated_TSPAN1+","P04_Unciliated_ADAM23+","P04_Unciliated_IG+BPIFA1+","P04_Unciliated_IG-BPIFA1+","P04_Unciliated_IG+",
  "P05_Ciliated_AGBL4+","P05_Ciliated_BPIFA1+","P05_Ciliated_IG+","P05_Ciliated_VWF+","P05_Ciliated_XKR4+","P05_Unciliated_BPIFA1+",
  "P06_Unciliated","P06_Unciliated_ABCA13+","P06_Unciliated_IG+","P06_Unciliated_MT+","P06_Unciliated_ZMAT4+",
  "P07_Unciliated","P07_Unciliated_MT+"
))

NPC.tmp@meta.data$order0 <- as.numeric(as.factor(NPC.tmp@meta.data$Annotation.Epi))
NPC.tmp@meta.data$order0 <- factor(NPC.tmp@meta.data$order0, levels = c(1:27))
NPC.tmp@meta.data$order1 <- paste(as.numeric(NPC.tmp@meta.data$order0),NPC.tmp@meta.data$Annotation.Epi,sep=': ')
NPC.tmp@meta.data$order1 <- factor(NPC.tmp@meta.data$order1, levels = c(
  "1: scNormal_Ciliated","2: scNormal_Unciliated",
  "3: snNormal_Ciliated","4: snNormal_Unciliated",
  "5: P01_Unciliated","6: P02_Unciliated",
  "7: P03_Unciliated_MT-","8: P03_Unciliated_MT+",
  "9: P04_Ciliated_AGBL4+","10: P04_Ciliated_TSPAN1+","11: P04_Unciliated_ADAM23+","12: P04_Unciliated_IG+BPIFA1+","13: P04_Unciliated_IG-BPIFA1+","14: P04_Unciliated_IG+",
  "15: P05_Ciliated_AGBL4+","16: P05_Ciliated_BPIFA1+","17: P05_Ciliated_IG+","18: P05_Ciliated_VWF+","19: P05_Ciliated_XKR4+","20: P05_Unciliated_BPIFA1+",
  "21: P06_Unciliated","22: P06_Unciliated_ABCA13+","23: P06_Unciliated_IG+","24: P06_Unciliated_MT+","25: P06_Unciliated_ZMAT4+",
  "26: P07_Unciliated","27: P07_Unciliated_MT+"
))


NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 5
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
#NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
#NPC.tmp <- FindClusters(NPC.tmp, resolution = 1)#resoulution越小分的越粗
DimPlot(NPC.tmp,group.by = "Annotation.Epi")
DimPlot(NPC.tmp, pt.size = 0.6, repel = T, group.by = "order0", reduction = "umap",label = T, cols= colors1)

p1 <- FeaturePlot(NPC.tmp,features = c("KRT5")) + NoAxes() + NoLegend() + theme(title = element_text(size= 20,color="black",face="plain"))
p2 <- FeaturePlot(NPC.tmp,features = c("KRT19")) + NoAxes() + NoLegend() + theme(title = element_text(size= 20,color="black",face="plain"))

ggsave("G:/manuscript/Figure/p1.pdf", egg::set_panel_size(p1, width=unit(165, "mm"), height=unit(165, "mm")), 
       width = 200, height = 200, units = 'mm', dpi = 300)
ggsave("G:/manuscript/Figure/p2.pdf", egg::set_panel_size(p2, width=unit(165, "mm"), height=unit(165, "mm")), 
       width = 200, height = 200, units = 'mm', dpi = 300)

# save as 8 inch
# resize to 50mm
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


#NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
#NPC.tmp@assays$RNA@var.features <- NPC.tmp@assays$RNA@var.features[!grepl("^RP[SL]", NPC.tmp@assays$RNA@var.features, ignore.case = F)]
#NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))



NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 5
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
#NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
#NPC.tmp <- FindClusters(NPC.tmp, resolution = 1)#resoulution越小分的越粗
DimPlot(NPC.tmp,group.by = "Annotation.Epi")
DimPlot(NPC.tmp, pt.size = 0.6, repel = T, group.by = "order0", reduction = "umap")

p1 <- DimPlot(NPC.tmp, label = F, cols= colors1, pt.size = 0.6, repel = T, group.by = "order0", reduction = "umap")+NoLegend()+ NoAxes()+labs(title = NULL)
p2 <- DimPlot(NPC.tmp, label = T, cols= colors1, pt.size = 0.6, repel = T, group.by = "order1", reduction = "umap")+
  theme(legend.text = element_text(size=5),
        legend.key.height = unit(0.2, "cm"))
p3 <- DimPlot(NPC.tmp, label = T, cols= colors1, pt.size = 0.6, repel = T, group.by = "order0", reduction = "umap")+NoLegend()+ NoAxes()+labs(title = NULL)
# save p1 as 8inch*8inch pdf
# Transfer p1 to semi-png figure
# Adjust p1 to 50mm*50mm
# Add labels and legend from p2/3 to p1
# Adjust word to 5pt, label gap 2mm, square-word left-left 2.5mm
# Add Axes


###############################################################
# violin

a <- c("CAPS","RSPH1", # Ciliated
       "COL1A1","COL1A2", #Fibroblasts
       "CDH1","EPCAM","KRT19","KRT5","VIM", # EMT
       "HIF1A","EPAS1","HIF3A", # Hypoxia
       "SREBF1","ACACA", # Lipid metabolism
       "GPI","LDHA","LDHB", # Glucose Metabolism
       "SOX2","CD44","CD163","ALDH1A1", #CSC
       "LYZ","ULK1","WIPI1","WIPI2","ATG5","MAP1LC3A" #Autophagy
)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
E_EMT <- c("COL1A1","KRT5","KRT19","CAPS","RSPH1","CDH1","LAMA3","MUC1","SDC1")
M_EMT <- c("ACTA2","FN1","CDH2","S100A4","VIM","MMP2","TWIST1")
#M_EMT <- c("ACTA2","FN1","CDH2","S100A4","SNAI1","SNAI2","VIM","ITGB6","MMP2","MMP3","MMP9","SOX10","TWIST1")

Hypoxia <- c("HIF1A","EPAS1","HIF3A",
             "SLC2A1","PDK1","HK1","HK2","HK3",
             "GAPDH","ALDOA","ENO1","PGK1",
             "PFKM","PFKP","PFKL","LDHA") # Aerobic Glycolysis
# Merge metabolism
Hypoxia <- c("HIF1A","EPAS1","HIF3A",
             "SLC2A1","PDK1","HK1","HK2","HK3",
             "GAPDH","ALDOA","ENO1","PGK1",
             "PFKM","PFKP","PFKL","LDHA","SLC27A3","SREBF1","ACLY","ACACA","ACACB","FASN",
             "SCD",
             "NANOG","MYC","CD36","CPT1A","CPT1B")

Lipid_Metabolism <- c("SLC27A3","SREBF1","ACLY","ACACA","ACACB","FASN",
                      "SCD",
                      "NANOG","MYC","CD36","CPT1A","CPT1B")

Mitochondrion <- grep(pattern = "^MT-", x = rownames(NPC.Epi@assays[["RNA"]]), value = TRUE)

CSC <- c("CD44","PLAUR","THY1","PROM1","ALCAM","EPCAM","ALDH2","NANOG","POU5F1",
         "CD24","ITGB1","ITGA6","ITGB3","CD70","CXCR4","LGR5","PROCR","BMI1","NOTCH1","SOX2",
         "LINGO2","LETM1","MSI2","AFP","SALL4","CLEC12A","HAVCR2",
         "IL2RA","DPP4","CD33","KIT","IL3RA","IL1RAP")
# https://www.frontiersin.org/articles/10.3389/fimmu.2020.01280/full
#CSC <- c("CD44","PLAUR","THY1","PROM1","ALCAM","EPCAM","ALDH2","NANOG","POU5F1")


Autophagy <- c("LYZ","ATG5","ATG7","ATG12","MAP1LC3A","MAP1LC3B","BECN1","HIF1A","FIS1")
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6460398/

MHC <- grep(pattern = "^HLA-", x = rownames(NPC.Epi@assays[["RNA"]]), value = TRUE)

NPC.tt <- NPC.tmp
for (i in 1:2){
  if (i==1) {NPC.tmp <- subset(NPC.tt, subset = orig.ident %in% c("N01","N02","P01","P02","P03"))}
  if (i==2) {NPC.tmp <- subset(NPC.tt, subset = orig.ident %in% c("N03","P04","P05","P06","P07"))}
  
  NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
  NPC.tmp@assays$RNA@var.features <- NPC.tmp@assays$RNA@var.features[!grepl("^RP[SL]", NPC.tmp@assays$RNA@var.features, ignore.case = F)]
  NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))
  
  ###############################################################
  # E-EMT
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(E_EMT),ctrl = 100,name = "score.E_EMT")
  FeaturePlot(NPC.tmp, features="score.E_EMT1")
  # M-EMT
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(M_EMT),ctrl = 100,name = "score.M_EMT")
  FeaturePlot(NPC.tmp, features="score.M_EMT1")
  
  NPC.tmp@meta.data$score.M_EMT1 <- range01(NPC.tmp@meta.data$score.M_EMT1)
  NPC.tmp@meta.data$score.E_EMT1 <- range01(NPC.tmp@meta.data$score.E_EMT1)
  NPC.tmp@meta.data$score.EMT <- (NPC.tmp@meta.data$score.E_EMT1-NPC.tmp@meta.data$score.M_EMT1)/sqrt(2)
  
  ###############################################################
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(Hypoxia),ctrl = 100,name = "score.Hypoxia")
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(Lipid_Metabolism),ctrl = 100,name = "score.Lipid_Metabolism")
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(Mitochondrion),ctrl = 100,name = "score.Mitochondrion")
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(CSC),ctrl = 100,name = "score.CSC")
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(Autophagy),ctrl = 100,name = "score.Autophagy")
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(MHC),ctrl = 100,name = "score.MHC")
  
  NPC.tmp@meta.data$stemness<- NPC.tmp@meta.data$CytoTRACE
  NPC.tmp@meta.data$E.EMT <- NPC.tmp@meta.data$score.E_EMT
  NPC.tmp@meta.data$M.EMT <- NPC.tmp@meta.data$score.M_EMT
  NPC.tmp@meta.data$EMT <- NPC.tmp@meta.data$score.EMT
  NPC.tmp@meta.data$Hypoxia <- NPC.tmp@meta.data$score.Hypoxia1
  NPC.tmp@meta.data$Lipid_Metabolism <- NPC.tmp@meta.data$score.Lipid_Metabolism1
  NPC.tmp@meta.data$Mitochondrion <- NPC.tmp@meta.data$score.Mitochondrion1
  NPC.tmp@meta.data$CSC <- NPC.tmp@meta.data$score.CSC1
  NPC.tmp@meta.data$Autophagy <- NPC.tmp@meta.data$score.Autophagy1
  NPC.tmp@meta.data$MHC <- NPC.tmp@meta.data$score.MHC1
  if (i==1) {NPC.sc <- NPC.tmp}
  if (i==2) {NPC.sn <- NPC.tmp}
}

##############################################################
# Figure
#install.packages("remotes")
#remotes::install_github("lyc-1995/MySeuratWrappers")#通过链接安装包
library(MySeuratWrappers)


markers <- c("EMT","Hypoxia","Mitochondrion","CSC","Autophagy","MHC")

##############################################################
# vertical
p.sc.vertical <- VlnPlot(NPC.sc, features = markers,
                         stacked=T,
                         pt.size=0,
                         cols = colors1[c(1,2,5,6,7,8)],#颜色
                         #direction = "vertical",
                         #direction = "horizontal",
                         group.by = "order0",
                         x.lab = "", y.lab = "") +
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(size=7, face="bold",color = "black"),#分面标签样式
        strip.background.y = element_rect(colour="white", fill="white",size = 0)
  )#不显示坐标刻度


p.sn.vertical <- VlnPlot(NPC.sn, features = markers,
                         stacked=T,
                         pt.size=0,
                         cols = colors1[c(3,4,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)],#颜色
                         #direction = "vertical",
                         #direction = "horizontal",
                         group.by = "order0",
                         x.lab = "", y.lab = "") +
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(size=7, face="bold",color = "black"),#分面标签样式
        strip.background.y = element_rect(colour="white", fill="white",size = 0)
  )#不显示坐标刻度

p.sc.vertical+p.sn.vertical


#gridExtra::grid.arrange(egg::set_panel_size(p=p.sc, width=unit(9, "cm"), height=unit(5, "cm")))

#6 21

ggsave("G:/manuscript/Figure/sc.pdf", egg::set_panel_size(p.sc.vertical, width=unit(0.6, "in"), height=unit(0.2, "in")), 
       width = 7, height = 9, units = 'in', dpi = 300)
ggsave("G:/manuscript/Figure/sn.pdf", egg::set_panel_size(p.sn.vertical, width=unit(2.1, "in"), height=unit(0.2, "in")), 
       width = 22, height = 9, units = 'in', dpi = 300)

##############################################################
# HLA-DR
VlnPlot(NPC.sn,features = c("HLA-DRA"),group.by = "order1")
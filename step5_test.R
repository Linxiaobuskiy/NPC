tmp <- NPC.tmp
tmp <- subset(NPC.tmp, subset = seurat_clusters %in% c(2,29))
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'seurat_clusters'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


tmp <- NPC.tmp
tmp <- subset(NPC.tmp, subset = order0 %in% c(8,26))
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'order0'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DimPlot(tmp,group.by = "Annotation.Epi")



tmp <- NPC.tmp
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'Annotation.Fib'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)




p1 <-FeaturePlot(NPC.tmp, features = c(
  "AGBL4","PDGFRA",
  "BPIFA1",#mCAF matrix
  "CAPS",#vCAF vascular
  "XKR4", # dCAFs developmental
  "PTPRC",#cCAF cycling
  "VWF" #iCAF
),min.cutoff = 0.5,raster=FALSE)
p2 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE)
p1+p2




DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

FeaturePlot(NPC.tmp, features = c(
  "LYZ"
),min.cutoff = 2,raster=FALSE)

VlnPlot(NPC.tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
NPC.tmp <- subset(NPC.tmp, subset = nFeature_RNA>500)


FeaturePlot(NPC.tmp, features = c(
  "ACTA2","PDGFRA",
  "FBLN1",#mCAF matrix
  "NID2","PECAM1",#vCAF vascular
  "SCRG1", # dCAFs developmental
  #"",#cCAF cycling
  "PDPN","IL6" #iCAF
),min.cutoff = 0.5,raster=FALSE)



FeaturePlot(NPC.tmp, features = c(
  "CAPS"
),min.cutoff = 1,raster=FALSE)


##########################################################
##########################  Go  ##########################
##########################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

tmp <- NPC.tmp
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'Annotation.CNV'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
up.genes.enrich <- subset(markers,subset = cluster %in% c("non-CNV"))
up.genes.enrich <- rownames(up.genes.enrich)
up.genes.enrich <- up.genes.enrich[!grepl("^RP[SL]", up.genes.enrich, ignore.case = F)]
################################################################################
# GO
up.genes.enrich <- bitr(up.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
#up.genes.enrich <- pull(up.genes.enrich,ENTREZID)  

up.genes.go <- enrichGO(gene          = up.genes.enrich$ENTREZID,
                        #universe     = row.names(dge.celltype),
                        OrgDb         = 'org.Hs.eg.db',
                        #keyType       = 'SYMBOL',
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.05)
#head(up.genes.go)
barplot(up.genes.go, showCategory=20)

############################################################
##########################  Score ##########################
############################################################

ImmunoSuppression <- c("TGFB1","IL6","IL11","LIF","CSF2","CXCL1","CXCL2","CXCL8","CXCL9","CXCL10","CXCL12","CCL2","CCL8","PTK2")
Mitochondrion <- grep(pattern = "^MT-", x = rownames(NPC.tmp@assays[["RNA"]]), value = TRUE)
MHC <- grep(pattern = "^HLA-", x = rownames(NPC.tmp@assays[["RNA"]]), value = TRUE)
MHC_I <- c("HLA-A","HLA-B","HLA-C")
MHC_II <- c("HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1")
ComplementRegulation = c("CFD", "C1QA","C1QB","C1QC")
Collagen <- grep(pattern = "^COL.A.", x = rownames(NPC.tmp@assays[["RNA"]]), value = TRUE)
# https://www.nature.com/articles/s41417-021-00318-4/tables/1
CAFs <- c(
  # Cytokines
  "IL6","IL11","IL32","CXCL12","CXCL14","CXCL16","CCL2","CXCL8","CCL5","CXCL8","GDF15","CSF2","KDM1A",
  # Growth Factors
  "TGFB1","VEGFA","VEGFC","HGF","EGF","CCN2","CNPY1","CNPY2","CNPY3","CNPY4",
  # ECM components
  "TNC","COL5A2","COL6A3","COL1A2","ASPN","POSTN","FN1","LOXL1","MMP9","ROCK1","ROCK2",
  # receptors and other membrane-bound proteins
  "ANXA11","PDGFRB","ITGB4","VCAM1","DDR2","FGFR1","FGFR2","FGFR3","FGFR4","PDPN","FAP","CAV1",
  # Cytoskeleton components and CPPs
  "VIM","ACTA2","S100A4","NTN1","GFPT2",
  # extracellular vesicles
  "EXOSC1","EXOSC2","EXOSC3","EXOSC4","EXOSC5","EXOSC6","EXOSC7","EXOSC8","EXOSC9","EXOSC10","POSTN","ANXA6","SNAI1","SNAI2"
)


NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(ImmunoSuppression),ctrl = 100,name = "score.ImmunoSuppression")
NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(Mitochondrion),ctrl = 100,name = "score.Mitochondrion")
NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(MHC),ctrl = 100,name = "score.MHC")
NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(ComplementRegulation),ctrl = 100,name = "score.ComplementRegulation")
NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(Collagen),ctrl = 100,name = "score.Collagen")

VlnPlot(NPC.tmp, features="score.ImmunoSuppression1")
VlnPlot(NPC.tmp, features="score.Mitochondrion1")
VlnPlot(NPC.tmp, features="score.MHC1")
VlnPlot(NPC.tmp, features="score.ComplementRegulation1")
VlnPlot(NPC.tmp, features="score.Collagen1")
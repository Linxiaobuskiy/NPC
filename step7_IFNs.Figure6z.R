rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)
library(scDblFinder)
library(scater)
library(patchwork)
library(stringr)

load("G:/scRNA data/DS0/step7.AnnotationImmune.RData")
NPC.tmp <- subset(NPC.tmp, subset = Annotation.Immune2 %in% c("Memory B","Memory B_FCRL4+",
                                                              "Naive B","NK_CD8-","T_CD4+","T_CD8+",
                                                              "T_MKI67+","T-reg"))
NPC.tmp@meta.data$source <- substr(NPC.tmp@meta.data$orig.ident,1,1)
NPC.tmp@meta.data$Annotation.Immune3 <- paste(NPC.tmp@meta.data$source,NPC.tmp@meta.data$Annotation.Immune2,sep="_")
celltypes <- unique(NPC.tmp@meta.data$Annotation.Immune2) %>% as.character %>% str_sort

## genes
# antiviral innate immune response GO:0140374
Anti_Viral <- c("AKAP1","CLPB","DHX15",
                "IFIT2","MARCHF2","MBL2",
                "NLRP1","NLRP6","NMB","NMBR",
                "OAS1","PHB2","RIGI","RNF135",
                "SENP7","TRAF3IP3",
                "ZDHHC1","ZDHHC11")

# positive regulation of immune response to tumor cell GO:0002839
Anti_Tumor <- c("CD24","CD40LG","CD160","CD226","CD274",
                "CRTAM","FBXO38","GSDME","HRG","HSPD1",
                "IL12A","IL12B","KLRK1","MR1",
                "NECTIN2","PVR","RAE1","SLC22A13","ULBP1","XCL1")


# positive regulation of interferon-alpha production, GO:0032727
# http://www.informatics.jax.org/go/term/GO:0032727
IFN_alpha <- c("CHUK","DDX3X","DHX9","DHX36","FLT3","HMGB1",
               "HSPD1","IFIH1","IRF3","IRF7","MAVS","MMP12",
               "NMB","NMBR","PTPN22","RIGI","SETD2","STAT1",
               "TBK1","TLR3","TLR4","TLR7","TLR8","TLR9",
               "TRAF3IP3","ZC3HAV1")

# positive regulation of type II interferon production, GO:0032729
# http://www.informatics.jax.org/go/term/GO:0032729
IFN_gamma <- c("ABL1","ABL2","APP","ARID5A","BCL3","CARMIL2",
               "CCR2","CD1D","CD2","CD3E","CD14","CD27","CD40LG",
               "CD160","CD226","CD244","CD276","CEBPG","CRTAM",
               "CYRIB","F2RL1","FADD","FLT3","FZD5","HAVCR2",
               "HRAS","HSPD1","IFNAR1","IL1B","IL1R1","IL2",
               "IL12A","IL12B","IL12RB1","IL12RB2","IL18","IL18R1",
               "IL21","IL23A","IL23R","IL27","IL27RA","IRF8","ISG15",
               "ISL1","JAK2","KLRK1","LTA","MIR155","NRAS",
               "PDE4B","PDE4D","PTPN22","PYCARD","RASGRP1","RIPK2",
               "RUNX1","SASH3","SCRIB","SLAMF1","SLAMF6","SLC7A5",
               "SLC11A1","TICAM2","TLR3","TLR4","TLR7","TLR8","TLR9",
               "TNF","TNFRSF13C","TNFSF4","TNFSF9","TRIM27","TXK",
               "TYK2","ULBP1","WNT5A","ZFPM1","ZP3")

tmp <- subset(NPC.tmp, subset = SeqTech %in% "Single-cell")
#tmp <- subset(NPC.tmp, subset = SeqTech %in% "Single-nucleus")

# single cell
tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
tmp <- ScaleData(tmp,  vars.to.regress=c("S.Score", "G2M.Score"))
tmp <- AddModuleScore(object = tmp,features = list(Anti_Viral),ctrl = 100,name = "score.Anti_Viral")
tmp <- AddModuleScore(object = tmp,features = list(Anti_Tumor),ctrl = 100,name = "score.Anti_Tumor")
tmp <- AddModuleScore(object = tmp,features = list(IFN_alpha),ctrl = 100,name = "score.IFN_alpha")
tmp <- AddModuleScore(object = tmp,features = list(IFN_gamma),ctrl = 100,name = "score.IFN_gamma")

test.size <- 5
p1<-VlnPlot(tmp,features = c("score.IFN_alpha1"),group.by = "Annotation.Immune2",pt.size=0, split.by = "source",cols=c("#1BB3B7","#EA7369"))+
  theme(axis.text.x = element_text(size=test.size,color="black",face="plain"),
        axis.text.y = element_text(size=test.size,color="black",face="plain"))+
  theme(legend.text = element_text(size = test.size),
        legend.title = element_text(size = test.size),
        plot.title = element_text(size = test.size))
p2<-VlnPlot(tmp,features = c("score.IFN_gamma1"),group.by = "Annotation.Immune2",pt.size=0, split.by = "source",cols=c("#1BB3B7","#EA7369"))+
  theme(axis.text.x = element_text(size=test.size,color="black",face="plain"),
        axis.text.y = element_text(size=test.size,color="black",face="plain"))+
  theme(legend.text = element_text(size = test.size),
        legend.title = element_text(size = test.size),
        plot.title = element_text(size = test.size))

p1+p2

ggsave("G:/manuscript/Figure/IFN1.pdf", egg::set_panel_size(p1, width=unit(50, "mm"), height=unit(30, "mm")), 
       width = 100, height = 100, units = 'mm', dpi = 300)
ggsave("G:/manuscript/Figure/IFN2.pdf", egg::set_panel_size(p2, width=unit(50, "mm"), height=unit(30, "mm")), 
       width = 100, height = 100, units = 'mm', dpi = 300)



# wilcox.test
datalist <- c("Memory B","Memory B_FCRL4+","Naive B","NK_CD8-","T_CD4+","T_CD8+","T_MKI67+","T-reg")

p.result <- data.frame(
  sample = datalist,
  p.1=rep(0,length(datalist)),
  p.2=rep(0,length(datalist)),
  mean.1.n=rep(0,length(datalist)),
  mean.2.n=rep(0,length(datalist)),
  mean.1.p=rep(0,length(datalist)),
  mean.2.p=rep(0,length(datalist))
)


for (i in datalist){
  
  xx <- subset(tmp, subset = source %in% c("N"))
  xx <- subset(xx, subset = Annotation.Immune2 %in% c(i))
  x.1 <- xx$score.IFN_alpha1
  x.2 <- xx$score.IFN_gamma1
  yy <- subset(tmp, subset = source %in% c("P"))
  yy <- subset(yy, subset = Annotation.Immune2 %in% c(i))
  y.1 <- yy$score.IFN_alpha1
  y.2 <- yy$score.IFN_gamma1
  tt.1 <- wilcox.test(x = x.1 , y = y.1, alternative = "two.sided")
  tt.2 <- wilcox.test(x = x.2 , y = y.2, alternative = "two.sided")
  
  p.result$p.1[match(i,datalist)] <- tt.1$p.value
  p.result$p.2[match(i,datalist)] <- tt.2$p.value
  p.result$mean.1.n[match(i,datalist)] <- mean(x.1)
  p.result$mean.2.n[match(i,datalist)] <- mean(x.2)
  p.result$mean.1.p[match(i,datalist)] <- mean(y.1)
  p.result$mean.2.p[match(i,datalist)] <- mean(y.2)
}

# * 0.05
# ** 0.01
# *** 0.001
# P01 M1 * M2 





a <- data.frame(score.IFN_alpha = tmp$score.IFN_alpha1, score.IFN_gamma = tmp$score.IFN_gamma1, group = tmp$source)
ggplot(a,aes(score.IFN_alpha,score.IFN_gamma))+ geom_point(aes(colour = factor(group)))
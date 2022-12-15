rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)
library(monocle3)
library(scater)
library(patchwork)
`%notin%` <- Negate(`%in%`)

load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
#load("/media/linxiaobuskiy/NPC Project/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
#load("/media/linxiaobuskiy/NPC Project/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse

#############################################################
# loading fibroblasts CNV annotation
NPC.Fib <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Fibroblast"))
load("G:/scRNA data/DS0/step5.FibroblastAnnotationCNV.RData")
#load("/media/linxiaobuskiy/NPC Project/scRNA data/DS0/step5.FibroblastAnnotationCNV.RData")
NPC.Fib@meta.data$Annotation.CNV <- Annotation.Fib

#############################################################
# loading Epithelial cells subtypes annotation
NPC.Epi <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Epithelial"))
load("G:/scRNA data/DS0/step4.AnnotatedEpithelial.RData")
#load("/media/linxiaobuskiy/NPC Project/scRNA data/DS0/step4.AnnotatedEpithelial.RData")
NPC.Epi@meta.data$Annotation.Epi <- Annotation.Epi


#############################################################
##########################   PXX   ##########################
#############################################################

for (tt in c("P04","P05","P06","P07")){
  NPC.tmp <- subset(NPC.combined, subset = (orig.ident %in% c(tt))&
                      (Annotation.Coarse %in% c("Fibroblast","Epithelial")))
  
  
  NPC.tmp@meta.data$Annotation.test <- "TEST"
  
  # Add fib annotation
  tmp.a <- NPC.Fib
  tmp.b <- NPC.tmp
  
  index1 <- match(colnames(tmp.a),colnames(tmp.b))
  tmp.dataframe <- data.frame(col1 <- index1, col2 <- 1:length(index1))
  tmp.dataframe <- tmp.dataframe[!is.na(tmp.dataframe$col1),]
  tmp.b@meta.data$Annotation.test[tmp.dataframe$col1] <- paste("Fibroblast", tmp.a@meta.data$Annotation.CNV[tmp.dataframe$col2],sep="_")
  
  NPC.tmp <- tmp.b
  
  # Add epi annotation
  tmp.a <- NPC.Epi
  tmp.b <- NPC.tmp
  
  index1 <- match(colnames(tmp.a),colnames(tmp.b))
  tmp.dataframe <- data.frame(col1 <- index1, col2 <- 1:length(index1))
  tmp.dataframe <- tmp.dataframe[!is.na(tmp.dataframe$col1),]
  tmp.b@meta.data$Annotation.test[tmp.dataframe$col1] <- tmp.a@meta.data$Annotation.Epi[tmp.dataframe$col2]
  
  NPC.tmp <- tmp.b
  
  
  # order
  level.P = list()
  level.P$P04 <- c("Fibroblast_non-CNV","Fibroblast_CNV",
                   "P04_Ciliated_AGBL4+","P04_Ciliated_TSPAN1+","P04_Unciliated_ADAM23+","P04_Unciliated_IG+BPIFA1+","P04_Unciliated_IG-BPIFA1+","P04_Unciliated_IG+")
  level.P$P05 <- c("Fibroblast_non-CNV","Fibroblast_CNV","P05_Ciliated_AGBL4+","P05_Ciliated_BPIFA1+","P05_Ciliated_IG+","P05_Ciliated_VWF+","P05_Ciliated_XKR4+","P05_Unciliated_BPIFA1+")
  level.P$P06 <- c("Fibroblast_non-CNV","Fibroblast_CNV","P06_Unciliated","P06_Unciliated_ABCA13+","P06_Unciliated_IG+","P06_Unciliated_MT+","P06_Unciliated_ZMAT4+")
  level.P$P07 <- c("Fibroblast_non-CNV","Fibroblast_CNV","P07_Unciliated","P07_Unciliated_MT+")
  
  NPC.tmp@meta.data$Annotation.test <- factor(NPC.tmp@meta.data$Annotation.test, levels = level.P[[tt]])
  
  
  
  
  # score
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  E_EMT <- c("KRT5","KRT19","CAPS","RSPH1","CDH1","LAMA3","MUC1","SDC1")
  M_EMT <- c("ACTA2","FN1","CDH2","S100A4","VIM","MMP2","TWIST1")
  MHC <- grep(pattern = "^HLA-", x = rownames(NPC.combined@assays[["RNA"]]), value = TRUE)
  MHC_I <- c("HLA-A","HLA-B","HLA-C","HLA-E", "HLA-F", "HLA-G")
  MHC_II <- c("HLA-DRB5", "HLA-DQA2", "HLA-DQB2", "HLA-DOA", "HLA-DOB", "HLA-DMA", "HLA-DMB")
  CSC <- c("CD44","PLAUR","THY1","PROM1","ALCAM","EPCAM","ALDH2","NANOG","POU5F1",
           "CD24","ITGB1","ITGA6","ITGB3","CD70","CXCR4","LGR5","PROCR","BMI1","NOTCH1","SOX2",
           "LINGO2","LETM1","MSI2","AFP","SALL4","CLEC12A","HAVCR2",
           "IL2RA","DPP4","CD33","KIT","IL3RA","IL1RAP")
  
  #M_EMT <- c("ACTA2","FN1","CDH2","S100A4","SNAI1","SNAI2","VIM","ITGB6","MMP2","MMP3","MMP9","SOX10","TWIST1")
  
  
  # E-EMT
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(E_EMT),ctrl = 100,name = "score.E_EMT")
  FeaturePlot(NPC.tmp, features="score.E_EMT1")
  # M-EMT
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(M_EMT),ctrl = 100,name = "score.M_EMT")
  FeaturePlot(NPC.tmp, features="score.M_EMT1")
  # MHC
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(MHC),ctrl = 100,name = "score.MHC")
  # MHC_I
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(MHC_I),ctrl = 100,name = "score.MHC_I")
  # MHC_II
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(MHC_II),ctrl = 100,name = "score.MHC_II")
  # CSC
  NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(CSC),ctrl = 100,name = "score.CSC")
  
  
  NPC.tmp@meta.data$score.M_EMT <- range01(NPC.tmp@meta.data$score.M_EMT1)
  NPC.tmp@meta.data$score.E_EMT <- range01(NPC.tmp@meta.data$score.E_EMT1)
  NPC.tmp@meta.data$score.EMT <- (NPC.tmp@meta.data$score.E_EMT-NPC.tmp@meta.data$score.M_EMT)/sqrt(2)
  NPC.tmp@meta.data$score.EMT <- range01(NPC.tmp@meta.data$score.EMT)
  NPC.tmp@meta.data$score.MHC <- range01(NPC.tmp@meta.data$score.MHC1)
  NPC.tmp@meta.data$score.MHC_I <- range01(NPC.tmp@meta.data$score.MHC_I1)
  NPC.tmp@meta.data$score.MHC_II <- range01(NPC.tmp@meta.data$score.MHC_II1)
  NPC.tmp@meta.data$score.CSC <- range01(NPC.tmp@meta.data$score.CSC1)
  # resize factor
  t <- 5
  
  NPC.tmp@meta.data$Annotation.label <- NPC.tmp@meta.data$Annotation.test
  NPC.tmp@meta.data$Annotation.test <- as.numeric(as.factor(NPC.tmp@meta.data$Annotation.test))
  
  
  p1 <- VlnPlot(NPC.tmp,features = c("score.EMT"),group.by = "Annotation.test",pt.size = 0)+labs(title="EMT.Score",x="")+NoLegend()+
    theme(plot.title = element_text(size=7*t),
          axis.text.x = element_text(size=5*t,color="black",face="plain"),
          axis.text.y = element_text(size=5*t,color="black",face="plain"),
          legend.text = element_text(size=5*t,color="black",face="plain"))
  p2 <- VlnPlot(NPC.tmp,features = c("score.M_EMT"),group.by = "Annotation.test",pt.size = 0)+labs(title="M.EMT.Score",x="")+NoLegend()+
    theme(plot.title = element_text(size=7*t),
          axis.text.x = element_text(size=5*t,color="black",face="plain"),
          axis.text.y = element_text(size=5*t,color="black",face="plain"),
          legend.text = element_text(size=5*t,color="black",face="plain"))
  p3 <- VlnPlot(NPC.tmp,features = c("score.E_EMT"),group.by = "Annotation.test",pt.size = 0)+labs(title="E.EMT.Score",x="")+NoLegend()+
    theme(plot.title = element_text(size=7*t),
          axis.text.x = element_text(size=5*t,color="black",face="plain"),
          axis.text.y = element_text(size=5*t,color="black",face="plain"),
          legend.text = element_text(size=5*t,color="black",face="plain"))
  p4 <- VlnPlot(NPC.tmp,features = c("score.MHC"),group.by = "Annotation.test",pt.size = 0)+labs(title="MHC.Score",x="")+NoLegend()+
    theme(plot.title = element_text(size=7*t),
          axis.text.x = element_text(size=5*t,color="black",face="plain"),
          axis.text.y = element_text(size=5*t,color="black",face="plain"),
          legend.text = element_text(size=5*t,color="black",face="plain"))
  p5 <- VlnPlot(NPC.tmp,features = c("score.MHC_I"),group.by = "Annotation.test",pt.size = 0)+labs(title="MHC_I.Score",x="")+NoLegend()+
    theme(plot.title = element_text(size=7*t),
          axis.text.x = element_text(size=5*t,color="black",face="plain"),
          axis.text.y = element_text(size=5*t,color="black",face="plain"),
          legend.text = element_text(size=5*t,color="black",face="plain"))
  p6 <- VlnPlot(NPC.tmp,features = c("score.MHC_II"),group.by = "Annotation.test",pt.size = 0)+labs(title="MHC_II.Score",x="")+NoLegend()+
    theme(plot.title = element_text(size=7*t),
          axis.text.x = element_text(size=5*t,color="black",face="plain"),
          axis.text.y = element_text(size=5*t,color="black",face="plain"),
          legend.text = element_text(size=5*t,color="black",face="plain"))
  p7 <- VlnPlot(NPC.tmp,features = c("score.CSC"),group.by = "Annotation.test",pt.size = 0)+labs(title="CSC.Score",x="")+NoLegend()+
    theme(plot.title = element_text(size=7*t),
          axis.text.x = element_text(size=5*t,color="black",face="plain"),
          axis.text.y = element_text(size=5*t,color="black",face="plain"),
          legend.text = element_text(size=5*t,color="black",face="plain"))
  p8 <- VlnPlot(NPC.tmp,features = c("score.MHC_II"),group.by = "Annotation.label",pt.size = 0)+labs(title="MHC_II.Score",x="")+
    theme(plot.title = element_text(size=7*t),
          axis.text.x = element_text(size=5*t,color="black",face="plain"),
          axis.text.y = element_text(size=5*t,color="black",face="plain"),
          legend.text = element_text(size=5*t,color="black",face="plain"))+
    theme(legend.key.size = unit(2*t, 'mm'))
  
  
  
  
  
  
  path1 <- paste("G:/manuscript/Figure/",tt,"p1.pdf",sep="")
  path2 <- paste("G:/manuscript/Figure/",tt,"p2.pdf",sep="")
  path3 <- paste("G:/manuscript/Figure/",tt,"p3.pdf",sep="")
  path4 <- paste("G:/manuscript/Figure/",tt,"p4.pdf",sep="")
  path5 <- paste("G:/manuscript/Figure/",tt,"p5.pdf",sep="")
  path6 <- paste("G:/manuscript/Figure/",tt,"p6.pdf",sep="")
  path7 <- paste("G:/manuscript/Figure/",tt,"p7.pdf",sep="")
  path8 <- paste("G:/manuscript/Figure/",tt,"p8.pdf",sep="")
  len = length(table(NPC.tmp@meta.data$Annotation.test))*4
  ggsave(path1, egg::set_panel_size(p1, width=unit(len*t, "mm"), height=unit(20*t, "mm")), width = 100*t, height = 100*t, units = 'mm', dpi = 300)
  ggsave(path2, egg::set_panel_size(p2, width=unit(len*t, "mm"), height=unit(20*t, "mm")), width = 100*t, height = 100*t, units = 'mm', dpi = 300)
  ggsave(path3, egg::set_panel_size(p3, width=unit(len*t, "mm"), height=unit(20*t, "mm")), width = 100*t, height = 100*t, units = 'mm', dpi = 300)
  ggsave(path4, egg::set_panel_size(p4, width=unit(len*t, "mm"), height=unit(20*t, "mm")), width = 100*t, height = 100*t, units = 'mm', dpi = 300)
  ggsave(path5, egg::set_panel_size(p5, width=unit(len*t, "mm"), height=unit(20*t, "mm")), width = 100*t, height = 100*t, units = 'mm', dpi = 300)
  ggsave(path6, egg::set_panel_size(p6, width=unit(len*t, "mm"), height=unit(20*t, "mm")), width = 100*t, height = 100*t, units = 'mm', dpi = 300)
  ggsave(path7, egg::set_panel_size(p7, width=unit(len*t, "mm"), height=unit(20*t, "mm")), width = 100*t, height = 100*t, units = 'mm', dpi = 300)
  ggsave(path8, egg::set_panel_size(p8, width=unit(len*t, "mm"), height=unit(20*t, "mm")), width = 100*t, height = 100*t, units = 'mm', dpi = 300)
  
  #resize to 1/t
}



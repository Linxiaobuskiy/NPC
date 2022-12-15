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
Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.1")
library(infercnv)

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
## infercnv reference
Reference <- subset(NPC.combined, subset = orig.ident %in% c("N03"))
Reference <- subset(Reference, subset = Annotation.Coarse %in% c("Epithelial","Fibroblast"))
Reference@meta.data$infercnv <- paste("N03",Reference@meta.data$Annotation.Coarse,sep="_")


#############################################################
##########################   PXX   ##########################
#############################################################
for (tt in c("P05")){
  #for (tt in c("P04","P05","P06","P07")){
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
  
  ## umap tmp
  NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
  NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))
  
  NPC.tmp <- RunPCA(NPC.tmp)
  ElbowPlot(NPC.tmp)
  dims_parameter <- 10
  NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
  NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
  NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.2)#resoulution越小分的越粗
  DimPlot(NPC.tmp, reduction = "umap",label = TRUE)
  DimPlot(NPC.tmp, reduction = "umap",label = TRUE, group.by = "Annotation.test")
  
  
  ## infercnv
  NPC.tmp@meta.data$infercnv <- NPC.tmp@meta.data$Annotation.test
  NPC.test <- merge(Reference, y = c(NPC.tmp),merge.data = TRUE)
  
  Ref <- c("N03_Epithelial","N03_Fibroblast")
  
  counts_matrix = as.data.frame(NPC.test@assays$RNA@counts)
  tmp_ann <- colnames(counts_matrix)
  file_ann <- matrix(nrow = length(tmp_ann),ncol=2)
  file_ann[1:length(tmp_ann),1] <- tmp_ann
  file_ann[1:length(tmp_ann),2] <- NPC.test@meta.data$infercnv
  
  ##############################################################
  # infercnv
  path.annotation = paste("G:/scRNA data/DS0/infercnv/",tt,"_Epi+Fib/annotations.txt",sep="")
  write.table(file_ann, file= path.annotation, row.names=FALSE, col.names=FALSE, quote =FALSE, sep ="\t")
  path.save = paste("G:/scRNA data/DS0/infercnv/",tt,"_Epi+Fib/counts_matrix",sep="")
  save(counts_matrix, file = path.save)
  
  #path.load = path.save
  #load(path.load)
  path.out <- paste("G:/scRNA data/DS0/infercnv/",tt,"_Epi+Fib",sep="")
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                      annotations_file = path.annotation,
                                      delim="\t",
                                      gene_order_file = "G:/scRNA data/DS0/infercnv/gencode_v19_gene_pos.txt",
                                      #gene_order_file = system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                      #ref_group_names=c("N01","N02","N03")
                                      ref_group_names= Ref,
                                      chr_exclude = c("chrY","chrM")
  )
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               denoise=FALSE,
                               num_threads=20,
                               HMM=FALSE,
                               cluster_by_groups=TRUE,
                               output_format = 'pdf',
                               out_dir=path.out)
  
}

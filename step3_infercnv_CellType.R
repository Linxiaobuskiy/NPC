rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)
library(scDblFinder)

Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.1")
library(infercnv)

load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse

###############################################################
########################## Reference ##########################
###############################################################
Reference.sc <- subset(NPC.combined, subset = orig.ident %in% c("N01","N02"))
#Reference.sc <- subset(Reference.sc, subset = CellType %in% c("Epithelial"))
Reference.sc@meta.data$infercnv <- Reference.sc@meta.data$orig.ident

Reference.sn <- subset(NPC.combined, subset = orig.ident %in% c("N03"))
#Reference.sn <- subset(Reference.sn, subset = CellType %in% c("Epithelial"))
Reference.sn@meta.data$infercnv <- Reference.sn@meta.data$orig.ident


##############################################################
########################## Patients ##########################
##############################################################
datalist <- c("P01","P02","P03","P04","P05","P06","P07")

for (i in datalist) {
  NPC.tmp <- subset(NPC.combined, subset = orig.ident %in% c(i))
  NPC.tmp@meta.data$infercnv <- NPC.tmp@meta.data$Annotation.Coarse
  
  ##############################################################
  # prepare data
  if (i %in% c("P01","P02","P03")) {
    NPC.test <- merge(Reference.sc, y = c(NPC.tmp),merge.data = TRUE)
    Ref <- c("N01","N02")
    }
  if (i %in% c("P04","P05","P06","P07")) {
    NPC.test <- merge(Reference.sn, y = c(NPC.tmp),merge.data = TRUE)
    Ref <- c("N03")
    }
  
  counts_matrix = as.data.frame(NPC.test@assays$RNA@counts)
  tmp_ann <- colnames(counts_matrix)
  file_ann <- matrix(nrow = length(tmp_ann),ncol=2)
  file_ann[1:length(tmp_ann),1] <- tmp_ann
  file_ann[1:length(tmp_ann),2] <- NPC.test@meta.data$infercnv
  
  ##############################################################
  # infercnv
  path.annotation = paste("G:/scRNA data/DS0/infercnv/",i,"/annotations.txt",sep="")
  write.table(file_ann, file= path.annotation, row.names=FALSE, col.names=FALSE, quote =FALSE, sep ="\t")
  path.save = paste("G:/scRNA data/DS0/infercnv/",i,"/counts_matrix",sep="")
  save(counts_matrix, file = path.save)
  
  #path.load = path.save
  #load(path.load)
  path.out <- paste("G:/scRNA data/DS0/infercnv/",i,sep="")
  
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

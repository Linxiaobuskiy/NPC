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


NPC.Fib <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Fibroblast"))
NPC.Fib@meta.data$number <- 1:length(NPC.Fib@meta.data$orig.ident)
load("G:/scRNA data/DS0/step5.FibroblastAnnotationCNV.RData")
NPC.Fib@meta.data$Annotation.CNV <- Annotation.Fib


###############################################################
########################## Reference ##########################
###############################################################
Reference.sn <- subset(NPC.Fib, subset = orig.ident %in% c("N03"))
#Reference.sn <- subset(Reference.sn, subset = CellType %in% c("Fibthelial"))
Reference.sn@meta.data$infercnv <- "N03"


##############################################################
########################## Patients ##########################
##############################################################
datalist <- c("P04","P05","P06","P07")

NPC.tmp <- subset(NPC.Fib, subset = orig.ident %in% datalist)
NPC.tmp@meta.data$infercnv <- paste(NPC.tmp@meta.data$orig.ident,NPC.tmp@meta.data$Annotation.CNV,sep=": ")

##############################################################
# prepare data

NPC.test <- merge(Reference.sn, y = c(NPC.tmp),merge.data = TRUE)
Ref <- c("N03")

counts_matrix = as.data.frame(NPC.test@assays$RNA@counts)
tmp_ann <- colnames(counts_matrix)
file_ann <- matrix(nrow = length(tmp_ann),ncol=2)
file_ann[1:length(tmp_ann),1] <- tmp_ann
file_ann[1:length(tmp_ann),2] <- NPC.test@meta.data$infercnv

##############################################################
# infercnv
path.annotation = paste("G:/scRNA data/DS0/infercnv/Fibroblast/annotations.txt",sep="")
write.table(file_ann, file= path.annotation, row.names=FALSE, col.names=FALSE, quote =FALSE, sep ="\t")
path.save = paste("G:/scRNA data/DS0/infercnv/Fibroblast/counts_matrix",sep="")
save(counts_matrix, file = path.save)

#path.load = path.save
#load(path.load)
path.out <- "G:/scRNA data/DS0/infercnv/Fibroblast"

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



rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)
library(scater)
library(patchwork)
library(infercnv)

load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse


NPC.Fib <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Fibroblast"))
NPC.Fib@meta.data$number <- 1:length(NPC.Fib@meta.data$orig.ident)
#DimPlot(NPC.Fib, reduction = "umap", group.by = "Annotation.Fib",label = TRUE)


###############################################################
########################## Reference ##########################
###############################################################
Reference.sn <- subset(NPC.Fib, subset = orig.ident %in% c("N03"))
#Reference.sn <- subset(Reference.sn, subset = CellType %in% c("Fibthelial"))
Reference.sn@meta.data$infercnv <- Reference.sn@meta.data$orig.ident


tt = "P07"
#############################################################
##########################   P0?   ##########################
#############################################################
# First order annotation
NPC.tmp <- subset(NPC.Fib, subset = orig.ident %in% c(tt))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 30
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.8)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE,group.by = "seurat_clusters")
FeaturePlot(NPC.tmp,features = c("KRT19"))

##############################################################
########################## Patients ##########################
##############################################################

NPC.tmp@meta.data$infercnv <- NPC.tmp@meta.data$seurat_clusters

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
path.annotation = paste("G:/scRNA data/DS0/infercnv/test_Fibroblast/annotations.txt",sep="")
write.table(file_ann, file= path.annotation, row.names=FALSE, col.names=FALSE, quote =FALSE, sep ="\t")
path.save = paste("G:/scRNA data/DS0/infercnv/test_Fibroblast/counts_matrix",sep="")
save(counts_matrix, file = path.save)

#path.load = path.save
#load(path.load)
path.out <- paste("G:/scRNA data/DS0/infercnv/test_Fibroblast",sep="")

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






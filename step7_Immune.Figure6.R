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

load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse
NPC.combined@meta.data$number <- 1:length(NPC.combined@meta.data$orig.ident)


NPC.combined@meta.data$Annotation.Immune2 <- "tmp"

###############################################################
##########################   PTPRC   ##########################
###############################################################
## first
NPC.tmp <- subset(NPC.combined, subset = Annotation.Coarse %in% c("B","NK/T"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 2)#resoulution越小分的越粗
p1 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE, group.by = "seurat_clusters")

p3 <- FeaturePlot(NPC.tmp, features = c(
  "PTPRC",
  "GNLY",#NK
  "CD79A",
  "CD4",#T-help
  "CD8A",#T-cytoxic
  "FOXP3",#Treg
  "LAG3","TIGIT",#T-exhausted
  "CD3D", # T
  "MKI67",
  "IGHM","IGHD",# Naive CD27-
  "FCRL4",
  "SERPINA9",
  "CD27" # Memory
))

#p3+p1


Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("T_CD4+","Memory B","T_CD8+","Memory B","Naive B",
                     "T_CD4+","T_CD8+","Memory B","T_CD4+","Memory B",
                     "T_CD4+","Naive B","T-reg","T-reg","T_CD8+",
                     "GCB_MKI67-","Naive B","T_CD4+","B/T_CD4+","GCB_MKI67+",
                     "T_CD4+","T_CD8+","T_CD4+","Memory B","B/T_CD8+",
                     "Memory B","T_CD8+","T_MKI67+","GCB_MKI67+","Naive B",
                     "T_CD8+","T_MKI67+","Naive B","Memory B")

names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Immune2 <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Immune2")
NPC.combined@meta.data$Annotation.Immune2[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Immune2)



## second NK
NPC.tmp <- subset(NPC.tmp, subset = Annotation.Immune2 %in% c("T_CD8+"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 10
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 1)#resoulution越小分的越粗

p1 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE, group.by = "seurat_clusters")
p2 <- FeaturePlot(NPC.tmp, features = c("GNLY","CD8A"))

#p2+p1


Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("T_CD8+","T_CD8+","T_CD8+","T_CD8+","T_CD8+",
                     "T_CD8+","T_CD8+","T_CD8+","T_CD8+","T_CD8+",
                     "T_CD8+","NK_CD8+","NK_CD8-","NK_CD8-","T_CD8+"
)

names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Immune2 <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Immune2")
NPC.combined@meta.data$Annotation.Immune2[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Immune2)



## third FCRL4+
NPC.tmp <- subset(NPC.combined, subset = Annotation.Immune2 %in% c("Memory B"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 2)#resoulution越小分的越粗

p1 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE, group.by = "seurat_clusters")
p2 <- FeaturePlot(NPC.tmp, features = c("FCRL4"))

#p2+p1


Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("Memory B","Memory B","Memory B","Memory B","Memory B","Memory B","Memory B","Memory B","Memory B",
                     "Memory B","Memory B","Memory B_FCRL4+","Memory B","Memory B","Memory B","Memory B","Memory B_FCRL4+","Memory B",
                     "Memory B","Memory B","Memory B","Memory B","Memory B","Memory B","Memory B","Memory B","Memory B")

names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Immune2 <- NPC.tmp@active.ident
DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Immune2")
NPC.combined@meta.data$Annotation.Immune2[NPC.tmp@meta.data$number] <- as.character(NPC.tmp@meta.data$Annotation.Immune2)


## FIGURE
NPC.tmp <- subset(NPC.combined, subset = Annotation.Coarse %in% c("B","NK/T"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
p1<-DimPlot(NPC.tmp, reduction = "umap",label = TRUE, group.by = "Annotation.Immune2")
p2<-FeaturePlot(NPC.tmp, features = c("GNLY"))
#p1+p2


save(NPC.tmp, file = "G:/scRNA data/DS0/step7.AnnotationImmune.RData")


################################################################
##########################   Figure   ##########################
################################################################
load("G:/scRNA data/DS0/step7.AnnotationImmune.RData")

## supplementary figures
NPC.sc <- subset(NPC.tmp,subset = SeqTech == "Single-cell")
VlnPlot(NPC.sc, features = c("nFeature_RNA","nCount_RNA"),group.by = "Annotation.Immune2",pt.size = 0)


## figure 6a
colors1 <- c(
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3",
           "#57C3F3", "#476D87", "#E95C59", "#91D0BE", "#AB3282",
           "#23452F", "#BD956A", "#8C549C")
           
p1 <- DimPlot(NPC.tmp, label = F, cols= colors1, pt.size = 0.6, repel = T, group.by = "Annotation.Immune2", reduction = "umap")+
  NoLegend()+ NoAxes()+labs(title = NULL)
p2 <- DimPlot(NPC.tmp, label = T, cols= colors1, pt.size = 0.6, repel = T, group.by = "Annotation.Immune2", reduction = "umap")+
  theme(legend.text = element_text(size=5),
        legend.key.height = unit(0.2, "cm"))

# save p1 as 8inch*8inch pdf
# Transfer p1 to semi-png figure
# Adjust p1 to 50mm*50mm
# Add labels and legend from p2 to p1
# Adjust word to 5pt, label gap 2mm, square-word left-left 2.5mm
# Add Axes

## figure 6b
p1 <- DimPlot(NPC.tmp, label = F, cols= colors1, pt.size = 0.6, repel = T, group.by = "orig.ident", reduction = "umap")+
  NoLegend()+ NoAxes()+labs(title = NULL)
p2 <- DimPlot(NPC.tmp, label = T, cols= colors1, pt.size = 0.6, repel = T, group.by = "orig.ident", reduction = "umap")+
  theme(legend.text = element_text(size=5),
        legend.key.height = unit(0.2, "cm"))

# save p1 as 8inch*8inch pdf
# Transfer p1 to semi-png figure
# Adjust p1 to 50mm*50mm
# Add labels and legend from p2 to p1
# Adjust word to 5pt, label gap 2mm, square-word left-left 2.5mm
# Add Axes

## figure 6c
markergenes <- c(
  "PTPRC",
  "CD79A","MS4A1",#B
  "CD3D","CD4","CD8A",#T
  "GNLY",#NK
  "FOXP3",#Treg
  "CD27", # Memory
  "IGHM","IGHD",# Naive CD27-
  "FCRL4",
  "SERPINA9",
  "LAG3","TIGIT",#T-exhausted
  "MKI67")
p3 <- FeaturePlot(NPC.tmp, features = markergenes)
for(i in 1:length(markergenes)) {
  p3[[i]] <- p3[[i]] + 
    NoLegend() + 
    NoAxes() + 
    theme(plot.title = element_text(size = 5))
}
p3


for(i in markergenes) {
  p <- FeaturePlot(NPC.tmp, features = i) + NoLegend() + NoAxes() + theme(plot.title = element_text(size = 10))
  ggsave(paste("G:/manuscript/Figure/BT_",i,".pdf",sep=""), 
         egg::set_panel_size(p, width=unit(40, "mm"), height=unit(40, "mm")), 
         width = 100, height = 100, units = 'mm', dpi = 300)
  #resize to 1/2
}

###########################################################################
######################### Figure Cell Composition #########################
###########################################################################

tmp <- subset(NPC.tmp, subset = Annotation.Immune2 %in% c("NK_CD8+","NK_CD8-"))

library(ggpubr)
pie <- table(tmp$orig.ident,tmp$Annotation.Immune2) %>% data.frame
df.CD8.neg <- data.frame(
  group = c("Normal", "Tumor"),
  value = c(222,382))
df.CD8.neg$feq <- round(df.CD8.neg$value/sum(df.CD8.neg$value)*100,2)
df.CD8.pos <- data.frame(
  group = c("Normal", "Tumor"),
  value = c(0,400)
)
df.CD8.pos$feq <- round(df.CD8.pos$value/sum(df.CD8.pos$value)*100,2)
labs.CD8.neg <- paste0(df.CD8.neg$group, " (", df.CD8.neg$feq, "%)")
labs.CD8.pos <- paste0(df.CD8.pos$group, " (", df.CD8.pos$feq, "%)")
p.neg <- ggpie(df.CD8.neg, 'feq',
               fill = "group",
               palette = c("#00AFBB", "#FC4E07"),
               label = labs.CD8.neg, lab.pos = 'in', lab.font = c(5, 'white'))
p.pos <- ggpie(df.CD8.pos, 'feq',
               fill = "group",
               palette = c("#00AFBB", "#FC4E07"),
               label = labs.CD8.pos, lab.pos = 'in', lab.font = c(5, 'white'))

p.neg+p.pos


############################################################
Idents(NPC.tmp) <- 'Annotation.Immune2'
NPC.tmp@meta.data$source <- substr(NPC.tmp@meta.data$orig.ident,1,1)
pt <- table(NPC.tmp$source,Idents(NPC.tmp))
pt <- as.data.frame(pt)
celltypes <- unique(NPC.tmp@meta.data$Annotation.Immune2) %>% as.character %>% str_sort
pt$Var1 <- factor(pt$Var1,levels=c('P','N'))
pt$Var2 <- factor(pt$Var2,levels=celltypes)

p1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  #theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("Sample ID") +
  ylab("Proportion") +
  #theme(legend.position="none") +
  #scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=5,color="black",face="plain"),
        axis.text.x = element_text(size=5,color="black",face="plain",angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=5,color="black",face="plain"))

path.save <- paste("G:/manuscript/Figure/CC.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p1, width=unit(30, "mm"), height=unit(30, "mm")), 
       width = 100, height = 100, units = 'mm', dpi = 300)

#################################################################
######################### Figure NK DEG #########################
#################################################################

## NK_CD8+ vs NK_CD8- in only Patient Group
#NPC.tmp@meta.data$source <- substr(NPC.tmp@meta.data$orig.ident,1,1)
#NPC.tmp <- subset(NPC.tmp, subset = source %in% c("P"))

## NK_CD8+ vs NK_CD8- in both Patient and Normal Group
tmp <- subset(NPC.tmp, subset = Annotation.Immune2 %in% c("NK_CD8+","NK_CD8-"))
tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
tmp <- ScaleData(tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

##heatmap
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'Annotation.Immune2'
#markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(tmp, features = top10$gene,group.by = "Annotation.Immune2")

## Go
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(stringr)

up.genes.enrich <- subset(markers,subset = cluster %in% c("NK_CD8+") )
up.genes.enrich <- rownames(up.genes.enrich)
up.genes.enrich <- up.genes.enrich[!grepl("^RP[SL]", up.genes.enrich, ignore.case = F)]
down.genes.enrich <- subset(markers,subset = cluster %in% c("NK_CD8-") )
down.genes.enrich <- rownames(down.genes.enrich)
down.genes.enrich <- down.genes.enrich[!grepl("^RP[SL]", down.genes.enrich, ignore.case = F)]

up.genes.enrich <- bitr(up.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
down.genes.enrich <- bitr(down.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
#up.genes.enrich <- pull(up.genes.enrich,ENTREZID)  

up.genes.go <- enrichGO(gene          = up.genes.enrich$SYMBOL,
                        #universe     = row.names(dge.celltype),
                        OrgDb         = 'org.Hs.eg.db',
                        keyType       = 'SYMBOL',
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.05)
down.genes.go <- enrichGO(gene          = down.genes.enrich$SYMBOL,
                          #universe     = row.names(dge.celltype),
                          OrgDb         = 'org.Hs.eg.db',
                          keyType       = 'SYMBOL',
                          ont           = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.1,
                          qvalueCutoff  = 0.05)
#head(up.genes.go)

barplot(up.genes.go, showCategory=20)




n.top = 15

up.genes.go@result$logg <- (-log(up.genes.go@result$p.adjust))
tmp <- up.genes.go %>% top_n(n = n.top, wt = logg)
p.up <- ggplot(tmp@result)+
  geom_bar(aes(x = Description,y = -log(p.adjust), fill = ONTOLOGY), stat='identity') + 
  coord_flip() +
  scale_x_discrete(limits=tmp@result$Description,position = "top") +
  #scale_x_discrete(limits=tmp@result$Description,position = "top",labels=function(x) str_wrap(x, width=30)) +
  theme_bw() + #theme(panel.grid=element_blank())+
  #scale_y_continuous(limits=c(0,150))+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"),
  )+NoLegend()



path.save <- paste("G:/manuscript/Figure/CD8+NK.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p.up, width=unit(0.8, "in"), height=unit(1.2, "in")), 
       width = 10, height = 10, units = 'in', dpi = 300)

#############################################################################
######################### Figure NP DEG (intersect) #########################
#############################################################################

NPC.tmp@meta.data$source <- substr(NPC.tmp@meta.data$orig.ident,1,1)

celltypes <- unique(NPC.tmp@meta.data$Annotation.Immune2) %>% as.character %>% str_sort

# T cells
tmp <- subset(NPC.tmp, subset = Annotation.Immune2 %in% celltypes[c(1,2,9,10,11,12,13)])

tmp <- subset(NPC.tmp, subset = Annotation.Immune2 %in% celltypes[c(11)])

DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'source'
tmp.sc <- subset(tmp, subset = SeqTech %in% c("Single-cell"))
tmp.sn <- subset(tmp, subset = SeqTech %in% c("Single-nucleus"))

markers.sc <- FindAllMarkers(object = tmp.sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
markers.sn <- FindAllMarkers(object = tmp.sn, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)

up.genes.enrich.sc <- subset(markers.sc,subset = cluster %in% c("P") )
up.genes.enrich.sc <- rownames(up.genes.enrich.sc)
up.genes.enrich.sc <- up.genes.enrich.sc[!grepl("^RP[SL]", up.genes.enrich.sc, ignore.case = F)]
down.genes.enrich.sc <- subset(markers.sc,subset = cluster %in% c("N") )
down.genes.enrich.sc <- rownames(down.genes.enrich.sc)
down.genes.enrich.sc <- down.genes.enrich.sc[!grepl("^RP[SL]", down.genes.enrich.sc, ignore.case = F)]

up.genes.enrich.sn <- subset(markers.sn,subset = cluster %in% c("P") )
up.genes.enrich.sn <- rownames(up.genes.enrich.sn)
up.genes.enrich.sn <- up.genes.enrich.sn[!grepl("^RP[SL]", up.genes.enrich.sn, ignore.case = F)]
down.genes.enrich.sn <- subset(markers.sn,subset = cluster %in% c("N") )
down.genes.enrich.sn <- rownames(down.genes.enrich.sn)
down.genes.enrich.sn <- down.genes.enrich.sn[!grepl("^RP[SL]", down.genes.enrich.sn, ignore.case = F)]
up.genes.enrich <- intersect(x=up.genes.enrich.sc, y = up.genes.enrich.sn)
down.genes.enrich <- intersect(x=down.genes.enrich.sc, y = down.genes.enrich.sn)

## Go
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(stringr)

up.genes.enrich <- bitr(up.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
down.genes.enrich <- bitr(down.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
#up.genes.enrich <- pull(up.genes.enrich,ENTREZID)  

up.genes.go <- enrichGO(gene          = up.genes.enrich$SYMBOL,
                        #universe     = row.names(dge.celltype),
                        OrgDb         = 'org.Hs.eg.db',
                        keyType       = 'SYMBOL',
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.05)
down.genes.go <- enrichGO(gene          = down.genes.enrich$SYMBOL,
                          #universe     = row.names(dge.celltype),
                          OrgDb         = 'org.Hs.eg.db',
                          keyType       = 'SYMBOL',
                          ont           = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.1,
                          qvalueCutoff  = 0.05)
#head(up.genes.go)

barplot(up.genes.go, showCategory=20)




n.top = 15

up.genes.go@result$logg <- (-log(up.genes.go@result$p.adjust))
tmp <- up.genes.go %>% top_n(n = n.top, wt = logg)
p.up <- ggplot(tmp@result)+
  geom_bar(aes(x = Description,y = -log(p.adjust), fill = ONTOLOGY), stat='identity') + 
  coord_flip() +
  scale_x_discrete(limits=tmp@result$Description,position = "top") +
  #scale_x_discrete(limits=tmp@result$Description,position = "top",labels=function(x) str_wrap(x, width=30)) +
  theme_bw() + #theme(panel.grid=element_blank())+
  #scale_y_continuous(limits=c(0,150))+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"),
  )+NoLegend()



path.save <- paste("G:/manuscript/Figure/CD8+NK.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p.up, width=unit(0.8, "in"), height=unit(1.2, "in")), 
       width = 10, height = 10, units = 'in', dpi = 300)


#################################################################################
######################### Figure NP DEG (not intersect) #########################
#################################################################################

NPC.tmp@meta.data$source <- substr(NPC.tmp@meta.data$orig.ident,1,1)

celltypes <- unique(NPC.tmp@meta.data$Annotation.Immune2) %>% as.character %>% str_sort

# T cells
tmp <- subset(NPC.tmp, subset = Annotation.Immune2 %in% celltypes[c(13)])


##heatmap
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'source'
#markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(tmp, features = top10$gene,group.by = "source")

## Go
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(stringr)

up.genes.enrich <- subset(markers,subset = cluster %in% c("P") )
up.genes.enrich <- rownames(up.genes.enrich)
up.genes.enrich <- up.genes.enrich[!grepl("^RP[SL]", up.genes.enrich, ignore.case = F)]
down.genes.enrich <- subset(markers,subset = cluster %in% c("N") )
down.genes.enrich <- rownames(down.genes.enrich)
down.genes.enrich <- down.genes.enrich[!grepl("^RP[SL]", down.genes.enrich, ignore.case = F)]

up.genes.enrich <- bitr(up.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')

up.genes.go <- enrichGO(gene          = up.genes.enrich$SYMBOL,
                        #universe     = row.names(dge.celltype),
                        OrgDb         = 'org.Hs.eg.db',
                        keyType       = 'SYMBOL',
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.05)

#head(up.genes.go)

barplot(up.genes.go, showCategory=20)




n.top = 15

up.genes.go@result$logg <- (-log(up.genes.go@result$p.adjust))
tmp <- up.genes.go %>% top_n(n = n.top, wt = logg)
p.up <- ggplot(tmp@result)+
  geom_bar(aes(x = Description,y = -log(p.adjust), fill = ONTOLOGY), stat='identity') + 
  coord_flip() +
  scale_x_discrete(limits=tmp@result$Description,position = "top") +
  #scale_x_discrete(limits=tmp@result$Description,position = "top",labels=function(x) str_wrap(x, width=30)) +
  theme_bw() + #theme(panel.grid=element_blank())+
  #scale_y_continuous(limits=c(0,150))+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain")
  )+NoLegend()



path.save <- paste("G:/manuscript/Figure/T-reg.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p.up, width=unit(0.8, "in"), height=unit(1.2, "in")), 
       width = 10, height = 10, units = 'in', dpi = 300)

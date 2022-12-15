rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)

library(scater)
library(patchwork)

# https://www.nature.com/articles/s41467-020-18916-5

load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse


NPC.Fib <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Fibroblast"))
NPC.Fib@meta.data$number <- 1:length(NPC.Fib@meta.data$orig.ident)
load("G:/scRNA data/DS0/step5.FibroblastAnnotationCNV.RData")
NPC.Fib@meta.data$Annotation.CNV <- Annotation.Fib
NPC.Fib@meta.data$Annotation.Fib <- "Fibroblast"
DimPlot(NPC.Fib, reduction = "umap", group.by = "Annotation.CNV",label = TRUE)

#############################################################
##########################   End   ##########################
#############################################################
NPC.Nor <- subset(NPC.Fib, subset = orig.ident %in% c("N03"))
NPC.tmp <- subset(NPC.Fib, subset = orig.ident %in% c("P04","P05","P06","P07"))

NPC.tmp <- NPC.Fib

NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp <- ScaleData(NPC.tmp)

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 10
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 0.05)#resoulution越小分的越粗

# Annotation
Idents(NPC.tmp) <- 'seurat_clusters'
new.cluster.ids <- c("C0","C1","C2","C3")
names(new.cluster.ids) <- levels(NPC.tmp)
NPC.tmp <- RenameIdents(NPC.tmp, new.cluster.ids)
NPC.tmp@meta.data$Annotation.Fib <- NPC.tmp@active.ident

# NPC.Fib@meta.data$Annotation.Fib[NPC.tmp@meta.data$number] <- paste(NPC.tmp@meta.data$orig.ident,"_",NPC.tmp@meta.data$Annotation.Fib,sep="")
NPC.Fib@meta.data$Annotation.Fib[match(colnames(NPC.tmp),colnames(NPC.Fib))] <- as.character(NPC.tmp@meta.data$Annotation.Fib)

# View
p1 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE,group.by = "Annotation.CNV")
p2 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE)
p3 <- DimPlot(NPC.tmp, reduction = "umap",label = TRUE,group.by = "orig.ident")
p1+p2+p3

# NPC.Fib <- NPC.tmp
# save(NPC.Fib, file = "G:/scRNA data/DS0/step5.FibroblastAnnotatedSubtype.RData")


# Figures 3 abc
p.cnv.1 <- DimPlot(NPC.tmp, reduction = "umap",group.by = "Annotation.CNV", pt.size = 2) + NoLegend() + NoAxes() + labs(title = NULL)
p.cnv.2 <- DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.CNV",label=TRUE)

p.ann.1 <- DimPlot(NPC.tmp, reduction = "umap",group.by = "Annotation.Fib", pt.size = 2) + NoLegend() + NoAxes() + labs(title = NULL)
p.ann.2 <- DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Fib",label=TRUE)

p.source.1 <- DimPlot(NPC.tmp, reduction = "umap",group.by = "orig.ident", pt.size = 2) + NoLegend() + NoAxes() + labs(title = NULL)
p.source.2 <- DimPlot(NPC.tmp, reduction = "umap", group.by = "orig.ident",label=TRUE)

# save as 10inch*10inch pdf
# label d=1.5mm gap=2.5mm


#############################################################
########################## Dotplot ##########################
#############################################################

# Figures 3 d
tmp <- NPC.tmp
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'Annotation.Fib'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

dot.plot <- DotPlot(NPC.tmp, features = top10$gene, dot.scale = 2)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,size=5,color="black",face="plain",hjust = 1,vjust=0.5),
        axis.text.y = element_text(size=5,color="black",face="plain")
  )+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))


ggsave("G:/manuscript/Figure/dot.plot.pdf", egg::set_panel_size(dot.plot, width=unit(10, "mm"), height=unit(100, "mm")), 
       width = 100, height = 250, units = 'mm', dpi = 300)

############################################################
########################## Scores ##########################
############################################################

ImmunoSuppression <- c("TGFB1","IL6","IL11","LIF","CSF2","CXCL1","CXCL2","CXCL8","CXCL9","CXCL10","CXCL12","CCL2","CCL8","PTK2")
Mitochondrion <- grep(pattern = "^MT-", x = rownames(NPC.tmp@assays[["RNA"]]), value = TRUE)
MHC <- grep(pattern = "^HLA-", x = rownames(NPC.combined@assays[["RNA"]]), value = TRUE)
MHC_I <- c("HLA-A","HLA-B","HLA-C","HLA-E", "HLA-F", "HLA-G")
MHC_II <- c("HLA-DRB5", "HLA-DQA2", "HLA-DQB2", "HLA-DOA", "HLA-DOB", "HLA-DMA", "HLA-DMB")
#MHC_II <- c("HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1")
ComplementRegulation = c("CFD", "C1QA","C1QB","C1QC")
Collagen <- grep(pattern = "^COL.A.", x = rownames(NPC.tmp@assays[["RNA"]]), value = TRUE)

NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(ImmunoSuppression),ctrl = 100,name = "score.ImmunoSuppression")
NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(MHC),ctrl = 100,name = "score.MHC")
NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(MHC_I),ctrl = 100,name = "score.MHC_I")
NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(MHC_II),ctrl = 100,name = "score.MHC_II")
NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(ComplementRegulation),ctrl = 100,name = "score.ComplementRegulation")
NPC.tmp <- AddModuleScore(object = NPC.tmp,features = list(Collagen),ctrl = 100,name = "score.Collagen")

VlnPlot(NPC.tmp, features=c("score.MHC1","score.MHC_I1","score.MHC_II1"),pt.size = 0)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
NPC.tmp@meta.data$score.MHC <- range01(NPC.tmp@meta.data$score.MHC1)
NPC.tmp@meta.data[["score.MHC-I"]] <- range01(NPC.tmp@meta.data$score.MHC_I1)
NPC.tmp@meta.data[["score.MHC-II"]] <- range01(NPC.tmp@meta.data$score.MHC_II1)



figure.MHC <- VlnPlot(NPC.tmp, features=c("score.MHC"),pt.size = 0)+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"))+
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.title = element_text(size = 5))
figure.MHC1 <- VlnPlot(NPC.tmp, features=c("score.MHC-I"),pt.size = 0)+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"))+
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.title = element_text(size = 5)) 
figure.MHC2 <- VlnPlot(NPC.tmp, features=c("score.MHC-II"),pt.size = 0)+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"))+
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.title = element_text(size = 5))
ggsave("G:/manuscript/Figure/MHC.pdf", egg::set_panel_size(figure.MHC, width=unit(20, "mm"), height=unit(20, "mm")), 
       width = 100, height = 100, units = 'mm', dpi = 300)
ggsave("G:/manuscript/Figure/MHC1.pdf", egg::set_panel_size(figure.MHC1, width=unit(20, "mm"), height=unit(20, "mm")), 
       width = 100, height = 100, units = 'mm', dpi = 300)
ggsave("G:/manuscript/Figure/MHC2.pdf", egg::set_panel_size(figure.MHC2, width=unit(20, "mm"), height=unit(20, "mm")), 
       width = 100, height = 100, units = 'mm', dpi = 300)

# Figures S
figure.PDGFRA <- VlnPlot(NPC.tmp, features=c("PDGFRA"),pt.size = 0)+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"))+
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.title = element_text(size = 5))
ggsave("G:/manuscript/Figure/PDGFRA.pdf", egg::set_panel_size(figure.PDGFRA, width=unit(20, "mm"), height=unit(20, "mm")), 
       width = 100, height = 100, units = 'mm', dpi = 300)

####################################################################
##########################  Go: Clusters  ##########################
####################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

tmp <- NPC.tmp
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'Annotation.Fib'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

new.cluster.ids <- c("C0","C1","C2","C3")

for (i in 1:4) {
  
  
  up.genes.enrich <- subset(markers,subset = cluster %in% new.cluster.ids[i])
  up.genes.enrich <- rownames(up.genes.enrich)
  up.genes.enrich <- up.genes.enrich[!grepl("^RP[SL]", up.genes.enrich, ignore.case = F)]
  ################################################################################
  # GO
  up.genes.enrich <- bitr(up.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
  #up.genes.enrich <- pull(up.genes.enrich,ENTREZID)  
  
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
  
  
  
  ################################################################################
  # figure Enrichment
  n.top = 4
  
  up.genes.go@result$logg <- (-log(up.genes.go@result$p.adjust))
  tmp <- up.genes.go %>% group_by(ONTOLOGY) %>% top_n(n = n.top, wt = logg)
  p.up <- ggplot(tmp@result)+
    geom_bar(aes(x = Description,y = -log(p.adjust), fill = ONTOLOGY), stat='identity') + 
    coord_flip() +
    scale_x_discrete(limits=tmp@result$Description,position = "top") +
    theme_bw() + #theme(panel.grid=element_blank())+
    #scale_y_continuous(limits=c(0,150))+
    theme(axis.text.x = element_text(size=5,color="black",face="plain"),
          axis.text.y = element_text(size=5,color="black",face="plain"),
    )
  
  p.up
  
  path.save <- paste("G:/manuscript/Figure/p.up.",i,".pdf",sep="")
  ggsave(path.save, egg::set_panel_size(p.up, width=unit(1, "in"), height=unit(1, "in")), 
         width = 10, height = 10, units = 'in', dpi = 300)
}

###############################################################
##########################  Go: CNV  ##########################
###############################################################
#载入所需的R包；
library(ggplot2)
library(ggsci)
library(sf)
library(ggVennDiagram)

color1 <- alpha("#BC3C29FF",1)
color2 <- alpha("#0072B5FF",0.5)

################################################################################
# DEGs

DefaultAssay(NPC.tmp) <- "RNA"
Idents(NPC.tmp) <- 'Annotation.CNV'




for (i in c("P04","P05","P06","P07")) {
  tmp <- subset(NPC.tmp, subset = orig.ident %in% i)
  markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  CNV.markers <- subset(markers, subset = cluster %in% c("CNV"))
  nonCNV.markers <- subset(markers, subset = cluster %in% c("non-CNV"))
  DEGs <- list(CNV.markers,nonCNV.markers)
  names(DEGs) <- c("CNV.markers","nonCNV.markers")
  assign(i,DEGs)
  rm(DEGs)
  rm(tmp)
}

################################################################################
################################################################################
# CNV up-regular genes
CNV.markers = list()
CNV.markers[["P04"]] <- P04$CNV.markers$gene
CNV.markers[["P05"]] <- P05$CNV.markers$gene
CNV.markers[["P06"]] <- P06$CNV.markers$gene
CNV.markers[["P07"]] <- P07$CNV.markers$gene


################################################################################
#venn

#label_alpha = 0去除文字标签底色；
#category.names参数用于设定样本名称；
p1 <- ggVennDiagram(CNV.markers, label_alpha=0, label_size =5, edge_size = 1,label = "count", edge_lty = "dashed") +
  scale_color_brewer(palette = "Paired")+
  #scale_fill_gradient(low="white",high = color1)
  scale_fill_distiller(palette = "Reds", direction = 1)+scale_color_brewer(palette = "Set1")

CNV.genes <- Reduce(intersect, CNV.markers)

CNV.genes.enrich <- CNV.genes

################################################################################
# GO
CNV.genes.enrich <- bitr(CNV.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
CNV.genes.go <- enrichGO(gene          = CNV.genes.enrich$ENTREZID,
                         #universe     = row.names(dge.celltype),
                         OrgDb         = 'org.Hs.eg.db',
                         #keyType       = 'SYMBOL',
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.1,
                         qvalueCutoff  = 0.05)

################################################################################
################################################################################
# nonCNV up-regular genes
nonCNV.markers = list()
nonCNV.markers[["P04"]] <- P04$nonCNV.markers$gene
nonCNV.markers[["P05"]] <- P05$nonCNV.markers$gene
nonCNV.markers[["P06"]] <- P06$nonCNV.markers$gene
nonCNV.markers[["P07"]] <- P07$nonCNV.markers$gene


################################################################################
#venn

#label_alpha = 0去除文字标签底色；
#category.names参数用于设定样本名称；
p2 <- ggVennDiagram(nonCNV.markers, label_alpha=0, label_size =5, edge_size = 1,label = "count", edge_lty = "dashed") +
  scale_color_brewer(palette = "Paired")+
  #scale_fill_gradient(low="white",high = color1)
  scale_fill_distiller(palette = "Blues", direction = 1)+scale_color_brewer(palette = "Set1")+
  theme(legend.text = element_text(size = 5))

nonCNV.genes <- Reduce(intersect, nonCNV.markers)

nonCNV.genes.enrich <- nonCNV.genes

################################################################################
# GO
nonCNV.genes.enrich <- bitr(nonCNV.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
nonCNV.genes.go <- enrichGO(gene          = nonCNV.genes.enrich$ENTREZID,
                            #universe     = row.names(dge.celltype),
                            OrgDb         = 'org.Hs.eg.db',
                            #keyType       = 'SYMBOL',
                            ont           = "ALL",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.1,
                            qvalueCutoff  = 0.05)

################################################################################
################################################################################
################################################################################
# figure Enrichment
n.top = 10


nonCNV.genes.go@result$logg <- (log(nonCNV.genes.go@result$p.adjust))

tmp <- nonCNV.genes.go %>% group_by(ONTOLOGY) %>% top_n(n = n.top, wt = logg)

p.nonCNV <- ggplot(tmp@result)+
  geom_bar(aes(x = Description,y = -log(p.adjust), fill = ONTOLOGY), stat='identity') + 
  coord_flip() + 
  scale_x_discrete(limits=tmp@result$Description,position = "bottom") +
  theme_bw() + #theme(panel.grid=element_blank())+
  #scale_y_continuous(limits=c(-150,0))+
  ylab("-log(p.adjust)")+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"),
  )


CNV.genes.go@result$logg <- (-log(CNV.genes.go@result$p.adjust))

tmp <- CNV.genes.go %>% group_by(ONTOLOGY) %>% top_n(n = n.top, wt = logg)

p.CNV <- ggplot(tmp@result)+
  geom_bar(aes(x = Description,y = -log(p.adjust), fill = ONTOLOGY), stat='identity') + 
  coord_flip() +
  scale_x_discrete(limits=tmp@result$Description,position = "top") +
  theme_bw() + #theme(panel.grid=element_blank())+
  #scale_y_continuous(limits=c(0,150))+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"),
  )

p.nonCNV+p.CNV


ggsave("G:/manuscript/Figure/venn.nonCNV.pdf", egg::set_panel_size(p2, width=unit(140, "mm"), height=unit(100, "mm")), 
       width = 10, height = 10, units = 'in', dpi = 300)
ggsave("G:/manuscript/Figure/venn.CNV.pdf", egg::set_panel_size(p1, width=unit(140, "mm"), height=unit(100, "mm")),  
       width = 10, height = 10, units = 'in', dpi = 300)
# default 14 pt
# 140 -> 50
ggsave("G:/manuscript/Figure/p.nonCNV.pdf", egg::set_panel_size(p.nonCNV, width=unit(1, "in"), height=unit(2, "in")), 
       width = 10, height = 10, units = 'in', dpi = 300)
ggsave("G:/manuscript/Figure/p.CNV.pdf", egg::set_panel_size(p.CNV, width=unit(2, "in"), height=unit(2, "in")), 
       width = 10, height = 10, units = 'in', dpi = 300)


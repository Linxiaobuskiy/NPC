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


load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse
NPC.combined@meta.data$Annotation.Immune <- NPC.combined@meta.data$Annotation.Coarse
NPC.combined@meta.data$number <- 1:length(NPC.combined@meta.data$orig.ident)


load("G:/scRNA data/DS0/step7.AnnotatedImmune.RData")
NPC.combined@meta.data$Annotation.Immune <- Annotation.Immune

#################################################################
##########################   Myeloid   ##########################
#################################################################
## Myeloid
NPC.tmp <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Myeloid"))

NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp@assays$RNA@var.features <- NPC.tmp@assays$RNA@var.features[!grepl("^RP[SL]", NPC.tmp@assays$RNA@var.features, ignore.case = F)]
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp@meta.data$Annotation.Immune <- factor(NPC.tmp@meta.data$Annotation.Immune,
                                              levels = c("DC_CCL5+","DC_CCR7+","DC_CD52+","DC_CD83+","DC_CLEC9A+",
                                                         "Mon/Mac_C1Q+","Mon/Mac_CLNK+","Mon/Mac_CTSB+",
                                                         "Mon/Mac_FCGBP+","Mon/Mac_FN1+","Mon/Mac_IL18+","Mon/Mac_PTPRS+",
                                                         "Neutrophil")
)

NPC.tmp <- RunPCA(NPC.tmp)
ElbowPlot(NPC.tmp)
dims_parameter <- 15
NPC.tmp <- RunUMAP(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindNeighbors(NPC.tmp, dims = 1:dims_parameter)
NPC.tmp <- FindClusters(NPC.tmp, resolution = 1)#resoulution越小分的越粗
DimPlot(NPC.tmp, reduction = "umap",label = TRUE)

DimPlot(NPC.tmp, reduction = "umap",label = TRUE,group.by = "orig.ident")

## Figures 4 ab
p.immune.1 <- DimPlot(NPC.tmp, reduction = "umap",group.by = "Annotation.Immune", pt.size = 2) + NoLegend() + NoAxes() + labs(title = NULL)
p.immune.2 <- DimPlot(NPC.tmp, reduction = "umap", group.by = "Annotation.Immune",label=TRUE)

p.id.1 <- DimPlot(NPC.tmp, reduction = "umap",group.by = "orig.ident", pt.size = 2) + NoLegend() + NoAxes() + labs(title = NULL)
p.id.2 <- DimPlot(NPC.tmp, reduction = "umap", group.by = "orig.ident",label=TRUE)

# save as 10inch*10inch pdf
# label d=1.5mm gap=2.5mm

#######################################################################
##########################  Figures 4b.merge ##########################
#######################################################################
##heatmap
## sc
NPC.sc <- subset(NPC.combined, subset = SeqTech %in% c("Single-cell"))
NPC.sc <- subset(NPC.sc, subset = Annotation.Coarse %in% c("Myeloid"))
NPC.sc <- FindVariableFeatures(NPC.sc, selection.method = "vst", nfeatures = 2000)
NPC.sc <- ScaleData(NPC.sc,  vars.to.regress=c("S.Score", "G2M.Score"))

## DEG
DefaultAssay(NPC.sc) <- "RNA"
Idents(NPC.sc) <- 'Annotation.Immune'
markers.sc <- FindAllMarkers(object = NPC.sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10.sc <- markers.sc %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# sn
NPC.sn <- subset(NPC.combined, subset = SeqTech %in% c("Single-nucleus"))
NPC.sn <- subset(NPC.sn, subset = Annotation.Coarse %in% c("Myeloid"))
NPC.sn <- FindVariableFeatures(NPC.sn, selection.method = "vst", nfeatures = 2000)
NPC.sn <- ScaleData(NPC.sn,  vars.to.regress=c("S.Score", "G2M.Score"))

## DEG
DefaultAssay(NPC.sn) <- "RNA"
Idents(NPC.sn) <- 'Annotation.Immune'
markers.sn <- FindAllMarkers(object = NPC.sn, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10.sn <- markers.sn %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

## merge
mmclusters = c("DC_CCL5+","DC_CCR7+","DC_CD52+","DC_CD83+","DC_CLEC9A+",
           "Mon/Mac_C1Q+","Mon/Mac_CLNK+","Mon/Mac_CTSB+",
           "Mon/Mac_FCGBP+","Mon/Mac_FN1+","Mon/Mac_IL18+","Mon/Mac_PTPRS+",
           "Neutrophil")
top10 <- rbind(top10.sc[which(top10.sc$cluster==mmclusters[1]),],
               top10.sc[which(top10.sc$cluster==mmclusters[2]),],
               top10.sc[which(top10.sc$cluster==mmclusters[3]),],
               top10.sn[which(top10.sn$cluster==mmclusters[4]),],
               top10.sc[which(top10.sc$cluster==mmclusters[5]),],
               top10.sc[which(top10.sc$cluster==mmclusters[6]),],
               top10.sn[which(top10.sn$cluster==mmclusters[7]),],
               top10.sn[which(top10.sn$cluster==mmclusters[8]),],
               top10.sn[which(top10.sn$cluster==mmclusters[9]),],
               top10.sc[which(top10.sc$cluster==mmclusters[10]),],
               top10.sn[which(top10.sn$cluster==mmclusters[11]),],
               top10.sn[which(top10.sn$cluster==mmclusters[12]),],
               top10.sc[which(top10.sc$cluster==mmclusters[13]),])

## NPC.tmp
NPC.tmp <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Myeloid"))
NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp@assays$RNA@var.features <- NPC.tmp@assays$RNA@var.features[!grepl("^RP[SL]", NPC.tmp@assays$RNA@var.features, ignore.case = F)]
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.tmp@meta.data$Annotation.Immune <- factor(NPC.tmp@meta.data$Annotation.Immune,
                                              levels = c("DC_CCL5+","DC_CCR7+","DC_CD52+","DC_CD83+","DC_CLEC9A+",
                                                         "Mon/Mac_C1Q+","Mon/Mac_CLNK+","Mon/Mac_CTSB+",
                                                         "Mon/Mac_FCGBP+","Mon/Mac_FN1+","Mon/Mac_IL18+","Mon/Mac_PTPRS+",
                                                         "Neutrophil"))
## DEG
DefaultAssay(NPC.tmp) <- "RNA"
Idents(NPC.tmp) <- 'Annotation.Immune'
markers <- FindAllMarkers(object = NPC.tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


avgexp = AverageExpression(NPC.tmp, return.seurat = T,group.by = "Annotation.Immune")

p<-DoHeatmap(avgexp, features = top10$gene, group.by = "orig.ident", 
             size = 3, angle = 0, combine = TRUE, draw.lines = FALSE) + 
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust=1,colour = "black"),
        axis.text.y = element_text(size = 2,colour = "black")) + 
  scale_fill_gradient2(low = '#1000ff', 
                       mid = "#ffffff", 
                       high = '#aa0101', 
                       space = "Lab", 
                       na.value = "#ffffff", 
                       midpoint = 0, 
                       guide = "colourbar", 
                       aesthetics = "fill") #+ coord_flip(expand = TRUE, clip = "on")

path.save <- paste("G:/manuscript/Figure/HeatmapMM.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p, width=unit(40, "mm"), height=unit(100, "mm")), 
       width = 500, height = 500, units = 'mm', dpi = 300)

#####################################################################
##########################  Figures 4b.sep ##########################
#####################################################################
##heatmap
NPC.tmp <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Myeloid"))
NPC.tmp@meta.data$Annotation.Immune <- factor(NPC.tmp@meta.data$Annotation.Immune,
                                              levels = c("DC_CCL5+","DC_CCR7+","DC_CD52+","DC_CD83+","DC_CLEC9A+",
                                                         "Mon/Mac_C1Q+","Mon/Mac_CLNK+","Mon/Mac_CTSB+",
                                                         "Mon/Mac_FCGBP+","Mon/Mac_FN1+","Mon/Mac_IL18+","Mon/Mac_PTPRS+",
                                                         "Neutrophil"))
## sc
NPC.sc <- subset(NPC.tmp, subset = SeqTech %in% c("Single-cell"))
NPC.sc <- FindVariableFeatures(NPC.sc, selection.method = "vst", nfeatures = 2000)
NPC.sc <- ScaleData(NPC.sc,  vars.to.regress=c("S.Score", "G2M.Score"))

## DEG
DefaultAssay(NPC.sc) <- "RNA"
Idents(NPC.sc) <- 'Annotation.Immune'
markers.sc <- FindAllMarkers(object = NPC.sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10.sc <- markers.sc %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

avgexp.sc = AverageExpression(NPC.sc, return.seurat = T,group.by = "Annotation.Immune")

p.sc<-DoHeatmap(avgexp.sc, features = top10.sc$gene, group.by = "orig.ident", 
             size = 3, angle = 0, combine = TRUE, draw.lines = FALSE) + 
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust=1,colour = "black"),
        axis.text.y = element_text(size = 3,colour = "black")) + 
  scale_fill_gradient2(low = '#1000ff', 
                       mid = "#ffffff", 
                       high = '#aa0101', 
                       space = "Lab", 
                       na.value = "#ffffff", 
                       midpoint = 0, 
                       guide = "colourbar", 
                       aesthetics = "fill") #+ coord_flip(expand = TRUE, clip = "on")
## sn
NPC.sn <- subset(NPC.tmp, subset = SeqTech %in% c("Single-nucleus"))
NPC.sn <- FindVariableFeatures(NPC.sn, selection.method = "vst", nfeatures = 2000)
NPC.sn <- ScaleData(NPC.sn,  vars.to.regress=c("S.Score", "G2M.Score"))

## DEG
DefaultAssay(NPC.sn) <- "RNA"
Idents(NPC.sn) <- 'Annotation.Immune'
markers.sn <- FindAllMarkers(object = NPC.sn, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10.sn <- markers.sn %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

avgexp.sn = AverageExpression(NPC.sn, return.seurat = T,group.by = "Annotation.Immune")

p.sn<-DoHeatmap(avgexp.sn, features = top10.sn$gene, group.by = "orig.ident", 
                size = 3, angle = 0, combine = TRUE, draw.lines = FALSE) + 
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust=1,colour = "black"),
        axis.text.y = element_text(size = 3,colour = "black")) + 
  scale_fill_gradient2(low = '#1000ff', 
                       mid = "#ffffff", 
                       high = '#aa0101', 
                       space = "Lab", 
                       na.value = "#ffffff", 
                       midpoint = 0, 
                       guide = "colourbar", 
                       aesthetics = "fill") #+ coord_flip(expand = TRUE, clip = "on")
## figure


path.save <- paste("G:/manuscript/Figure/HeatmapMM.sc.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p.sc, width=unit(21, "mm"), height=unit(105, "mm")), 
       width = 500, height = 500, units = 'mm', dpi = 300)
path.save <- paste("G:/manuscript/Figure/HeatmapMM.sn.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p.sn, width=unit(18, "mm"), height=unit(90, "mm")), 
       width = 500, height = 500, units = 'mm', dpi = 300)







###########################################################
##########################  DEGs ##########################
###########################################################
NPC.combined@meta.data$source <- substr(NPC.combined@meta.data$orig.ident,1,1)
NPC.tmp <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Myeloid"))
NPC.tmp <- subset(NPC.combined, subset = Annotation.Immune %in% c("Mon/Mac_C1Q+","Mon/Mac_CLNK+","Mon/Mac_CTSB+",
                                                                  "Mon/Mac_FCGBP+","Mon/Mac_FN1+","Mon/Mac_IL18+","Mon/Mac_PTPRS+"))


NPC.sc <- subset(NPC.tmp, subset = SeqTech %in% c("Single-cell"))
NPC.sc <- FindVariableFeatures(NPC.sc, selection.method = "vst", nfeatures = 2000)
NPC.sc <- ScaleData(NPC.sc,  vars.to.regress=c("S.Score", "G2M.Score"))
DefaultAssay(NPC.sc) <- "RNA"
Idents(NPC.sc) <- 'source'
markers.sc <- FindAllMarkers(object = NPC.sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

up.genes.enrich.sc <- subset(markers.sc,subset = cluster %in% c("P") )
up.genes.enrich.sc <- rownames(up.genes.enrich.sc)
up.genes.enrich.sc <- up.genes.enrich.sc[!grepl("^RP[SL]", up.genes.enrich.sc, ignore.case = F)]
down.genes.enrich.sc <- subset(markers.sc,subset = cluster %in% c("N") )
down.genes.enrich.sc <- rownames(down.genes.enrich.sc)
down.genes.enrich.sc <- down.genes.enrich.sc[!grepl("^RP[SL]", down.genes.enrich.sc, ignore.case = F)]

NPC.sn <- subset(NPC.tmp, subset = SeqTech %in% c("Single-nucleus"))
NPC.sn <- FindVariableFeatures(NPC.sn, selection.method = "vst", nfeatures = 2000)
NPC.sn <- ScaleData(NPC.sn,  vars.to.regress=c("S.Score", "G2M.Score"))
DefaultAssay(NPC.sn) <- "RNA"
Idents(NPC.sn) <- 'source'
markers.sn <- FindAllMarkers(object = NPC.sn, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

up.genes.enrich.sn <- subset(markers.sn,subset = cluster %in% c("P") )
up.genes.enrich.sn <- rownames(up.genes.enrich.sn)
up.genes.enrich.sn <- up.genes.enrich.sn[!grepl("^RP[SL]", up.genes.enrich.sn, ignore.case = F)]
down.genes.enrich.sn <- subset(markers.sn,subset = cluster %in% c("N") )
down.genes.enrich.sn <- rownames(down.genes.enrich.sn)
down.genes.enrich.sn <- down.genes.enrich.sn[!grepl("^RP[SL]", down.genes.enrich.sn, ignore.case = F)]


up.genes.enrich <- intersect(up.genes.enrich.sn,up.genes.enrich.sc)
down.genes.enrich <- intersect(down.genes.enrich.sn,down.genes.enrich.sc)

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
barplot(down.genes.go, showCategory=20)

################################################################################
# figure Enrichment
n.top = 5

up.genes.go@result$logg <- (-log(up.genes.go@result$p.adjust))
tmp <- up.genes.go %>% group_by(ONTOLOGY) %>% top_n(n = n.top, wt = logg)
p.up <- ggplot(tmp@result)+
  geom_bar(aes(x = Description,y = -log(p.adjust), fill = ONTOLOGY), stat='identity') + 
  coord_flip() +
  scale_x_discrete(limits=tmp@result$Description,position = "top") +
  #scale_x_discrete(limits=tmp@result$Description,position = "top",labels=function(x) str_wrap(x, width=30)) +
  theme_bw() + #theme(panel.grid=element_blank())+
  #scale_y_continuous(limits=c(0,150))+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"),
  )
down.genes.go@result$logg <- (-log(down.genes.go@result$p.adjust))
tmp <- down.genes.go %>% group_by(ONTOLOGY) %>% top_n(n = n.top, wt = logg)
p.down <- ggplot(tmp@result)+
  geom_bar(aes(x = Description,y = -log(p.adjust), fill = ONTOLOGY), stat='identity') + 
  coord_flip() +
  scale_x_discrete(limits=tmp@result$Description,position = "top") +
  #scale_x_discrete(limits=tmp@result$Description,position = "top",labels=function(x) str_wrap(x, width=30)) +
  theme_bw() + #theme(panel.grid=element_blank())+
  #scale_y_continuous(limits=c(0,150))+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"),
  )

p.up
p.down

path.save <- paste("G:/manuscript/Figure/p.down.Myeloid.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p.down, width=unit(1.2, "in"), height=unit(1.2, "in")), 
       width = 10, height = 10, units = 'in', dpi = 300)
path.save <- paste("G:/manuscript/Figure/p.up.Myeloid.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p.up, width=unit(1.2, "in"), height=unit(1.2, "in")), 
       width = 10, height = 10, units = 'in', dpi = 300)



################################################################################
# VENN

library(ggplot2)
library(ggsci)
library(sf)
library(ggVennDiagram)

color1 <- alpha("#BC3C29FF",1)
color2 <- alpha("#0072B5FF",0.5)


################################################################################
################################################################################
# CNV up-regular genes
up.markers = list()
up.markers[["sc"]] <- up.genes.enrich.sc
up.markers[["sn"]] <- up.genes.enrich.sn

################################################################################
#venn

#label_alpha = 0去除文字标签底色；
#category.names参数用于设定样本名称；
p1 <- ggVennDiagram(up.markers, label_alpha=0, label_size =5, edge_size = 1,label = "count", edge_lty = "dashed") +
  scale_color_brewer(palette = "Paired")+
  #scale_fill_gradient(low="white",high = color1)
  scale_fill_distiller(palette = "Reds", direction = 1)+scale_color_brewer(palette = "Set1")



ggsave("G:/manuscript/Figure/venn.TAM.pdf", egg::set_panel_size(p1, width=unit(140, "mm"), height=unit(100, "mm")), 
       width = 700, height = 700, units = 'mm', dpi = 300)
# resize to 700 to 250


##################################################################################
#########################  Figure 5 M1/M2 Normal vs NPC  #########################
##################################################################################

# https://www.sciencedirect.com/science/article/pii/S1471490602023025?via%3Dihub
M1 <- c("CD80","CD86",
        "FCGR1A","FCGR2B","FCGR3B",#Fcr-RI,II,III
        "TLR2","TLR4",
        "TNF","IL1A","IL1B","IL6","IL12A","IL12B",
        "IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA21",#Type I IFN
        "IL1R1",
        "CXCL9","CXCL10","CXCL11","CCL2","CCL3","CCL4","CCL5","CXCL8",
        "CCR7"
)
M2 <- c("MRC1","CD163","CD68",
        "MSR1",#Scavenger receptor A
        "SCARB1","SCARB2","CD36",#Scavenger receptor B 
        "CD14",
        "FCER2",#FceRII,CD23
        "IL1RN","IL10",
        "IL1R2",
        "CCL17","CCL22","CCL24","CCL18","CCL16",
        "CCR2","CXCR1","CXCR2",
        "TGFB1","TGFB2","TGFB3","IL4"
)
# https://www.nature.com/articles/nrc.2016.54


NPC.tmp <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Myeloid"))
NPC.tmp@meta.data$source <- substr(NPC.tmp@meta.data$orig.ident,1,1)
NPC.tmp@meta.data$source[which(NPC.tmp@meta.data$source=="N")] = "Normal"
NPC.tmp@meta.data$source[which(NPC.tmp@meta.data$source=="P")] = "NPC"

NPC.TAM <- subset(NPC.tmp, subset = Annotation.Immune %in% c("Mon/Mac_C1Q+","Mon/Mac_CLNK+","Mon/Mac_CTSB+",
                                                             "Mon/Mac_FCGBP+","Mon/Mac_FN1+","Mon/Mac_IL18+","Mon/Mac_PTPRS+"))

NPC.sc <- subset(NPC.TAM, subset = orig.ident %in% c("N01","N02","P01","P02","P03"))
NPC.sn <- subset(NPC.TAM, subset = orig.ident %in% c("N03","P04","P05","P06","P07"))

NPC.sc <- FindVariableFeatures(NPC.sc, selection.method = "vst", nfeatures = 2000)
NPC.sc <- ScaleData(NPC.sc,  vars.to.regress=c("S.Score", "G2M.Score"))
NPC.sn <- FindVariableFeatures(NPC.sn, selection.method = "vst", nfeatures = 2000)
NPC.sn <- ScaleData(NPC.sn,  vars.to.regress=c("S.Score", "G2M.Score"))

NPC.sc <- AddModuleScore(object = NPC.sc,features = list(M1),ctrl = 100,name = "score.M1")
NPC.sc <- AddModuleScore(object = NPC.sc,features = list(M2),ctrl = 100,name = "score.M2")
NPC.sn <- AddModuleScore(object = NPC.sn,features = list(M1),ctrl = 100,name = "score.M1")
NPC.sn <- AddModuleScore(object = NPC.sn,features = list(M2),ctrl = 100,name = "score.M2")

NPC.sc$score.M1 <- NPC.sc$score.M11
NPC.sc$score.M2 <- NPC.sc$score.M21
NPC.sn$score.M1 <- NPC.sn$score.M11
NPC.sn$score.M2 <- NPC.sn$score.M21

test.size = 5*5

p1 <- VlnPlot(NPC.sc, features="score.M1",group.by = "source", pt.size = 0)+NoLegend()+ggtitle("M1 Score(scRNA)")+
  theme(axis.text.x = element_text(size=test.size,color="black",face="plain"),
        axis.text.y = element_text(size=test.size,color="black",face="plain"))+
  theme(legend.text = element_text(size = test.size),
        legend.title = element_text(size = test.size),
        plot.title = element_text(size = test.size))+
  geom_boxplot(width=0.2,fill="white",outlier.size=-1)
p2 <- VlnPlot(NPC.sc, features="score.M2",group.by = "source", pt.size = 0)+NoLegend()+ggtitle("M2 Score(scRNA)")+
  theme(axis.text.x = element_text(size=test.size,color="black",face="plain"),
        axis.text.y = element_text(size=test.size,color="black",face="plain"))+
  theme(legend.text = element_text(size = test.size),
        legend.title = element_text(size = test.size),
        plot.title = element_text(size = test.size))+
  geom_boxplot(width=0.2,fill="white",outlier.size=-1)
p3 <- VlnPlot(NPC.sn, features="score.M1",group.by = "source", pt.size = 0)+NoLegend()+ggtitle("M1 Score(snRNA)")+
  theme(axis.text.x = element_text(size=test.size,color="black",face="plain"),
        axis.text.y = element_text(size=test.size,color="black",face="plain"))+
  theme(legend.text = element_text(size = test.size),
        legend.title = element_text(size = test.size),
        plot.title = element_text(size = test.size))+
  geom_boxplot(width=0.2,fill="white",outlier.size=-1)
p4 <- VlnPlot(NPC.sn, features="score.M2",group.by = "source", pt.size = 0)+NoLegend()+ggtitle("M2 Score(snRNA)")+
  theme(axis.text.x = element_text(size=test.size,color="black",face="plain"),
        axis.text.y = element_text(size=test.size,color="black",face="plain"))+
  theme(legend.text = element_text(size = test.size),
        legend.title = element_text(size = test.size),
        plot.title = element_text(size = test.size))+
  geom_boxplot(width=0.2,fill="white",outlier.size=-1)

ggsave("G:/manuscript/Figure/sc.M1.pdf", egg::set_panel_size(p1, width=unit(10*5, "mm"), height=unit(25*5, "mm")), 
       width = 500, height = 500, units = 'mm', dpi = 300)
ggsave("G:/manuscript/Figure/sc.M2.pdf", egg::set_panel_size(p2, width=unit(10*5, "mm"), height=unit(25*5, "mm")), 
       width = 500, height = 500, units = 'mm', dpi = 300)
ggsave("G:/manuscript/Figure/sn.M1.pdf", egg::set_panel_size(p3, width=unit(10*5, "mm"), height=unit(25*5, "mm")), 
       width = 500, height = 500, units = 'mm', dpi = 300)
ggsave("G:/manuscript/Figure/sn.M2.pdf", egg::set_panel_size(p4, width=unit(10*5, "mm"), height=unit(25*5, "mm")), 
       width = 500, height = 500, units = 'mm', dpi = 300)
# resize to 1/5


# wilcox.test
datalist <- c("single cell","single nucleus")
p.result <- data.frame(
  sample = datalist,
  p.M1=c(0,0),
  p.M2=c(0,0),
  mean.M1.n=c(0,0),
  mean.M2.n=c(0,0),
  mean.M1.p=c(0,0),
  mean.M2.p=c(0,0)
)
for (i in datalist){
  if (i %in% c("single cell")) {
    xx <- subset(NPC.sc, subset = source %in% c("Normal"))
    tmp <- NPC.sc
  } else {
    xx <- subset(NPC.sn, subset = source %in% c("Normal"))
    tmp <- NPC.sn
  }
  
  x.1 <- xx$score.M1
  x.2 <- xx$score.M2
  yy <- subset(tmp, subset = source %in% "NPC")
  y.1 <- yy$score.M1
  y.2 <- yy$score.M2
  tt.1 <- wilcox.test(x = x.1 , y = y.1, alternative = "two.sided")
  tt.2 <- wilcox.test(x = x.2 , y = y.2, alternative = "two.sided")
  
  p.result$p.M1[match(i,datalist)] <- tt.1$p.value
  p.result$p.M2[match(i,datalist)] <- tt.2$p.value
  p.result$mean.M1.n[match(i,datalist)] <- mean(x.1)
  p.result$mean.M2.n[match(i,datalist)] <- mean(x.2)
  p.result$mean.M1.p[match(i,datalist)] <- mean(y.1)
  p.result$mean.M2.p[match(i,datalist)] <- mean(y.2)
}

# * 0.05
# ** 0.01
# *** 0.001
# P01 M1 * M2 


# x-y plot
tmp.sc <- data.frame(M1=NPC.sc$score.M1,M2=NPC.sc$score.M2,group = NPC.sc$orig.ident)
tmp.sn <- data.frame(M1=NPC.sn$score.M1,M2=NPC.sn$score.M2,group = NPC.sn$orig.ident)

p1 <- ggplot(data=tmp.sc, aes(x= M1, y= M2,color=group))+ geom_point()
p2 <- ggplot(data=tmp.sn, aes(x= M1, y= M2,color=group))+ geom_point()
p1+p2


##heatmap
tmp <- NPC.tmp
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'Annotation.Immune'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(NPC.tmp, features = top10$gene,group.by = "Annotation.Immune")

######################################################################
#########################  Figure CD86/CD68  #########################
######################################################################
## Myeloid
NPC.tmp <- subset(NPC.combined, subset = Annotation.Coarse %in% c("Myeloid"))
NPC.tmp <- subset(NPC.tmp, subset = Annotation.Immune %in% c("Mon/Mac_C1Q+","Mon/Mac_CLNK+","Mon/Mac_CTSB+",
                                                             "Mon/Mac_FCGBP+","Mon/Mac_FN1+","Mon/Mac_IL18+","Mon/Mac_PTPRS+"))

NPC.tmp <- FindVariableFeatures(NPC.tmp, selection.method = "vst", nfeatures = 2000)
NPC.tmp@assays$RNA@var.features <- NPC.tmp@assays$RNA@var.features[!grepl("^RP[SL]", NPC.tmp@assays$RNA@var.features, ignore.case = F)]
NPC.tmp <- ScaleData(NPC.tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

genes <- c("CD86","CD68")
tt = 5

for (i in genes){
p <- VlnPlot(NPC.tmp, features=i,group.by = "Annotation.Immune", pt.size = 0)+NoLegend()+
  theme(axis.text.x = element_text(size=tt*5,color="black",face="plain"),
        axis.text.y = element_text(size=tt*5,color="black",face="plain"))+
  theme(legend.text = element_text(size = tt*5),
        legend.title = element_text(size = tt*5),
        plot.title = element_text(size = tt*5))

ppath <- paste("G:/manuscript/Figure/",i,".pdf",sep="")
ggsave(ppath, egg::set_panel_size(p, width=unit(25*tt, "mm"), height=unit(25*tt, "mm")), 
       width = 500, height = 500, units = 'mm', dpi = 300)

}
# resize to 1/tt

############################################################################
#########################  Figure GO of cell type  #########################
############################################################################
##heatmap
## sc
NPC.sc <- subset(NPC.combined, subset = SeqTech %in% c("Single-cell"))
NPC.sc <- subset(NPC.sc, subset = Annotation.Coarse %in% c("Myeloid"))
NPC.sc <- FindVariableFeatures(NPC.sc, selection.method = "vst", nfeatures = 2000)
NPC.sc <- ScaleData(NPC.sc,  vars.to.regress=c("S.Score", "G2M.Score"))

## DEG
DefaultAssay(NPC.sc) <- "RNA"
Idents(NPC.sc) <- 'Annotation.Immune'
markers.sc <- FindAllMarkers(object = NPC.sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)

genes <- markers.sc[which(markers.sc$cluster=="DC_CCL5+"),]

## Go
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(stringr)

up.genes.enrich <- genes$gene
up.genes.enrich <- bitr(up.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
up.genes.go <- enrichGO(gene          = up.genes.enrich$SYMBOL,
                        #universe     = row.names(dge.celltype),
                        OrgDb         = 'org.Hs.eg.db',
                        keyType       = 'SYMBOL',
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.05)
barplot(up.genes.go, showCategory=20)

################################################################################
# figure Enrichment
n.top = 5

up.genes.go@result$logg <- (-log(up.genes.go@result$p.adjust))
tmp <- up.genes.go %>% group_by(ONTOLOGY) %>% top_n(n = n.top, wt = logg)
p.up <- ggplot(tmp@result)+
  geom_bar(aes(x = Description,y = -log(p.adjust), fill = ONTOLOGY), stat='identity') + 
  coord_flip() +
  scale_x_discrete(limits=tmp@result$Description,position = "top") +
  #scale_x_discrete(limits=tmp@result$Description,position = "top",labels=function(x) str_wrap(x, width=30)) +
  theme_bw() + #theme(panel.grid=element_blank())+
  #scale_y_continuous(limits=c(0,150))+
  theme(axis.text.x = element_text(size=5,color="black",face="plain"),
        axis.text.y = element_text(size=5,color="black",face="plain"),
  )

p.up

path.save <- paste("G:/manuscript/Figure/p.up.xxMyeloid.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p.up, width=unit(1.2, "in"), height=unit(1.2, "in")), 
       width = 10, height = 10, units = 'in', dpi = 300)


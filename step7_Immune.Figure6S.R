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
load("G:/scRNA data/DS0/step7.AnnotationImmune.RData")

############################################################
######################### Figure a #########################
############################################################

p1 <- VlnPlot(NPC.tmp, features = c("nFeature_RNA"), group.by = "Annotation.Immune2",pt.size = 0) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 5, colour = "black"),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size=0),
    axis.title.y = element_text(size=6)
  )
p2 <- VlnPlot(NPC.tmp, features = c("nCount_RNA"), group.by = "Annotation.Immune2",pt.size = 0) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 5, colour = "black"),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size=0),
    axis.title.y = element_text(size=6)
  )
p3 <- VlnPlot(NPC.tmp, features = c("percent.mt"), group.by = "Annotation.Immune2",pt.size = 0) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 5, colour = "black"),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size=0),
    axis.title.y = element_text(size=6)
  )
path.save <- paste("G:/manuscript/Figure/BT_Violin1.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p1, width=unit(40, "mm"), height=unit(30, "mm")), 
       width = 100, height = 100, units = 'mm', dpi = 300)
path.save <- paste("G:/manuscript/Figure/BT_Violin2.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p2, width=unit(40, "mm"), height=unit(30, "mm")), 
       width = 100, height = 100, units = 'mm', dpi = 300)
path.save <- paste("G:/manuscript/Figure/BT_Violin3.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p3, width=unit(40, "mm"), height=unit(40, "mm")), 
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

DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'Annotation.Immune2'
#markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
markers.1 <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)



## NK_CD8+ vs T_CD8+ in both Patient and Normal Group
tmp <- subset(NPC.tmp, subset = Annotation.Immune2 %in% c("T_CD8+","NK_CD8+"))
tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
tmp <- ScaleData(tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'Annotation.Immune2'
#markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
markers.2 <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)

markers <- intersect(markers.1$gene,markers.2$gene)


## Go
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(stringr)

up.genes.enrich <- markers
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




n.top = 10

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


#############################################################
######################### huoshantu #########################
#############################################################

VolcanoPlot=function(dif, log2FC=log2(1.5), padj=0.05, 
                     label.symbols=NULL, label.max=30,
                     cols=c("#497aa2", "#ae3137"), title=""){
  if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(dif) )){
    stop("Colnames must include: log2FoldChange, padj, symbol")
  }
  rownames(dif)=dif$symbol
  
  # (1) define up and down
  dif$threshold="ns";
  dif[which(dif$log2FoldChange > log2FC & dif$padj <padj),]$threshold="up";
  dif[which(dif$log2FoldChange < (-log2FC) & dif$padj < padj),]$threshold="down";
  dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
  #head(dif)
  #
  tb2=table(dif$threshold); print(tb2)
  library(ggplot2)
  # (2) plot
  g1 = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
    geom_point(alpha=0.8, size=0.8) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    labs(title= ifelse(""==title, "", paste("DEG:", title)))+
    xlab(bquote(Log[2]*FoldChange))+
    ylab(bquote(-Log[10]*italic(P.adj)) )+
    theme_classic(base_size = 14) +
    theme(legend.box = "horizontal",
          legend.position="top",
          legend.spacing.x = unit(0, 'pt'),
          legend.text = element_text( margin = margin(r = 20) ),
          legend.margin=margin(b= -10, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size=11)
    ) +
    scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns',
                                   paste0("up(",tb2[[3]],')' )),
                       values=c(cols[1], "grey", cols[2]) )+
    guides(color=guide_legend(override.aes = list(size=3, alpha=1))); g1;
  # (3)label genes
  if(is.null(label.symbols)){
    dif.sig=dif[which(dif$threshold != "ns" ), ]
    len=nrow(dif.sig)
    if(len<label.max){
      label.symbols=rownames(dif.sig)
    }else{
      dif.sig=dif.sig[order(dif.sig$log2FoldChange), ]
      dif.sig= rbind(dif.sig[1:(label.max/2),], dif.sig[(len-label.max/2):len,])
      label.symbols=rownames(dif.sig)
    }
  }
  dd_text = dif[label.symbols, ]
  print((dd_text))
  # add text
  library(ggrepel)
  g1 + geom_text_repel(data=dd_text,
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)),
                       #size=2.5, 
                       colour="black",alpha=1)
}


######################################################################
######################### NK_CD8- vs NK_CD8+ #########################
######################################################################

tmp <- subset(NPC.tmp, subset = Annotation.Immune2 %in% c("T_CD8+","NK_CD8+"))
tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
tmp <- ScaleData(tmp,  vars.to.regress=c("S.Score", "G2M.Score"))

##DEG
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'Annotation.Immune2'

deg <- FindMarkers(tmp,
                   ident.1 = 'NK_CD8+',
                   ident.2 = 'T_CD8+',
                   group.by="Annotation.Immune2")

dim(deg)
head(deg)

## 中间数据
dif=data.frame(
  symbol=rownames(deg),
  log2FoldChange=deg$avg_log2FC,
  padj=deg$p_val_adj
)

## FIGURE
# 可以指定要标记的DEG数量，选出FC最大和最小的基因标记
p1 <- VolcanoPlot(dif, padj=0.05, title="T_CD8+ vs NK_CD8+", label.max = 50)
# 自定义颜色
VolcanoPlot(dif, padj=0.05, title="T_CD8+ vs NK_CD8+", label.max = 50, cols=c("blue", "red"))


# 也可以指定要标记的基因名字
VolcanoPlot(dif, padj=1e-10, title="T_CD8+ vs NK_CD8+", 
            label.symbols=dif[ ((abs(dif$log2FoldChange) > 2) & (dif$padj < 1e-50) ) | 
                                 abs(dif$log2FoldChange) > 4,]$symbol )



p1 <- VolcanoPlot(dif, padj=0.05, title="T_CD8+ vs NK_CD8+", label.max = 20)
path.save <- paste("G:/manuscript/Figure/hst.pdf",sep="")
ggsave(path.save, egg::set_panel_size(p1, width=unit(66, "mm"), height=unit(66, "mm")), 
       width = 110, height = 110, units = 'mm', dpi = 300)

#resize to 5/11


rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(rhdf5)

library(scater)
library(patchwork)


load("G:/scRNA data/DS0/step1.NPC.Merge.RData")
load("G:/scRNA data/DS0/step2.Annotation.RData")
NPC.combined@meta.data$Annotation.Coarse <- Annotation.Coarse

############################################################
#########################Figure 1b##########################
############################################################

p1 <- DimPlot(NPC.combined, reduction = "umap", group.by = 'Annotation.Coarse', raster=FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
p2 <- DimPlot(NPC.combined, reduction = "umap", group.by = 'Annotation.Coarse',label=TRUE)

# save p1 as 8inch*8inch pdf
# save p2 as 8inch*8inch pdf
# Transfer p1 to semi-png figure
# Adjust p1 to 50mm*50mm
# Add labels and legend from p2 to p1
# Adjust word to 5pt, label gap 2mm, square-word left-left 2.5mm
# Add Axes

############################################################
#########################Figure 1c##########################
############################################################

p1 <- DimPlot(NPC.combined, reduction = "umap", group.by = 'SeqTech', raster=FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
p2 <- DimPlot(NPC.combined, reduction = "umap", group.by = 'SeqTech',label=TRUE)

# save p1 as 8inch*8inch pdf
# save p2 as 8inch*8inch pdf
# Transfer p1 to semi-png figure
# Adjust p1 to 50mm*50mm
# Add labels and legend from p2 to p1
# Adjust word to 5pt, label gap 2mm, square-word left-left 2.5mm
# Add Axes

############################################################
#########################Figure 1d##########################
############################################################
Idents(NPC.combined) <- 'Annotation.Coarse'
pt <- table(Idents(NPC.combined), NPC.combined$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

p1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  #theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("Sample ID") +
  ylab("Proportion") +
  theme(legend.position="none") +
  #scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank()) # + coord_flip()


# save p1 as 5inch*5inch pdf
# Adjust p1 to 40mm*40mm


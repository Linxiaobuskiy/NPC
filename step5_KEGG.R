# kegg
tmp <- NPC.tmp
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'Annotation.Fib'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
new.cluster.ids <- c("PDGFRA+/NOVA1+","PDGFRA-/SAA+","PDGFRA+/NOVA1-","PDGFRA-/SAA-")


up.genes.enrich <- subset(markers,subset = cluster %in% new.cluster.ids[4])
up.genes.enrich <- rownames(up.genes.enrich)
up.genes.enrich <- up.genes.enrich[!grepl("^RP[SL]", up.genes.enrich, ignore.case = F)]
up.genes.enrich <- bitr(up.genes.enrich, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
up.genes.kegg <- enrichKEGG(gene = up.genes.enrich$ENTREZID,
                            organism = 'hsa',
                            pvalueCutoff = 0.05)
head(up.genes.kegg)

p1.up.kegg <- barplot(up.genes.kegg, showCategory=20)
p2.up.kegg <- dotplot(up.genes.kegg, showCategory=20)
plotc.up.kegg = p1.up.kegg/p2.up.kegg
plotc.up.kegg
library(ggpubr)


# https://www.biostars.org/p/458261/

## Follow the Seurat tutorial until UMAP at: https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

## create fake 'Reponse' variable

#pbmc@meta.data$Response <- rep(c("CTRL", "TRT"), times = c(1500, 2638-1500))

#vp_case1 <- function(gene_signature, test_sign){
vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(NPC.tmp, features = signature,
            pt.size = 0, 
            group.by = "Annotation.Fib", 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 14, height = 8)
}

# gene_sig <- c("PDGFRA","PECAM1")
# comparisons <- list(c("mCAF", "iCAF"),c("mCAF", "vCAF"), c("iCAF", "vCAF"))
# vp_case1(gene_signature = gene_sig,  file_name ="G:/manuscript/Figure/a", test_sign = comparisons)

gene_sig <- c("PDGFRA","PECAM1")
comparisons <- list(c("mCAF", "iCAF"), c("mCAF", "vCAF"))
vp_case1(gene_signature = gene_sig,  file_name ="G:/manuscript/Figure/mCAF", test_sign = comparisons)

comparisons <- list(c("iCAF", "mCAF"), c("iCAF", "vCAF"))
vp_case1(gene_signature = gene_sig,  file_name ="G:/manuscript/Figure/iCAF", test_sign = comparisons)

comparisons <- list(c("vCAF", "mCAF"), c("vCAF", "iCAF"))
vp_case1(gene_signature = gene_sig,  file_name ="G:/manuscript/Figure/vCAF", test_sign = comparisons)
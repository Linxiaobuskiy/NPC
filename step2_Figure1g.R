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

a <- table(NPC.combined@meta.data$Annotation.Coarse,NPC.combined@meta.data$orig.ident)

IMC <- c(0.561330876,
         0.525869704,
         0.541477671,
         0.276974491,
         0.360451357,
         0.436475339,
         0.906196575)
IMC <- IMC * 100

sc <- c(0.003242242,
        0.009082533,
        0.059461166)
sc <- sc * 100

sn <- c(0.895272041,
        0.953079534,
        0.323795668,
        0.482403787)
sn <- sn * 100

sc1 <- c(0.827244023,
         0.019041366,
         0.048899756,
         0.2875,
         0.092290988,
         0.104860486,
         0,
         0,
         0,
         0,
         0.762842796,
         0.558516196,
         0.594690967,
         0.085702843,
         0.010195975)
sc1 <- sc1 * 100

sc2 <- c(0.049859769,
         0.015996827,
         0.00308642,
         0.006747086,
         0.001941936,
         0.019647778,
         0.006415668,
         0.075375329,
         0.032997615,
         0.153950768)
sc2 <- sc2 * 100

IMC <- as.data.frame(IMC)
IMC$group <- "IMC"
colnames(IMC) <- c('Epi%','Dataset')
sc <- as.data.frame(sc)
sc$group <- "Single-cell"
colnames(sc) <- c('Epi%','Dataset')
sn <- as.data.frame(sn)
sn$group <- "Single-nucleus"
colnames(sn) <- c('Epi%','Dataset')
sc1 <- as.data.frame(sc1)
sc1$group <- "Yu-Pei Chen etc."
colnames(sc1) <- c('Epi%','Dataset')
sc2 <- as.data.frame(sc2)
sc2$group <- "Yang Liu etc."
colnames(sc2) <- c('Epi%','Dataset')




library(ggplot2)

test <- rbind(IMC,sc,sn,sc1,sc2)
test$Dataset <- factor(test$Dataset, levels = c("IMC","Single-cell","Single-nucleus","Yu-Pei Chen etc.","Yang Liu etc."))
test$pair <- 1:39
test$pair[8:14] <- 1:7

p <-
test %>%
  ggplot(aes(Dataset,`Epi%`, fill=Dataset)) +
  geom_violin(alpha = 1) +
  geom_line(aes(group=pair), position = position_dodge(0.2),color='gray', alpha=0.6) +
  geom_point(aes(fill=Dataset,group=pair), position = position_dodge(0.2),size = 0.2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 5, colour = "black"),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    )

ggsave("G:/manuscript/Figure/SeqTechCompare.pdf", egg::set_panel_size(p, width=unit(80, "mm"), height=unit(20, "mm")), 
       width = 120, height = 120, units = 'mm', dpi = 300)


# wilcox.test
tt.1 <- wilcox.test(x = sc$`Epi%` , y = IMC$`Epi%`, alternative = "less")
tt.2 <- wilcox.test(x = sn$`Epi%` , y = IMC$`Epi%`, alternative = "less")
tt.3 <- wilcox.test(x = sc1$`Epi%` , y = IMC$`Epi%`, alternative = "less")
tt.4 <- wilcox.test(x = sc2$`Epi%` , y = IMC$`Epi%`, alternative = "less")



wilcox.test(x = c(1.03,1.01,1.1,1.05,0.95,1,1.11,1.12,1.112,1.13) , y = c(5,6,7,8,9,10), alternative = "two.sided")







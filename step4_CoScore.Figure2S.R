library("ggpubr")

NPC.tmp <- subset(NPC.sc, subset = orig.ident %in% c("P01","P02","P03"))
x <- NPC.tmp@meta.data$Mitochondrion
y <- NPC.tmp@meta.data$Hypoxia
cor(x, y, method = c("pearson", "kendall", "spearman"))

df <- data.frame(NPC.tmp@meta.data$Mitochondrion,NPC.tmp@meta.data$Hypoxia,NPC.tmp@meta.data$order1)
ggplot(df, aes(x=NPC.tmp.meta.data.Mitochondrion, y=NPC.tmp.meta.data.Hypoxia)) + geom_point()




NPC.tmp <- subset(NPC.sn, subset = orig.ident %in% c("P04","P05","P06","P07"))
x <- NPC.tmp@meta.data$Mitochondrion
y <- NPC.tmp@meta.data$Hypoxia
cor(x, y, method = c("pearson", "kendall", "spearman"))


df <- data.frame(NPC.tmp@meta.data$Mitochondrion,NPC.tmp@meta.data$Hypoxia,NPC.tmp@meta.data$order1)
ggplot(df, aes(x=NPC.tmp.meta.data.Mitochondrion, y=NPC.tmp.meta.data.Hypoxia)) + geom_point()










NPC.tmp <- subset(NPC.figure, subset = orig.ident %in% c("P01","P02","P03","P04","P05","P06","P07"))
x <- NPC.tmp@meta.data$Mitochondrion
y <- NPC.tmp@meta.data$G2M.Score
cor(x, y, method = c("pearson", "kendall", "spearman"))


df <- data.frame(
  NPC.tmp@meta.data$Mitochondrion,
  NPC.tmp@meta.data$Hypoxia,
  NPC.tmp@meta.data$order1,
  NPC.tmp@meta.data$S.Score,
  NPC.tmp@meta.data$G2M.Score)
ggplot(df, aes(x=NPC.tmp.meta.data.Mitochondrion, y=NPC.tmp.meta.data.Hypoxia)) + geom_point()
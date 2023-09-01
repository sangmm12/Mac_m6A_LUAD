setwd("~/projects/m6A_LUAD/result/cor1/DE")


seurat_obj_path <-
  "~/projects/m6A_LUAD/result/cor1/seurat_all.rds"

seurat_obj <- readRDS(seurat_obj_path)



gene_DE <- c("BCL2A1","DUSP2","FTL","HLA-E","HMGB2","ID2","JUNB","LGMN","LYZ","MALAT1","NFKB1","NR4A2","NR4A3","PIM3","PLA2G7","ZFP36","ZNF331")




data <- seurat_obj@assays$RNA@data[match(gene_DE,rownames(seurat_obj)),]
data <- as.matrix(data)

data1 <- reshape2::melt(data)
data1 <- as.data.frame(data1)
data1$Type <- as.character(data1$Var2)

data1$Type[grep("LUNG_N",data1$Type)] <-  "Normal"
data1$Type[grep("LUNG_T",data1$Type)] <- "Tumor"



#### 加载包 ----
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
consistentcolors = colors <- c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", 
                                        "#005AC8", "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", 
                                        "#FA7850", "#14D2DC", "#FA78FA")


pdf("DE_gene.pdf", width = 12, 
    height = 5)
p <- ggboxplot(data1, x = "Var1", y = "value", color = "Type", 
          add = "mean_se") + stat_compare_means(aes(group = Type), label = "p.signif", 
                                                              method = "wilcox.test", size = 2, label.y = 8) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])

p <- p+xlab("Gene")+ylab("Gene Expression")
p
dev.off()







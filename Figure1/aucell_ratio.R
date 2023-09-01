setwd("~/projects/m6A_LUAD/result/ratio")

source("~/projects/R_code/custom_function.R")

seurat_obj_path <-
  "~/projects/m6A_LUAD/result/cor2/seurat_all.rds"

selected_genes <-  c('METTL3',"METTL14","METTL16","WTAP","VIRMA","ZC3H13","CBLL1","RBM15","RBM15B","KIAA1429",
                     "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","HNRNPC","HNRNPG","RBMX","FMR1","HNRNPA2B1","IGFBP1","IGFBP2","IGFBP3","PRRC2A",
                     "FTO","ALKBH5")
# Ready to analyze cell type?
selected_cell_type <- "Macrophage"


k <- 50
# Output path?
outdir <- "./"


#### Preprocessing ----

seurat_obj <- seurat_obj_path %>% readRDS()
metacell_obj <- seurat_obj %>% hdWGCNA::GetMetacellObject()
metacell_obj %>% dim()
cell_types <- names(table(metacell_obj$cell_type))
selected_genes <-
  selected_genes[selected_genes %in% rownames(metacell_obj)]


#### Correlations between genes ----
all_genes <- rownames(metacell_obj)


genesets <- selected_genes


cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data) 

cells_AUC <- AUCell_calcAUC(genesets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

seurat_obj$m6A <- as.numeric(getAUC(cells_AUC))



data <- data.frame(aucell=seurat_obj$m6A,
                   Type=seurat_obj$Type)


library(ggplot2)
library(pacman)
pacman::p_load(tidyverse,ggpubr,rstatix,ggsci,ggsignif,reshape2)


pdf("aucell_Ratio_sample.pdf",width=4,height=5)
ggplot(data,aes(Type,aucell,fill=Type)) + 
  geom_boxplot()+
  scale_fill_jco()+
  geom_jitter(shape=16,size=2,position=position_jitter(0.2))+
  geom_signif(comparisons = list(c("Normal", "Tumor")),
              map_signif_level=T,
              textsize=6,test=wilcox.test,step_increase=0.2)+
  guides(fill=F)+xlab("Type")+ ylab("Aucell")+ theme_classic()  # , title = "Genus"          


dev.off()


pdf("aucell_Ratio_sample.pdf",width=4,height=5)
p <- ggboxplot(data, x = "Type", y = "aucell",
               color = "Type", palette = "jco", 
               add = "mean_se") # palette可以按照期刊选择相应的配色，如"npg"等
p <- p+xlab("Type")+ylab("Percentage")
p <- p + stat_compare_means(aes(group = Type), label = "p.format")
p+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

dev.off()




library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
consistentcolors = colors <- c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", 
                                        "#005AC8", "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", 
                                        "#FA7850", "#14D2DC", "#FA78FA")



ggboxplot(data, x = "Type", y = "aucell", color = "Type", 
          add = "mean_se") + ylim(0, 0.4) + stat_compare_means(aes(group = Type), label = "p.format", 
                                                              method = "wilcox.test", size = 2, label.y = 0.55) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])









setwd("~/projects/m6A_LUAD/result/ratio")


library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratObject)

data <- readRDS("~/projects/m6A_LUAD/object.integrated.rds")




Cell_type <- c("0"="Others",
               "1"="Macrophage",
               "2"="Others",
               "3"="Others",
               "4"="Others",
               "5"="Others",
               "6"="Others",
               "7"="Others",
               "8"="Others",
               "9"="Others",
               "10"="Others",
               "11"="Macrophage",
               "12"="Others",
               "13"="Others")


data[['cell_type']] <- unname(Cell_type[data@meta.data$seurat_clusters])




Idents(data) <- "cell_type"


table(data$cell_type)


#table(Idents(data),data_Mac$cell_type)


data$Type <- data$orig.ident
data$Type[grep("LUNG_N",data$orig.ident)] <- "Normal"
data$Type[grep("LUNG_T",data$orig.ident)] <- "Tumor"




Cellratio <- prop.table(table(Idents(data),data$orig.ident),margin=2)


Cellratio <- as.data.frame(Cellratio[1,])



colnames(Cellratio)[1] <- c("celltype")
colourCount <- length(unique(Cellratio$Var1))


Cellratio$Type <- rownames(Cellratio)
Cellratio$Type[grep("LUNG_N",Cellratio$Type)] <- "Normal"
Cellratio$Type[grep("LUNG_T",Cellratio$Type)] <- "Tumor"

colnames(Cellratio)[1] <- "Value"


library(ggplot2)
library(pacman)
pacman::p_load(tidyverse,ggpubr,rstatix,ggsci,ggsignif,reshape2)


pdf("cell_Ratio_sample.pdf",width=4,height=5)
ggplot(Cellratio,aes(Type,Value,fill=Type)) + 
  geom_boxplot()+
  scale_fill_jco()+
  geom_jitter(shape=16,size=2,position=position_jitter(0.2))+
  geom_signif(comparisons = list(c("Normal", "Tumor")),
              map_signif_level=T,
              textsize=6,test=wilcox.test,step_increase=0.2)+
  guides(fill=F)+xlab("Type")+ ylab("Percentage")+ theme_classic()  # , title = "Genus"          


dev.off()







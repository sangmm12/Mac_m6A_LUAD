#Figure 1


library(Seurat)

data <- readRDS("GSE131907_Lung_Cancer_raw_UMI_matrix.rds")

library(data.table)
cell_anno <- fread("GSE131907_Lung_Cancer_cell_annotation.txt.gz")


Sample_name <- unique(cell_anno$Sample)


for(i in 1:length(Sample_name)){
  cell_name <- cell_anno$Index[which(cell_anno$Sample==Sample_name[i])]
  n1 <- match(cell_name,colnames(data))
  dat <- data[,n1]
  name1 <- paste0("data/",Sample_name[i],".rds")
  
  saveRDS(dat,file=name1)
}


#########################



ana_Sample <- Sample_name[grep("LUNG",Sample_name)]

scRNAlist <- list()
for(i in 1:length(ana_Sample)){
  name <- paste0("data/",ana_Sample[i],".rds")
  counts <- readRDS(name)
  
  scRNAlist[[i]] <-
    CreateSeuratObject(
      counts,
      project = ana_Sample[i],
      min.cells = 5,
      min.features = 200,
    )
  scRNAlist[[i]] <- PercentageFeatureSet( scRNAlist[[i]], pattern = "^mt-", col.name = "percent.mt")
}


obj <- merge(scRNAlist[[1]],scRNAlist[-1])


obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5 )


sce.all.list <- SplitObject(obj, split.by = "orig.ident")

sce.all.list <- lapply(X = sce.all.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1500)
})





select.integrate.feature=1500
pc.num=15
resolution=0.12

remove(data_all)
object.features <- SelectIntegrationFeatures(object.list = sce.all.list, nfeatures = select.integrate.feature)

object.anchors <- FindIntegrationAnchors(object.list = sce.all.list,dims=1:15,
                                         anchor.features = object.features)
object.integrated <- IntegrateData(anchorset = object.anchors,dims=1:15)


remove(sce.all.list)
remove(object.anchors)
DefaultAssay(object.integrated) <- "integrated"


object.integrated <- ScaleData(object.integrated, verbose = FALSE)
object.integrated <- RunPCA(object.integrated, npcs = 15, verbose = FALSE)
object.integrated <- RunUMAP(object.integrated, reduction = "pca", dims = 1:15)
object.integrated <- RunTSNE(object.integrated, dims = 1:15, verbose = FALSE)
object.integrated <- FindNeighbors(object.integrated, reduction = "pca", dims = 1:15)
object.integrated<- FindClusters(object.integrated, resolution = 0.12)

saveRDS(object.integrated, file = "object.integrated.rds")

p1 <- DimPlot(object.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object.integrated, reduction = "tsne", group.by = "orig.ident")
p3 <- DimPlot(object.integrated, label = TRUE, reduction="umap") + NoLegend()
p4 <- DimPlot(object.integrated, label = TRUE, reduction="tsne") + NoLegend()

#merged clusters with label
p5 <- DimPlot(object.integrated, split.by="orig.ident",reduction="umap",ncol=2)
p6 <- DimPlot(object.integrated, split.by="orig.ident",reduction="tsne",ncol=2)
#merged clusters without label
p7 <- DimPlot(object.integrated, split.by="orig.ident",reduction="umap",label = TRUE ,ncol=2) + NoLegend()
p8 <- DimPlot(object.integrated, split.by="orig.ident",reduction="tsne",label = TRUE, ncol=2) + NoLegend()


################## Figure1 A ######################

Idents(data) <- c("Cell_subtype")
pdf("Cell_subtype.pdf",height = 10,width = 10)
DimPlot(data, reduction = "tsne",raster=FALSE)
dev.off()




#######################



library(Seurat)
library(tidyverse)



cat_construct_metacells <-
  function(seurat_obj,
           k,
           name,
           min_cells = 50,
           max_shared = 10,
           target_metacells = 1000,
           assay = "RNA") {
    if ("cell_type" %in% colnames(seurat_obj@meta.data) &
        "orig.ident" %in% colnames(seurat_obj@meta.data)) {
      set.seed(717)
      #### Set up Seurat object for WGCNA ----
      seurat_obj <- hdWGCNA::SetupForWGCNA(
        seurat_obj,
        gene_select = "fraction",
        # the gene selection approach
        fraction = 0.05,
        # fraction of cells that a gene needs to be expressed in order to be included
        wgcna_name = name # the name of the hdWGCNA experiment
      )
      
      #### Construct metacells ----
      # construct metacells  in each group
      seurat_obj <- hdWGCNA::MetacellsByGroups(
        seurat_obj = seurat_obj,
        group.by = c("orig.ident","cell_type" ),
        # specify the columns in seurat_obj@meta.data to group by
        k = k,
        min_cells = min_cells,
        max_shared = max_shared,
        mode = "sum",
        assay = assay, target_metacells = target_metacells,
        # nearest-neighbors parameter
        ident.group = "cell_type" # set the Idents of the metacell seurat object
      )
      # normalize metacell expression matrix:
      seurat_obj <- hdWGCNA::NormalizeMetacells(seurat_obj)
      return(seurat_obj)
    }
    stop(
      "The column name of cell_type or donor does not exist in the meta.data of your Seurat object, please add it!"
    )
  }

library(hdWGCNA)
selected_cell_type <- "Macrophage"

seurat_obj <- data_Mac <- readRDS("~/projects/m6A_LUAD/data_Mac.rds")

sub_seurat_obj <- seurat_obj %>%
  subset(cell_type == selected_cell_type)


k <- 50

sub_seurat_obj <- sub_seurat_obj %>%
  cat_construct_metacells(k = k, name = selected_cell_type)

sub_seurat_obj %>% saveRDS(file.path(outdir, "seurat_all.rds"))


#############################




source("~/projects/R_code/custom_function.R")

seurat_obj_path <-
  "~/projects/m6A_LUAD/result/cor2/seurat_all.rds"

selected_genes <-  c('METTL3',"METTL14","METTL16","WTAP","VIRMA","ZC3H13","CBLL1","RBM15","RBM15B","KIAA1429",
                     "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","HNRNPC","HNRNPG","RBMX","FMR1","HNRNPA2B1","IGF2BP1","IGF2BP2","IGF2BP3","PRRC2A",
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




#############################  Figure1  B & C #############



data <- readRDS("object.integrated.rds")



DefaultAssay(data) <- "RNA"



pdf("Mac-1.pdf",height=8,width=6)
plots <- VlnPlot(data, features = c("TREM2","CD163","MARCO"), group.by = "integrated_snn_res.0.12",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()


pdf("Mac.pdf",height=8,width=8)
p9 <- FeaturePlot(data,
                  reduction = "umap",
                  features = c("TREM2","CD163","MARCO"),
                  label = TRUE)
p9

dev.off()


data_Mac <- subset(data,subset=integrated_snn_res.0.12==1 | integrated_snn_res.0.12==11)


pdf("fig1.pdf",height=4,width=4)
p3 <- DimPlot(object.integrated, label = TRUE, reduction="umap") + NoLegend()
p3
dev.off()



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




pdf("figure-celltype.pdf",height=4,width=4)
DimPlot(data,
        reduction = "umap",
        group.by = "cell_type",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()

dev.off()



data$Type <- data$orig.ident
data$Type[grep("LUNG_N",data$orig.ident)] <- "Normal"
data$Type[grep("LUNG_T",data$orig.ident)] <- "Tumor"



data_Mac <- subset(data,subset=cell_type=="Macrophage")



pdf("figure4.pdf",height=4,width=4)
p1 <- DimPlot(data_Mac, reduction = "umap", group.by = "Type")
p1
dev.off()

saveRDS(data_Mac,file="data_Mac.rds")


Idents(data_Mac) <- "Type"
df_Mac <- FindMarkers(data_Mac, ident.1 = "Tumor", ident.2 = "Normal", verbose = FALSE)



data_Mac_Normal <- subset(data_Mac,subset=Type=="Normal")
cells_rankings <- AUCell_buildRankings(scRNA@assays$RNA@data) 

data_Mac_Tumor <- subset(data_Mac,subset=Type=="Tumor")



saveRDS(data_Mac_Normal, file = "data_Mac_Normal.rds")
saveRDS(data_Mac_Tumor, file = "data_Mac_Tumor.rds")




#################### Figure1 D


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


#################### Figure1 E


load("gene.cor.Normal.RData")
load("gene.cor.Tumor.RData")

load("df_Mac.RData")

cor_Normal <- gene.cor.Normal
cor_Tumor <- gene.cor.Tumor
DE <- df_Mac

gene_Normal <- cor_Normal$gene[which(abs(cor_Normal$estimate) > 0.1)]
gene_Tumor  <- cor_Tumor$gene[which(abs(cor_Tumor$estimate) > 0.1)]

gene_DE <- DE[DE$p_val < 0.05 &  DE$avg_log2FC > 0.25 ,]
gene_DE <- rownames(gene_DE)



which(table(c(gene_Normal,gene_Tumor,gene_DE))==3)


which(table(c(gene_Tumor,gene_DE))==2)

#BCL2A1      DUSP2     FTL   HLA-E   HMGB2     ID2       JUNB    LGMN     LYZ  MALAT1   NFKB1    NR4A2   NR4A3    PIM3  PLA2G7  

#   ZFP36  ZNF331 

############################


library(VennDiagram)

venn.diagram(list(cor.Normal=gene_Normal,cor.Tumor=gene_Tumor,DEGs=gene_DE), 
             fill=c("#729ECE","#FF9E4A","#67BF5C"), #,"#67BF5C","#ED665D","#AD8BC9"
             alpha=c(0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C"),
             cex=1,
             # cat.pos = 10,
             cat.dist = 0.01,
             cat.fontface=4, 
             # fontfamily=1,
             filename = "Venn.tiff",
             height = 1600, 
             width = 1600, resolution = 500)



#################### Figure1 F


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



seurat_obj <- subset(seurat_obj,subset=Type=="Tumor")

gene_DE <- c("BCL2A1","DUSP2","FTL","HLA-E","HMGB2","ID2","JUNB","LGMN","LYZ","MALAT1","NFKB1","NR4A2","NR4A3","PIM3","PLA2G7","ZFP36","ZNF331")


for(i in 1:length(gene_DE)){
  g1 <- seurat_obj$m6A
  g2 <- seurat_obj@assays$RNA@data[match(gene_DE[i],rownames(seurat_obj)),]
  
  data <- data.frame(g1=g1,g2=g2)
  
  
  x.name <- c("m6A")
  y.name <- c("Gene Expression")
  
  name1 <- paste0(gene_DE[i],".pdf")
  pdf(name1)
  p1 <- ggplot(data, aes(g1, g2)) + geom_point()
  p1 <- p1 +geom_smooth(method = 'lm', formula = y ~ x,colour="orange2",size=2)+ theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12, vjust = 3)) +labs(x=x.name,y=y.name) + stat_cor(method = "pearson",label.x =0.05, label.y =4.5)
  print(p1)
  dev.off()
}



#################### Figure1 G


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







################################


load("gene.cor.Tumor.RData")

gene_DE <- c("BCL2A1","DUSP2","FTL","HLA-E","HMGB2","ID2","JUNB","LGMN","LYZ","MALAT1","NFKB1","NR4A2","NR4A3","PIM3","PLA2G7","ZFP36","ZNF331")


cor1 <- gene.cor.Tumor
data <- cor1[match(gene_DE,cor1$gene),]

write.csv(data,file="sig_gene.csv",quote=F)



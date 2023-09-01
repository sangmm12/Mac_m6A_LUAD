setwd("~/projects/m6A_LUAD/result/cor1/cor_plot")

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











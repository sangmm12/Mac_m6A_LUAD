                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
setwd("~/projects/m6A_LUAD/result/cor1")

source("~/projects/R_code/custom_function.R")

seurat_obj_path <-
  "~/projects/m6A_LUAD/result/cor1/seurat_Tumor.rds"
# Ready to analyze genes?
selected_genes <-  c('METTL3',"METTL14","METTL16","WTAP","VIRMA","ZC3H13","CBLL1","RBM15","RBM15B","KIAA1429",
                     "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","HNRNPC","HNRNPG","RBMX","FMR1","HNRNPA2B1","IGFBP1","IGFBP2","IGFBP3","PRRC2A",
                     "FTO","ALKBH5")
# Ready to analyze cell type?
selected_cell_type <- "Macrophage"
# Number of nearest neighbors to aggregate?
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

cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE,assign=TRUE) 



p_value <- c()
estimate <- c()
for(i in 1:dim(seurat_obj)[1]){
  g1 <- seurat_obj@assays$RNA@data[i,]
  cor.test.res <- cor.test(seurat_obj$m6A,g1)
   p.value <- cor.test.res$p.value
   estimate1 <- cor.test.res$estimate
   p_value <- c( p_value,p.value )
   estimate <- c( estimate,estimate1)
}

 
gene.cor.Tumor <- data.frame(gene=all_genes,
                             p_value=p_value,
                             estimate=estimate)


save(gene.cor.Tumor,file="gene.cor.Tumor.RData")



















setwd("~/projects/m6A_LUAD")

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











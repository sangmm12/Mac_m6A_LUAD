setwd("~/projects/m6A_LUAD/result/cor1")
load("gene.cor.Tumor.RData")

gene_DE <- c("BCL2A1","DUSP2","FTL","HLA-E","HMGB2","ID2","JUNB","LGMN","LYZ","MALAT1","NFKB1","NR4A2","NR4A3","PIM3","PLA2G7","ZFP36","ZNF331")


cor1 <- gene.cor.Tumor
data <- cor1[match(gene_DE,cor1$gene),]

write.csv(data,file="sig_gene.csv",quote=F)


setwd("~/projects/m6A_LUAD/result/cor1")

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







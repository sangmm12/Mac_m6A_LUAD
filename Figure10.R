
##############Figure 10

###############Figure 10 A

cancer <- c("LUAD")

sur <- read.csv("D:/R/TCGA_Fancancer/lasso/sur_all.csv",header = F)

n1 <- grep(cancer,sur$V1)

sur <- sur[n1,]

sur<- sur[,-1]
rownames(sur) <- sur[,1]
sur <- sur[,-1]


library(data.table)
dat <- fread("D:/R/Mac-m6A_LUAD/lasso/lasso.csv",header = T)
dat <- as.data.frame(dat)
rownames(dat) <- dat[,1]



cluster <- read.csv(paste("D:/R/Mac-m6A_LUAD/xcell cluster/xcell cluster_hc.csv",sep=''),header = F)
rownames(cluster) <- cluster[,1]

all_name <- names(which(table(c(rownames(dat),rownames(sur),rownames(cluster) ))==3))


sur1 <- sur[match(all_name,rownames(sur)),]

dat1 <- dat[match(all_name,rownames(dat)),]

cluster1 <- cluster[match(all_name,rownames(cluster)),]


data <- data.frame(Time=sur1$V3,Status=sur1$V4,lasso = dat1[,2], cluster=cluster1[,2])

rownames(data) <- all_name

write.csv(data,file=paste("sur_lasso_xcellcluster.csv",sep=''),quote=F)

dat<- data
#########################################################################
high_high <- subset(dat,dat$lasso=="high" & dat$cluster=="cold")
high_low <- subset(dat,dat$lasso=="high" & dat$cluster=="hot")
low_high <- subset(dat,dat$lasso=="low" & dat$cluster=="cold")
low_low <- subset(dat,dat$lasso=="low" & dat$cluster=="hot")

class1 <- rep("high-cold", times=length(rownames(high_high)))
class2 <- rep("high-hot", times=length(rownames(high_low)))
class3 <- rep("low-cold", times=length(rownames(low_high)))
class4 <- rep("low-hot", times=length(rownames(low_low)))

length(rownames(high_high))
length(rownames(high_low))
length(rownames(low_high))
length(rownames(low_low))

group<- c(class1,class2,class3,class4)


KM4 <- rbind(high_high,high_low,low_high,low_low)

KM4 <- cbind(group,KM4)

write.csv(KM4,file=paste("sur_lasso_xcellcluster_KM4",".csv",sep=''),quote=F)


###############Figure 10 B



library(GSVA)
library(GSEABase)
library(limma)
library(Seurat)
library(msigdbr)

#all_gene_sets = msigdbr(species = "Mus musculus")
human <- msigdbr(species = "Homo sapiens")
#hallmark
human_GO_bp = msigdbr(species = "human",
                      category = "H") %>% 
  dplyr::select(gs_name,gene_symbol)


human_GO_bp_Set = human_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)


s.sets <- human_GO_bp_Set


###########

library(data.table)
data <- read.table("D:/R/Mac-m6A_LUAD/LUAD_Count/convert_exp.txt",row.names = 1)
colnames(data) <- data[1,]
data <- data[-1,]



cluster <- read.csv("D:/R/Mac-m6A_LUAD/lasso&xcell/KM4/KM2.csv",header = F)


# c1 <- cluster[which(cluster$V2=="hot"),]
# c2 <- cluster[which(cluster$V2=="cold"),]

c1 <- cluster[which(cluster$V2=="low-hot"),]
c2 <- cluster[which(cluster$V2=="high-cold"),]



data_c1 <- as.matrix(data)[,na.omit(match(c1[,1],colnames(data)))]
data_c2 <- as.matrix(data)[,na.omit(match(c2[,1],colnames(data)))]


data_merge <- cbind(data_c1,data_c2)
data_merge=apply(data_merge,2,as.numeric)


rownames(data_merge) <- rownames(data)

##########

expr1 <- as.matrix(data_merge)


expr <- expr1



es.matrix <- gsva(
  expr,
  s.sets,
  min.sz = 10,
  max.sz = Inf,
  tau = 1,
  method = "gsva",
  abs.ranking = FALSE,
  verbose = TRUE,
  parallel.sz = 1
)

saveRDS(es.matrix,file="hallmark.rds")





n1 <- 1:dim(data_c1)[2]
#grep("male",seurat_obj@meta.data$sex)

n2 <- dim(data_c1)[2]+1:dim(data_c2)[2]

#grep("female",seurat_obj@meta.data$sex)

es.matrix.1 <-
  as.data.frame(es.matrix[, n1],
                row.names = row.names(es.matrix))
es.matrix.2 <-
  as.data.frame(es.matrix[, n2],
                row.names = row.names(es.matrix))


es.matrix.f <- cbind(es.matrix.1, es.matrix.2)

grouP <-
  c(rep("case", dim(es.matrix.1)[2]),
    rep("control", dim(es.matrix.2)[2]))

grouP <- as.factor( grouP)
design <- model.matrix(~ grouP + 0)



row.names(design)<-c(colnames(es.matrix.1), colnames(es.matrix.2))

comparE <-
  makeContrasts(grouPcase - grouPcontrol, levels = design)

fit <- lmFit(es.matrix, design)
fit2 <- contrasts.fit(fit, comparE)
fit3 <- eBayes(fit2)


diff <- topTable(fit3, coef = 1, number = dim(es.matrix)[1])

t_results <-
  as.data.frame(diff$t, row.names = rownames(es.matrix))
head(t_results)
colnames(t_results) <- c("t_value")




saveRDS(t_results,file="t_results_hallmark.rds")
write.csv(t_results,file="t_results_hallmark.csv")

#t_results <- readRDS("t_results_c5.go.bp_Elongating.rds")

library(ggplot2)

library(pheatmap)
rownames(t_results) <- gsub("HALLMARK_", "", rownames(t_results))
rownames(t_results) <- gsub("_", " ", rownames(t_results))
focus.cluster <- "t_value"
sub_t_results <- as.data.frame(t_results[, focus.cluster],
                               row.names = rownames(t_results))
sub_t_results$hallmark <- rownames(sub_t_results)
colnames(sub_t_results) <- c("t", "hallmark")

sub_t_results$hallmark = with(sub_t_results, reorder(hallmark, t))
sub_t_results$fill <- ""
sub_t_results[sub_t_results$t >= 2.58,]$fill <-
  "up"
sub_t_results[sub_t_results$t <= -2.58,]$fill <-
  "down"
sub_t_results[abs(sub_t_results$t) < 2.58,]$fill <-
  "no"
sub_t_results$color <- ""
sub_t_results[abs(sub_t_results$t) < 2.58,]$color <-
  "n"
sub_t_results[abs(sub_t_results$t) >= 2.58,]$color <-
  "y"


sub_t_results <- sub_t_results[c(1:50),]

p <-
  ggplot(sub_t_results, aes(x = hallmark, y = t, fill = fill)) +
  geom_bar(stat = "identity", width = .8) +
  scale_fill_manual(
    values = c(
      "down" = "#36648b",
      "up" = "#e94644",
      "no" = "#cccccc"
    ),
    guide = F
  ) + ylim(-10,14)+
  geom_hline(
    yintercept = c(-2.58, 2.58),
    color = "white",
    linetype = "dotted",
    size = 0.5
  ) +
  coord_flip() +
  xlab("") +
  geom_text(
    data = subset(sub_t_results, t < 0),
    aes(
      x = hallmark,
      y = 0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 1.8,
    hjust = "inward"
  ) +geom_text(
    data = subset(sub_t_results, t > 0),
    aes(
      x = hallmark,
      y = -0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 1.8,
    hjust = "outward"
  ) +
  scale_colour_manual(values = c("y" = "black", "n" = "#cccccc"),
                      guide = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    ),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),
  )
#p
ggsave(
  filename = "hallmark.pdf",
  plot = p,
  height = 6.5,
  width =5
)



###############Figure 10 C D



library(data.table)
data <-  fread(paste("D:/R/Mac-m6A_LUAD/xCELL_MAC.csv",sep=''))
data <- as.data.frame(data)

n1 <- grep("LUAD",data$CODE)

data <- data[n1,]

data<- data[,-1]
rownames(data) <- data[,5]
data <- data[,-5]

dat_gene <- data[,-1]


data_cell <- readRDS(paste("D:/R/Mac-m6A_LUAD/lasso&xcell/KM4/GSVA_hallmark/hallmark.rds",sep=''))

DEG_BLCA <- read.csv("D:/R/Mac-m6A_LUAD/lasso&xcell/KM4/GSVA_hallmark/t_results_hallmark.csv",header = T)

GSVAdown <- DEG_BLCA$X[which(DEG_BLCA$t_value<= -2.58)]

dat_GSVAdown <- data_cell[match(GSVAdown,rownames(data_cell)),]

dat_im <- t(dat_GSVAdown)

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene)))==2))

dat_im <- dat_im[match(all_name,rownames(dat_im)),]

dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]

colnames_dat_im <- DEG_BLCA$A[which(DEG_BLCA$t_value<= -2.58)]

colnames(dat_im) <- colnames_dat_im 


# for(i in 1:length(colnames(dat_gene)))
# {
#   dat_gene[,i] = as.numeric(unlist(dat_gene[,i]))
# }
# i=1
# nrow(dat_gene)
# ncol(dat_gene)
# 

colSums(dat_im)



library(psych)
data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

paste("data.r_",cancer,".csv",sep='')
write.csv(data.r,file= paste("data.r_xcell.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_xcell.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_KM2_xcell.pdf",sep=''),width =4.8,height = 8)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()

















#############################################################################




library(data.table)
data <-  fread(paste("D:/R/Mac-m6A_LUAD/LUAD_TPM/LUAD_convert_exp.txt",sep=''))
data <- as.data.frame(data)

rownames(data) <- data[,1]
data <- data[,-1]

genelist <- c("NR4A2","HMGB2","PIM3","ZNF331","ID2","ZFP36","NR4A3","MALAT1","LYZ","FTL",
              "LGMN","NFKB1","JUNB","HLA-E","BCL2A1","PLA2G7","DUSP2")

dat_gene <- data[match(genelist,rownames(data)),]

dat_gene=as.data.frame(lapply(dat_gene,as.numeric),check.names=F)

rownames(dat_gene) <- genelist

dat_gene <- t(dat_gene)






data_cell <- readRDS(paste("D:/R/Mac-m6A_LUAD/lasso&xcell/KM4/GSVA_hallmark/hallmark.rds",sep=''))

DEG_BLCA <- read.csv("D:/R/Mac-m6A_LUAD/lasso&xcell/KM4/GSVA_hallmark/t_results_hallmark.csv",header = T)

GSVAdown <- DEG_BLCA$X[which(DEG_BLCA$t_value<= -2.58)]

dat_GSVAdown <- data_cell[match(GSVAdown,rownames(data_cell)),]

dat_im <- t(dat_GSVAdown)

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene)))==2))

dat_im <- dat_im[match(all_name,rownames(dat_im)),]

dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]

colnames_dat_im <- DEG_BLCA$A[which(DEG_BLCA$t_value<= -2.58)]

colnames(dat_im) <- colnames_dat_im 

# for(i in 1:length(colnames(dat_gene)))
# {
#   dat_gene[,i] = as.numeric(unlist(dat_gene[,i]))
# }
# i=1
# nrow(dat_gene)
# ncol(dat_gene)
# 

colSums(dat_im)



library(psych)
data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

paste("data.r_",cancer,".csv",sep='')
write.csv(data.r,file= paste("data.r_gene.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_gene.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_KM2_gene.pdf",sep=''),width =9,height = 8)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()










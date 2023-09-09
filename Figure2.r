#Figure2A
library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, cyto.exclude<-c("chrX", "chrY"),tracks.inside=10, tracks.outside=0 )
out.file<-"RCircosDemoHumanGenome.pdf"#output file name
pdf(file=out.file,height=8,width=8,compress=TRUE)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
data(RCircos.Gene.Label.Data) 
name.col <- 4 
side <- "in" 
track.num <- 1 
data = read.csv('mac-m6a.csv')
head(data)
RCircos.Gene.Connector.Plot(data,
                              + track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(data,
                         + name.col,track.num, side)#gene name
dev.off()

#Figure2B
setwd("~/cor")

my_data = read.csv('Mac-m6A-LUAD.csv',header=T,row.names=1,check.names=F)
res <- cor(my_data)


library(corrplot)
pdf("out.pdf",width=30,height=30)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
dev.off()

#Figure2D
setwd("~")
rm(list = ls())
library(dplyr)

tumor_id = "LUAD"	#Tumor type  

dat_input <- read.csv("D:/R/home/cor/exp/Mac-m6a-LUAD.csv",check.names = F)

#dat_input <- dat_input[grep(tumor_id,dat_input$CODE),]
#dat_input <- dat_input[,-1]




rownames(dat_input) <- make.unique(sub('\n','',dat_input$SampleName))

dat_input <- dat_input[,-1]

data_all_tumor <- read.csv(paste("m6a_gene.csv",sep=''),check.names = F)

dat_one_tumor <- data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]

rownames(dat_one_tumor) <- make.unique(sub('\n','',dat_one_tumor$SampleName))

colscluster = length(rownames(dat_one_tumor))/3

dat_one_tumor <- dat_one_tumor[,!duplicated(colnames(dat_one_tumor))]

dat_one_tumor <- dat_one_tumor[,3:length(colnames(dat_one_tumor))]

one_tumor_sample <- unlist(rownames(dat_one_tumor))

all_name <- names(which(table(c(rownames(dat_input),one_tumor_sample))==2))

dat_gene <- dat_one_tumor[match(all_name,rownames(dat_one_tumor)),]

write.csv(t(dat_one_tumor),file=paste("M6A",".csv",sep=''),quote=F)

dat_im <- dat_input[match(all_name,rownames(dat_input)),]

library(psych)

data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")

data.r <- data.corr$r
data.p <- data.corr$p



write.csv(data.r,file=paste(tumor_id,"_r.csv",sep=''),quote=F)
write.csv(data.p,file=paste(tumor_id,"_p.csv",sep=''),quote=F)
library(pheatmap)


getSig <- function(dc) {
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
myColor <- colorRampPalette(c("#FF8500", "white", "#024C68"))(paletteLength)



test <- data.r
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))


pdf(paste(tumor_id,".pdf",sep=''),width = length(colnames(dat_input))/1.5,height = length(colnames(data_all_tumor))/1.5)
if(colscluster > 1){
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, display_numbers=sig.mat,fontsize=length(colnames(dat_input)))
}else{
   pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,cluster_cols=F, display_numbers=sig.mat)
}
dev.off()


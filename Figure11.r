#Figure9A
setwd("~")
rm(list = ls())
library(dplyr)

tumor_id = "LUAD"


dat_input <- read.csv("Mac-m6A-LUAD.csv",check.names = F)

#带CODE种类列表用
dat_input <- dat_input[grep(tumor_id,dat_input$CODE),]
dat_input <- dat_input[,-1]




rownames(dat_input) <- make.unique(sub('\n','',dat_input$SampleName))

dat_input <- dat_input[,-1]

dat_one_tumor <- read.csv("drug.csv",check.names = F)

rownames(dat_one_tumor) <- make.unique(sub('\n','',dat_one_tumor$SampleName))

dat_one_tumor <- dat_one_tumor[,!duplicated(colnames(dat_one_tumor))]

dat_one_tumor <- dat_one_tumor[,-1]

one_tumor_sample <- unlist(rownames(dat_one_tumor))

all_name <- names(which(table(c(rownames(dat_input),one_tumor_sample))==2))

dat_gene <- dat_one_tumor[match(all_name,rownames(dat_one_tumor)),]

#write.csv(t(dat_one_tumor),file=paste("M6A",".csv",sep=''),quote=F)

dat_im <- dat_input[match(all_name,rownames(dat_input)),]

colscluster = length(colnames(dat_im))/3

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


pdf(paste(tumor_id,".pdf",sep=''),width = length(colnames(dat_input))/2,height = length(colnames(dat_one_tumor))/2)
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

#Figure9B
setwd("~")
rm(list = ls())
library(dplyr)

tumor_id = "LUAD"


dat_input <- read.csv("Mac-m6A-LUAD.csv",check.names = F)

#带CODE种类列表用
dat_input <- dat_input[grep(tumor_id,dat_input$CODE),]
dat_input <- dat_input[,-1]




rownames(dat_input) <- make.unique(sub('\n','',dat_input$SampleName))

dat_input <- dat_input[,-1]

dat_one_tumor <- read.csv("drug.csv",check.names = F)

rownames(dat_one_tumor) <- make.unique(sub('\n','',dat_one_tumor$SampleName))

dat_one_tumor <- dat_one_tumor[,!duplicated(colnames(dat_one_tumor))]

dat_one_tumor <- dat_one_tumor[,-1]

one_tumor_sample <- unlist(rownames(dat_one_tumor))

all_name <- names(which(table(c(rownames(dat_input),one_tumor_sample))==2))

dat_gene <- dat_one_tumor[match(all_name,rownames(dat_one_tumor)),]

#write.csv(t(dat_one_tumor),file=paste("M6A",".csv",sep=''),quote=F)

dat_im <- dat_input[match(all_name,rownames(dat_input)),]

colscluster = length(colnames(dat_im))/3

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


pdf(paste(tumor_id,".pdf",sep=''),width = length(colnames(dat_input))/2,height = length(colnames(dat_one_tumor))/2)
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

#Figure9C

data_all_tumor <- read.csv(paste("data/m5c-laml/mirna.csv",sep=''),check.names = F)
data_all_tumor<- data_all_tumor[,colnames(data_all_tumor)%in%c('SampleName',mirna)]

dat_input <- read.csv("~/diff/data/EPIC.csv",check.names = F,header=T)
dat_input<-dat_input[,colnames(dat_input)%in%c('SampleName',sample)]

out<- merge(dat_input,data_all_tumor, by = "SampleName")
out <- subset(out,select=-c(SampleName))
library(ggplot2)

fit <- lm(out[,2] ~ out[,1])

summary_fit <- summary(fit)

r_squared <- summary_fit$r.squared
p_value <- summary_fit$coefficients[2, 4]

p <- ggplot(out,aes(x=out[,1],y=out[,2])) + geom_point(shape=19) + xlab(sample) + ylab(mirna)+geom_smooth(method = lm)+theme_classic()+theme(axis.text = element_text(size = 25),axis.title = element_text(size = 25))+
  annotate('text',size=8, y=max(out[,2])+max(out[,2])/5,x = max(out[,1])/4,label = paste("R=",format(sqrt(round(r_squared, 3)),digits = 2), " P=", format(p_value, scientific=TRUE, digits=2),sep=''),colour="black")

#Figure9D
data_all_tumor <- read.csv(paste("data/m5c-laml/mirna.csv",sep=''),check.names = F)
data_all_tumor<- data_all_tumor[,colnames(data_all_tumor)%in%c('SampleName',mirna)]

dat_input <- read.csv("~/diff/data/EPIC.csv",check.names = F,header=T)
dat_input<-dat_input[,colnames(dat_input)%in%c('SampleName',sample)]

out<- merge(dat_input,data_all_tumor, by = "SampleName")
out <- subset(out,select=-c(SampleName))
library(ggplot2)

fit <- lm(out[,2] ~ out[,1])

summary_fit <- summary(fit)

r_squared <- summary_fit$r.squared
p_value <- summary_fit$coefficients[2, 4]

p <- ggplot(out,aes(x=out[,1],y=out[,2])) + geom_point(shape=19) + xlab(sample) + ylab(mirna)+geom_smooth(method = lm)+theme_classic()+theme(axis.text = element_text(size = 25),axis.title = element_text(size = 25))+
  annotate('text',size=8, y=max(out[,2])+max(out[,2])/5,x = max(out[,1])/4,label = paste("R=",format(sqrt(round(r_squared, 3)),digits = 2), " P=", format(p_value, scientific=TRUE, digits=2),sep=''),colour="black")

#Figure9E
rm(list=ls())
file_dir = "~/diff/mac-m6a"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)
need_list = c('T','N','M','STAGE','AGE','GENDER')	#This is the group document

for (file_name in need_list)
{

  tumor_list = c('LUAD')	#Use this to change what type of tumor you want to analyse
  
  final_tumor_list = c()
  
  my_data = read.csv('tcgaexp',header=T,check.names=F)	#Expression matrix document
  
  other_data = read.csv(paste(file_name,".csv",sep=''),header=T,check.names=F)
  
  p_value_csv <- data.frame()
  
  for(name in tumor_list)
  {
  
    dat <- data.frame(check.names = F)
    
    
    other_file <- other_data[grep(name,other_data$CODE),]
    
    if (length(rownames(other_file))==0)
    {
        next
    }
    
    other_file <- other_file[!duplicated(other_file$SampleName),]
    rownames(other_file) = gsub('\n','',other_file$SampleName)
    other_file <- subset(other_file,select=-c(SampleName,CODE))
    exp_file <- my_data[grep(name,my_data$CODE),]
    exp_file <- exp_file[grep("Tumor",exp_file$Group),]
    exp_file <- exp_file[!duplicated(exp_file$SampleName),]
    rownames(exp_file) = gsub('\n','',exp_file$SampleName)
    exp_file <- subset(exp_file,select=-c(Group,CODE,SampleName))
    gene_list <- colnames(exp_file)
    all_name <- names(which(table(c(rownames(other_file),rownames(exp_file)))==2))
    if (length(all_name)==0)
    {
        next
    }
    
    for(gene_name in gene_list)
    {
        for(i in all_name)
        {
            dat <- rbind(dat,c(gene_name,other_file[match(i,rownames(other_file)),],exp_file[c(gene_name)][match(i,rownames(exp_file)),]))
        }
    
    }
    colnames(dat) <- c("Gene","Group","value")
    dat[,3] = as.numeric(dat[,3])
    dat[,3] = log(dat[,3]+1)
    dat <- na.omit(dat)
    dat <- dat[dat[,1]%in%temp_name,]
    final_tumor_list <- append(final_tumor_list,name)
    print(name)
    pdf(paste(file_name,"-",name,".pdf",sep=''),width=length(unique(dat[,1])),height = 8)
    p <- ggboxplot(dat, x = "Gene", y = "value",
                   color = "Group", palette = 'jama',
                   add = "jitter",x.text.angle=60)
    p <- p + xlab("")+ylab("Gene Expression(log(x+1))")
    p <- p + theme(axis.text = element_text(size = 15),axis.title=element_text(size=30))
    print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
    dev.off()
  }
}

#Figure9F
rm(list=ls())
file_dir = "~/diff/mac-m6a"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)
need_list = c('T','N','M','STAGE','AGE','GENDER')	#This is the group document

for (file_name in need_list)
{

  tumor_list = c('LUAD')	#Use this to change what type of tumor you want to analyse
  
  final_tumor_list = c()
  
  my_data = read.csv('tcgaexp',header=T,check.names=F)	#Expression matrix document
  
  other_data = read.csv(paste(file_name,".csv",sep=''),header=T,check.names=F)
  
  p_value_csv <- data.frame()
  
  for(name in tumor_list)
  {
  
    dat <- data.frame(check.names = F)
    
    
    other_file <- other_data[grep(name,other_data$CODE),]
    
    if (length(rownames(other_file))==0)
    {
        next
    }
    
    other_file <- other_file[!duplicated(other_file$SampleName),]
    rownames(other_file) = gsub('\n','',other_file$SampleName)
    other_file <- subset(other_file,select=-c(SampleName,CODE))
    exp_file <- my_data[grep(name,my_data$CODE),]
    exp_file <- exp_file[grep("Tumor",exp_file$Group),]
    exp_file <- exp_file[!duplicated(exp_file$SampleName),]
    rownames(exp_file) = gsub('\n','',exp_file$SampleName)
    exp_file <- subset(exp_file,select=-c(Group,CODE,SampleName))
    gene_list <- colnames(exp_file)
    all_name <- names(which(table(c(rownames(other_file),rownames(exp_file)))==2))
    if (length(all_name)==0)
    {
        next
    }
    
    for(gene_name in gene_list)
    {
        for(i in all_name)
        {
            dat <- rbind(dat,c(gene_name,other_file[match(i,rownames(other_file)),],exp_file[c(gene_name)][match(i,rownames(exp_file)),]))
        }
    
    }
    colnames(dat) <- c("Gene","Group","value")
    dat[,3] = as.numeric(dat[,3])
    dat[,3] = log(dat[,3]+1)
    dat <- na.omit(dat)
    dat <- dat[dat[,1]%in%temp_name,]
    final_tumor_list <- append(final_tumor_list,name)
    print(name)
    pdf(paste(file_name,"-",name,".pdf",sep=''),width=length(unique(dat[,1])),height = 8)
    p <- ggboxplot(dat, x = "Gene", y = "value",
                   color = "Group", palette = 'jama',
                   add = "jitter",x.text.angle=60)
    p <- p + xlab("")+ylab("Gene Expression(log(x+1))")
    p <- p + theme(axis.text = element_text(size = 15),axis.title=element_text(size=30))
    print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
    dev.off()
  }
}
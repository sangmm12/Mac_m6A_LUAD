#Figure5A,B,C,D,E,F
tumor_id = 'LUAD'
name_fold = 'GSE68465'
main_fold = 'mac-m6a'
setwd(paste('~/',main_fold,'/',name_fold,sep=''))
all_sur_data <- read.table(paste(tumor_id,'_sur.txt',sep=''),header=T,check.names=F)

gene_exp <- as.data.frame(t(read.table(paste(name_fold,'_exp.txt',sep=''),header=T,check.names=F)))
colnames(gene_exp) <- gene_exp[1,]
gene_exp <- gene_exp[-1,]

gene_exp$SampleName = rownames(gene_exp)

erged_df <- merge(gene_exp, all_sur_data, by = "SampleName")

erged_df$Time <- erged_df$Time/365
erged_df <- na.omit(erged_df)

survival_time <- erged_df$Time
survival_status <- erged_df$Status
patient_id <- erged_df$SampleName
erged_df <-subset(erged_df, select = -Status)
erged_df <-subset(erged_df, select = -SampleName)
erged_df <-subset(erged_df, select = -Time)
erged_df=as.data.frame(lapply(erged_df,as.numeric))
rownames(erged_df) <- patient_id
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean,"-")
  rv <- sweep(rv, 1, rowsd, "/")
  return(rv)
}

#########################################

library(survival)
library(glmnet)
seed <- sample(1:10000,size=1)

set.seed(seed)
surv_obj <- Surv(event = survival_status,time = survival_time)

fit <- glmnet(as.matrix(erged_df), surv_obj, family = "cox",maxit=1000)

cvfit <- cv.glmnet(as.matrix(erged_df), surv_obj, family = "cox",maxit=1000)

erged_df$Time <- survival_time
erged_df$Status <- survival_status

if(file.exists(paste('~/',main_fold,'/coef.Rdata',sep='')))
{
  load(paste('~/',main_fold,'/coef.Rdata',sep=''))
  print('load....')
} else {
  coef=coef(fit, s = cvfit$lambda.min)
  if(length(coef[which(coef != 0),])==0)
  {
    stop('coef = 0')
  }
  save(coef,file = paste('~/',main_fold,'/coef.Rdata',sep=''))
  print('save....')
}

index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

trainFinalGeneExp=erged_df[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("Time","Status",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
num_out <- standardize(as.data.frame(lapply(erged_df[,lassoGene],as.numeric)))
rownames(num_out) <- rownames(erged_df)
outTab=cbind(erged_df[,c("Time","Status")],num_out,riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file=paste(name_fold,"Risk.txt",sep=''),sep="\t",quote=F,row.names=F)

################################################

library(timeROC)

bioROC=function(inputFile=null,rocFile=null){

  rt=read.table(inputFile,header=T,sep="\t")
  ROC_rt=timeROC(T=rt$Time,delta=rt$Status,
                 marker=rt$riskScore,cause=1,
                 weighting='aalen',
                 times=c(1,2,3),ROC=TRUE)
  pdf(file=rocFile,width=5,height=5)
  plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC_rt,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green",'blue','red'),lwd=2,bty = 'n')
  dev.off()
  return(max(c(ROC_rt$AUC[1],ROC_rt$AUC[2],ROC_rt$AUC[3])))
}
roc <- bioROC(inputFile=paste(name_fold,"Risk.txt",sep=''),rocFile=paste(name_fold,".ROC.pdf",sep=''))

######################################

library(survminer)
bioSurvival=function(inputFile=null,outFile=null){
  rt=read.table(inputFile,header=T,sep="\t")                  
  diff=survdiff(Surv(Time, Status) ~ risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
  surPlot=ggsurvplot(fit, 
                     conf.int = TRUE,
                     data=rt,
                     pval=paste0("p=",pValue),
                     pval.size=5,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     risk.table=F,
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 5,height =4.5)
  print(surPlot)
  dev.off()
  return(pValue)
}
p_value <- bioSurvival(inputFile=paste(name_fold,"Risk.txt",sep=''),outFile=paste(name_fold,".survival.pdf",sep=''))

#######################################################################3

library(pheatmap)
bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null){
  rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)  
  rt=rt[order(rt$riskScore),]    
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=riskScoreFile,width = 8,height = 6)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("lightblue",lowLength),rep("red",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
  dev.off()
  color=as.vector(rt$Status)
  color[color==1]="red"
  color[color==0]="lightblue"
  pdf(file=survStatFile,width = 8,height = 6)
  plot(rt$Time, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
}

bioRiskPlot(inputFile=paste(name_fold,"Risk.txt",sep=''),riskScoreFile=paste(name_fold,".riskScore.pdf",sep=''),survStatFile=paste(name_fold,".survStat.pdf",sep=''))

##########################################################

risk.out <- read.table(paste(name_fold,'Risk.txt',sep=''),row.names = 1,header=T)
risk.out <- risk.out[order(risk.out$risk),]
risk.out <-subset(risk.out, select = -c(Time,Status,riskScore))
Group <- risk.out$risk
risk.out <-subset(risk.out, select = -risk)
Group <- data.frame(Risk = Group)
rownames(Group) = rownames(risk.out)
library(pheatmap)
pdf('phep.pdf',width=6,height=4)
mycolors <- colorRampPalette(c("darkblue","white", "darkred"))(15)
tmp=as.data.frame(lapply(risk.out,as.numeric))
tmp <- as.matrix(t(tmp))
colnames(tmp) <- rownames(Group)
pheatmap(scale(tmp),annotation_col = Group,color = mycolors,cluster_rows=F,show_colnames=F,cluster_cols=F,cellheight = 10,cellwidth = 200/length(colnames(tmp)),border=F)
dev.off()
print(paste('p:',p_value,' roc:',roc,sep=''))

#Figure5H
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

#Figure5G
setwd(dir = "~/rate_histogram/txt")
library("stringr")
library("do")
temp = list.files(pattern="*.txt")
for (names in temp)
{
name = strsplit(names,'.',1)[[1]][1]
print(name)
tables = read.table(str_c(name,".txt"),sep="\t",header=T)
tables$Status = Replace(tables$Status,"1", "Dead")
tables$Status = Replace(tables$Status,"0", "Alive")
library(tidyverse)
tables <- count(tables,risk,Status)
print(tables)
row_num <- length(row(tables))/length(tables)
temp_pct <-matrix(tables$n, nrow = row_num)

for(i in 1:(row_num/2))
  {
    temp_sum = temp_pct[i*2-1,1] + temp_pct[i*2,1]
    temp_pct[i*2-1,1] <- temp_pct[i*2-1,1]/temp_sum
    temp_pct[i*2,1] <- temp_pct[i*2,1]/temp_sum
  }

pdf(str_c("~/rate_histogram/",name,".pdf"),width=length(unique(tables$risk))*2,height = 8)

library("ggplot2")
p <- ggplot(tables,aes(fill=Status, y= temp_pct, x = risk)) + 
  geom_col()+
  geom_text(aes(label = scales::percent(temp_pct)))+
  scale_y_continuous(labels = scales::percent) + 
  labs(y = "Rate",x=name) +
  theme_gray()
print(p)
dev.off()
}


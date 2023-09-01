file_dir = "~"
setwd(dir = file_dir)
temp = list.files(pattern="*.csv")
library(ggplot2)
library(ggpubr)


namess <- c()
for(i in temp){

now_name <- regmatches(i,regexpr('_([A-Z]*)', i))
namess <- append(namess,sub('_','', now_name))

name = strsplit(i,'.',1)
name = name[[1]][1]
print(name)

dat <- read.csv(paste(file_dir,name,".csv",sep=''))
#dat[,3] = log2(dat[,3]+1)
name_class = dat[,2][!duplicated(dat[,2])]
xxxx <- matrix(unique(dat[,1])) #colname
temp_length <- length(xxxx)*2

#colnames(xxx) <- namess
colnames(dat) <- c("Gene","Group","value")

# xx <-compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
# p_value <- as.matrix(xx$p)

# p_value[is.na(p_value)] <- 0.5

#xxx<-cbind(xxx,as.matrix(xx$p))

# write.csv(xxx,file=paste(file_dir,"p/p.csv",sep=''),quote=F)

# for(j in 1:length(p_value))
# {
#     if (p_value[j] > 0.05)
#     {
#         dat <- dat[dat[,"Gene"]!=xxxx[j],]
#         temp_length <- temp_length - 2
#     }
# }
pdf(paste(name,".pdf",sep=''),width=temp_length/3+2,height = 5)
p <- ggboxplot(dat, x = "Group", y = "value",
                 palette = "npg", 
                fill="Group",x.text.angle=60,x.text.size=15)

p <- p + xlab("Cluster")+ylab(paste(name,"sensitivity (IC50)"))

p <- p + theme(axis.text = element_text(size = 15),axis.title=element_text(size=15))
my_comparisons <- list(  name_class )
p <- p + stat_compare_means(comparisons = my_comparisons)
print(p)
dev.off()
}
print(namess)

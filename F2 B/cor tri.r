setwd("~/cor")

my_data = read.csv('Mac-m6A-LUAD.csv',header=T,row.names=1,check.names=F)
res <- cor(my_data)


library(corrplot)
pdf("out.pdf",width=30,height=30)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
dev.off()
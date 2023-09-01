


Idents(data) <- c("Cell_subtype")
pdf("Cell_subtype.pdf",height = 10,width = 10)
DimPlot(data, reduction = "tsne",raster=FALSE)
dev.off()



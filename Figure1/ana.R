
setwd("~/projects/m6A_LUAD")



ana_Sample <- Sample_name[grep("LUNG",Sample_name)]

scRNAlist <- list()
for(i in 1:length(ana_Sample)){
  name <- paste0("data/",ana_Sample[i],".rds")
  counts <- readRDS(name)
  
  scRNAlist[[i]] <-
    CreateSeuratObject(
      counts,
      project = ana_Sample[i],
      min.cells = 5,
      min.features = 200,
    )
  scRNAlist[[i]] <- PercentageFeatureSet( scRNAlist[[i]], pattern = "^mt-", col.name = "percent.mt")
}


obj <- merge(scRNAlist[[1]],scRNAlist[-1])


obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5 )


sce.all.list <- SplitObject(obj, split.by = "orig.ident")

sce.all.list <- lapply(X = sce.all.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1500)
})





select.integrate.feature=1500
pc.num=15
resolution=0.12

remove(data_all)
object.features <- SelectIntegrationFeatures(object.list = sce.all.list, nfeatures = select.integrate.feature)

object.anchors <- FindIntegrationAnchors(object.list = sce.all.list,dims=1:15,
                                         anchor.features = object.features)
object.integrated <- IntegrateData(anchorset = object.anchors,dims=1:15)


remove(sce.all.list)
remove(object.anchors)
DefaultAssay(object.integrated) <- "integrated"


object.integrated <- ScaleData(object.integrated, verbose = FALSE)
object.integrated <- RunPCA(object.integrated, npcs = 15, verbose = FALSE)
object.integrated <- RunUMAP(object.integrated, reduction = "pca", dims = 1:15)
object.integrated <- RunTSNE(object.integrated, dims = 1:15, verbose = FALSE)
object.integrated <- FindNeighbors(object.integrated, reduction = "pca", dims = 1:15)
object.integrated<- FindClusters(object.integrated, resolution = 0.12)

saveRDS(object.integrated, file = "object.integrated.rds")

p1 <- DimPlot(object.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object.integrated, reduction = "tsne", group.by = "orig.ident")
p3 <- DimPlot(object.integrated, label = TRUE, reduction="umap") + NoLegend()
p4 <- DimPlot(object.integrated, label = TRUE, reduction="tsne") + NoLegend()

#merged clusters with label
p5 <- DimPlot(object.integrated, split.by="orig.ident",reduction="umap",ncol=2)
p6 <- DimPlot(object.integrated, split.by="orig.ident",reduction="tsne",ncol=2)
#merged clusters without label
p7 <- DimPlot(object.integrated, split.by="orig.ident",reduction="umap",label = TRUE ,ncol=2) + NoLegend()
p8 <- DimPlot(object.integrated, split.by="orig.ident",reduction="tsne",label = TRUE, ncol=2) + NoLegend()


#QC plot
p9 <- VlnPlot(object.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "orig.ident")

p <- plot_grid(p1, p2,p3,p4,ncol=2)
pdf("figure-merge-1.pdf",height=8,width=8)
print(p)
dev.off()

pdf("figure-merge-2.pdf",height=8,width=8)
print(p5)
dev.off()

pdf("figure-merge-3.pdf",height=8,width=8)
print(p6)
dev.off()

pdf("figure-merge-4.pdf",height=8,width=8)
print(p7)
dev.off()

pdf("figure-merge-5.pdf",height=8,width=8)
print(p8)
dev.off()

pdf("figure-merge-6.pdf",height=6,width=9)
print(p9)
dev.off()






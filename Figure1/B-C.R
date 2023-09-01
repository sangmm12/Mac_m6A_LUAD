

data <- readRDS("object.integrated.rds")



DefaultAssay(data) <- "RNA"



pdf("Mac-1.pdf",height=8,width=6)
plots <- VlnPlot(data, features = c("TREM2","CD163","MARCO"), group.by = "integrated_snn_res.0.12",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()


pdf("Mac.pdf",height=8,width=8)
p9 <- FeaturePlot(data,
                  reduction = "umap",
                  features = c("TREM2","CD163","MARCO"),
                  label = TRUE)
p9

dev.off()


data_Mac <- subset(data,subset=integrated_snn_res.0.12==1 | integrated_snn_res.0.12==11)


pdf("fig1.pdf",height=4,width=4)
p3 <- DimPlot(object.integrated, label = TRUE, reduction="umap") + NoLegend()
p3
dev.off()



Cell_type <- c("0"="Others",
               "1"="Macrophage",
               "2"="Others",
               "3"="Others",
               "4"="Others",
               "5"="Others",
               "6"="Others",
               "7"="Others",
               "8"="Others",
               "9"="Others",
               "10"="Others",
               "11"="Macrophage",
               "12"="Others",
               "13"="Others")


data[['cell_type']] <- unname(Cell_type[data@meta.data$seurat_clusters])




pdf("figure-celltype.pdf",height=4,width=4)
DimPlot(data,
        reduction = "umap",
        group.by = "cell_type",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()

dev.off()



data$Type <- data$orig.ident
data$Type[grep("LUNG_N",data$orig.ident)] <- "Normal"
data$Type[grep("LUNG_T",data$orig.ident)] <- "Tumor"



data_Mac <- subset(data,subset=cell_type=="Macrophage")



pdf("figure4.pdf",height=4,width=4)
p1 <- DimPlot(data_Mac, reduction = "umap", group.by = "Type")
p1
dev.off()

saveRDS(data_Mac,file="data_Mac.rds")


Idents(data_Mac) <- "Type"
df_Mac <- FindMarkers(data_Mac, ident.1 = "Tumor", ident.2 = "Normal", verbose = FALSE)



data_Mac_Normal <- subset(data_Mac,subset=Type=="Normal")
cells_rankings <- AUCell_buildRankings(scRNA@assays$RNA@data) 

data_Mac_Tumor <- subset(data_Mac,subset=Type=="Tumor")



saveRDS(data_Mac_Normal, file = "data_Mac_Normal.rds")
saveRDS(data_Mac_Tumor, file = "data_Mac_Tumor.rds")



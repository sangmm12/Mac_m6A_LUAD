
setwd("~/projects/m6A_LUAD/result/cor1")

library(Seurat)
library(tidyverse)



cat_construct_metacells <-
  function(seurat_obj,
           k,
           name,
           min_cells = 50,
           max_shared = 10,
           target_metacells = 1000,
           assay = "RNA") {
    if ("cell_type" %in% colnames(seurat_obj@meta.data) &
        "orig.ident" %in% colnames(seurat_obj@meta.data)) {
      set.seed(717)
      #### Set up Seurat object for WGCNA ----
      seurat_obj <- hdWGCNA::SetupForWGCNA(
        seurat_obj,
        gene_select = "fraction",
        # the gene selection approach
        fraction = 0.05,
        # fraction of cells that a gene needs to be expressed in order to be included
        wgcna_name = name # the name of the hdWGCNA experiment
      )
      
      #### Construct metacells ----
      # construct metacells  in each group
      seurat_obj <- hdWGCNA::MetacellsByGroups(
        seurat_obj = seurat_obj,
        group.by = c("orig.ident","cell_type" ),
        # specify the columns in seurat_obj@meta.data to group by
        k = k,
        min_cells = min_cells,
        max_shared = max_shared,
        mode = "sum",
        assay = assay, target_metacells = target_metacells,
        # nearest-neighbors parameter
        ident.group = "cell_type" # set the Idents of the metacell seurat object
      )
      # normalize metacell expression matrix:
      seurat_obj <- hdWGCNA::NormalizeMetacells(seurat_obj)
      return(seurat_obj)
    }
    stop(
      "The column name of cell_type or donor does not exist in the meta.data of your Seurat object, please add it!"
    )
  }

library(hdWGCNA)
selected_cell_type <- "Macrophage"

seurat_obj <- subset(data_Mac,subset=Type=="Tumor")


sub_seurat_obj <- seurat_obj %>%
  subset(cell_type == selected_cell_type)


k <- 50

sub_seurat_obj <- sub_seurat_obj %>%
  cat_construct_metacells(k = k, name = selected_cell_type)

sub_seurat_obj %>% saveRDS(file.path(outdir, "seurat_Tumor.rds"))


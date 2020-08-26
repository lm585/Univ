#
suppressPackageStartupMessages(library(package = "knitr"))
suppressPackageStartupMessages(library(package = "Seurat"))
suppressPackageStartupMessages(library(package = "MAST"))
suppressPackageStartupMessages(library(package = "parallel"))
suppressPackageStartupMessages(library(package = "ComplexHeatmap"))
suppressPackageStartupMessages(library(package = "GSEABase"))
suppressPackageStartupMessages(library(package = "tidyverse"))

setwd("/Users/linyongmao/Documents/antiIL10 monkey/dir-IL10-covid19/dir-belg-scRNAseq/ugent.covid-master/run")
workDir <- dirname(getwd())
opts_chunk$set(tidy = FALSE, fig.path = "../figure/")
options(stringsAsFactors       = FALSE,
        dplyr.summarise.inform = FALSE,
	mc.cores               = 1)
# result will be written in directory called advanced
gseaDir <- file.path(workDir, "advanced")
if (!file.exists(gseaDir)) {
  flag <- dir.create(path = gseaDir)
}

load(file = file.path(workDir, "output/ugent.seurat.RData"))
bComb <- seuratComb[, seuratComb$Sample..type %in% "BAL" &
                         ( seuratComb$COVID19.status == "POS")]
dim(bComb)
# [1] 60688 75016
bComb <- AddMetaData(object   = bComb,
		     metadata = colnames(bComb),
		     col.name = "wellKey")
bAssay <- FromMatrix(exprsArray = list(et = bComb$RNA@data),
		     cData      = bComb@meta.data)
cl <- read.table("BAL_barcodes_meta.subset.POS.bcellsONLY.txt", header = T, sep = "\t") ###6437 cells
for ( i in rownames(cl)) {bComb$HTO_classification[i] = cl[i,1] }
 dim(bComb)
# [1] 60688 75016
unique(bComb$HTO_classification)
bCombBcell <- bComb[, grepl(pattern = "cell", bComb$HTO_classification)]
dim(bCombBcell)
#[1] 60688  1542 cells with b cell annot from Ruth
bCombBcell$RNA@counts <- rbind(bCombBcell$RNA@counts, bCombBcell$ADT@counts)
bCombBcell$RNA@data   <- rbind(bCombBcell$RNA@data, bCombBcell$ADT@data)
 dim(bCombBcell)
splitLS <- SplitObject(object = bCombBcell, split.by = "orig.ident")
splitLS <- lapply(X = splitLS, FUN = function(x) {
    x <- NormalizeData(x, verbose = T)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = T)
})
intFeats <- SelectIntegrationFeatures(splitLS)
intFeats <- setdiff(intFeats, row.names(bCombBcell@assays$ADT))
splitAnchors <- FindIntegrationAnchors(object.list = splitLS,
                                         dims = 1:30, k.filter = 35, anchor.features = c(intFeats, row.names(bCombBcell@assays$ADT)))
#####Retained 115 anchors#####
bCombBcell.integrated <- IntegrateData(anchorset = splitAnchors, dims = 1:30)
bCombBcell.integrated <- ScaleData(bCombBcell.integrated, verbose = T)
bCombBcell.integrated <- RunPCA(bCombBcell.integrated, npcs = 30, verbose = T)
bCombBcell.integrated <- RunUMAP(bCombBcell.integrated, reduction = "pca", dims = 1:30)
DimPlot(bCombBcell.integrated, reduction = "umap", group.by = "Severity.Score.worse")
DimPlot(bCombBcell.integrated, reduction = "umap", group.by = "orig.ident")
bCombBcell.integrated <- FindNeighbors(object = bCombBcell.integrated)
bCombBcell.integrated <- FindClusters(object = bCombBcell.integrated)
#####Number of nodes: 1542, Number of edges: 48583, Number of communities: 10
DimPlot(object = bCombBcell.integrated, reduction = "umap", label = TRUE, repel = T ) + NoLegend()
DimPlot(object = bCombBcell.integrated, reduction = "umap", label = TRUE, repel = T, cols = c("purple", "yellow", "orange", "green", "black", "blue", "red", "yellow", "red", "blue") )
DefaultAssay(bCombBcell.integrated) <- "RNA"

a1 <- c("GRCh38.99-------------MZB1", "GRCh38.99-------------TNFSF13B", "GRCh38.99-------------XBP1", "GRCh38.99-------------IGHM", "GRCh38.99-------------IL10RA", "CD5")
a.emory <- c("CD19", "CD27", "CD37", "CD38", "CD40", "CD44", "GRCh38.99-------------BACH2", "GRCh38.99-------------BCL6", "GRCh38.99-------------IRF4", "GRCh38.99-------------MS4A1", "GRCh38.99-------------PAX5", "GRCh38.99-------------PRDM1", "GRCh38.99-------------SDC1", "GRCh38.99-------------TNFRSF17", "GRCh38.99-------------XBP1")
a2 <- c("GRCh38.99-------------CD27", "GRCh38.99-------------CD37", "GRCh38.99-------------CD38")
a.proposal <- c("CD11c", "CD138-A0055", "CD138-A0831", "CD19", "CD21", "CD27-A0154", "CD27-A0191", "CD38-A0389", "CD38-A0410", "IgA", "IgD", "IgG", "IgM")
a.TF <- c("GRCh38.99-------------PRDM1", "GRCh38.99-------------IRF4", "GRCh38.99-------------XBP1", "GRCh38.99-------------PAX5", "GRCh38.99-------------BCL6", "GRCh38.99-------------BACH2")
FeaturePlot(bCombBcell.integrated, features = c(a.proposal), slot = "data", pt.size = 0.05, label=T, label.size = 3)
FeaturePlot(bCombBcell.integrated, features = c(a1, a.emory), slot = "data", pt.size = 0.05, label=T, label.size = 3)
 FeaturePlot(bCombBcell.integrated, features = c("CD20"), slot = "data")

VlnPlot(bCombBcell.integrated, features = c(a.TF), slot = "data", pt.size = 0.05, group.by = "seurat_clusters", size.title.use = 10)
patch <- VlnPlot(bCombBcell.integrated, features = c(a.TF), slot = "data", pt.size = 0.05, group.by = "seurat_clusters")
### install.packages('patchwork')
library(patchwork)
patch <- VlnPlot(bCombBcell.integrated, features = c(a.TF), slot = "data", pt.size = 0.05, group.by = "seurat_clusters", combine = F)
p1 <- list() 
for (i in seq_along(patch)){
    #Change x and y tick label font size.
    p1[[i]] = patch[[i]] + theme(title = element_text(size = 8) ) +  NoLegend() 
}
library(cowplot)
plot_grid(plotlist = p1, ncol = 3)


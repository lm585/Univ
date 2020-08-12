#
load(file = file.path(workDir, "output/ugent.seurat.RData"))
bComb <- seuratComb[, seuratComb$Sample..type %in% "BAL" &
                         ( seuratComb$COVID19.status == "POS")]
cl <- read.table("BAL_barcodes_meta.subset.POS.bcellsONLY.txt", header = T, sep = "\t") ###6437 cells
for ( i in rownames(cl)) {bComb$HTO_classification[i] = cl[i,1] }
 dim(bComb)
# [1] 60688 75016
unique(bComb$HTO_classification)

bCombBcell <- bComb[, grepl(pattern = "cell", bComb$HTO_classification)]
dim(bCombBcell)

splitLS <- SplitObject(object = bCombBcell, split.by = "orig.ident")
splitLS <- lapply(X = splitLS, FUN = function(x) {
    x <- NormalizeData(x, verbose = T)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = T)
})
intFeats <- SelectIntegrationFeatures(splitLS)
splitAnchors <- FindIntegrationAnchors(object.list = splitLS,
                                         dims = 1:30, k.filter = 35, anchor.features = c(intFeats))
####Retained 101 anchors
bCombBcell.integrated <- IntegrateData(anchorset = splitAnchors, dims = 1:30)
bCombBcell.integrated <- ScaleData(bCombBcell.integrated, verbose = T)
bCombBcell.integrated <- RunPCA(bCombBcell.integrated, npcs = 30, verbose = T)
bCombBcell.integrated <- RunUMAP(bCombBcell.integrated, reduction = "pca", dims = 1:30)
DimPlot(bCombBcell.integrated, reduction = "umap", group.by = "HTO_classification")

scaledMat <- as.matrix(bCombBcell.integrated$integrated@data)
rownames(scaledMat) <- gsub(pattern = ".+---", replacement = "",
                            rownames(scaledMat))

hpca <- HumanPrimaryCellAtlasData()
hpcaB <- hpca[, hpca@colData$label.main %in% "B_cell" &
		     hpca@colData$label.fine %in% c("B_cell",
"B_cell:Germinal_center", "B_cell:Plasma_cell",        "B_cell:Naive",             
 "B_cell:Memory",             "B_cell:CXCR4+_centroblast", "B_cell:CXCR4-_centrocyte",  "B_cell:immature" 
						    )]

predSubsetB <- SingleR(test    = scaledMat,
		       ref     = hpcaB,
		       labels  = hpcaB$label.fine,
		       )
bCombBcell.integrated <- AddMetaData(object = bCombBcell.integrated,
			 metadata = predSubsetB$labels,
			 col.name = "bcellSubset")
DimPlot(bCombBcell.integrated, reduction = "umap", group.by = "bcellSubset")


suppressPackageStartupMessages(library(package = "knitr"))
suppressPackageStartupMessages(library(package = "Seurat"))
suppressPackageStartupMessages(library(package = "MAST"))
suppressPackageStartupMessages(library(package = "parallel"))
suppressPackageStartupMessages(library(package = "ComplexHeatmap"))
suppressPackageStartupMessages(library(package = "GSEABase"))
suppressPackageStartupMessages(library(package = "tidyverse"))

setwd("/Users/linyongmao/Documents/antiIL10 monkey/dir-IL10-covid19/dir-belg-scRNAseq/ugent.covid-master/run")

pbmc <-readRDS("/Users/linyongmao/Documents/antiIL10 monkey/dir-IL10-covid19//dir-belg-scRNAseq/ugent.covid-master/output/seuratObj_PBMC_CZI_meta.rds")
addAnot <- read.table("PBMC_CZI.meta-annot.comb.txt", header = T, sep = "\t")
####Age Age.num Sex Sex.num Race Race.num Ethnicity Ethnicity.num BMI BMI.num
pbmc$Age <- addAnot$orig.ident.1 ###COV003_AAACCCAAGACTTCGT-1
pbmc$Age.num <- addAnot$notApplicable.2 ###covid status, POS
pbmc$Sex <- addAnot$notApplicable.3 ####time point, T1
pbmc$Sex.num <- addAnot$notApplicable.21 ####Critical Moderate Severe
pbmc$Ethnicity <- addAnot$notApplicable.17 ####Monaco main label
pbmc$Ethnicity.num <- addAnot$notApplicable.18 ###Monaco fine label
for(i in seq(100000,300000,20000)) {print(colnames(pbmc)[i] == pbmc$Age[i])}

pbmc.pos <- pbmc[, pbmc$Age.num == "POS" ]  ###24 cells in COV029
dim(pbmc.pos)
df <- as.data.frame(pbmc.pos$ADT@counts)
rownames(df) <- paste("adt--", rownames(df), sep = "")
pbmc.pos[["newADT"]] <- CreateAssayObject(df)

pbmc.pos$RNA@counts <- rbind(pbmc.pos$RNA@counts, pbmc.pos$newADT@counts)
pbmc.pos$RNA@data   <- rbind(pbmc.pos$RNA@data, pbmc.pos$newADT@data) 
dim(pbmc.pos)
### 23373 11048
splitLS <- SplitObject(object = pbmc.pos, split.by = "orig.ident")
splitLS <- lapply(X = splitLS, FUN = function(x) {
    x <- NormalizeData(x, verbose = T)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = T)
})
intFeats <- SelectIntegrationFeatures(splitLS)
intFeats <- setdiff(intFeats, row.names(pbmc.pos@assays$newADT))
splitAnchors <- FindIntegrationAnchors(object.list = splitLS,
                                         dims = 1:30, k.filter = 200, anchor.features = c(intFeats, row.names(pbmc.pos@assays$newADT)))
###Retained 3750 anchors
pbmc.pos.integrated <- IntegrateData(anchorset = splitAnchors, dims = 1:30)
pbmc.pos.integrated <- ScaleData(pbmc.pos.integrated, verbose = T)
pbmc.pos.integrated <- RunPCA(pbmc.pos.integrated, npcs = 30, verbose = T)
pbmc.pos.integrated <- RunUMAP(pbmc.pos.integrated, reduction = "pca", dims = 1:30)
## cols = c("purple", "yellow", "orange", "green", "black", "blue", "red", "yellow", "red", "blue")
## cols = c("red","orange","yellow","green","rosybrown1", "blue", "purple","black", "skyblue1","slategray1")

write.table(pbmc.pos.integrated$Ethnicity, file = "temp.xls", sep = "\t")
cut -f 1 temp.xls > temp.cellID.xls
./scrnaseq-metaannot-combine temp.cellID.xls PBMC_barcodes_meta.subset.annot.Ruth.txt out
cat  out | cut -f 1,17,19 > temp ###atlas_main monaco_main
paste temp.xls temp | awk 'BEGIN {FS = "\t"; OFS = "\t"} $2 != $5' | wc -l
paste temp.xls temp |cut -f 2,4 | sort | uniq -c | sort -n | grep "B cells"
paste temp.xls temp |cut -f 2,4 > temp-339
vi temp-339
cut -f 1 temp-339 > monaco-atlas-T-B-cell.txt


atlasTcell.bcell <- read.table("monaco-atlas-T-B-cell.txt", header = T, sep = "\t")
pbmc.pos.integrated$BMI <- atlasTcell.bcell$x
DimPlot(pbmc.pos.integrated, reduction = "umap", group.by = "BMI", cols = c("red", "brown", "blueviolet", "grey", "orange","yellow","green","rosybrown1", "blue", "purple","black", "skyblue1","slategray1"))
DimPlot(pbmc.pos.integrated, reduction = "umap", group.by = "BMI", cols = c("blue", "green", "purple","grey", "grey", "grey", "grey","grey", "grey","grey", "grey", "grey", "grey","grey"))
fea = c( "adt--CD19", "adt--CD20")
fea = c("adt--CD3-A0034", "adt--CD3-A0049", "adt--CD4-A0045", "adt--CD4-A0072", "adt--CD4-A0922", "adt--CD8")
fea = c("adt--IgG", "adt--IgD", "adt--CD21", "adt--CD40", "adt--CD44", "adt--CD44-mh") 
fea = c("adt--CD14", "adt--CD14-A0051", "adt--CD14-A0081", "adt--CD16")
VlnPlot(pbmc.pos.integrated, features=fea, slot = "data", pt.size=0, group.by="BMI")

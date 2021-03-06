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
###[1] 60688 75016
bComb <- AddMetaData(object   = bComb,
		     metadata = colnames(bComb),
		     col.name = "wellKey")
###bAssay <- FromMatrix(exprsArray = list(et = bComb$RNA@data),
###		     cData      = bComb@meta.data)

cl <- read.table("BAL_barcodes_meta.subset.POS.bcellsONLY.txt", header = T, sep = "\t") ###6437 cells
for ( i in rownames(cl)) {bComb$HTO_classification[i] = cl[i,1] }
 dim(bComb)
# [1] 60688 75016
unique(bComb$HTO_classification)
bCombBcell <- bComb[, grepl(pattern = "cell", bComb$HTO_classification)]
dim(bCombBcell)
#[1] 60688  1542 cells with b cell annot from Ruth

bCombBcell <- NormalizeData(object = bCombBcell, assay = "ADT", normalization.method = "LogNormalize",verbose = TRUE)

 dim(bCombBcell)
###[1] 60988  1542
clustNum = 9
clustNumStr = as.character(clustNum)

cl <- read.table("bal.bcell.covidPos.clust.startBelg.xls", header = T, sep = "\t") ### cells
for ( i in rownames(cl)) {bCombBcell$HTO_classification[i] = cl[i,1] }
dim(bCombBcell)
###[1] 60988  1542
unique(bCombBcell$HTO_classification)

bCombBcell <- AddMetaData(object   = bCombBcell,
                     metadata = colnames(bCombBcell),
                     col.name = "wellKey")

bCombBcellplasmablast  <- bCombBcell[, bCombBcell$HTO_classification == 1 | bCombBcell$HTO_classification == 3 | bCombBcell$HTO_classification == 5 | bCombBcell$HTO_classification == 8]
dim(bCombBcellplasmablast)
#############################################
sev <- c("Moderate", "Severe", "Critical")
for(index in 1:3) {
 for(j in 1:3) {
  if(index != j) {
   bCombBcell2clusts <- bCombBcellplasmablast[,bCombBcellplasmablast$Severity.Score.worse == sev[index] | bCombBcellplasmablast$Severity.Score.worse == sev[j] ]
   for ( i in colnames(bCombBcell2clusts)) {
   if(bCombBcell2clusts$Severity.Score.worse[i] == sev[index]) {bCombBcell2clusts$Severity.Score.worse[i] = "targetClust"}
   else {bCombBcell2clusts$Severity.Score.worse[i] = "others"}
   }
   print(sum(bCombBcell2clusts$Severity.Score.worse == "targetClust"))
   print(sum(bCombBcell2clusts$Severity.Score.worse == "others"))
   bAssay <- FromMatrix(exprsArray = list(et = bCombBcell2clusts$RNA@data),
                      cData      = bCombBcell2clusts@meta.data) ######, check_sanity = FALSE)
  fit <- zlm(formula  = ~Severity.Score.worse+orig.ident,
           sca      = bAssay,
           parallel = FALSE)

  summaryCond <- MAST::summary(fit, doLRT = "Severity.Score.worsetargetClust")

  top <- summaryCond$datatable %>%
  filter(contrast %in% "Severity.Score.worsetargetClust" &
           # hurdle test
           component %in% c("H", "logFC")) %>%
  pivot_longer(cols = c("Pr(>Chisq)", "coef")) %>%
  filter(!is.na(value)) %>%
  select(primerid, contrast, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  dplyr::rename(`logFC` = `coef`) %>%
  mutate(adjp = p.adjust(`Pr(>Chisq)`, method = "BH"))

 top %>%
  mutate(primerid = gsub(pattern = ".+---", replacement = "", primerid)) %>%
  top_n(n = 20, wt = -adjp) %>%
  print()

 top50 <- top_n(top, n = 98765, wt = -adjp)
 file = paste("bal.bcell.startBelg.plasmablastSev", sev[index], "vs", sev[j], ".", sev[index], "-Pos.txt", sep = "")
 print(file)
 write.table(top50, file, sep="\t")

####ADT log nomaliized already
 bAssay <- FromMatrix(exprsArray = list(et = bCombBcell2clusts$ADT@data),
                      cData      = bCombBcell2clusts@meta.data) ####ADT log nomaliize
 fit <- zlm(formula  = ~Severity.Score.worse+orig.ident,
           sca      = bAssay,
           parallel = FALSE)

  summaryCond <- MAST::summary(fit, doLRT = "Severity.Score.worsetargetClust")

  top <- summaryCond$datatable %>%
  filter(contrast %in% "Severity.Score.worsetargetClust" &
           # hurdle test
           component %in% c("H", "logFC")) %>%
  pivot_longer(cols = c("Pr(>Chisq)", "coef")) %>%
  filter(!is.na(value)) %>%
  select(primerid, contrast, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  dplyr::rename(`logFC` = `coef`) %>%
  mutate(adjp = p.adjust(`Pr(>Chisq)`, method = "BH"))

 top %>%
  mutate(primerid = gsub(pattern = "hello---world", replacement = "", primerid)) %>%
  top_n(n = 20, wt = -adjp) %>%
  print()

 top50 <- top_n(top, n = 98765, wt = -adjp)
 file = paste("bal.bcell.startBelg.plasmablastSev", sev[index], "vs", sev[j], ".", sev[index], "-Pos.ADT.txt", sep = "")

 print(file)
 write.table(top50, file, sep="\t")
  }
 }
} ##end of for(index in 1:3)



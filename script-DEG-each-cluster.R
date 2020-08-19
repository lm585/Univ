suppressPackageStartupMessages(library(package = "knitr"))
suppressPackageStartupMessages(library(package = "Seurat"))
suppressPackageStartupMessages(library(package = "MAST"))
suppressPackageStartupMessages(library(package = "parallel"))
suppressPackageStartupMessages(library(package = "ComplexHeatmap"))
suppressPackageStartupMessages(library(package = "GSEABase"))
suppressPackageStartupMessages(library(package = "tidyverse"))

setwd("/Users/linyongmao/Documents/antiIL10 monkey/dir-IL10-covid19/dir-belg-scRNAseq/ugent.covid-master/run2")
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

clustNum = 6
load(file = file.path(workDir, "output/ugent.seurat.RData"))
ls()
bComb <- seuratComb[, seuratComb$Sample..type %in% "BAL" &
                         ( seuratComb$COVID19.status == "POS")]
cl <- read.table("bal.bcell.covidPos.clust.startBelg.xls", header = T, sep = "\t") ### cells
for ( i in rownames(cl)) {bComb$HTO_classification[i] = cl[i,1] }
 dim(bComb)
# [1] 60688 75016
unique(bComb$HTO_classification)

bCombBcell <- bComb[,  bComb$HTO_classification == clustNum] 
dim(bCombBcell)

bCombBcell <- bComb[,  bComb$HTO_classification == clustNum &
                       (bComb$Severity.Score.worse == "Severe" | bComb$Severity.Score.worse == "Moderate")]
bCombBcell <- AddMetaData(object   = bCombBcell,
                     metadata = colnames(bCombBcell),
                     col.name = "wellKey")
bAssay <- FromMatrix(exprsArray = list(et = bCombBcell$RNA@data),
                     cData      = bCombBcell@meta.data)
dim(bCombBcell)

fit <- zlm(formula  = ~Severity.Score.worse+orig.ident,
           sca      = bAssay,
           parallel = FALSE)
summaryCond <- MAST::summary(fit, doLRT = "Severity.Score.worseSevere")

top <- summaryCond$datatable %>%
  filter(contrast %in% "Severity.Score.worseSevere" &
           # hurdle test
           component %in% c("H", "logFC")) %>%
  pivot_longer(cols = c("Pr(>Chisq)", "coef")) %>%
  filter(!is.na(value)) %>%
  select(primerid, contrast, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  dplyr::rename(`logFC` = `coef`) %>%
  mutate(adjp = p.adjust(`Pr(>Chisq)`, method = "BH"))

# print top 10 DEG
top %>%
  mutate(primerid = gsub(pattern = ".+---", replacement = "", primerid)) %>%
  top_n(n = 20, wt = -adjp) %>%
  print()

top50 <- top_n(top, n = 98765, wt = -adjp)
write.table(top50, "bal.bcell.startBelg.cluster6.diff2-1.SvM.Spos.txt", sep="\t")
##################################################################################################################

bCombBcell <- bComb[,  bComb$HTO_classification == clustNum & 
                       (bComb$Severity.Score.worse == "Critical" | bComb$Severity.Score.worse == "Moderate")]
bCombBcell <- AddMetaData(object   = bCombBcell,
		     metadata = colnames(bCombBcell),
		     col.name = "wellKey")
bAssay <- FromMatrix(exprsArray = list(et = bCombBcell$RNA@data),
		     cData      = bCombBcell@meta.data)
dim(bCombBcell)

fit <- zlm(formula  = ~Severity.Score.worse+orig.ident,
           sca      = bAssay,
           parallel = FALSE)
summaryCond <- MAST::summary(fit, doLRT = "Severity.Score.worseModerate")

top <- summaryCond$datatable %>%
  filter(contrast %in% "Severity.Score.worseModerate" &
           # hurdle test
           component %in% c("H", "logFC")) %>%
  pivot_longer(cols = c("Pr(>Chisq)", "coef")) %>%
  filter(!is.na(value)) %>%
  select(primerid, contrast, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  dplyr::rename(`logFC` = `coef`) %>%
  mutate(adjp = p.adjust(`Pr(>Chisq)`, method = "BH"))

# print top 10 DEG
top %>%
  mutate(primerid = gsub(pattern = ".+---", replacement = "", primerid)) %>%
  top_n(n = 20, wt = -adjp) %>%
  print()

top50 <- top_n(top, n = 98765, wt = -adjp)
write.table(top50, "bal.bcell.startBelg.cluster6.diff2-1.CvM.Mpos.txt", sep="\t")

##################################################################################################################

bCombBcell <- bComb[,  bComb$HTO_classification == clustNum &
                       (bComb$Severity.Score.worse == "Critical" | bComb$Severity.Score.worse == "Severe")]
bCombBcell <- AddMetaData(object   = bCombBcell,
                     metadata = colnames(bCombBcell),
                     col.name = "wellKey")
bAssay <- FromMatrix(exprsArray = list(et = bCombBcell$RNA@data),
                     cData      = bCombBcell@meta.data)
dim(bCombBcell)

fit <- zlm(formula  = ~Severity.Score.worse+orig.ident,
           sca      = bAssay,
           parallel = FALSE)
summaryCond <- MAST::summary(fit, doLRT = "Severity.Score.worseSevere")

top <- summaryCond$datatable %>%
  filter(contrast %in% "Severity.Score.worseSevere" &
           # hurdle test
           component %in% c("H", "logFC")) %>%
  pivot_longer(cols = c("Pr(>Chisq)", "coef")) %>%
  filter(!is.na(value)) %>%
  select(primerid, contrast, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  dplyr::rename(`logFC` = `coef`) %>%
  mutate(adjp = p.adjust(`Pr(>Chisq)`, method = "BH"))

# print top 10 DEG
top %>%
  mutate(primerid = gsub(pattern = ".+---", replacement = "", primerid)) %>%
  top_n(n = 20, wt = -adjp) %>%
  print()

top50 <- top_n(top, n = 98765, wt = -adjp)
write.table(top50, "bal.bcell.startBelg.cluster6.diff2-1.CvS.Spos.txt", sep="\t")


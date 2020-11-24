#
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
pbmc$atlas.main <- addAnot$notApplicable.15 ###cell atlas main
for(i in seq(100000,300000,20000)) {print(colnames(pbmc)[i] == pbmc$Age[i])}

pbmc.bcell.pos <- pbmc[, pbmc$Age.num == "POS" & pbmc$Ethnicity == "B cells" & pbmc$orig.ident != "COV029" & grepl(pattern = "B_cell", pbmc$atlas.main) ]  ###24 cells in COV029
dim(pbmc.bcell.pos)
df <- as.data.frame(pbmc.bcell.pos$ADT@counts)
rownames(df) <- paste("adt--", rownames(df), sep = "")
pbmc.bcell.pos[["newADT"]] <- CreateAssayObject(df)
pbmc.bcell.pos <- NormalizeData(object = pbmc.bcell.pos, assay = "newADT",  normalization.method = "LogNormalize",verbose = TRUE)
dim(pbmc.bcell.pos)
clustNum = 4
clustNumStr = as.character(clustNum)

cl <- read.table("pbmc.bcell.pos.integ.10132cells.cellID-clusters.txt", header = T, sep = "\t") ### cells
for ( i in rownames(cl)) {pbmc.bcell.pos$Smoking[i] = cl[i,1] }
dim(pbmc.bcell.pos)
###[1] 60988  1542
unique(pbmc.bcell.pos$Smoking)

pbmc.bcell.pos <- AddMetaData(object   = pbmc.bcell.pos,
                     metadata = colnames(pbmc.bcell.pos),
                     col.name = "wellKey")

######start for loop
for(index in 13:13) {
 if (index != 4) {
  query = as.character(index)
  pbmcBcell2clusts <- pbmc.bcell.pos[, pbmc.bcell.pos$Smoking == clustNum | pbmc.bcell.pos$Smoking == index]
  for ( i in colnames(pbmcBcell2clusts)) {
  if(pbmcBcell2clusts$Smoking[i] ==  clustNum) {pbmcBcell2clusts$Smoking[i] = "targetClust"}
  else {pbmcBcell2clusts$Smoking[i] = "others"}
  }
  print(sum(pbmcBcell2clusts$Smoking == "targetClust"))
  print(sum(pbmcBcell2clusts$Smoking == "others"))
  bAssay <- FromMatrix(exprsArray = list(et = pbmcBcell2clusts$RNA@data),
                      cData      = pbmcBcell2clusts@meta.data) ######, check_sanity = FALSE)
  fit <- zlm(formula  = ~Smoking+orig.ident,
           sca      = bAssay,
           parallel = FALSE)

  summaryCond <- MAST::summary(fit, doLRT = "SmokingtargetClust")

  top <- summaryCond$datatable %>%
  filter(contrast %in% "SmokingtargetClust" &
           # hurdle test
           component %in% c("H", "logFC")) %>%
  pivot_longer(cols = c("Pr(>Chisq)", "coef")) %>%
  filter(!is.na(value)) %>%
  select(primerid, contrast, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  dplyr::rename(`logFC` = `coef`) %>%
  mutate(adjp = p.adjust(`Pr(>Chisq)`, method = "BH"))

 top %>%
  mutate(primerid = gsub(pattern = "helloWorld", replacement = "", primerid)) %>%
  top_n(n = 20, wt = -adjp) %>%
  print()

 top50 <- top_n(top, n = 98765, wt = -adjp)
 file = paste("pbmc.bcell.startBelg.cluster", clustNum, "vs", index, ".clust", clustNum, "Pos.txt", sep = "")
 print(file)
 write.table(top50, file, sep="\t")

####ADT log nomaliized already
 bAssay <- FromMatrix(exprsArray = list(et = pbmcBcell2clusts$newADT@data),
                      cData      = pbmcBcell2clusts@meta.data) ####newADT log nomaliize
 fit <- zlm(formula  = ~Smoking+orig.ident,
           sca      = bAssay,
           parallel = FALSE)
  
  summaryCond <- MAST::summary(fit, doLRT = "SmokingtargetClust")
  
  top <- summaryCond$datatable %>%
  filter(contrast %in% "SmokingtargetClust" &
           # hurdle test
           component %in% c("H", "logFC")) %>% 
  pivot_longer(cols = c("Pr(>Chisq)", "coef")) %>%
  filter(!is.na(value)) %>%
  select(primerid, contrast, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  dplyr::rename(`logFC` = `coef`) %>%
  mutate(adjp = p.adjust(`Pr(>Chisq)`, method = "BH"))
 
 top %>%
  mutate(primerid = gsub(pattern = "helloWorld", replacement = "", primerid)) %>%
  top_n(n = 20, wt = -adjp) %>%
  print()
 
 top50 <- top_n(top, n = 98765, wt = -adjp)
 file = paste("pbmc.bcell.startBelg.cluster", clustNum, "vs", index, ".clust", clustNum, "Pos.newADT.txt", sep = "")
 print(file)
 write.table(top50, file, sep="\t")


 } ####if (index != 4)
} ###end for loop



####newADT log nomaliize
#######heatmap

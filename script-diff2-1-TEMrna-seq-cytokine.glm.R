#
library(edgeR)

setwd("/Users/linyongmao/Documents/dir-star-hiv")

h <- read.delim("temp-cyto-p2-trans-with-rna-seq", header = T)
for(index in 2:dim(h)[2]) {
 hiv <- log10(h[,index])
 x <- read.delim("temp-gene-count-tem-4samp-star568-10.txt", row.names=1)
 y <- DGEList(counts=x)
 y <- calcNormFactors(y)
 keep <-rowSums(cpm(y)>0) >= 3
 y<-y[keep,]
 print(dim(y))
 design<-model.matrix(~hiv) 
 y <- estimateGLMCommonDisp(y,design) 
 y <- estimateGLMTrendedDisp(y,design) 
 y <- estimateGLMTagwiseDisp(y,design) 
 fit<-glmFit(y,design) 

 lrt.2vs1 <- glmLRT(fit, coef = "hiv")
 top2v1 <- topTags(lrt.2vs1, n=150000) 
 # file = paste("pbmc.bcell.startBelg.cluster", clustNum, "vs", index, ".clust", clustNum, "Pos.txt", sep = "")
 file = paste("star-CPM.tem.GLM.", colnames(h)[index], ".diff2-1.txt", sep = "")
 print(file)
 write.table(top2v1, file, sep="\t")
} ###end of index loop

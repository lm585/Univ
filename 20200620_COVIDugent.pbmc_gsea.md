---
title: Geneset analysis of B cells in PBMCs
author: Slim Fourati
date: "27 June, 2020"
output: github_documents
---

Load packages

```r
suppressPackageStartupMessages(library(package = "knitr"))
suppressPackageStartupMessages(library(package = "Seurat"))
suppressPackageStartupMessages(library(package = "MAST"))
suppressPackageStartupMessages(library(package = "ComplexHeatmap"))
suppressPackageStartupMessages(library(package = "GSEABase"))
suppressPackageStartupMessages(library(package = "tidyverse"))
```

Set session options

```r
workDir <- dirname(getwd())
opts_chunk$set(tidy = FALSE, fig.path = "../figure/")
options(stringsAsFactors       = FALSE,                                                   
        dplyr.summarise.inform = FALSE)
# result will be written in directory called advanced
gseaDir <- file.path(workDir, "advanced")
if (!file.exists(gseaDir)) {
  flag <- dir.create(path = gseaDir)
}
```

Load scRNAseq data

```r
load(file = file.path(workDir, "output/ugent.seurat.RData"))
```

Use MAST do differential expression (adjust for sample of origin)

```r
# subset on PBMCs B-cell at T1 (where we have COVID19 pos vs neg)
bComb <- seuratComb[, seuratComb$singler.subset %in% "B_cell" &
			seuratComb$Sample..type %in% "PBMC" &
			grepl(pattern = "T1", seuratComb$Patient.number)]
bAssay <- FromMatrix(exprsArray = list(et = bComb$RNA@counts),
		     cData      = bComb@meta.data,
		     check_sanity = FALSE)
fit <- zlm(formula = ~COVID19.status+orig.ident, sca = bAssay)
summaryCond <- MAST::summary(fit, doLRT = "COVID19.statusPOS")
top <- summaryCond$datatable %>%
  filter(contrast %in% "COVID19.statusPOS" &
	   # hurdle test
	   component %in% c("H", "logFC")) %>% 
  pivot_longer(cols = c("Pr(>Chisq)", "coef")) %>%
  filter(!is.na(value)) %>%
  select(primerid, contrast, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  dplyr::rename(`logFC` = `coef`) %>%
  mutate(adjp = p.adjust(`Pr(>Chisq)`, method = "BH"))
```


```r
# print top 10 DEG
top %>%
  mutate(primerid = gsub(pattern = ".+---", replacement = "", primerid)) %>%
  top_n(n = 10, wt = -adjp) %>%
  print()
```

```
## # A tibble: 10 x 5
##    primerid   contrast          `Pr(>Chisq)`   logFC     adjp
##    <chr>      <fct>                    <dbl>   <dbl>    <dbl>
##  1 ACADVL     COVID19.statusPOS     2.15e- 9  1.80   1.31e- 5
##  2 AL450405.1 COVID19.statusPOS     3.83e-11 -2.38   4.65e- 7
##  3 HLA-A      COVID19.statusPOS     3.48e-10 28.2    2.64e- 6
##  4 MYDGF      COVID19.statusPOS     1.66e-11 14.5    2.51e- 7
##  5 MZB1       COVID19.statusPOS     8.25e-11 50.6    8.34e- 7
##  6 RGPD5      COVID19.statusPOS     1.13e-19 -0.0601 6.85e-15
##  7 RPS4Y1     COVID19.statusPOS     2.05e-16 -2.54   6.21e-12
##  8 SELENOS    COVID19.statusPOS     1.22e-10  5.48   1.05e- 6
##  9 XBP1       COVID19.statusPOS     7.73e-10  9.25   5.21e- 6
## 10 XIST       COVID19.statusPOS     5.37e-15  2.92   1.09e-10
```

Print heatmap based on top50 DEGs

```r
top50 <- top_n(top, n = 50, wt = -adjp)

mat <- as.matrix(bComb$RNA@counts[top50$primerid, ])
flag <- rowMeans(mat > 0) >= 0.2
mat <- mat[flag, ]
top50 <- top50[flag, ]
mat <- log2(mat + 0.25)
rownames(mat) <- gsub(pattern = ".+---", replacement = "", rownames(mat))
ha <- HeatmapAnnotation(ident = bComb@meta.data$orig.ident)
Heatmap(matrix = mat,
	column_split = bComb@meta.data$COVID19.status,
	top_annotation = ha,
	row_split = c("DN", "UP")[(top50$logFC > 0) + 1],
	show_column_names = FALSE,
	name = "log2(counts)")
```

![plot of chunk pbmc-b-heatmap-top50](../figure/pbmc-b-heatmap-top50-1.png)

Run GSEA

```r
unzip(zipfile = file.path(workDir, "utils/GSEA_4.0.3.zip"),
      exdir = gseaDir)
gseaCli <- file.path(gseaDir, "GSEA_4.0.3/gsea-cli.sh")
gmtFile <- file.path(workDir, "utils/msigdb.v7.1.symbols.gmt")
```


```r
modelName <- "PBMC_BCELL_T1"
contrast <- "COVID19.statusPOS"
ranked <- top %>%
  filter(!is.na(logFC)) %>%
  mutate(score = sign(logFC) * -log10(`Pr(>Chisq)`),
	 gene_name = gsub(pattern = ".+---",
			  replacement = "",
			  primerid)) %>%
  select(gene_name, score)
rnkFile <- paste0("gsea_", modelName, "_", contrast,  ".rnk")
rnkFile <- make.names(rnkFile)
rnkFile <- file.path(workDir, "advanced", rnkFile)
write(paste(c("#", colnames(ranked)), collapse = " "), file = rnkFile)
write_tsv(ranked,
	  path      = rnkFile,
	  append    = TRUE,
	  col_names = FALSE)
rnkList <- data.frame(modelName = modelName,
		      contrast  = contrast,
		      rnk       = rnkFile)
```


```r
gseaRnk <- rnkList$rnk
logFileName <- gsub(pattern = "rnk$", replacement = "log", gseaRnk)
gseaRpt <- paste(gsub(pattern     = "rnk",
                          replacement = "",
                          basename(gseaRnk)),
                     gsub(pattern     = "^([^\\.]+).+$",
                          replacement = "\\1",
                          basename(gmtFile)),
                     sep = ".")
gseaCall <- paste("bash",
		  gseaCli,
		  "GSEAPreranked",
		  "-gmx",
		  gmtFile,
		  "-collapse No_Collapse",
		  "-mode Max_probe",
		  "-norm None",
		  "-nperm 1000",
		  "-rnk",
		  gseaRnk,
		  "-scoring_scheme weighted",
		  "-rpt_label",
		  gseaRpt,
		  "-create_svgs false",
		  "-include_only_symbols true",
		  "-make_sets true",
		  "-plot_top_x 1",
		  "-rnd_seed 111",
		  "-set_max 2000",
		  "-set_min 15",
		  "-zip_report false",
		  "-out",
		  gseaDir,
		  ">",
		  logFileName)
gseaIntern <- system(command       = gseaCall,
                       intern        = TRUE,
		     ignore.stderr = TRUE)

gseaIndex <- data.frame(rnk = gseaRnk, rpt = file.path(gseaDir, gseaRpt))
gseaIndex <- merge(rnkList, gseaIndex, by = "rnk")
```


```r
# remove previous gsea run from the advanced directory
dirLS <- list.dirs(path = file.path(workDir, "advanced"), recursive = FALSE)
dirLS <- cbind(directory = dirLS,
               rpt       = gsub(pattern = ".GseaPreranked.+$",
                   replacement = "",
                   dirLS))
gseaIndex <- merge(gseaIndex, dirLS, by = "rpt")
```


```r
# delete temporary directory create during gsea run
dirName <- tolower(format(Sys.Date(), "%b%d"))
flag <- file.remove(dirName)
```


```r
# read gsea output directories
gseaOutput <- apply(gseaIndex, MARGIN = 1, FUN = function(gseaRun) {
  print(gseaRun[["modelName"]])
  gseaDir <- gseaRun[["directory"]]
  # read rpt file in gsea output directory
  rptFile <- list.files(path = gseaDir, pattern = "rpt", full.names = TRUE)
  rpt <- read_tsv(file      = rptFile,
                  col_names = c("type", "name", "value"))
  # read gmt file
  gmxFile <- rpt$value[rpt$name %in% "gmx"]
  cNames <- count_fields(file = gmxFile, tokenizer = tokenizer_tsv()) %>%
    max() %>%
    seq(from = 1) %>%
    as.character()
  gmx <- read_tsv(file = gmxFile, col_names = cNames)
  # remove geneset name and description column
  gsNames <- toupper(gmx$"1")
  gmx <- apply(select(gmx, -(1:2)), MARGIN = 1, FUN = function(x) {
    return(value = setdiff(unname(x), NA))
  })
  names(gmx) <- gsNames
  # read result files
  resFile <- grep(pattern = "gsea_report.*xls",
                  dir(path = gseaDir, full.names = TRUE),
                  value   = TRUE)
  resOut <- lapply(resFile, FUN = function(fileName) {
    resTable <- read_tsv(file = fileName)
  })
  resOut <- do.call(what = rbind, args = resOut)
  # extract leading edge genes
  rnk <- read_tsv(file      = gseaRun[["rnk"]],
                  skip      = 1,
                  col_names = c("SYMBOL", "score")) %>%
         arrange(desc(score))
  leGenes <- group_by(resOut, NAME) %>%
             do(LEADING_EDGE = ifelse(test = sign(.$NES) %in% 1,
                    yes = paste(intersect(rnk$SYMBOL[seq(from = 1,
                        to = .$"RANK AT MAX" + 1)],
                        gmx[[.$NAME]]), collapse = ","),
                    no  = paste(intersect(rnk$SYMBOL[seq(from = nrow(rnk) -
                        .$"RANK AT MAX",
                        to = nrow(rnk))],
                        gmx[[.$NAME]]), collapse = ","))) %>%
             ungroup() %>%
             mutate(LEADING_EDGE = unlist(LEADING_EDGE))
  resOut <- merge(resOut, leGenes, by = "NAME")
  # append directory name
  resOut <- mutate(resOut, directory = gseaDir)
  return(value = resOut)
})
gseaOutput <- do.call(what = rbind, args = gseaOutput)
gseaOutput <- merge(gseaOutput, gseaIndex, by = "directory")

save(gseaOutput, file = file.path(workDir, "output/ugent.pbmc.bcell.gseaOutput.RData"))
```

remove advanced directory

```r
unlink(gseaDir, recursive = TRUE)
```



Print B cell activation pathway enriched between COVID19 POS vs NEG

```r
# read msigdb
xmlFile <- file.path(workDir, "utils/msigdb_v7.1.xml")
msig <- getBroadSets(uri = xmlFile)
descDF <- sapply(msig, FUN = description) %>%
  data.frame(DESC = .) %>%
  mutate(NAME = names(msig))

# Identify genesets of interest based on keywords
# keywords: B-CELL and activation
gsNames <- descDF %>%
  filter(grepl(pattern = "B.?CELL", NAME) &
	   grepl(pattern = "activation", DESC, ignore.case = TRUE)) %>%
  filter(!grepl(pattern = "NEGATIVE_REGULATION", NAME))

# print b cell activation genesets
print(gsNames$NAME)
```

```
## [1] "BIOCARTA_ASBCELL_PATHWAY"                           
## [2] "REACTOME_ACTIVATION_OF_NF_KAPPAB_IN_B_CELLS"        
## [3] "GO_B_CELL_PROLIFERATION"                            
## [4] "GO_REGULATION_OF_B_CELL_ACTIVATION"                 
## [5] "GO_POSITIVE_REGULATION_OF_B_CELL_ACTIVATION"        
## [6] "BIOCARTA_BBCELL_PATHWAY"                            
## [7] "GO_B_CELL_PROLIFERATION_INVOLVED_IN_IMMUNE_RESPONSE"
## [8] "REACTOME_ACTIVATION_OF_RAS_IN_B_CELLS"
```

```r
# print b cell activation genesets meeting size restriction (min=15, max=2000)
print(intersect(gsNames$NAME, gseaOutput$NAME))
```

```
## [1] "GO_B_CELL_PROLIFERATION"                    
## [2] "GO_REGULATION_OF_B_CELL_ACTIVATION"         
## [3] "GO_POSITIVE_REGULATION_OF_B_CELL_ACTIVATION"
```

```r
# print number significantly associated with COVID19
filter(gseaOutput, NAME %in% gsNames$NAME & `NOM p-val` <= 0.05) %>%
  select(modelName, contrast, NAME, NES, `NOM p-val`, `FDR q-val`)
```

```
##       modelName          contrast                               NAME       NES
## 1 PBMC_BCELL_T1 COVID19.statusPOS GO_REGULATION_OF_B_CELL_ACTIVATION 0.4881553
##   NOM p-val FDR q-val
## 1     0.038 0.3707107
```
Only one out of three of the genesets is significantly associated with COVID19

Print leading edge heatmap

```r
sigDF <- filter(gseaOutput, NAME %in% gsNames$NAME & `NOM p-val` <= 0.05)

leLS <- strsplit(sigDF$LEADING_EDGE, split = ",") %>%
  unlist()

rNames <- gsub(pattern = ".+---", replacement = "", rownames(bComb))
mat <- as.matrix(bComb$RNA@counts[match(leLS, table = rNames), ])
rownames(mat) <- leLS
flag <- rowMeans(mat > 0) >= 0.2
mat <- mat[flag, ]
mat <- log2(mat + 0.25)
Heatmap(matrix = mat,
	column_split = bComb@meta.data$COVID19.status,
	show_column_names = FALSE,
	name = "log2(counts)")
```

![plot of chunk le-heat](../figure/le-heat-1.png)

Print session info

```r
sessionInfo()
```

```
## R version 4.0.0 (2020-04-24)
## Platform: x86_64-apple-darwin19.4.0 (64-bit)
## Running under: macOS Catalina 10.15.5
## 
## Matrix products: default
## BLAS/LAPACK: /usr/local/Cellar/openblas/0.3.9/lib/libopenblasp-r0.3.9.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] forcats_0.5.0               stringr_1.4.0              
##  [3] dplyr_1.0.0                 purrr_0.3.4                
##  [5] readr_1.3.1                 tidyr_1.1.0                
##  [7] tibble_3.0.1                ggplot2_3.3.1              
##  [9] tidyverse_1.3.0             GSEABase_1.50.1            
## [11] graph_1.66.0                annotate_1.66.0            
## [13] XML_3.99-0.3                AnnotationDbi_1.50.0       
## [15] ComplexHeatmap_2.4.2        MAST_1.14.0                
## [17] SingleCellExperiment_1.10.1 SummarizedExperiment_1.18.1
## [19] DelayedArray_0.14.0         matrixStats_0.56.0         
## [21] Biobase_2.48.0              GenomicRanges_1.40.0       
## [23] GenomeInfoDb_1.24.0         IRanges_2.22.2             
## [25] S4Vectors_0.26.1            BiocGenerics_0.34.0        
## [27] Seurat_3.1.5                knitr_1.28                 
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.3.1           backports_1.1.7        circlize_0.4.9        
##   [4] plyr_1.8.6             igraph_1.2.5           lazyeval_0.2.2        
##   [7] splines_4.0.0          listenv_0.8.0          digest_0.6.25         
##  [10] htmltools_0.4.0        fansi_0.4.1            magrittr_1.5          
##  [13] memoise_1.1.0          cluster_2.1.0          ROCR_1.0-11           
##  [16] globals_0.12.5         modelr_0.1.8           prettyunits_1.1.1     
##  [19] colorspace_1.4-1       blob_1.2.1             rvest_0.3.5           
##  [22] ggrepel_0.8.2          haven_2.3.1            xfun_0.14             
##  [25] crayon_1.3.4           RCurl_1.98-1.2         jsonlite_1.6.1        
##  [28] survival_3.1-12        zoo_1.8-8              ape_5.4               
##  [31] glue_1.4.1             gtable_0.3.0           zlibbioc_1.34.0       
##  [34] XVector_0.28.0         leiden_0.3.3           GetoptLong_0.1.8      
##  [37] future.apply_1.5.0     shape_1.4.4            abind_1.4-5           
##  [40] scales_1.1.1           DBI_1.1.0              Rcpp_1.0.4.6          
##  [43] progress_1.2.2         viridisLite_0.3.0      xtable_1.8-4          
##  [46] clue_0.3-57            reticulate_1.16        bit_1.1-15.2          
##  [49] rsvd_1.0.3             tsne_0.1-3             htmlwidgets_1.5.1     
##  [52] httr_1.4.1             RColorBrewer_1.1-2     ellipsis_0.3.1        
##  [55] ica_1.0-2              pkgconfig_2.0.3        uwot_0.1.8            
##  [58] dbplyr_1.4.4           utf8_1.1.4             tidyselect_1.1.0      
##  [61] rlang_0.4.6            reshape2_1.4.4         munsell_0.5.0         
##  [64] cellranger_1.1.0       tools_4.0.0            cli_2.0.2             
##  [67] generics_0.0.2         RSQLite_2.2.0          broom_0.5.6           
##  [70] ggridges_0.5.2         evaluate_0.14          bit64_0.9-7           
##  [73] fs_1.4.1               fitdistrplus_1.1-1     RANN_2.6.1            
##  [76] pbapply_1.4-2          future_1.17.0          nlme_3.1-148          
##  [79] xml2_1.3.2             compiler_4.0.0         rstudioapi_0.11       
##  [82] plotly_4.9.2.1         png_0.1-7              reprex_0.3.0          
##  [85] stringi_1.4.6          highr_0.8              lattice_0.20-41       
##  [88] Matrix_1.2-18          vctrs_0.3.1            pillar_1.4.4          
##  [91] lifecycle_0.2.0        lmtest_0.9-37          GlobalOptions_0.1.1   
##  [94] RcppAnnoy_0.0.16       data.table_1.12.8      cowplot_1.0.0         
##  [97] bitops_1.0-6           irlba_2.3.3            patchwork_1.0.0       
## [100] R6_2.4.1               KernSmooth_2.23-17     gridExtra_2.3         
## [103] codetools_0.2-16       MASS_7.3-51.6          assertthat_0.2.1      
## [106] rjson_0.2.20           withr_2.2.0            sctransform_0.2.1     
## [109] GenomeInfoDbData_1.2.3 hms_0.5.3              Rtsne_0.15            
## [112] lubridate_1.7.8
```


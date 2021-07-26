#

c <- read.delim("temp.plasma.cytokine.txt", row.names=1)
r <- read.delim("temp-path-avg-exp-2", row.names=1)
for(i in 1:dim(c)[1]) {
 for(j in 1:dim(r)[1]) {
 t <- dim(c)[2] - rowSums(c[i,] == "NaN")
 if(t > dim(c)[2] / 2) {
  res <- cor.test(as.numeric(c[i,]), as.numeric(r[j,]), method="spearman")
  str <- paste(rownames(c)[i], rownames(r)[j], t, res$estimate, res$p.value, sep = "\t")
  cat(str);
  cat("\n")
 }
 }
}


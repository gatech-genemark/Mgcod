install.packages("glmnet", dependencies=TRUE)
install.packages("Rcpp", dependencies=TRUE)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install("qvalue")
install.packages("~/Downloads/VirFinder_1.1.tar.gz", type="source", repos=NULL)

library(VirFinder)
predResult <- VF.pred("/home/admin-aaron/genetic_code/all_human_gut_contigs.fna")
predResult[order(predResult$pvalue),]
write.table(predResult, file='/home/admin-aaron/genetic_code/VirFinder_human_gut_contigs.tab', quote=FALSE, sep='\t', col.names = TRUE, row.names=TRUE)
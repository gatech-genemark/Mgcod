#!/usr/bin/env Rscript

install.packages("glmnet", dependencies=TRUE)
install.packages("Rcpp", dependencies=TRUE)
install.packages("optparse", dependencies=TRUE)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("qvalue")

library("optparse")
option_list = list(
  make_option(c("--virfinder_tar"), type="character", default="~/Downloads/VirFinder_1.1.tar.gz",
              help="path to VirFinder tar, which is used to install the package [default= %default]", metavar="character"),
  make_option(c("-i", "--input"), type="character",
              help="input file name", metavar="character"),
  make_option(c("-o", "--output"), type="character",
              help="output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
install.packages(opt$virfinder_tar, type="source", repos=NULL)

library(VirFinder)
predResult <- VF.pred(opt$input)
predResult[order(predResult$pvalue),]
write.table(predResult, file=opt$output, quote=FALSE, sep='\t', col.names = TRUE, row.names=TRUE)
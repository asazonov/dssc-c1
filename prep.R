source("https://bioconductor.org/biocLite.R")

bertie = read.csv("~/Dropbox/Work/Cambridge/UCLhack/bertie.csv")
bertie[1:5, 1:5]

library(devtools)
install_github("satijalab/seurat")
library(Seurat)

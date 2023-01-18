# install.packages("devtools")
library(devtools)
# install_version("BiocManager", version = "1.30.18")
# BiocManager::install("SingleR")

library(SingleR)
library(Seurat)
library(zellkonverter)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(arrow)

annotation <- "nowakowski"
ref_in <- "D:/Data/setd1a/4-singler/in/nowakowski_bx_v1.cpm.fth"
mat_in <- "D:/Data/setd1a/4-singler/in/setd1a_b2_bx_v1.cpm.fth"
labels_in <- "D:/Data/setd1a/4-singler/in/nowakowski.label_med.csv"
prefix <- "setd1a_b2_bx_29_12_2022"
setwd("~/Coding/R/scrnaseq")

# obtain annotations 
ref <- as.data.frame(read_feather(ref_in))
rownames(ref) <- ref[,1]
ref <- ref[,-1]

mat <- as.data.frame(read_feather(mat_in))
rownames(mat) <- mat[,1]
mat <- mat[,-1]

labels <- read.table(
  sep=',',
  header=TRUE, 
  check.names=FALSE,
  labels_in,
  row.names=1,
)

# find common celltypes
if (annotation == 'nowakowski') {
  cells <- intersect(rownames(ref), rownames(labels))
  ref <- ref[cells,]
}

# log1p regularisation
if (annotation == 'nowakowski') {
  ref <- log1p(t(ref))
}
mat <- log1p(t(mat))

if (annotation == 'nowakowski') {
  labels.full <- labels[colnames(ref),]
}

# prediction
# BiocManager::install("scran")
library(scran)
pred <- SingleR(
  test=mat,
  ref=ref,
  labels=labels.full,
  de.method='wilcox',
)

write.table(
  pred,
  file=paste(prefix,'_', annotation, '.singler.csv', sep=''),
  sep=',',
  quote=FALSE,
)

pred <- read.csv("~/Coding/R/scrnaseq/setd1a_b2_bx_19.09.2022_nowakowski.singler.csv")

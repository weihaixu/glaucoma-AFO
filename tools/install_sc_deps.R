#!/usr/bin/env Rscript
# 把 monocle3 与 CellChat 的缺失依赖从国内镜像装齐（不含 GitHub-only 的本体）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"),
        BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

deps <- unique(c(
  # monocle3 DESCRIPTION（Depends/Imports）
  "BPCells", "Biobase", "BiocGenerics", "BiocParallel", "DelayedArray",
  "DelayedMatrixStats", "HDF5Array", "MASS", "Matrix", "RANN", "RColorBrewer",
  "Rcpp", "RcppAnnoy", "RcppHNSW", "RhpcBLASctl", "Rtsne", "S4Vectors",
  "SingleCellExperiment", "SummarizedExperiment", "assertthat", "batchelor",
  "digest", "dplyr", "future", "ggdist", "ggforce", "ggplot2", "ggrastr",
  "ggrepel", "igraph", "irlba", "leidenbase", "limma", "lme4", "lmtest",
  "methods", "proxyC", "ps", "reshape2", "rlang", "rtracklayer", "scales",
  "scater", "scran", "shiny", "slingshot", "tibble", "tidyr", "uwot",
  "viridis", "zoo", "pryr",
  # CellChat Imports（对照 git clone 的 DESCRIPTION）
  "circlize", "clusterProfiler", "ComplexHeatmap", "NMF", "plyr",
  "reticulate", "matrixStats", "future.apply", "ggraph", "Rcpp"
))
deps <- setdiff(unique(deps), rownames(installed.packages()))
cat("missing:", length(deps), ":", paste(deps, collapse = ", "), "\n")

if (length(deps)) {
  # 先 CRAN
  cran <- setdiff(deps, rownames(installed.packages()))
  try(install.packages(unique(cran), type = "win.binary"), silent = TRUE)
  # 剩下的走 Bioc
  still <- setdiff(deps, rownames(installed.packages()))
  if (length(still)) {
    cat("after CRAN, still missing:", paste(still, collapse = ", "), "\n")
    try(BiocManager::install(unique(still), ask = FALSE, update = FALSE), silent = TRUE)
  }
}

left <- setdiff(deps, rownames(installed.packages()))
cat("still missing at end:", if (length(left)) paste(left, collapse = ", ") else "NONE", "\n")

#!/usr/bin/env Rscript
# 从 tuna 的 Bioc release 源码索引找 4 个包的版本号
ct <- readLines(gzfile("downloads/bioc_packages.gz"))
want <- c("batchelor", "BiocNeighbors", "DelayedMatrixStats", "BPCells")
pkg <- sub("^Package: ", "", grep("^Package: ", ct, value = TRUE))
ver <- sub("^Version: ", "", grep("^Version: ", ct, value = TRUE))
idx <- match(want, pkg)
for (i in seq_along(want)) {
  cat(want[i], if (!is.na(idx[i])) ver[idx[i]] else "NOT_FOUND", "\n")
}

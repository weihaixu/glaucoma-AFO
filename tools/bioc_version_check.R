#!/usr/bin/env Rscript
# 查 Bioc release 索引里若干包的版本 vs 本机已装版本
ct <- readLines(gzfile("downloads/bioc_packages.gz"))
pkgs <- c("assorthead", "knncolle", "BiocNeighbors", "batchelor", "DelayedMatrixStats",
          "beachmat", "MAGEE", "Sputnik")  # 后两个仅验证索引解析
for (p in pkgs) {
  i <- grep(paste0("^Package: ", p, "$"), ct)
  if (length(i)) {
    v <- sub("^Version: ", "", grep("^Version:", ct[i[1]:(i[1] + 2)], value = TRUE)[1])
    have <- tryCatch(as.character(packageVersion(p)), error = function(e) "无")
    cat(sprintf("%-20s Bioc: %-10s 已装: %s\n", p, v, have))
  } else {
    cat(sprintf("%-20s 不在 Bioc bioc 索引（已装: %s）\n", p,
                tryCatch(as.character(packageVersion(p)), error = function(e) "无")))
  }
}

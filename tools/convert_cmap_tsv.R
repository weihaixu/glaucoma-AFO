#!/usr/bin/env Rscript
# 把 cmap_connectivity.src.csv 转为账本读取器(read.delim)可解析的制表符 TSV
d <- read.csv("results/tables/cmap_connectivity.src.csv", stringsAsFactors = FALSE)
stopifnot(identical(names(d), c("drug", "cs")))
write.table(d, "results/tables/cmap_connectivity.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(d), "rows;", sum(d$cs < 0), "negative\n")

# 从 CELLxGENE h5ad 导出的 10x 风格目录组装 Seurat 参考对象
# 用法：Rscript tools/build_seurat_ref.R <dir> <out_rds>
#   <dir> 必须含 matrix.mtx / genes.tsv / barcodes.tsv / obs.csv
# cell_type / donor 等 meta 从 obs.csv 取回；没有 cell_type 的场合交给 06 的
# marker 打分兜底。

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
dir_in <- args[1]
out_rds <- args[2]

suppressMessages(library(Seurat))
suppressMessages(library(Matrix))

m <- readMM(file.path(dir_in, "matrix.mtx"))
genes <- read.delim(file.path(dir_in, "genes.tsv"), header = FALSE,
                    stringsAsFactors = FALSE)[, 1]
bc <- readLines(file.path(dir_in, "barcodes.tsv"))
rownames(m) <- make.unique(genes)
colnames(m) <- bc

# 去重基因：细胞图谱里偶尔有同名 symbol，保留表达量高的
if (anyDuplicated(genes)) {
  keep <- tapply(rowSums(m), genes, which.max)
  m <- m[as.integer(keep), ]
}

obs <- read.csv(file.path(dir_in, "obs.csv"), stringsAsFactors = FALSE,
                check.names = FALSE, row.names = 1)
obs <- obs[bc, , drop = FALSE]

sc <- CreateSeuratObject(m, project = "scRNA_ref", min.cells = 3, min.features = 200)
for (cn in intersect(c("cell_type", "donor_id", "sample", "condition",
                       "tissue", "assay", "development_stage"), colnames(obs))) {
  sc[[cn]] <- obs[[cn]]
}
if (!"sample" %in% colnames(sc@meta.data) && "donor_id" %in% colnames(sc@meta.data)) {
  sc$sample <- obs$donor_id
}

saveRDS(sc, out_rds)
cat("saved:", out_rds, "cells:", ncol(sc), "genes:", nrow(sc),
    "cell_type n =", if ("cell_type" %in% colnames(sc@meta.data)) length(unique(sc$cell_type)) else 0L, "\n")

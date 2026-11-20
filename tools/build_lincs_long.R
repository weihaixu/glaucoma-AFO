# tools/build_lincs_long.R —— 把 Trans-Phar 的 QCed L1000 宽表合并成 (drug, gene, zscore) 长表
# 每个文件一个细胞系（phase），列名 {cell}_{drug}_{dose}_{time}；按 (drug, gene) 取均值。
suppressMessages(library(data.table))

src <- "data/raw/lincs/Cmap_QCeddata"
out <- "data/raw/lincs/lincs_l1000_long.tsv.gz"
files <- list.files(src, pattern = "\\.gz$", full.names = TRUE)
files <- files[file.size(files) > 5000]
cat("files:", length(files), "\n")

agg <- list()
for (f in files) {
  d <- fread(f, showProgress = FALSE, encoding = "Latin-1")
  gene_col <- "Gene"
  if (!gene_col %in% names(d)) { cat("SKIP (no Gene col):", basename(f), "\n"); next }
  m <- melt(d, id.vars = gene_col, variable.name = "profile", value.name = "zscore")
  m <- m[is.finite(zscore)]
  # 列名: {cell}_{drug}_{dose}_{time}h -> drug = 第 2 段；源文件是 latin-1，转 UTF-8
  m[, profile := iconv(as.character(profile), from = "Latin1", to = "UTF-8")]
  m[, drug := sub("^[^_]+_([^_]+)_.*$", "\\1", profile)]
  a <- m[, .(zscore = mean(zscore)), by = .(drug, gene = get(gene_col))]
  agg[[basename(f)]] <- a
  cat(sprintf("%-42s %5d profiles  -> %5d (drug,gene) pairs\n", basename(f),
              length(unique(m$profile)), nrow(a)))
}
all <- rbindlist(agg)
final <- all[, .(zscore = mean(zscore)), by = .(drug, gene)]
setorder(final, drug, gene)
fwrite(final, out)
cat("long table:", nrow(final), "rows,", uniqueN(final$drug), "drugs,",
    uniqueN(final$gene), "genes ->", out, "\n")
cat("total size:", file.size(out) / 1e6, "MB\n")

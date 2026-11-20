# tools/cmap_fill_file.R -- v4 断点续跑：只处理一个面板文件，追加到桶文件
# 用法: Rscript tools/cmap_fill_file.R <k>   (k = panel 中 1-based 序号)
# 字典构建逻辑与 cmap_pipeline_v4.R 完全一致（同序 union -> int 编码一致）。
suppressMessages(library(data.table))
setDTthreads(4)

k        <- as.integer(commandArgs(trailingOnly = TRUE)[1])
src      <- "data/raw/lincs/Cmap_QCeddata"
work     <- "data/raw/lincs/cmap_work"
NBUCKET  <- 12L
CHUNK    <- 200L

panel <- c(
  "lung_A549_ph1_QCed.txt.gz", "lung_A549_ph2_QCed.txt.gz",
  "skin_A375_ph1_QCed.txt.gz", "skin_A375_ph2_QCed.txt.gz",
  "liver_HEPG2_ph1_QCed.txt.gz", "liver_HEPG2_ph2_QCed.txt.gz",
  "lung_HCC515_ph1_QCed.txt.gz", "lung_HCC515_ph2_QCed.txt.gz",
  "digestive_HT29_ph1_QCed.txt.gz", "digestive_HT29_ph2_QCed.txt.gz",
  "breast_MCF7_ph1_QCed.txt.gz", "breast_MCF7_ph2_QCed.txt.gz",
  "prostate_PC3_ph1_QCed.txt.gz", "prostate_PC3_ph2_QCed.txt.gz",
  "prostate_VCAP_ph1_QCed.txt.gz", "skin_SKMEL28_ph1_QCed.txt.gz"
)
bfiles <- file.path(work, sprintf("part_%02d.tsv", seq_len(NBUCKET) - 1L))
stopifnot(k >= 1L, k <= length(panel))
f <- file.path(src, panel[k])

# ---- 字典（与 v4 完全同序的全局构建，保证 int 编码跨文件一致）----
drug_levels <- character(0)
file_info <- list()
for (g in panel) {
  h <- fread(file.path(src, g), nrows = 0, encoding = "Latin-1", showProgress = FALSE)
  cs <- names(h)[names(h) != "Gene"]
  pf <- iconv(cs, from = "Latin1", to = "UTF-8")
  d2 <- sub("^[^_]+_([^_]+)_.*$", "\\1", pf)
  file_info[[basename(g)]] <- list(cols = cs, dr = d2)
  drug_levels <- union(drug_levels, d2)
}
G <- character(0)
for (g in panel) {
  d3 <- fread(file.path(src, g), select = "Gene", encoding = "Latin-1", showProgress = FALSE)
  G <- union(G, d3$Gene)
  rm(d3); invisible(gc())
}
drug_id <- setNames(seq_along(drug_levels), drug_levels)
cat("global drugs:", length(drug_levels), " genes:", length(G), "\n")

# ---- 单文件处理：chunk 聚合 -> 12 桶累积 -> 追加落盘 ----
t0 <- Sys.time()
info <- file_info[[basename(f)]]
d <- fread(f, encoding = "Latin-1", showProgress = FALSE)
gene_idx <- match(d$Gene, G)
ok_g <- !is.na(gene_idx)
gene_idx <- as.integer(gene_idx[ok_g])
ng <- length(gene_idx)
cols_all <- intersect(info$cols, names(d)[names(d) != "Gene"])
dr_idx <- as.integer(drug_id[info$dr[match(cols_all, info$cols)]])
tot <- vector("list", NBUCKET)
for (st in seq(1L, length(cols_all), by = CHUNK)) {
  en <- min(st + CHUNK - 1L, length(cols_all))
  idx <- st:en
  sub <- d[, cols_all[idx], with = FALSE]
  zm <- do.call(cbind, lapply(sub, function(x) as.numeric(as.character(x))))
  dim(zm) <- c(ng, length(idx))
  z <- as.vector(zm)
  ok <- is.finite(z)
  if (!any(ok)) { rm(sub, zm); next }
  drv <- rep(dr_idx[idx], each = ng)[ok]
  giv <- rep(gene_idx, length(idx))[ok]
  zv <- z[ok]
  a <- data.table(dr = drv, gi = giv, z = zv)
  a <- a[, .(zsum = sum(z), n = .N), by = .(dr, gi)]
  bsel <- ((a$dr - 1L) %% NBUCKET) + 1L
  for (bb in unique(bsel)) {
    s <- a[bsel == bb]
    tot[[bb]] <- rbindlist(list(tot[[bb]], s))
    tot[[bb]] <- tot[[bb]][, .(zsum = sum(zsum), n = sum(n)), by = .(dr, gi)]
  }
  rm(sub, zm, a, s); invisible(gc(full = FALSE))
}
rm(d); invisible(gc())
for (bb in seq_len(NBUCKET)) {
  t <- tot[[bb]]
  if (is.null(t) || !nrow(t)) next
  t <- t[, .(zsum = sum(zsum), n = sum(n)), by = .(dr, gi)]
  setorder(t, dr, gi)
  fwrite(t, bfiles[bb], sep = "\t", quote = FALSE,
         append = file.exists(bfiles[bb]) && file.size(bfiles[bb]) > 0)
}
cat(sprintf("done %2d/16 %s elapsed=%s\n", k, panel[k], format(Sys.time() - t0)))

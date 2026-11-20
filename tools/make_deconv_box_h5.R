# tools/make_deconv_box_h5.R -- 从 vanZyl h5ad 分块流式读，构建签名 + NNLS + 箱线图
# 背景：scRNA_ref.rds 整体读入需要 ~2.7GB 峰值，本机提交内存不够；
# h5ad 是 CSR 稀疏（X/data 整数 raw counts），按细胞分块读、只累积各类型
# 基因总 counts 与文库总大小，内存增量 <200MB。
# 类型口径与 scRNA_ref.rds 的 cell_type 一致（author_cell_type，TM_Fibro/Schlemm_Endo
# 分别对应 TM/SC）；签名平均表达口径 log1p(mean_raw / mean_lib * 1e4)，
# 与 NormalizeData 后 rowMeans 的做法在签名选择用途下等价。
source("R/helpers.R")
cfg <- read_config()
log_init("deconv_h5", cfg)
load_pkgs(c("ggplot2", "reshape2", "Matrix", "hdf5r"))
set.seed(cfg_get(cfg, "project.seed", 1))

gse   <- if (isTRUE(cfg_get(cfg, "run.demo", FALSE))) "DEMO" else cfg_get(cfg, "bulk.train.gse", "GSE27276")
expr  <- readRDS(proj_path("data", "processed", paste0(gse, "_expr.rds")))
pheno <- readRDS(proj_path("data", "processed", paste0(gse, "_pheno.rds")))

h5_path <- "data/raw/scRNA_cellxgene/vanZyl_TM_antseg.h5ad"
f   <- h5file(h5_path, mode = "r")
ptr <- f[["X/indptr"]][]
# var/_index 是 Ensembl ID；基因符号在 feature_name（categorical 存储，对齐 bulk 行名）
fn_cats <- as.character(f[["var/feature_name/categories"]][])
fn_codes <- as.integer(f[["var/feature_name/codes"]][] + 1L)
genes_h5 <- fn_cats[fn_codes]
# NA / 重复符号行在聚合后剔除（保留首个）
g_ok <- !is.na(genes_h5) & nzchar(genes_h5) & !duplicated(genes_h5)
cats <- f[["obs/author_cell_type/categories"]][]
codes <- as.integer(f[["obs/author_cell_type/codes"]][] + 1L)  # 0-based -> 1-based
ct_all <- cats[codes]
n_cells <- length(ptr) - 1L
n_genes <- length(genes_h5)
cat("cells:", n_cells, " genes:", n_genes, " nnz:", ptr[n_cells + 1], "\n")

types <- names(sort(table(ct_all), decreasing = TRUE))
types <- types[table(ct_all)[types] >= 150]
cat("signature types (>=150 cells):", length(types), "->", paste(types, collapse = ", "), "\n")

# ---- 分块累积：sum_raw[genes x types], sum_lib[type] ----
sum_raw <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 0, dims = c(n_genes, length(types)))
sum_lib_type <- setNames(numeric(length(types)), types)
ds_data <- f[["X/data"]]; ds_idx <- f[["X/indices"]]
CHUNK <- 800L
n_chunk <- ceiling(n_cells / CHUNK)
t0 <- Sys.time()
for (i1 in seq_len(n_chunk)) {
  lo <- (i1 - 1L) * CHUNK + 1L
  hi <- min(lo + CHUNK - 1L, n_cells)
  p1 <- ptr[lo] + 1L; p2 <- ptr[hi + 1L]
  dv <- as.numeric(ds_data[p1:p2])
  iv <- as.integer(ds_idx[p1:p2]) + 1L
  cnt_per_cell <- diff(ptr[lo:(hi + 1L)])
  jv <- rep.int(seq_len(hi - lo + 1L), cnt_per_cell)
  M <- Matrix::sparseMatrix(i = iv, j = jv, x = dv,
                            dims = c(n_genes, hi - lo + 1L))
  rm(dv, iv, jv); invisible(gc(full = FALSE))
  lib <- Matrix::colSums(M)
  ct_chunk <- ct_all[lo:hi]
  valid <- which(ct_chunk %in% types)          # 排除 NA 与 <150 细胞的小类型
  Ik <- Matrix::sparseMatrix(i = valid, j = match(ct_chunk[valid], types), x = 1,
                             dims = c(hi - lo + 1L, length(types)))
  sum_raw <- sum_raw + M %*% Ik   # genes x cells %*% cells x types -> genes x types
  rm(M, Ik); invisible(gc(full = FALSE))
  lib_by_type <- tapply(lib, ct_chunk, sum)
  for (k in names(lib_by_type)) if (k %in% types) sum_lib_type[k] <- sum_lib_type[k] + lib_by_type[[k]]
  cat(sprintf("  chunk %d/%d done (%s elapsed)\n", i1, n_chunk, format(Sys.time() - t0)))
}
f$close_all()
n_k <- as.numeric(table(ct_all)[types])
sum_raw <- as.matrix(sum_raw)
mean_raw <- sweep(sum_raw, 2, n_k, "/")
mean_lib <- sum_lib_type / n_k
means <- log1p(sweep(mean_raw, 2, mean_lib, "/") * 1e4)
colnames(means) <- types
rm(sum_raw); invisible(gc(full = TRUE))

# ---- 签名：bulk 交集基因（仅有效、去重后的符号），每类型特异性 top 300 ----
keep_idx <- which(g_ok & genes_h5 %in% rownames(expr))
cat("bulk intersect genes:", length(keep_idx), "/", n_genes, "\n")
means <- means[keep_idx, , drop = FALSE]
rownames(means) <- genes_h5[keep_idx]
others <- apply(means, 1, function(r) sort(r, decreasing = TRUE)[2])
sig <- list()
for (k in types) {
  spec <- means[, k] - others
  cand <- names(sort(spec[means[, k] > 0.5], decreasing = TRUE))
  sig[[k]] <- head(cand, cfg_get(cfg, "bulk.deconv.max_genes_per_celltype", 300))
}
g_sel <- unique(unlist(sig))
S <- means[g_sel, , drop = FALSE]
S[S < 0] <- 0
S <- S[rowSums(S) > 0, , drop = FALSE]
cat("signature matrix:", nrow(S), "genes x", ncol(S), "types\n")
if (nrow(S) < 50) stop("signature 基因过少（", nrow(S), "），疑似基因符号未对齐，中止")

# ---- NNLS ----
common <- rownames(S)
X <- expr[common, , drop = FALSE]
frac <- matrix(NA_real_, ncol(X), ncol(S), dimnames = list(colnames(X), colnames(S)))
for (j in seq_len(ncol(X))) {
  fit <- nnls::nnls(S, X[, j])$x
  if (sum(fit) > 0) frac[j, ] <- fit / sum(fit)
}
cat("NNLS done:", nrow(frac), "samples x", ncol(frac), "types\n")
frac <- as.data.frame(frac)   # 必须先转 data.frame，否则 $<- 会把矩阵摊平成逐元素列表

frac$sample <- rownames(frac)
frac <- merge(pheno[, c("sample", "group")], frac, by = "sample")
write_tsv(frac, "deconv_fraction", cfg)

cell_cols <- setdiff(names(frac), c("sample", "group"))
ddiff <- do.call(rbind, lapply(cell_cols, function(cc) {
  x <- frac[[cc]][frac$group == "case"]; y <- frac[[cc]][frac$group == "control"]
  if (length(x) < 3 || length(y) < 3) return(NULL)
  wt <- stats::wilcox.test(x, y)
  data.frame(cell_type = cc, mean_case = mean(x), mean_control = mean(y),
             diff = mean(x) - mean(y), p = wt$p.value)
}))
ddiff$fdr <- stats::p.adjust(ddiff$p, "BH")
ddiff <- ddiff[order(ddiff$p), ]
write_tsv(ddiff, "deconv_group_diff", cfg)
print(ddiff)

long <- reshape2::melt(frac[, c("sample", "group", cell_cols)], id.vars = c("sample", "group"))
g_box <- ggplot2::ggplot(long, ggplot2::aes(group, value, fill = group)) +
  ggplot2::geom_boxplot(outlier.size = 0.4, width = 0.6) +
  ggplot2::facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  ggplot2::scale_fill_manual(values = c(control = "#1F6FB4", case = "#C0392B"), guide = "none") +
  ggplot2::labs(x = NULL, y = "estimated fraction") +
  theme_pub()
save_fig(g_box, "deconv_box", cfg, 9, 6)
cat("written: deconv_fraction.tsv, deconv_group_diff.tsv, deconv_box\n")

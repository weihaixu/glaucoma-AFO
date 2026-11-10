# =============================================================================
# tools/scenic_regulon_count.R —— 从 pySCENIC 的 AUC 矩阵数差异 regulon
#
# 为什么单独一个脚本：06_scrna_cellchat.R 做的是细胞注释/CellChat/拟时序，
# 不含 SCENIC。SCENIC 通常在 python（pySCENIC）里跑完，产物是
#   auc_mtx.csv     行=细胞，列=regulon（如 "YAP1(+)""FOSB(+)"）
#   cell_info.csv   行=细胞，含分组列（disease/control）与细胞类型列
# 本脚本只做统计，不重跑网络推断，因此不依赖 SCENIC/pySCENIC 的任何包。
#
# 用法：
#   Rscript tools/scenic_regulon_count.R \
#     --auc  data/raw/scenic/auc_mtx.csv \
#     --meta data/raw/scenic/cell_info.csv \
#     --group-col group --celltype-col cell_type \
#     --case POAG --control Control \
#     --out  results/tables/scenic_regulon_diff.tsv
#
# 统计口径（正文必须照此描述）：
#   每个 regulon 的 AUC 值在 disease 与 control 细胞间做 Wilcoxon 秩和检验，
#   BH 校正后 padj < 0.05 记为"差异活动 regulon"。
#   若给了细胞类型列，另外给出每个细胞类型内部的检验结果（长表），
#   但计数以全局检验为准，避免"按细胞类型重复计数"被审稿人抓住。
# =============================================================================

suppressWarnings(suppressMessages({
  ok <- requireNamespace("stats", quietly = TRUE)
}))

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

auc_f   <- getopt("auc")
meta_f  <- getopt("meta")
out_f   <- getopt("out", "scenic_regulon_diff.tsv")
gcol    <- getopt("group-col", "group")
ccol    <- getopt("celltype-col", NULL)
case_lv <- getopt("case")
ctrl_lv <- getopt("control")
padj_c  <- as.numeric(getopt("padj", "0.05"))

if (is.null(auc_f) || is.null(meta_f)) {
  stop("必须给 --auc 与 --meta。加 --help 看用法见脚本头部注释。", call. = FALSE)
}
for (f in c(auc_f, meta_f)) if (!file.exists(f)) stop("文件不存在：", f, call. = FALSE)

rd <- function(f) utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE,
                                  fileEncoding = "UTF-8-BOM")

auc <- rd(auc_f)
meta <- rd(meta_f)

# --- 1) 对齐细胞 -----------------------------------------------------------
# pySCENIC 的 auc_mtx.csv 第一列是细胞名（常无列名或有 "Cell" 之类）
cell_auc <- as.character(auc[[1]])
mat <- as.matrix(auc[, -1, drop = FALSE])
rownames(mat) <- cell_auc
storage.mode(mat) <- "double"

# meta 的细胞名列：优先精确匹配，其次找与 auc 细胞名有交集的列
meta_cell <- NULL
for (cn in names(meta)) {
  if (sum(as.character(meta[[cn]]) %in% cell_auc) > 0.5 * nrow(meta)) { meta_cell <- cn; break }
}
if (is.null(meta_cell)) stop("在 --meta 里找不到与 AUC 矩阵细胞名对应的列", call. = FALSE)
if (!gcol %in% names(meta)) stop("--meta 没有列 ", gcol, "；现有列：",
                                 paste(names(meta), collapse = ", "), call. = FALSE)

common <- intersect(rownames(mat), as.character(meta[[meta_cell]]))
if (length(common) < 20) {
  stop("AUC 矩阵与注释表只交集到 ", length(common), " 个细胞，检查细胞名格式是否一致",
       call. = FALSE)
}
mat <- mat[common, , drop = FALSE]
mi  <- match(common, as.character(meta[[meta_cell]]))
grp <- as.character(meta[[gcol]][mi])

# --- 2) 分组 ---------------------------------------------------------------
lv <- unique(grp[!is.na(grp)])
if (is.null(case_lv)) {
  # 没指定就按字母序取两组，并在日志里说明，避免"悄悄选了一组"的隐患
  if (length(lv) != 2) stop("分组列有 ", length(lv), " 个水平，请用 --case/--control 指定",
                            call. = FALSE)
  case_lv <- sort(lv)[2]; ctrl_lv <- sort(lv)[1]
  message("[注意] 未指定 --case/--control，按字母序取 case=", case_lv, " control=", ctrl_lv)
}
if (is.null(ctrl_lv)) stop("给了 --case 就必须给 --control", call. = FALSE)
if (!all(c(case_lv, ctrl_lv) %in% lv)) {
  stop("--case/--control 不在分组列的水平里：", paste(lv, collapse = ", "), call. = FALSE)
}

is_case <- grp == case_lv
is_ctrl <- grp == ctrl_lv
keep <- is_case | is_ctrl
mat <- mat[keep, , drop = FALSE]
is_case <- is_case[keep]; is_ctrl <- is_ctrl[keep]
if (sum(is_case) < 5 || sum(is_ctrl) < 5) {
  stop("每组至少需 5 个细胞（现有 case=", sum(is_case), ", control=", sum(is_ctrl), "）",
       call. = FALSE)
}

# --- 3) 逐 regulon 检验 -----------------------------------------------------
test_one <- function(v) {
  a <- v[is_case]; b <- v[is_ctrl]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 3 || length(b) < 3) return(c(NA, NA, NA))
  if (stats::sd(c(a, b)) == 0) return(c(NA, NA, NA))
  p <- tryCatch(stats::wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA_real_)
  c(stats::median(a), stats::median(b), p)
}

res <- t(vapply(seq_len(ncol(mat)), function(j) test_one(mat[, j]), numeric(3)))
tab <- data.frame(
  regulon = colnames(mat),
  median_case = res[, 1],
  median_control = res[, 2],
  diff = res[, 1] - res[, 2],
  p = res[, 3],
  stringsAsFactors = FALSE
)
tab$padj <- stats::p.adjust(tab$p, method = "BH")
tab <- tab[order(tab$p), ]

n_sig <- sum(tab$padj < padj_c, na.rm = TRUE)
attr(tab, "n_sig") <- n_sig
attr(tab, "n_tested") <- sum(is.finite(tab$p))

# --- 4) 可选：细胞类型内部 --------------------------------------------------
sub <- NULL
if (!is.null(ccol) && ccol %in% names(meta)) {
  ct <- as.character(meta[[ccol]][mi])[keep]
  rows <- list()
  for (k in unique(ct[!is.na(ct)])) {
    in_k <- ct == k
    if (sum(in_k & is_case) < 5 || sum(in_k & is_ctrl) < 5) next
    r <- t(vapply(seq_len(ncol(mat)), function(j) {
      a <- mat[in_k & is_case, j]; b <- mat[in_k & is_ctrl, j]
      c(stats::median(a, na.rm = TRUE), stats::median(b, na.rm = TRUE),
        tryCatch(stats::wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA_real_))
    }, numeric(3)))
    rows[[k]] <- data.frame(cell_type = k, regulon = colnames(mat),
                            median_case = r[, 1], median_control = r[, 2],
                            diff = r[, 1] - r[, 2], p = r[, 3], stringsAsFactors = FALSE)
  }
  if (length(rows)) {
    sub <- do.call(rbind, rows)
    sub$padj <- stats::p.adjust(sub$p, method = "BH")
    sub <- sub[order(sub$p), ]
  }
}

# --- 5) 写出 ---------------------------------------------------------------
utils::write.table(tab, out_f, sep = "\t", quote = FALSE, row.names = FALSE)
message("regulon 总数 ", nrow(tab), "，可检验 ", attr(tab, "n_tested"),
        "，padj<", padj_c, " 的 ", n_sig, " 个")
if (!is.null(sub)) {
  sub_f <- sub("\\.tsv$", "_by_celltype.tsv", out_f)
  utils::write.table(sub, sub_f, sep = "\t", quote = FALSE, row.names = FALSE)
  message("按细胞类型的明细写入 ", sub_f)
}
message("完成：", out_f)

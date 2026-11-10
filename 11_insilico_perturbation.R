# =============================================================================
# 11_insilico_perturbation.R —— 课题三（干实验版）：计算扰动
#
# 目的：在没有湿实验的条件下，用单细胞转录组做"数字敲除"，回答
#       "如果把 LOX/YAP1 敲掉，TM 细胞的疾病相关基因程序会不会被推回去"。
#
# 方法（简化网络传播，不是因果推断，措辞必须写清）：
#   1. 取目标细胞类型（TM/JCT）的 log 归一化表达，选 HVG
#   2. 基因-基因相关网络，每个基因保留 top-k 正相关邻居，行归一化
#   3. 把目标基因的初值设为 -1（敲低），沿网络传播 5 步（阻尼 alpha）
#   4. 传播得到的 delta 与 bulk 疾病特征（deg_consensus 的 logFC）比较：
#        concordance = 方向一致基因比例
#        score       = mean(z(delta) * sign(logFC))   >0 表示扰动把程序推向疾病态
#   5. 置换检验：随机取度数相近的基因重复，得到经验 P
#
# 输入：config$scrna.rds（含 cell_type）、results/tables/deg_consensus.tsv
#       config$drylab$perturbation
# 输出：results/tables/perturb_summary.tsv
#       results/tables/perturb_<gene>_<celltype>.tsv
#       results/figures/perturb_<gene>_<celltype>.pdf
#
# 定位：产出的是"可被实验检验的预测"，不能写成"证明了因果关系"。
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("11_insilico_perturbation", cfg)
load_pkgs(c("ggplot2", "Matrix", "data.table"))
set.seed(cfg_get(cfg, "project.seed", 1))

pcfg <- cfg_get(cfg, "drylab.perturbation", list())
obj_path <- proj_path(cfg_get(cfg, "scrna.rds", "data/processed/scRNA_ref.rds"))
require_file(obj_path, "先准备单细胞对象（见操作手册 3.4），或在 config 改 scrna.rds")
if (!requireNamespace("Seurat", quietly = TRUE)) stop("需要 Seurat 包", call. = FALSE)

obj <- readRDS(obj_path)
if (!"cell_type" %in% names(obj@meta.data)) {
  stop("单细胞对象缺 cell_type 列：先用 06 注释，或把原文注释写进 meta.data", call. = FALSE)
}

get_mat <- function(o, what = "data") {
  tryCatch(Seurat::GetAssayData(o, assay = "RNA", layer = what),
           error = function(e) Seurat::GetAssayData(o, assay = "RNA", slot = what))
}

# ------------------------------------------------- 1) 疾病特征 --------------
sig_file <- res_path(cfg, "tables", "deg_consensus.tsv")
require_file(sig_file, "先跑 02_deg_limma.R 得到 deg_consensus.tsv")
sig <- utils::read.delim(sig_file, stringsAsFactors = FALSE)
lf_col <- grep("^logFC", names(sig), value = TRUE)[1]
if (is.na(lf_col)) stop("共识基因表里找不到 logFC 列", call. = FALSE)
sig_map <- stats::setNames(as.numeric(sig[[lf_col]]), sig$gene)
sig_map <- sig_map[is.finite(sig_map) & sig_map != 0]

cts <- intersect(pcfg$celltypes %||% c("TM", "JCT"), unique(as.character(obj$cell_type)))
if (!length(cts)) stop("cell_type 中没有目标细胞类型；现有：",
                       paste(unique(obj$cell_type), collapse = ", "), call. = FALSE)
targets <- pcfg$target_genes %||% c("LOX", "LOXL1", "YAP1")
log_msg("目标细胞类型：", paste(cts, collapse = "/"), "；目标基因：", paste(targets, collapse = "/"))

# ------------------------------------------------- 2) 工具函数 --------------

build_net <- function(mat, topk = 20) {
  cm <- stats::cor(t(mat))
  cm[!is.finite(cm)] <- 0
  diag(cm) <- 0
  cm[cm < 0] <- 0                       # 只保留正相关，避免"抑制关系"被当成激活传播
  n <- nrow(cm)
  W <- matrix(0, n, n, dimnames = dimnames(cm))
  for (i in seq_len(n)) {
    v <- cm[i, ]
    k <- min(topk, sum(v > 0))
    if (k < 1) next
    idx <- order(v, decreasing = TRUE)[seq_len(k)]
    W[i, idx] <- v[idx]
  }
  W <- pmax(W, t(W))                    # 对称化
  rs <- rowSums(W); rs[rs == 0] <- 1
  Matrix::Matrix(sweep(W, 1, rs, "/"), sparse = TRUE)
}

propagate <- function(W, gene, alpha = 0.7, steps = 5) {
  if (!gene %in% rownames(W)) return(NULL)
  d <- rep(0, nrow(W)); names(d) <- rownames(W)
  d[gene] <- -1                          # 敲低：该基因活性归零，等价于 -1 的扰动
  acc <- d
  cur <- d
  for (s in seq_len(steps)) {
    cur <- as.numeric(alpha * (W %*% cur))
    acc <- acc + cur
  }
  names(acc) <- rownames(W)
  acc
}

score_delta <- function(delta, ref) {
  g <- intersect(names(delta), names(ref))
  if (length(g) < 20) return(NULL)
  x <- delta[g]; y <- ref[g]
  list(n = length(g),
       concordance = mean(sign(x) == sign(y)),
       score = mean(as.numeric(scale(x)) * sign(y), na.rm = TRUE))
}

# ------------------------------------------------- 3) 逐细胞类型计算 --------
summ <- list()
for (ct in cts) {
  sub <- tryCatch(subset(obj, cell_type == ct), error = function(e) NULL)
  if (is.null(sub) || ncol(sub) < 30) { log_msg(ct, " 细胞数 <30，跳过"); next }
  mat <- get_mat(sub, "data")
  # 建网只需要基因间相关，细胞数上限控制内存（2000 基因 × 2 万细胞约 320 MB）
  max_cells <- pcfg$max_cells %||% 20000
  if (ncol(mat) > max_cells) {
    set.seed(cfg_get(cfg, "project.seed", 1))
    mat <- mat[, sample(ncol(mat), max_cells)]
    log_msg(ct, "：细胞数超过 ", max_cells, "，随机抽样建网")
  }
  hvg_n <- min(pcfg$hvg %||% 2000, nrow(mat))
  vf <- if (length(Seurat::VariableFeatures(sub))) Seurat::VariableFeatures(sub) else character(0)
  vf <- intersect(vf, rownames(mat))
  if (length(vf) < hvg_n) {
    mu <- Matrix::rowMeans(mat)
    vr <- Matrix::rowMeans(mat^2) - mu^2          # 稀疏矩阵上直接算方差，避免整块转稠密
    vf <- rownames(mat)[order(vr, decreasing = TRUE)]
  }
  vf <- vf[seq_len(min(hvg_n, length(vf)))]
  mat <- as.matrix(mat[vf, , drop = FALSE])
  # 只保留在疾病特征里有方向的基因做评分，其余基因仍参与建网（提供传播路径）
  log_msg(ct, "：", ncol(mat), " 细胞 × ", nrow(mat), " 基因建网")

  W <- build_net(mat, topk = pcfg$topk %||% 20)
  deg <- as.numeric(Matrix::rowSums(W != 0))
  names(deg) <- rownames(W)

  for (gn in targets) {
    if (!gn %in% rownames(W)) { log_msg(gn, " 不在 ", ct, " 的网络中，跳过"); next }
    delta <- propagate(W, gn, alpha = pcfg$alpha %||% 0.7, steps = pcfg$steps %||% 5)
    sc <- score_delta(delta, sig_map)
    if (is.null(sc)) { log_msg(gn, "/", ct, "：重叠基因不足，跳过"); next }

    # 度数匹配的置换检验
    n_perm <- pcfg$n_perm %||% 1000
    cand <- names(deg)[abs(deg - deg[gn]) <= max(5, 0.25 * deg[gn])]
    cand <- setdiff(cand, gn)
    perm_score <- rep(NA_real_, n_perm)
    for (b in seq_len(n_perm)) {
      g2 <- sample(cand, 1)
      d2 <- propagate(W, g2, alpha = pcfg$alpha %||% 0.7, steps = pcfg$steps %||% 5)
      s2 <- score_delta(d2, sig_map)
      perm_score[b] <- if (is.null(s2)) NA_real_ else s2$score
    }
    perm_score <- perm_score[is.finite(perm_score)]
    p_two <- if (length(perm_score)) {
      (1 + sum(abs(perm_score) >= abs(sc$score))) / (length(perm_score) + 1)
    } else NA_real_

    # 下游基因表：按 |delta| 排序，附疾病方向与是否一致
    dt <- data.frame(gene = names(delta), delta = as.numeric(delta),
                     disease_logFC = sig_map[names(delta)],
                     stringsAsFactors = FALSE)
    dt <- dt[order(-abs(dt$delta)), ]
    dt$direction_consistent <- sign(dt$delta) == sign(dt$disease_logFC)
    write_tsv(head(dt, 500), paste0("perturb_", gn, "_", ct), cfg)

    summ[[paste(gn, ct)]] <- data.frame(
      gene = gn, celltype = ct, n_cells = ncol(mat), n_genes_net = nrow(W),
      degree = deg[gn], n_overlap = sc$n, concordance = round(sc$concordance, 3),
      score = round(sc$score, 4), p_perm = p_two,
      n_perm_used = length(perm_score),
      top_downstream = paste(head(dt$gene[abs(dt$delta) > 0], 5), collapse = ";"),
      stringsAsFactors = FALSE)

    # 图：扰动效应 vs 疾病方向
    g <- ggplot2::ggplot(dt[is.finite(dt$disease_logFC), ],
                         ggplot2::aes(disease_logFC, delta)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
      ggplot2::geom_point(ggplot2::aes(colour = direction_consistent), size = 1.2, alpha = 0.7) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "#185FA5", `FALSE` = "#B4B2A9"),
                                   name = "direction consistent") +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                           colour = "#993C1D", linewidth = 0.6) +
      ggplot2::labs(x = "disease logFC (bulk consensus)", y = paste0("perturbation delta (", gn, " KD)"),
                    title = paste0(gn, " in silico knockdown - ", ct),
                    subtitle = paste0("concordance = ", round(sc$concordance, 2),
                                      "，score = ", round(sc$score, 3),
                                      "，P_perm = ", signif(p_two, 3))) +
      theme_pub()
    save_fig(g, paste0("perturb_", gn, "_", ct), cfg, 5.2, 4.6)
    log_msg(gn, "/", ct, "：concordance ", round(sc$concordance, 3),
            "，score ", round(sc$score, 3), "，P=", signif(p_two, 3))
  }
}

if (!length(summ)) stop("没有任何基因完成扰动分析，检查 cell_type 与基因名是否匹配", call. = FALSE)
out <- do.call(rbind, summ)
out <- out[order(out$p_perm), ]
write_tsv(out, "perturb_summary", cfg)

log_msg("判读提示：score > 0 表示敲低该基因会把转录组推向疾病方向（该基因是保护性的）；")
log_msg("          score < 0 表示敲低会逆转疾病程序（该基因是驱动性的，可作为干预靶点）。")
log_msg("11 完成")

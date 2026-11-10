# =============================================================================
# 13_robustness_cmap.R —— 课题三（干实验版）：稳健性、证伪与药物连通性
#
# 没有湿实验时，审稿人最集中的两个质疑是：
#   "这个 AUC/DEG 是不是调参调出来的？"
#   "这套特征有没有超出随机水平？"
# 本脚本用三件事回答：
#   A. 标签置换零分布：把分组打乱 1000 次，看真实 DEG 数与 AUC 落在什么位置
#   B. 多管线（multiverse）：换标准化、换基因集、换分类器，报 AUC 的范围
#      —— 单一数字不可信，范围才可信
#   C. 药物连通性（可选）：用 LINCS/CMap 类差异表达签名，算"哪些化合物能逆转
#      疾病特征"，输出候选药物清单（属于预测，不是疗效证据）
#
# 输入：data/processed/<GSE>_expr.rds / _pheno.rds、results/tables/deg_consensus.tsv
#       config$drylab$robustness（含可选的 lincs_path）
# 输出：results/tables/robustness_permutation.tsv / robustness_multiverse.tsv
#       results/tables/cmap_connectivity.tsv（有 LINCS 时）
#       results/figures/robustness_null_auc.pdf
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("13_robustness_cmap", cfg)
load_pkgs(c("ggplot2", "data.table"))
set.seed(cfg_get(cfg, "project.seed", 1))

rcfg <- cfg_get(cfg, "drylab.robustness", list())
gse <- if (isTRUE(cfg_get(cfg, "run.demo", FALSE))) "DEMO" else cfg_get(cfg, "bulk.train.gse", "GSE27276")

expr <- readRDS(proj_path("data", "processed", paste0(gse, "_expr.rds")))
pheno <- readRDS(proj_path("data", "processed", paste0(gse, "_pheno.rds")))
if (!"group" %in% names(pheno)) stop("pheno 缺 group 列", call. = FALSE)
y <- pheno$group == "case"
if (length(unique(y)) < 2) stop("分组只有一类，无法做置换检验", call. = FALSE)
x <- as.matrix(expr)
x <- x[stats::complete.cases(x), , drop = FALSE]
log_msg("训练集：", ncol(x), " 样本（", sum(y), " 例 / ", sum(!y), " 对照），", nrow(x), " 基因")

# ------------------------------------------------- 工具 ---------------------
row_tt <- function(m, g) {
  n1 <- sum(g); n2 <- sum(!g)
  m1 <- rowMeans(m[, g, drop = FALSE]); m2 <- rowMeans(m[, !g, drop = FALSE])
  v1 <- rowSums((m[, g, drop = FALSE] - m1)^2) / max(1, n1 - 1)
  v2 <- rowSums((m[, !g, drop = FALSE] - m2)^2) / max(1, n2 - 1)
  se <- sqrt(v1 / n1 + v2 / n2)
  t <- (m1 - m2) / se
  df <- (v1 / n1 + v2 / n2)^2 /
        ((v1 / n1)^2 / max(1, n1 - 1) + (v2 / n2)^2 / max(1, n2 - 1))
  p <- 2 * stats::pt(-abs(t), df = df)
  list(t = t, p = p, diff = m1 - m2)
}

auc_w <- function(score, g) {
  if (length(unique(g)) < 2) return(NA_real_)
  w <- stats::wilcox.test(score[g], score[!g], exact = FALSE)
  as.numeric(w$statistic) / (sum(g) * sum(!g))
}

score_cells <- function(xz, t, n = 10) {
  o <- order(-abs(t))[seq_len(min(n, length(t)))]
  as.numeric(colSums(xz[o, , drop = FALSE] * sign(t[o])))
}

deg_metrics <- function(m, g, xz = NULL, fdr = 0.05, lfc = 1, top_n = 10) {
  tt <- row_tt(m, g)
  padj <- stats::p.adjust(tt$p, method = "BH")
  n_sig <- sum(padj < fdr & abs(tt$diff) > lfc)
  z <- if (is.null(xz)) t(scale(t(m))) else xz     # 基因内 z 与分组无关，置换时复用
  s <- score_cells(z, tt$t, top_n)
  list(n_sig = n_sig, auc = auc_w(s, g), t = tt$t, p = tt$p)
}

xz_cache <- t(scale(t(x)))
obs <- deg_metrics(x, y, xz_cache)
log_msg("观测：DEG ", obs$n_sig, " 个，top-", rcfg$top_n %||% 10,
        " 评分 AUC = ", round(obs$auc, 3))

# ------------------------------------------------- A. 标签置换 --------------
n_perm <- rcfg$n_perm %||% 1000
perm_n <- integer(n_perm); perm_auc <- numeric(n_perm)
for (b in seq_len(n_perm)) {
  g2 <- sample(y)
  r <- deg_metrics(x, g2, xz_cache)
  perm_n[b] <- r$n_sig
  perm_auc[b] <- r$auc
}
p_deg <- (1 + sum(perm_n >= obs$n_sig)) / (n_perm + 1)
p_auc <- (1 + sum(perm_auc >= obs$auc)) / (n_perm + 1)
summ <- data.frame(
  metric = c("n_DEG", "AUC_topN"),
  observed = c(obs$n_sig, round(obs$auc, 3)),
  null_mean = c(round(mean(perm_n), 2), round(mean(perm_auc, na.rm = TRUE), 3)),
  null_q95 = c(stats::quantile(perm_n, 0.95), round(stats::quantile(perm_auc, 0.95, na.rm = TRUE), 3)),
  p_empirical = c(p_deg, p_auc), n_perm = n_perm, stringsAsFactors = FALSE)
write_tsv(summ, "robustness_permutation", cfg)
log_msg("置换检验：DEG P=", signif(p_deg, 3), "；AUC P=", signif(p_auc, 3))

d_null <- data.frame(auc = perm_auc[is.finite(perm_auc)])
if (nrow(d_null) > 10) {
  g_null <- ggplot2::ggplot(d_null, ggplot2::aes(auc)) +
    ggplot2::geom_histogram(bins = 40, fill = "#D3D1C7", colour = "white", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = obs$auc, colour = "#C0392B", linewidth = 0.8) +
    ggplot2::annotate("text", x = obs$auc, y = Inf, vjust = 1.6, hjust = -0.05,
                      label = paste0("observed = ", round(obs$auc, 3),
                                     "\nP = ", signif(p_auc, 2)),
                      size = 3, colour = "#C0392B") +
    ggplot2::labs(x = "AUC under label permutation", y = "count",
                  title = paste0("Null distribution (", n_perm, " permutations)")) +
    theme_pub()
  save_fig(g_null, "robustness_null_auc", cfg, 5.4, 4)
}

# ------------------------------------------------- B. 多管线 ----------------
sig_file <- res_path(cfg, "tables", "deg_consensus.tsv")
cons <- if (file.exists(sig_file)) utils::read.delim(sig_file, stringsAsFactors = FALSE)$gene else character(0)
xz_all <- t(scale(t(x)))

fit_cv <- function(m, g, model = "glm", folds = 5, seed = 1) {
  set.seed(seed)
  idx <- sample(rep(seq_len(folds), length.out = length(g)))
  pred <- rep(NA_real_, length(g))
  for (k in seq_len(folds)) {
    tr <- idx != k
    if (sum(tr) < 6 || length(unique(g[tr])) < 2) next
    d <- data.frame(t(m[, tr, drop = FALSE]), y = as.numeric(g[tr]))
    nd <- data.frame(t(m[, !tr, drop = FALSE]))
    names(nd) <- names(d)[seq_len(ncol(nd))]
    pr <- tryCatch({
      if (model == "glm") {
        f <- stats::glm(y ~ ., data = d, family = stats::binomial())
        as.numeric(stats::predict(f, newdata = nd, type = "response"))
      } else {
        if (!requireNamespace("e1071", quietly = TRUE)) return(rep(NA_real_, sum(!tr)))
        yfac <- factor(d$y, levels = c(0, 1), labels = c("ctrl", "case"))
        sv <- e1071::svm(x = as.matrix(d[, -ncol(d)]), y = yfac,
                         kernel = "linear", cost = 1, scale = TRUE, probability = TRUE)
        pr <- attr(stats::predict(sv, as.matrix(nd), probability = TRUE), "probabilities")
        as.numeric(pr[, "case"])
      }
    }, error = function(e) rep(NA_real_, sum(!tr)))
    pred[!tr] <- pr
  }
  if (all(is.na(pred))) return(NA_real_)
  auc_w(pred, g)
}

grid <- expand.grid(zscore = c(TRUE, FALSE), gene_set = c("consensus", "all"),
                    model = c("glm", "svm"), stringsAsFactors = FALSE)
multi <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  gset <- if (grid$gene_set[i] == "consensus" && length(cons) >= 5) intersect(cons, rownames(x)) else rownames(x)
  m <- if (grid$zscore[i]) xz_all[gset, , drop = FALSE] else x[gset, , drop = FALSE]
  # glm 在 p >> n 时不收敛：限制特征数；svm 可放宽
  cap <- if (grid$model[i] == "glm") 15 else 200
  if (nrow(m) > cap) m <- m[order(apply(m, 1, stats::var), decreasing = TRUE)[seq_len(cap)], , drop = FALSE]
  aucs <- vapply(1:5, function(s) fit_cv(m, y, grid$model[i], seed = s), numeric(1))
  multi[[i]] <- data.frame(grid[i, ], n_gene = nrow(m),
                           auc_mean = round(mean(aucs, na.rm = TRUE), 3),
                           auc_sd = round(stats::sd(aucs, na.rm = TRUE), 3),
                           stringsAsFactors = FALSE)
  log_msg("multiverse ", i, "/", nrow(grid), "：", grid$gene_set[i], " + ",
          grid$model[i], " + zscore=", grid$zscore[i], " → AUC ",
          multi[[i]]$auc_mean)
}
multi_df <- do.call(rbind, multi)
write_tsv(multi_df, "robustness_multiverse", cfg)
log_msg("多管线 AUC 范围：", round(min(multi_df$auc_mean, na.rm = TRUE), 3), "–",
        round(max(multi_df$auc_mean, na.rm = TRUE), 3),
        "（写文章时报这个范围，不要只报最好看的那个）")

# Fig9C：多管线 AUC 区间。图的作用是让人一眼看出「换管线后结论是否还站得住」，
# 所以画区间而不是只画均值；两条同色线是同一组基因集 × 分类器下的 z-score / 原始尺度。
if (nrow(multi_df) >= 2) {
  multi_df$label <- factor(
    paste0(multi_df$gene_set, " genes / ", multi_df$model),
    levels = rev(unique(paste0(multi_df$gene_set, " genes / ", multi_df$model))))
  multi_df$scaling <- ifelse(multi_df$zscore, "z-scored", "raw")
  g_mv <- ggplot2::ggplot(
    multi_df, ggplot2::aes(auc_mean, label, colour = scaling)) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = pmax(0, auc_mean - auc_sd), xmax = pmin(1, auc_mean + auc_sd)),
      height = 0.16, linewidth = 0.5,
      position = ggplot2::position_dodge(width = 0.45)) +
    ggplot2::geom_point(size = 2.1, position = ggplot2::position_dodge(width = 0.45)) +
    ggplot2::scale_colour_manual(values = c(`z-scored` = "#1F6FB4", raw = "#C0392B"),
                                 name = NULL) +
    ggplot2::scale_x_continuous(limits = c(min(0.4, min(multi_df$auc_mean - multi_df$auc_sd, na.rm = TRUE) - 0.02),
                                           max(1.0, max(multi_df$auc_mean + multi_df$auc_sd, na.rm = TRUE) + 0.02))) +
    ggplot2::labs(x = "Cross-validated AUC (mean +/- SD over 5 repeats)", y = NULL,
                  title = "Multiverse of modelling pipelines") +
    theme_pub() +
    ggplot2::theme(legend.position = "top")
  save_fig(g_mv, "robustness_multiverse", cfg, 6, 3.6)
}

# ------------------------------------------------- C. 药物连通性（可选） ----
lincs <- rcfg$lincs_path
if (!is.null(lincs) && nzchar(lincs) && file.exists(proj_path(lincs))) {
  lc <- data.table::fread(proj_path(lincs), showProgress = FALSE)
  need <- c("drug", "gene", "zscore")
  if (!all(need %in% names(lc))) {
    log_msg("LINCS 文件缺少列：", paste(setdiff(need, names(lc)), collapse = ", "), "，跳过 C 部分")
  } else {
    # 疾病特征：上调/下调基因集（取自 02 的共识基因 logFC）
    if (!file.exists(sig_file)) {
      log_msg("缺 deg_consensus.tsv，C 部分只输出空表")
      lc <- data.table::data.table(drug = character(0), gene = character(0), zscore = numeric(0))
    }
    cons_lfc <- if (file.exists(sig_file)) {
      d <- utils::read.delim(sig_file, stringsAsFactors = FALSE)
      stats::setNames(as.numeric(d[[grep("^logFC", names(d), value = TRUE)[1]]]), d$gene)
    } else numeric(0)
    up <- names(cons_lfc)[cons_lfc > 0]
    dn <- names(cons_lfc)[cons_lfc < 0]

    ks_cs <- function(d) {
      d <- d[order(-d$zscore), ]
      r <- rank(-d$zscore)
      nr <- nrow(d)
      hit_up <- r[match(up, d$gene)]; hit_up <- hit_up[!is.na(hit_up)]
      hit_dn <- r[match(dn, d$gene)]; hit_dn <- hit_dn[!is.na(hit_dn)]
      if (length(hit_up) < 5 || length(hit_dn) < 5) return(NA_real_)
      ks_up <- max(seq_along(hit_up) / length(hit_up) - hit_up / nr)
      ks_dn <- max(seq_along(hit_dn) / length(hit_dn) - hit_dn / nr)
      if (sign(ks_up) == sign(ks_dn)) 0 else ks_up - ks_dn
    }
    cs <- lc[, list(cs = ks_cs(.SD)), by = drug]
    cs <- cs[is.finite(cs)][order(cs)]
    write_tsv(as.data.frame(cs), "cmap_connectivity", cfg)
    log_msg("药物连通性：", nrow(cs), " 个化合物；cs 最负的 5 个（可能逆转疾病特征）：",
            paste(head(cs$drug, 5), collapse = ", "))
    log_msg("措辞提醒：这是基于转录组签名的预测，不是疗效证据，必须写成 candidate compounds")
  }
} else {
  log_msg("未配置 drylab.robustness.lincs_path，跳过药物连通性（可选模块）")
}

log_msg("13 完成")

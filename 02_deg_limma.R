# =============================================================================
# 02_deg_limma.R —— 课题一（2）：差异表达分析
#
# 输入：data/processed/<GSE>_expr.rds / _pheno.rds（由 01 生成）
# 输出：results/tables/deg_<GSE>.tsv          全基因结果（含 logFC / P / adj.P）
#       results/tables/deg_consensus.tsv      训练集与验证集方向一致的稳健差异基因
#       results/figures/volcano_<GSE>.pdf
#       results/figures/heatmap_top50_<GSE>.pdf
#
# 注意：
#   - limma 的 logFC 是 log2 尺度；芯片数据 FDR + |logFC| 双阈值是惯例
#   - 验证集不必单独卡阈值，但"方向一致"是后续外部验证成立的前提
#   - 想加协变量（年龄/性别）就改 design，不要事后回归残差
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("02_deg_limma", cfg)
load_pkgs(c("ggplot2", "limma", "pheatmap"))
set.seed(cfg_get(cfg, "project.seed", 1))

fdr_cut  <- cfg_get(cfg, "bulk.deg.fdr", 0.05)
lfc_cut  <- cfg_get(cfg, "bulk.deg.logfc", 1)
demo     <- isTRUE(cfg_get(cfg, "run.demo", FALSE))

load_pair <- function(gse) {
  e <- readRDS(proj_path("data", "processed", paste0(gse, "_expr.rds")))
  ph <- readRDS(proj_path("data", "processed", paste0(gse, "_pheno.rds")))
  stopifnot(identical(colnames(e), as.character(ph$sample)))
  list(expr = e, pheno = ph)
}

run_limma <- function(gse, cfg) {
  d <- load_pair(gse)
  expr <- d$expr; pheno <- d$pheno
  pheno$group <- factor(pheno$group, levels = c("control", "case"))

  # 设计矩阵：分组为主，批次（若存在且非全同）作为协变量
  rhs <- "group"
  if ("batch" %in% names(pheno) && length(unique(stats::na.omit(pheno$batch))) > 1) rhs <- "~ 0 + group + batch" else rhs <- "~ 0 + group"
  design <- stats::model.matrix(stats::as.formula(rhs), data = pheno)
  colnames(design) <- sub("^group", "", colnames(design))
  fit <- limma::lmFit(expr, design)
  cm <- limma::makeContrasts(case - control, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
  tt <- limma::topTable(fit2, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  tt$dataset <- gse
  tt$sig <- ifelse(tt$adj.P.Val < fdr_cut & abs(tt$logFC) > lfc_cut,
                   ifelse(tt$logFC > 0, "Up", "Down"), "NS")
  list(tt = tt, expr = expr, pheno = pheno)
}

plot_volcano <- function(tt, gse) {
  tt$neglogp <- -log10(pmax(tt$P.Value, 1e-300))
  pal <- c(Up = "#C0392B", Down = "#1F6FB4", NS = "grey75")
  g <- ggplot2::ggplot(tt, ggplot2::aes(logFC, neglogp, colour = sig)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.75) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", linewidth = 0.3) +
    ggplot2::labs(title = paste0(gse, ": POAG vs control"),
                  x = expression(log[2]~"fold change"), y = expression(-log[10]~"P")) +
    theme_pub()
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    lab <- subset(tt, sig != "NS")
    lab <- lab[order(lab$P.Value), ][seq_len(min(15, nrow(lab))), ]
    g <- g + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene),
                                      size = 2.6, max.overlaps = 30, show.legend = FALSE)
  }
  g
}

plot_heatmap <- function(expr, pheno, tt, gse, n = 50) {
  top <- tt[tt$sig != "NS", ]
  top <- top[order(top$P.Value), ]
  genes <- head(top$gene, n)
  if (length(genes) < 5) {
    log_msg(gse, " 显著基因不足 5 个，跳过热图")
    return(invisible(NULL))
  }
  m <- expr[genes, , drop = FALSE]
  m <- t(scale(t(m)))
  m[m > 2.5] <- 2.5; m[m < -2.5] <- -2.5
  ann <- data.frame(group = pheno$group, row.names = pheno$sample)
  pal <- if (nlevels(pheno$group) == 2) c("#1F6FB4", "#C0392B") else NULL
  f <- proj_path(cfg_get(cfg, "project.outdir", "results"), "figures",
                 paste0("heatmap_top", n, "_", gse, ".pdf"))
  # 先画一次拿到 gtable，再分别灌进 pdf / png。
  # 原来只 grDevices::pdf() 一发，投稿拼图（16 用 magick 拼 panel）时这张就没有 png 可用。
  ph <- pheatmap::pheatmap(m, annotation_col = ann,
                           annotation_colors = if (is.null(pal)) NULL else list(group = setNames(pal, levels(pheno$group))),
                           show_colnames = FALSE, cluster_cols = TRUE, cluster_rows = TRUE,
                           fontsize_row = 6, border_color = NA,
                           color = grDevices::colorRampPalette(c("#1F6FB4", "white", "#C0392B"))(101),
                           main = paste0(gse, " top ", length(genes), " DEGs"),
                           silent = TRUE)
  base <- sub("\\.pdf$", "", f)
  grDevices::pdf(paste0(base, ".pdf"), width = 7, height = 8)
  grid::grid.newpage(); grid::grid.draw(ph$gtable)
  grDevices::dev.off()
  grDevices::png(paste0(base, ".png"), width = 7, height = 8, units = "in",
                 res = cfg_get(cfg, "project.dpi", 300))
  grid::grid.newpage(); grid::grid.draw(ph$gtable)
  grDevices::dev.off()
  invisible(f)
}

# ------------------------------------------------------------------ 主流程 ----

train_gse <- cfg_get(cfg, "bulk.train.gse", "GSE27276")
if (demo) train_gse <- "DEMO"

res_train <- run_limma(train_gse, cfg)
tt_train <- res_train$tt
write_tsv(tt_train, paste0("deg_", train_gse), cfg)
save_fig(plot_volcano(tt_train, train_gse), paste0("volcano_", train_gse), cfg, 6, 5)
plot_heatmap(res_train$expr, res_train$pheno, tt_train, train_gse)

log_msg("训练集 ", train_gse, "：上调 ", sum(tt_train$sig == "Up"),
        "，下调 ", sum(tt_train$sig == "Down"))

# 验证集：只算差异，不重设阈值；用于筛"方向一致"的稳健基因
valid_gse <- if (demo) c("DEMO_V") else vapply(cfg_get(cfg, "bulk.valid", list()), function(x) x$gse, character(1))
valid_res <- list()
for (v in valid_gse) {
  f <- proj_path("data", "processed", paste0(v, "_expr.rds"))
  if (!file.exists(f)) { log_msg("跳过缺失验证集：", v); next }
  r <- run_limma(v, cfg)
  write_tsv(r$tt, paste0("deg_", v), cfg)
  save_fig(plot_volcano(r$tt, v), paste0("volcano_", v), cfg, 6, 5)
  valid_res[[v]] <- r$tt
  log_msg("验证集 ", v, "：上调 ", sum(r$tt$sig == "Up"), "，下调 ", sum(r$tt$sig == "Down"))
}

# 共识基因：训练集显著 + 在所有验证集都测得 + 至少在其中一个验证集同向复现。
# 原实现写的是 `is.na(lf) | (...)`，把"该验证集根本没测到这个基因"当成通过；
# 两个平台基因数不同（13,118 vs 13,468），凭空多算出一批"共识"，而 Venn 中心是 0——
# 图和表当场打架。缺失只能算"不可验证"，不能算"通过"。
if (length(valid_res)) {
  sig_tr <- tt_train[tt_train$sig != "NS", c("gene", "logFC", "adj.P.Val")]
  names(sig_tr)[2:3] <- c("logFC_train", "fdr_train")
  for (v in names(valid_res)) {
    vt <- valid_res[[v]][, c("gene", "logFC", "P.Value")]
    names(vt) <- c("gene", paste0(c("logFC_", "p_"), v))
    sig_tr <- merge(sig_tr, vt, by = "gene", all.x = TRUE)
  }
  measured <- rep(TRUE, nrow(sig_tr))
  replicated <- rep(FALSE, nrow(sig_tr))
  for (v in names(valid_res)) {
    lf <- sig_tr[[paste0("logFC_", v)]]
    pv <- sig_tr[[paste0("p_", v)]]
    m <- !is.na(lf)
    measured <- measured & m
    hit <- m & !is.na(pv) & pv < 0.05 & sign(lf) == sign(sig_tr$logFC_train)
    replicated <- replicated | hit
    log_msg("  ", v, "：同向复现 ", sum(hit), " / ", nrow(sig_tr),
            "（该平台未测到 ", sum(!m), " 个）")
  }
  consensus <- sig_tr[measured & replicated, ]
  consensus <- consensus[order(consensus$fdr_train), ]
  write_tsv(consensus, "deg_consensus", cfg)
  log_msg("共识基因：", nrow(consensus), " 个（各验证集均测得，且至少一个验证集同向复现）")
  if (!nrow(consensus)) {
    log_msg("[WARN] 共识集为空：04/05 会自动回退到训练集 DEG，",
            "但正文凡引用 {{N_CONSENSUS}} 的句子都要按真实结果改写")
  }

  # --- Fig2D：发现集与可复现集的关系 -------------------------------------------
  # 复现集必然是发现集的子集；画成两个交叠的圆，读者只会以为半径画错了。
  # 改成嵌套图：外圈 = 发现集，内圈 = 可复现集，面积按比例。
  n_a <- nrow(sig_tr); n_b <- nrow(consensus)
  circ <- function(cx, cy, r, n = 240) {
    th <- seq(0, 2 * pi, length.out = n)
    data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
  }
  r_out <- 0.95
  r_in <- max(r_out * sqrt(max(n_b, 1) / max(n_a, 1)), 0.17)
  g_venn <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = circ(1.6, 0, r_out), ggplot2::aes(x, y),
                          fill = "#1F6FB4", alpha = 0.12,
                          colour = "#1F6FB4", linewidth = 0.4) +
    ggplot2::geom_polygon(data = circ(1.6, 0, r_in), ggplot2::aes(x, y),
                          fill = "#C0392B", alpha = 0.22,
                          colour = "#C0392B", linewidth = 0.4) +
    ggplot2::annotate("text", x = 1.6, y = 0, label = n_b, size = 5.6, colour = "grey10") +
    ggplot2::annotate("text", x = 1.6, y = r_out + 0.10,
                      label = paste0("Discovery DEGs (n = ", n_a, ")"),
                      size = 3.1, colour = "grey30") +
    ggplot2::annotate("text", x = 1.6, y = -r_out - 0.10,
                      label = paste0("Replicated in >=1 cohort (n = ", n_b, ")"),
                      size = 3.1, colour = "grey30") +
    ggplot2::coord_equal(xlim = c(0.55, 2.65), ylim = c(-1.30, 1.30), expand = FALSE) +
    ggplot2::labs(title = "Replication of discovery DEGs",
                  subtitle = paste("inner set: nominal p < 0.05 with concordant",
                                   "direction, measured in every cohort")) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11),
                   plot.subtitle = ggplot2::element_text(size = 7.4, colour = "grey40"),
                   plot.margin = ggplot2::margin(6, 6, 6, 6))
  save_fig(g_venn, "venn_consensus", cfg, 5.2, 4.4)
}

log_msg("02 完成")

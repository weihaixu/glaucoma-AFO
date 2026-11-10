# =============================================================================
# 14_cross_species.R —— 课题一（干实验版）：跨物种保守性
#
# 目的：没有动物实验时，用「人 / 猴 / 鼠」三套单细胞数据回答同一个问题——
#       在人小梁网里找到的疾病相关表达程序，在猴 POAG 模型和鼠高眼压模型里
#       是否以同一方向重现。这是湿实验里"小鼠模型验证"的计算等价物，
#       但它证明的是「保守性」，不是「因果性」，措辞必须区分。
#
# 三步：
#   A. 伪批量差异：每个物种、每个目标细胞类型内部做 case vs control，
#      只在「人共识基因」这一小撮基因上做检验（省时间，也避免多重检验爆炸）
#   B. 一致性统计：方向一致比例 + 二项检验；logFC 的 Spearman 相关
#   C. 模块打分：用人的上调/下调基因集给每个细胞打分（AddModuleScore），
#      比较 case 与 control 的 AUC（0.5 = 无差异；>0.5 表示疾病组更高）
#
# 输入：
#   drylab.cross_species.objects      物种 -> Seurat 对象路径
#   results/tables/deg_consensus.tsv  人共识基因（02 的产物）
#   drylab.cross_species.ortholog_table（可选）两列 tsv：human_symbol / species_symbol
# 输出：
#   results/tables/xspecies_concordance.tsv   每个 物种×细胞类型 的一致性
#   results/tables/xspecies_de.tsv            逐基因的跨物种 logFC 对照
#   results/tables/xspecies_modulescore.tsv   模块打分的 AUC
#   results/figures/xspecies_*.pdf
#
# 重要局限（写论文时必须写进 Limitations）：
#   没有 ortholog 表时用 toupper() 做近似映射，会漏掉命名不一致的一对一同源基因；
#   伪批量只用两个均值相减，不给细胞水平的置信区间，不能当差异分析的正式结果。
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("14_cross_species", cfg)
load_pkgs(c("ggplot2", "Matrix", "data.table"))
set.seed(cfg_get(cfg, "project.seed", 1))

xs <- cfg_get(cfg, "drylab.cross_species", list())
if (!length(xs$objects)) {
  stop("config 里 drylab.cross_species.objects 为空：\n",
       "  至少给一个物种，例如 macaque: \"data/processed/scRNA_macaque.rds\"",
       call. = FALSE)
}
if (!requireNamespace("Seurat", quietly = TRUE)) stop("需要 Seurat 包", call. = FALSE)

# ---------------------------------------------------------------- 工具函数 ----

get_mat <- function(o, what = "data") {
  tryCatch(Seurat::GetAssayData(o, assay = "RNA", layer = what),
           error = function(e) Seurat::GetAssayData(o, assay = "RNA", slot = what))
}

# 分组列：优先用户指定列；否则在前几列 meta 里按 group_regex 推断
resolve_condition <- function(obj, col, regex) {
  md <- obj@meta.data
  if (!is.null(col) && col %in% names(md)) {
    v <- as.character(md[[col]])
    if (length(unique(v[!is.na(v)])) >= 2) return(v)
  }
  cols <- names(md)[seq_len(min(6, ncol(md)))]
  for (cc in cols) {
    v <- as.character(md[[cc]])
    if (length(unique(v)) < 2 || length(unique(v)) > 10) next
    hit <- grepl(regex, v, ignore.case = TRUE)
    if (sum(hit) > 0 && sum(!hit) > 0) {
      log_msg("  分组列推断为 meta$", cc, "（pattern: ", regex, "）")
      return(ifelse(hit, "case", "control"))
    }
  }
  NULL
}

# 表达矩阵确认是 log 归一化：max > 50 视为原始 counts，补做 LogNormalize
ensure_lognorm <- function(obj) {
  m <- get_mat(obj, "data")
  if (length(m) && max(m, na.rm = TRUE) > 50) {
    log_msg("  数值范围像原始 counts，执行 LogNormalize")
    obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize",
                                 scale.factor = 1e4, verbose = FALSE)
  }
  obj
}

# 同源映射：有表用表，没表用 toupper 近似
make_mapper <- function(species, orth_path) {
  if (!is.null(orth_path) && nzchar(orth_path) && file.exists(proj_path(orth_path))) {
    ot <- utils::read.delim(proj_path(orth_path), stringsAsFactors = FALSE)
    h <- grep("^human", names(ot), value = TRUE, ignore.case = TRUE)[1]
    s <- grep(species, names(ot), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(h) || is.na(s)) {
      log_msg("  ortholog_table 需要列名含 human 与 ", species, "，退回 toupper 映射")
    } else {
      map <- stats::setNames(toupper(ot[[h]]), toupper(ot[[s]]))
      map <- map[!duplicated(names(map))]
      log_msg("  使用同源表：", length(map), " 条映射")
      return(function(x) unname(map[toupper(x)]))
    }
  }
  function(x) {
    y <- toupper(x)
    y[is.na(y) | y == ""] <- NA_character_
    y
  }
}

# 伪批量：细胞类型 × 分组 的均值矩阵（基因 × 组）
pseudobulk <- function(obj, cond, celltypes, min_cells) {
  m <- get_mat(obj, "data")
  out <- list()
  keep_genes <- rownames(m)
  for (ct in celltypes) {
    idx <- which(as.character(obj$cell_type) == ct)
    if (length(idx) < min_cells) {
      log_msg("  ", ct, "：细胞数 ", length(idx), " < min_cells，跳过")
      next
    }
    for (grp in c("control", "case")) {
      sel <- idx[cond[idx] == grp]
      if (length(sel) < min_cells) {
        log_msg("  ", ct, "/", grp, "：细胞数 ", length(sel), " 不足，跳过该组")
        next
      }
      v <- Matrix::rowMeans(m[keep_genes, sel, drop = FALSE])
      key <- paste(ct, grp, sep = "|")
      out[[key]] <- list(vec = v, n = length(sel))
    }
  }
  out
}

# ---------------------------------------------------------------- 1) 人共识 ----
sig_file <- res_path(cfg, "tables", "deg_consensus.tsv")
require_file(sig_file, "先跑 02_deg_limma.R 得到 deg_consensus.tsv")
sig <- utils::read.delim(sig_file, stringsAsFactors = FALSE)
lf_col <- grep("^logFC", names(sig), value = TRUE)[1]
if (is.na(lf_col)) stop("共识基因表里找不到 logFC 列", call. = FALSE)
human_lfc <- stats::setNames(as.numeric(sig[[lf_col]]), sig$gene)
human_lfc <- human_lfc[is.finite(human_lfc) & human_lfc != 0]
log_msg("人共识基因：", length(human_lfc), " 个（上调 ", sum(human_lfc > 0),
        "，下调 ", sum(human_lfc < 0), "）")

min_cells <- xs$min_cells %||% 20
cond_col  <- xs$condition_col %||% "condition"
grp_regex <- xs$group_regex %||% "POAG|glaucoma|DEX|treated|chronic|disease"
orth_path <- xs$ortholog_table
ct_cfg    <- xs$celltypes %||% list()

# 人自身也跑一遍作为阳性对照：如果人 vs 人都不一致，说明流程有问题
if (isTRUE(xs$include_self) && file.exists(proj_path(cfg_get(cfg, "scrna.rds", "")))) {
  xs$objects <- c(list(human = cfg_get(cfg, "scrna.rds")), xs$objects)
}

# ---------------------------------------------------------------- 2) 主循环 ----
conc_rows <- list(); de_rows <- list(); ms_rows <- list()

for (sp in names(xs$objects)) {
  path <- proj_path(xs$objects[[sp]])
  if (!file.exists(path)) {
    log_msg("跳过 ", sp, "：找不到 ", path)
    next
  }
  log_msg("===== 处理物种：", sp, " =====")
  obj <- readRDS(path)
  if (!"cell_type" %in% names(obj@meta.data)) {
    log_msg("跳过 ", sp, "：对象缺 cell_type 列（先用 06 注释或写入 meta.data）")
    next
  }
  obj <- ensure_lognorm(obj)

  cts <- ct_cfg[[sp]] %||% xs$celltypes_default %||% c("TM", "JCT")
  cts <- intersect(cts, unique(as.character(obj$cell_type)))
  if (!length(cts)) {
    log_msg("跳过 ", sp, "：cell_type 中没有 ", paste(cts, collapse = "/"),
            "；现有：", paste(unique(obj$cell_type), collapse = ", "))
    next
  }

  cond <- resolve_condition(obj, cond_col, grp_regex)
  if (is.null(cond)) {
    log_msg("跳过 ", sp, " 的差异部分：无法确定 case/control 分组（在 config 指定 condition_col）")
  } else {
    obj$condition <- cond
  }

  to_h <- make_mapper(sp, orth_path)
  shared_h <- NULL

  # ---- A. 伪批量差异 ----
  if (!is.null(cond)) {
    pb <- pseudobulk(obj, cond, cts, min_cells)
    for (ct in cts) {
      k1 <- paste(ct, "control", sep = "|"); k2 <- paste(ct, "case", sep = "|")
      if (is.null(pb[[k1]]) || is.null(pb[[k2]])) next
      genes_sp <- rownames(get_mat(obj, "data"))
      gh <- to_h(genes_sp)
      ok <- !is.na(gh) & gh %in% names(human_lfc)
      if (sum(ok) < 50) {
        log_msg("  ", ct, "：同源映射后只剩 ", sum(ok), " 个基因，跳过（检查基因名格式）")
        next
      }
      lfc <- pb[[k2]]$vec - pb[[k1]]$vec     # log 空间均值差 ≈ logFC
      lfc_h <- stats::setNames(lfc[ok], gh[ok])
      lfc_h <- lfc_h[!duplicated(names(lfc_h))]

      # 细胞水平 Wilcoxon：只对共识基因做，控制计算量
      gsel <- genes_sp[ok]
      gsel <- gsel[!duplicated(gh[ok])]
      cell_ctrl <- which(obj$cell_type == ct & obj$condition == "control")
      cell_case <- which(obj$cell_type == ct & obj$condition == "case")
      m <- get_mat(obj, "data")
      pv <- rep(NA_real_, length(gsel))
      names(pv) <- to_h(gsel)
      for (g in gsel) {
        x <- as.numeric(m[g, cell_ctrl]); y <- as.numeric(m[g, cell_case])
        if (sum(x > 0) + sum(y > 0) < 10) next
        pv[to_h(g)] <- tryCatch(suppressWarnings(stats::wilcox.test(y, x)$p.value),
                                error = function(e) NA_real_)
      }
      padj <- stats::p.adjust(pv, method = "BH")
      shared <- intersect(names(lfc_h), names(pv))
      if (!length(shared)) next
      lfc_h <- lfc_h[shared]; padj <- padj[shared]

      hl <- human_lfc[shared]
      sl <- lfc_h[shared]
      match_dir <- sign(hl) == sign(sl)
      conc <- mean(match_dir, na.rm = TRUE)
      binom_p <- tryCatch(stats::binom.test(sum(match_dir, na.rm = TRUE),
                                            sum(!is.na(match_dir)), p = 0.5)$p.value,
                          error = function(e) NA_real_)
      sp_r <- suppressWarnings(stats::cor.test(hl, sl, method = "spearman"))

      conc_rows[[length(conc_rows) + 1]] <- data.frame(
        species = sp, celltype = ct,
        n_cells_case = length(cell_case), n_cells_control = length(cell_ctrl),
        n_shared = length(shared),
        n_up_shared = sum(hl > 0 & sl > 0), n_down_shared = sum(hl < 0 & sl < 0),
        concordance = round(conc, 3), binom_p = signif(binom_p, 3),
        spearman_rho = round(unname(sp_r$estimate), 3),
        spearman_p = signif(sp_r$p.value, 3),
        n_sig_padj05 = sum(padj < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      de_rows[[length(de_rows) + 1]] <- data.frame(
        species = sp, celltype = ct, gene = shared,
        human_logFC = round(hl, 4), species_logFC = round(sl, 4),
        direction_match = match_dir, padj_species = signif(padj, 3),
        stringsAsFactors = FALSE
      )
      log_msg("  ", ct, "：shared=", length(shared), "，方向一致 ",
              sprintf("%.1f%%", 100 * conc), "，rho=", round(unname(sp_r$estimate), 2))
      shared_h <- shared
    }
  }

  # ---- C. 模块打分（不依赖分组也能跑，只要一个对象）----
  if (!is.null(cond) && requireNamespace("Seurat", quietly = TRUE)) {
    up_h <- names(human_lfc)[human_lfc > 0]
    dn_h <- names(human_lfc)[human_lfc < 0]
    genes_sp <- rownames(get_mat(obj, "data"))
    gh <- to_h(genes_sp)
    map_up <- genes_sp[!is.na(gh) & gh %in% up_h]
    map_dn <- genes_sp[!is.na(gh) & gh %in% dn_h]
    for (set_name in c("up", "down")) {
      gs <- if (set_name == "up") map_up else map_dn
      if (length(gs) < 5) next
      obj2 <- tryCatch(Seurat::AddModuleScore(obj, features = list(gs), name = "xsm",
                                              seed = 1, nbin = 12),
                       error = function(e) NULL)
      if (is.null(obj2)) next
      sc <- obj2$xs1
      for (ct in cts) {
        i1 <- which(obj2$cell_type == ct & obj2$condition == "control")
        i2 <- which(obj2$cell_type == ct & obj2$condition == "case")
        if (length(i1) < min_cells || length(i2) < min_cells) next
        w <- suppressWarnings(stats::wilcox.test(sc[i2], sc[i1]))
        auc <- unname(w$statistic) / (length(i1) * length(i2))
        ms_rows[[length(ms_rows) + 1]] <- data.frame(
          species = sp, celltype = ct, gene_set = set_name,
          n_genes = length(gs), mean_control = round(mean(sc[i1]), 3),
          mean_case = round(mean(sc[i2]), 3),
          delta = round(mean(sc[i2]) - mean(sc[i1]), 3),
          auc = round(auc, 3), p = signif(w$p.value, 3),
          stringsAsFactors = FALSE
        )
        log_msg("  模块[", set_name, "] ", ct, "：AUC=", round(auc, 3),
                "，delta=", round(mean(sc[i2]) - mean(sc[i1]), 2))
      }
    }
  }
}

# ---------------------------------------------------------------- 3) 输出 ----

# 空结果也要写出带列名的空表，避免下游/审稿人以为脚本没跑
empty_df <- function(cols) {
  d <- as.data.frame(stats::setNames(replicate(length(cols), character(0), simplify = FALSE), cols))
  d
}
conc <- if (length(conc_rows)) data.table::rbindlist(conc_rows, fill = TRUE) else
  empty_df(c("species", "celltype", "n_cells_case", "n_cells_control", "n_shared",
             "n_up_shared", "n_down_shared", "concordance", "binom_p",
             "spearman_rho", "spearman_p", "n_sig_padj05"))
de <- if (length(de_rows)) data.table::rbindlist(de_rows, fill = TRUE) else
  empty_df(c("species", "celltype", "gene", "human_logFC", "species_logFC",
             "direction_match", "padj_species"))
ms <- if (length(ms_rows)) data.table::rbindlist(ms_rows, fill = TRUE) else
  empty_df(c("species", "celltype", "gene_set", "n_genes", "mean_control",
             "mean_case", "delta", "auc", "p"))

write_tsv(conc, "xspecies_concordance", cfg)
write_tsv(de, "xspecies_de", cfg)
write_tsv(ms, "xspecies_modulescore", cfg)

# ---------------------------------------------------------------- 4) 出图 ----

if (nrow(conc)) {
  conc$label <- sprintf("%.0f%%\n(n=%d)", 100 * conc$concordance, conc$n_shared)
  g1 <- ggplot2::ggplot(conc, ggplot2::aes(x = celltype, y = species, fill = concordance)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3, lineheight = 0.9,
                       colour = ifelse(conc$concordance > 0.75, "white", "grey15")) +
    # 图内文字必须 ASCII：本机 R 在 C locale 下跑，UTF-8 中文进 PNG 设备会变乱码
    ggplot2::scale_fill_gradient2(low = "#2C6FA6", mid = "#F2F2F2", high = "#B3322C",
                                  midpoint = 0.5, limits = c(0.2, 0.8),
                                  name = "Concordant\ndirection") +
    ggplot2::labs(x = NULL, y = NULL, title = "Cross-species concordance of disease programs") +
    theme_pub() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  save_fig(g1, "xspecies_concordance_heatmap", cfg, 6, 2.2 + 0.5 * length(unique(conc$species)))
}

if (nrow(de) && length(unique(de$species)) >= 1) {
  de$panel <- paste(de$species, de$celltype, sep = " / ")
  g2 <- ggplot2::ggplot(de, ggplot2::aes(human_logFC, species_logFC)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(colour = direction_match), size = 1.2, alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#B3322C", `FALSE` = "#7F8C99"),
                                 labels = c("Discordant", "Concordant"), name = NULL) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         linewidth = 0.4, colour = "grey30") +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::labs(x = "Human TM logFC", y = "Cross-species logFC",
                  title = "Consensus genes: effect-size comparison across species") +
    theme_pub()
  save_fig(g2, "xspecies_logfc_scatter", cfg, 3.2 * min(3, length(unique(de$panel))), 3)
}

if (nrow(ms)) {
  ms$panel <- paste(ms$species, ms$celltype, sep = " / ")
  g3 <- ggplot2::ggplot(ms, ggplot2::aes(x = gene_set, y = auc)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3,
                        colour = "grey50") +
    ggplot2::geom_col(ggplot2::aes(fill = gene_set), width = 0.6, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c(up = "#B3322C", down = "#2C6FA6")) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", auc)), vjust = -0.4, size = 2.8) +
    ggplot2::facet_wrap(~panel) +
    ggplot2::coord_cartesian(ylim = c(0.3, 0.75)) +
    ggplot2::labs(x = NULL, y = "AUC (case vs control)",
                  title = "Cross-species module scoring of human disease genes") +
    theme_pub()
  save_fig(g3, "xspecies_modulescore", cfg, 3 * min(3, length(unique(ms$panel))), 3)
}

# ---------------------------------------------------------------- 5) 结论提示 ----
if (nrow(conc)) {
  ok <- conc[conc$concordance >= 0.7 & conc$binom_p < 0.05, ]
  log_msg("跨物种保守性判据（concordance ≥ 0.70 且 binom_p < 0.05）：",
          if (nrow(ok)) paste(unique(paste(ok$species, ok$celltype)), collapse = "; ")
          else "无。此时不要写『跨物种保守』，改写为『人源数据中的发现未能在其他物种重现』")
  log_msg("提示：toupper 映射会漏一对多同源基因，正式分析请提供 ortholog_table")
} else {
  log_msg("没有可比较的数据对：确认各物种对象含 TM/JCT 标签且存在 case/control 分组")
}
log_msg("14_cross_species 完成")

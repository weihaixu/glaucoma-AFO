# =============================================================================
# 01_geo_prepare.R —— 课题一（1）：数据下载、探针重注释、标准化、批次校正、QC
#
# 输入：GEO series matrix（自动下载到 data/raw/）
#       GPL 注释表（自动下载到 data/annotation/，或你手动放一份）
# 输出：data/processed/<GSE>_expr.rds   基因 × 样本 表达矩阵（log2 尺度）
#       data/processed/<GSE>_pheno.rds  样本表型（含 group / batch / 临床变量）
#       results/tables/qc_summary.tsv
#       results/figures/pca_<GSE>_{raw,clean}.pdf
#
# 关键点（审稿人真正会查的）：
#   1. 不用 GEO 自带的 GDS 矩阵，必须回到 GPL 注释重新映射探针
#   2. 同基因多探针取平均表达最高者，不是取平均（取平均会系统性低估差异倍数）
#   3. 批次校正必须在「样本 × 基因」矩阵上做，且 mod 里保留分组变量
#   4. 离群样本先标记、人工看一下再删，不要脚本自动删
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("01_geo_prepare", cfg)
load_pkgs(c("ggplot2", "limma", "sva"))
set.seed(cfg_get(cfg, "project.seed", 1))

demo <- isTRUE(cfg_get(cfg, "run.demo", FALSE))

# ------------------------------------------------------------------ 工具 ----

# 从 series matrix 头部读平台号
series_platform <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  head_lines <- readLines(con, n = 300, warn = FALSE)
  hit <- grep("^!Series_platform_id", head_lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(gsub('"', "", strsplit(hit[1], "\t")[[1]][2]))
}

# 从表型表里按正则推分组：命中即病例，其余对照
make_group <- function(pheno, pattern) {
  txt <- apply(pheno, 1, function(r) paste(r, collapse = " | "))
  g <- ifelse(grepl(pattern, txt, ignore.case = TRUE), "case", "control")
  if (length(unique(g)) < 2) {
    warning("按 group_regex='", pattern, "' 未分出两组，请检查表型字段：\n  ",
            paste(head(txt, 3), collapse = "\n  "))
    return(factor(g, levels = c("control", "case")))
  }
  factor(g, levels = c("control", "case"))
}

# 读 pheno 里的临床数值列（IOP / MD / CCT / age），字段名宽松匹配
extract_numeric_traits <- function(pheno) {
  out <- list()
  keys <- list(IOP = "iop|眼压|intraocular", MD = "\\bmd\\b|mean deviation|视野",
               CCT = "cct|corneal thickness|角膜厚度", age = "^age|年龄")
  all_txt <- apply(pheno, 1, function(r) paste(r, collapse = " | "))
  for (nm in names(keys)) {
    hit_col <- grep(keys[[nm]], names(pheno), ignore.case = TRUE, value = TRUE)
    if (length(hit_col)) {
      v <- suppressWarnings(as.numeric(gsub("[^0-9.\\-]", "", pheno[[hit_col[1]]])))
      if (sum(!is.na(v)) >= max(5, 0.5 * nrow(pheno))) out[[nm]] <- v
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

pca_plot <- function(expr, pheno, title, ntop = 2000) {
  v <- apply(expr, 1, stats::var, na.rm = TRUE)
  top <- names(sort(v, decreasing = TRUE))[seq_len(min(ntop, length(v)))]
  pc <- stats::prcomp(t(expr[top, , drop = FALSE]), scale. = TRUE)
  ve <- round(100 * summary(pc)$importance[2, 1:2], 1)
  df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                   group = pheno$group, batch = pheno$batch,
                   sample = rownames(pheno))
  ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, fill = group, shape = batch)) +
    ggplot2::geom_point(size = 2.6, colour = "grey20", stroke = 0.3) +
    ggplot2::stat_ellipse(ggplot2::aes(group = group), type = "norm",
                          linewidth = 0.4, show.legend = FALSE) +
    ggplot2::labs(title = title,
                  x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) +
    theme_pub()
}

# --------------------------------------------------------------- 单数据集 ----

prepare_one <- function(gse, gpl, group_regex, cfg, role = "train") {
  log_msg("---- 处理 ", gse, " (", role, ") ----")
  f <- proj_path("data", "raw", paste0(gse, "_series_matrix.txt.gz"))
  if (!file.exists(f)) {
    download_if_missing(
      paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", geo_series_dir(gse),
             "/", gse, "/matrix/", gse, "_series_matrix.txt.gz"), f)
  }
  sm <- read_series_matrix(f)
  if (is.null(gpl) || is.na(gpl) || gpl == "") gpl <- series_platform(f)

  annot <- read_gpl_annotation(gpl, cfg)
  expr <- collapse_probes(sm$expr, annot)

  # 1) log2 判定与转换
  if (!is_log_scale(expr, cfg_get(cfg, "bulk.qc.log2_if_max_below", 50))) {
    expr[expr < 0] <- 0
    expr <- log2(expr + 1)
    log_msg("判定原始尺度，已 log2(x+1)")
  }

  # 2) 分位数标准化（跨样本可比）
  expr <- limma::normalizeBetweenArrays(expr,
            method = cfg_get(cfg, "bulk.qc.normalize_method", "quantile"))

  # 3) 低表达过滤：至少一半样本高于 25% 分位
  min_frac <- cfg_get(cfg, "bulk.deg.min_samples_expr", 0.5)
  thr <- apply(expr, 2, stats::quantile, probs = cfg_get(cfg, "bulk.qc.min_expr_quantile", 0.25), na.rm = TRUE)
  keep <- rowSums(sweep(expr, 2, thr, ">"), na.rm = TRUE) >= min_frac * ncol(expr)
  expr <- expr[keep, , drop = FALSE]
  log_msg("过滤后保留 ", nrow(expr), " 个基因")

  # 4) 缺失值：少量用行均值补，基因缺失过多直接剔除
  na_rate <- rowMeans(is.na(expr))
  expr <- expr[na_rate <= 0.2, , drop = FALSE]
  if (anyNA(expr)) {
    for (j in seq_len(ncol(expr))) {
      idx <- is.na(expr[, j])
      if (any(idx)) expr[idx, j] <- rowMeans(expr, na.rm = TRUE)[idx]
    }
  }

  # 5) 表型
  pheno <- sm$pheno
  pheno$sample <- rownames(pheno)
  pheno$group <- make_group(pheno, group_regex)
  pheno$dataset <- gse
  pheno$batch <- if (!is.null(cfg_get(cfg, paste0("bulk.", role, ".batch")))) {
    as.character(pheno[[cfg_get(cfg, paste0("bulk.", role, ".batch"))]])
  } else NA_character_
  pheno <- cbind(pheno, extract_numeric_traits(pheno))

  # 样本顺序与矩阵列对齐
  expr <- expr[, pheno$sample, drop = FALSE]

  save_fig(pca_plot(expr, pheno, paste0(gse, " — before batch correction")),
           paste0("pca_", gse, "_raw"), cfg, 5.5, 4.5)

  # 6) 批次校正
  if (isTRUE(cfg_get(cfg, "bulk.qc.combat", TRUE)) && length(unique(stats::na.omit(pheno$batch))) > 1) {
    mod <- stats::model.matrix(~ group, data = pheno)
    expr <- sva::ComBat(dat = expr, batch = pheno$batch, mod = mod, par.prior = TRUE)
    log_msg("已完成 ComBat 批次校正（", paste(unique(pheno$batch), collapse = "/"), "）")
  } else {
    log_msg("跳过 ComBat：未提供批次列或批次唯一")
  }

  # 7) 离群样本标记（稳健距离，不自动删除）
  v <- apply(expr, 1, stats::var, na.rm = TRUE)
  top <- names(sort(v, decreasing = TRUE))[seq_len(min(2000, length(v)))]
  pc <- stats::prcomp(t(expr[top, , drop = FALSE]), scale. = TRUE)$x[, 1:2]
  ctr <- apply(pc, 2, stats::median)
  mad_d <- apply(pc, 2, stats::mad)
  dist_z <- sqrt(rowSums(sweep(pc, 2, ctr, "-")^2 / (mad_d^2 + 1e-8)))
  sd_cut <- cfg_get(cfg, "bulk.qc.outlier_sd", 3)
  pheno$pc1 <- pc[, 1]; pheno$pc2 <- pc[, 2]
  pheno$outlier <- dist_z > sd_cut
  if (any(pheno$outlier)) {
    log_msg("标记离群样本：", paste(pheno$sample[pheno$outlier], collapse = ", "),
            "（当前不删除，确认后把 auto_remove_outliers 设为 true 重跑）")
  }
  if (isTRUE(cfg_get(cfg, "bulk.qc.auto_remove_outliers", FALSE)) && any(pheno$outlier)) {
    keep_s <- !pheno$outlier
    expr <- expr[, keep_s, drop = FALSE]; pheno <- pheno[keep_s, , drop = FALSE]
    log_msg("已按配置删除离群样本，剩 ", ncol(expr), " 例")
  }

  save_fig(pca_plot(expr, pheno, paste0(gse, " — after QC")),
           paste0("pca_", gse, "_clean"), cfg, 5.5, 4.5)

  out_dir <- ensure_dir("data", "processed")
  saveRDS(expr, file.path(out_dir, paste0(gse, "_expr.rds")))
  saveRDS(pheno, file.path(out_dir, paste0(gse, "_pheno.rds")))

  data.frame(dataset = gse, role = role, platform = gpl,
             n_sample = ncol(expr), n_case = sum(pheno$group == "case"),
             n_control = sum(pheno$group == "control"),
             n_gene = nrow(expr), n_batch = length(unique(stats::na.omit(pheno$batch))),
             n_outlier = sum(pheno$outlier), stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------ 主流程 ----

if (demo) {
  log_msg("demo 模式：使用 99_demo_data.R 生成的合成数据，跳过 GEO 下载")
  qc_tab <- readRDS(proj_path("data", "processed", "demo_qc.rds"))
  write_tsv(qc_tab, "qc_summary", cfg)   # demo 也要落盘，否则 15 会把 01 记成缺产物
  # demo 也出一张 PCA：真实流程里这由 prepare_one 触发，demo 下补画是为了让
  # pca_plot 这条路径也被走到（Fig2A），同时暴露数据有没有被造歪
  de <- readRDS(proj_path("data", "processed", "DEMO_expr.rds"))
  dp <- readRDS(proj_path("data", "processed", "DEMO_pheno.rds"))
  save_fig(pca_plot(de, dp, "DEMO — synthetic data (before QC)"),
           "pca_DEMO_raw", cfg, 5.5, 4.5)
} else {
  jobs <- list(list(gse = cfg_get(cfg, "bulk.train.gse"), gpl = cfg_get(cfg, "bulk.train.gpl"),
                    rx = cfg_get(cfg, "bulk.train.group_regex"), role = "train"))
  for (v in cfg_get(cfg, "bulk.valid", list())) {
    jobs[[length(jobs) + 1]] <- list(gse = v$gse, gpl = v$gpl %||% NA,
                                     rx = v$group_regex %||% "glaucoma|POAG", role = "valid")
  }
  qc_tab <- do.call(rbind, lapply(jobs, function(j)
    prepare_one(j$gse, j$gpl, j$rx, cfg, j$role)))
  write_tsv(qc_tab, "qc_summary", cfg)
}

print(qc_tab)
log_msg("01 完成")

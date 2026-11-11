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

# 分组与选样只能看"描述样本身份"的列。把整行拼起来会把 !Sample_data_processing
# 里"glaucoma tissue was compared with..."这类全局描述也算进来，于是所有样本
# 都被判成病例——GSE4316 上一版就是这样跑出 8 病例 0 对照，而且不报错。
informative_text <- function(pheno) {
  cols <- grep("^characteristics|^title$|^source_name_ch1$|^description$",
               names(pheno), value = TRUE)
  if (!length(cols)) cols <- names(pheno)
  apply(pheno[, cols, drop = FALSE], 1, function(r) paste(r, collapse = " | "))
}

# 从表型表里按正则推分组：命中即病例，其余对照
make_group <- function(pheno, pattern) {
  txt <- informative_text(pheno)
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
  # 一个临床列都没匹配上时要返回"0 列但 36 行"的空表；直接
  # as.data.frame(list()) 得到的是 0 行，cbind(36 行, 0 行) 会报
  # "arguments imply differing number of rows"。GSE27276 只登记了
  # disease state / collection / 取材延迟，没有 IOP 与视野，必然走到这里。
  if (!length(out)) return(data.frame(row.names = rownames(pheno)))
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
  # batch 全 NA 时不能把 shape 映射上去：shape 取到 NA，ggplot 会把【所有点】
  # 丢掉，只留下 stat_ellipse 的空椭圆——图看着"画出来了"，其实一个样本都没有。
  # 训练集没有批次列时必然踩到，所以按有没有第二个批次水平决定映射。
  has_batch <- length(unique(stats::na.omit(df$batch))) > 1
  # 无批次时用 shape=21（空心圆带 fill）：默认的 19 号实心圆不吃 fill，
  # 两组会画成一模一样的黑点，颜色图例也失去意义
  pts <- if (has_batch) {
    ggplot2::geom_point(size = 2.6, colour = "grey20", stroke = 0.3)
  } else {
    ggplot2::geom_point(size = 2.6, colour = "grey20", stroke = 0.3, shape = 21)
  }
  p <- if (has_batch) {
    ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, fill = group, shape = batch))
  } else {
    ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, fill = group))
  }
  # 某组不足 3 例时 stat_ellipse 只会刷 "Too few points to calculate an ellipse"
  # 然后什么都不画，所以干脆不加这一层
  ell <- if (min(table(df$group)) >= 3) {
    ggplot2::stat_ellipse(ggplot2::aes(group = group), type = "norm",
                          linewidth = 0.4, show.legend = FALSE)
  }
  p + pts + ell +
    ggplot2::labs(title = title,
                  x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) +
    theme_pub()
}

# --------------------------------------------------------------- 单数据集 ----

prepare_one <- function(gse, gpl, group_regex, cfg, role = "train",
                        exclude_regex = NA_character_) {
  log_msg("---- 处理 ", gse, " (", role, ") ----")
  f <- proj_path("data", "raw", paste0(gse, "_series_matrix.txt.gz"))
  if (!file.exists(f)) {
    download_if_missing(
      paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", geo_series_dir(gse),
             "/", gse, "/matrix/", gse, "_series_matrix.txt.gz"), f)
  }
  sm <- read_series_matrix(f)
  sm_plat <- series_platform(f)
  if (is.null(gpl) || is.na(gpl) || gpl == "") gpl <- sm_plat
  # 平台号对不上必须以错误终止：探针 ID 是平台专属的，拿错注释表不会报错，
  # 只会把大部分探针映射不上（剩几百个基因）或映射到错误基因上，
  # 而下游每一步都"正常"跑完，最后得到一篇错的稿子。
  if (!is.na(sm_plat) && !identical(toupper(gpl), toupper(sm_plat))) {
    stop(sprintf("%s 的配置平台是 %s，但 series matrix 自己写的是 %s。%s",
                 gse, gpl, sm_plat,
                 "请把 config/config.yml 里的 gpl 改成实测值，或留空让脚本自动读取。"),
         call. = FALSE)
  }

  # 选样：exclude_regex 命中的整条丢掉。GSE4316 里有 3 例是原代培养细胞
  # （不是组织），GSE45570 里有 6 例是高眼压（OHT，不是正常对照）——
  # 这类样本按"其余即对照"的默认逻辑会被当成对照，把效应量稀释掉，
  # 而且正文里没法交代"对照里有 6 例高眼压"。
  if (!is.na(exclude_regex) && nzchar(exclude_regex)) {
    drop <- grepl(exclude_regex, informative_text(sm$pheno), ignore.case = TRUE)
    if (any(drop)) {
      log_msg("按 exclude_regex 剔除 ", sum(drop), " 例：",
              paste(sm$pheno$geo_accession[drop], collapse = ", "))
      sm$pheno <- sm$pheno[!drop, , drop = FALSE]
      sm$expr <- sm$expr[, !drop, drop = FALSE]
    } else {
      log_msg("提示：exclude_regex 没匹配到样本，请核对写法")
    }
  }

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
  # 样本太少时 MAD 完全不可靠：n=5 的 GSE4316 上一版标出 3 例"离群"，
  # 而 MAD 是拿这 5 个点自己估的，等于用自己判自己。少于 10 例直接跳过，
  # 免得把 n_outlier 这个数字写进正文还煞有介事地解释。
  v <- apply(expr, 1, stats::var, na.rm = TRUE)
  top <- names(sort(v, decreasing = TRUE))[seq_len(min(2000, length(v)))]
  pc <- stats::prcomp(t(expr[top, , drop = FALSE]), scale. = TRUE)$x[, 1:2]
  pheno$pc1 <- pc[, 1]; pheno$pc2 <- pc[, 2]
  pheno$outlier <- FALSE
  if (ncol(expr) < 10) {
    log_msg("样本仅 ", ncol(expr), " 例，跳过离群判定（MAD 在这个量级上没有意义）")
  } else {
    ctr <- apply(pc, 2, stats::median)
    mad_d <- apply(pc, 2, stats::mad)
    dist_z <- sqrt(rowSums(sweep(pc, 2, ctr, "-")^2 / (mad_d^2 + 1e-8)))
    sd_cut <- cfg_get(cfg, "bulk.qc.outlier_sd", 3)
    pheno$outlier <- dist_z > sd_cut
    if (any(pheno$outlier)) {
      log_msg("标记离群样本：", paste(pheno$sample[pheno$outlier], collapse = ", "),
              "（当前不删除，确认后把 auto_remove_outliers 设为 true 重跑）")
    }
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
                    rx = cfg_get(cfg, "bulk.train.group_regex"), role = "train",
                    ex = cfg_get(cfg, "bulk.train.exclude_regex", NA)))
  for (v in cfg_get(cfg, "bulk.valid", list())) {
    jobs[[length(jobs) + 1]] <- list(gse = v$gse, gpl = v$gpl %||% NA,
                                     rx = v$group_regex %||% "glaucoma|POAG", role = "valid",
                                     ex = v$exclude_regex %||% NA_character_)
  }
  qc_tab <- do.call(rbind, lapply(jobs, function(j)
    prepare_one(j$gse, j$gpl, j$rx, cfg, j$role, j$ex)))
  write_tsv(qc_tab, "qc_summary", cfg)
}

print(qc_tab)
log_msg("01 完成")

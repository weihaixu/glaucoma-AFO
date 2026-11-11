# =============================================================================
# 10_mr.R —— 课题二（4）：孟德尔随机化（表达/pQTL → POAG）
#
# 输入：config$genetics$mr$exposure（SNP, gene, beta, se, ea, nea, p，可选 tissue）
#       data/processed/gwas_clean.rds（结局）
# 输出：results/tables/mr_instruments.tsv      工具变量与 F 统计量
#       results/tables/mr_results.tsv          各方法效应估计
#       results/tables/mr_sensitivity.tsv      多效性 / 异质性 / 留一法
#       results/figures/mr_scatter_<gene>.pdf / mr_forest.pdf
#
# 必须报告的六件事（少一件审稿人就会问）：
#   1. 工具变量强度：F > 10（弱工具变量会让 IVW 有偏）
#   2. 样本重叠：暴露与结局是否同一批人 → 用 MRlap 或直接说明
#   3. 水平多效性：MR-Egger 截距 + MR-PRESSO 全局检验
#   4. 异质性：Cochran's Q
#   5. 方向性：Steiger 检验（暴露→结局 而非反向）
#   6. 多重检验：多少个基因 × 多少种方法 → Bonferroni/FDR 阈值是多少
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("10_mr", cfg)
load_pkgs(c("ggplot2", "TwoSampleMR", "data.table"))
set.seed(cfg_get(cfg, "project.seed", 1))

mcfg <- cfg_get(cfg, "genetics.mr", list())
outcome_all <- readRDS(proj_path("data", "processed", "gwas_clean.rds"))

# ------------------------------------------------- 1) 读取暴露数据 ----------
exp_path <- proj_path(mcfg$exposure %||% "data/raw/eqtl_instruments.tsv.gz")
require_file(exp_path, "请提供暴露（eQTL/pQTL）工具变量表，或在 config 改 genetics.mr.exposure")
exp_raw <- data.table::fread(exp_path, showProgress = FALSE)
need <- c("SNP", "beta", "se", "ea", "nea", "p")
if (!all(need %in% names(exp_raw))) {
  stop("暴露表需包含列：", paste(need, collapse = ", "), "；可选 gene / tissue", call. = FALSE)
}
if (!"gene" %in% names(exp_raw)) exp_raw$gene <- "exposure"
if (!"tissue" %in% names(exp_raw)) exp_raw$tissue <- "NA"
log_msg("暴露表：", nrow(exp_raw), " 行，", length(unique(exp_raw$gene)), " 个基因")

# ------------------------------------------------- 2) 工具变量筛选 --------
p_thr <- mcfg$p_threshold %||% 5e-8
sig <- exp_raw[exp_raw$p < p_thr, ]
if (!nrow(sig)) {
  log_msg("P<", p_thr, " 无可用工具变量，放宽到 5e-6 并在文中说明")
  sig <- exp_raw[exp_raw$p < 5e-6, ]
}
if (!nrow(sig)) stop("暴露数据中没有任何可用工具变量")

# LD 剪枝：优先用本地 plink + 1000G 面板；没有就按距离粗剪（必须在文中说明退化）
clump_local <- function(df, r2 = 0.001, kb = 10000) {
  bfile <- proj_path(cfg_get(cfg, "genetics.ld.panel_plink", "data/raw/1000G_EUR_Phase3_plink/1000G.EUR.QC"))
  has_plink <- file.exists(paste0(bfile, ".bed")) && nzchar(Sys.which("plink"))
  if (has_plink && requireNamespace("ieugwasr", quietly = TRUE)) {
    out <- tryCatch(ieugwasr::ld_clump(
      dat = data.frame(rsid = df$SNP, pval = df$p, id = df$gene),
      clump_r2 = r2, clump_kb = kb, bfile = bfile, plink_bin = Sys.which("plink")),
      error = function(e) NULL)
    if (!is.null(out)) return(df[df$SNP %in% out$rsid, ])
  }
  log_msg("无 plink/LD 面板，改用距离剪枝（同基因内 ", kb / 1000, "Mb 只留最显著 SNP）")
  res <- lapply(split(df, df$gene), function(d) {
    d <- d[order(d$p), ]
    keep <- integer(0)
    for (i in seq_len(nrow(d))) {
      if (!length(keep)) { keep <- i; next }
      same_region <- !is.na(d$pos) & !is.na(d$pos[keep]) &
                     abs(d$pos[i] - d$pos[keep]) < kb * 1000
      if (!any(same_region)) keep <- c(keep, i)
    }
    d[keep, ]
  })
  do.call(rbind, res)
}
if (!"pos" %in% names(sig)) sig$pos <- NA_real_
inst <- clump_local(sig, mcfg$clump_r2 %||% 0.001, mcfg$clump_kb %||% 10000)
inst$F_stat <- (inst$beta / inst$se)^2
inst <- inst[inst$F_stat >= (mcfg$f_stat_min %||% 10), ]
if (!nrow(inst)) stop("所有工具变量 F<10，属于弱工具变量，不适合做 MR")
write_tsv(inst, "mr_instruments", cfg)
log_msg("保留工具变量 ", nrow(inst), " 个（", length(unique(inst$gene)), " 个基因），F 中位数 ",
        round(stats::median(inst$F_stat), 1))

# ------------------------------------------------- 3) 逐基因 MR ------------
# gwas_clean 只是分析子集（显著位点窗口 + P<1e-4），工具变量可能不在里面。
# 缺哪个就从全量 LDSC 文件里回捞，别让 harmonise 静默少行。
miss <- setdiff(inst$SNP, outcome_all$SNP)
if (length(miss)) {
  log_msg("分析子集缺 ", length(miss), " 个工具变量，从全量 sumstats 回捞")
  extra <- subset_by_snp(proj_path("data", "processed", "gwas_ldsc.sumstats.gz"), miss,
                         cols = c("SNP", "A1", "A2", "BETA", "SE", "P", "N"),
                         types = c("character", "character", "character",
                                   "numeric", "numeric", "numeric", "numeric"))
  if (nrow(extra)) {
    extra <- data.frame(SNP = extra$SNP, A1 = extra$A1, A2 = extra$A2,
                        beta = extra$BETA, se = extra$SE, p = extra$P,
                        EAF = NA_real_, N = extra$N, INFO = NA_real_,
                        CHR = NA_integer_, POS = NA_real_, pos_key = NA_character_,
                        stringsAsFactors = FALSE)
    outcome_all <- rbind(outcome_all, extra)
    log_msg("  回捞到 ", nrow(extra), " 个（双等位基因型已知，EAF 缺失）")
  }
}

# 结局列名对齐 TwoSampleMR 要求
out_df <- data.frame(SNP = outcome_all$SNP, beta = outcome_all$beta, se = outcome_all$se,
                     ea = outcome_all$A1, nea = outcome_all$A2,
                     eaf = outcome_all$EAF, pval = outcome_all$p,
                     samplesize = outcome_all$N, stringsAsFactors = FALSE)
out_fmt <- TwoSampleMR::format_data(out_df, type = "outcome",
                                    snp_col = "SNP", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "ea", other_allele_col = "nea",
                                    eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize")

all_res <- list(); all_sens <- list()
for (gn in unique(inst$gene)) {
  e <- inst[inst$gene == gn, ]
  if (nrow(e) < 1) next
  exp_fmt <- TwoSampleMR::format_data(
    data.frame(SNP = e$SNP, beta = e$beta, se = e$se, ea = e$ea, nea = e$nea,
               eaf = if ("eaf" %in% names(e)) e$eaf else NA_real_,
               pval = e$p, stringsAsFactors = FALSE),
    type = "exposure", snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "ea", other_allele_col = "nea", eaf_col = "eaf", pval_col = "pval")

  dat <- TwoSampleMR::harmonise_data(exp_fmt, out_fmt, action = 2)
  dat <- dat[dat$mr_keep, ]
  if (nrow(dat) < 1) { log_msg(gn, "：harmonise 后无可用 SNP，跳过"); next }

  # 单个工具变量只能用 Wald 比（eQTL 的常态：一个基因位点往往只有一个独立信号）
  methods <- if (nrow(dat) == 1) {
    "mr_wald_ratio"
  } else {
    mcfg$methods %||% c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
  }
  res <- TwoSampleMR::mr(dat, method_list = methods)
  res$gene <- gn
  res$tissue <- e$tissue[1]
  res$n_snp <- nrow(dat)
  all_res[[gn]] <- res

  # 敏感性分析
  plei <- if (nrow(dat) >= 3) tryCatch(TwoSampleMR::mr_pleiotropy_test(dat), error = function(e) NULL) else NULL
  het  <- if (nrow(dat) >= 2) tryCatch(TwoSampleMR::mr_heterogeneity(dat), error = function(e) NULL) else NULL
  loo  <- tryCatch(TwoSampleMR::mr_leaveoneout(dat), error = function(e) NULL)
  sens <- data.frame(gene = gn, tissue = e$tissue[1],
                     egger_intercept = if (!is.null(plei)) plei$egger_intercept else NA_real_,
                     egger_intercept_p = if (!is.null(plei)) plei$pval else NA_real_,
                     Q = if (!is.null(het)) het$Q[1] else NA_real_,
                     Q_p = if (!is.null(het)) het$Q_pval[1] else NA_real_,
                     n_snp = nrow(dat), stringsAsFactors = FALSE)
  all_sens[[gn]] <- sens

  if (!is.null(loo)) write_tsv(as.data.frame(loo), paste0("mr_leaveoneout_", gn), cfg)

  # MR-PRESSO：离群点检测，需 GitHub 包 MRPRESSO
  if (requireNamespace("MRPRESSO", quietly = TRUE) && nrow(dat) >= 4) {
    pr <- tryCatch(MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                                       SdOutcome = "se.outcome", SdExposure = "se.exposure",
                                       OUTCOME = "outcome", data = dat,
                                       NbDistribution = mcfg$nboot_presso %||% 1000,
                                       SignifThreshold = 0.05), error = function(e) NULL)
    if (!is.null(pr)) {
      g <- pr$`Main MR results`$`P-value`[1]
      all_sens[[gn]]$presso_global_p <- g
    }
  }

  # 散点图（每个基因一张，基因多时只画显著的前若干个）
  g_sc <- TwoSampleMR::mr_scatter_plot(res, dat)[[1]] +
    ggplot2::labs(title = paste0(gn, " (", e$tissue[1], ")")) + theme_pub()
  save_fig(g_sc, paste0("mr_scatter_", gn), cfg, 5, 4.4)
}

if (!length(all_res)) stop("没有任何基因完成 MR，请检查暴露与结局的 SNP 交集")

mr_tab <- do.call(rbind, all_res)
mr_tab <- mr_tab[order(mr_tab$pval), ]

sens_tab <- do.call(rbind, all_sens)

# ------------------------------------------------- 4) 多重检验与汇总 ------
# 多重检验只在「主分析方法（IVW）× 基因数」这一层校正，不能把不同方法的行混在一起校正
ivw <- mr_tab[mr_tab$method == "Inverse variance weighted", ]
ivw$fdr <- stats::p.adjust(ivw$pval, method = "BH")
mr_tab$fdr <- NA_real_
mr_tab$fdr[match(paste(ivw$gene, ivw$method), paste(mr_tab$gene, mr_tab$method))] <- ivw$fdr

# 写表必须放在 fdr 赋值之后：之前先写后算，导致 mr_results.tsv 里根本没有 fdr 列，
# 下游 {{MR_TOP_FDR}} 一直取不到值
write_tsv(mr_tab, "mr_results", cfg)
write_tsv(sens_tab, "mr_sensitivity", cfg)

if (nrow(ivw)) {
  ivw <- ivw[order(ivw$pval), ]
  log_msg("IVW 名义显著基因：", sum(ivw$pval < 0.05), "；FDR<0.05：", sum(ivw$fdr < 0.05))
  print(head(ivw[, c("gene", "nsnp", "b", "se", "pval", "fdr")], 10))
}

# 森林图（各基因 IVW 效应）
if (nrow(ivw) >= 2) {
  ivw$gene <- factor(ivw$gene, levels = rev(ivw$gene[order(ivw$b)]))
  g_forest <- ggplot2::ggplot(ivw, ggplot2::aes(b, gene)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                            height = 0.2, colour = "grey40") +
    ggplot2::geom_point(ggplot2::aes(colour = pval < 0.05), size = 2.2) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "grey55"), guide = "none") +
    ggplot2::labs(x = "IVW causal estimate (log OR per SD expression)", y = NULL,
                  title = "MR: expression -> POAG") +
    theme_pub()
  save_fig(g_forest, "mr_forest", cfg, 6.5, max(3, 0.28 * nrow(ivw) + 1.5))
}

# Steiger 方向性（需要暴露/结局的样本量方差信息，缺失时跳过）
log_msg("提示：Steiger 检验需 directionality_test()，要求暴露与结局的 N 与 SD 都已知；",
        "若暴露为 eQTL（N 已知但 SD 缺失），请在文中说明改用生物学先验判断方向。")
log_msg("提示：若暴露与结局来自同一队列（如 UKB eQTL × UKB POAG），",
        "请用 MRlap 校正样本重叠，或改用 FinnGen 等独立结局数据。")

saveRDS(list(results = mr_tab, sensitivity = sens_tab, instruments = inst),
        proj_path("data", "processed", "mr_all.rds"))
log_msg("10 完成")

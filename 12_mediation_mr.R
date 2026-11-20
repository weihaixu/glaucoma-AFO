# =============================================================================
# 12_mediation_mr.R —— 课题二（扩展）：两阶段 MR 与中介分析
#
# 要回答的问题：候选基因是通过"升高眼压"致病，还是独立于眼压起作用？
# 这是无湿实验方案里最接近"机制必要性"的证据——用遗传工具把通路拆成两段。
#
# 设计：
#   总效应   ：基因表达 --IVW--> POAG            （10_mr.R 已算）
#   第一段   ：基因表达 --IVW--> IOP             （本脚本，UKB IOP GWAS）
#   第二段   ：IOP      --IVW--> POAG            （本脚本，IOP 独立工具变量）
#   间接效应 ：b1 × b2，SE 用 delta 法
#   直接效应 ：多变量 MR（MVMR），同时放表达与 IOP
#   中介比例 ：indirect / total
#
# 输入：data/processed/mr_all.rds（10 的输出）
#       config$drylab$mediation$iop_gwas（眼压 GWAS 汇总统计）
#       data/processed/gwas_clean.rds（POAG 结局）
# 输出：results/tables/mr_mediation.tsv、mr_mediation_iop_stage2.tsv
#       results/figures/mr_mediation_forest.pdf
#
# 必须交代的局限：
#   1. 两阶段 MR 要求两段样本独立或近似独立，否则乘积法有偏
#   2. IOP 工具变量若与表达工具变量重叠，需用 MVMR 的条件 F 检验弱工具
#   3. 中介比例 >1 或 <0 时不要硬解释，报告为"方向不一致"
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("12_mediation_mr", cfg)
load_pkgs(c("ggplot2", "TwoSampleMR", "data.table"))
set.seed(cfg_get(cfg, "project.seed", 1))

mcfg <- cfg_get(cfg, "genetics.mr", list())
dcfg <- cfg_get(cfg, "drylab.mediation", list())

mr_rds <- proj_path("data", "processed", "mr_all.rds")
require_file(mr_rds, "先跑 10_mr.R 得到 mr_all.rds（含工具变量与总效应）")
mr <- readRDS(mr_rds)
inst <- mr$instruments
tot <- mr$results
# 10_mr.R 的每个工具变量只有一个 SNP，结果全是 Wald ratio，没有 IVW。
# 没有 IVW 时回退到 Wald ratio 作为总效应（每个基因保留一条）。
keep <- grepl("Inverse variance|IVW", tot$method, ignore.case = TRUE)
if (!any(keep)) keep <- rep(TRUE, nrow(tot))
tot <- tot[keep, c("gene", "b", "se", "pval")]
tot <- tot[!duplicated(tot$gene), , drop = FALSE]
names(tot) <- c("gene", "beta_total", "se_total", "p_total")
if (!nrow(tot)) stop("mr_results 里没有可用的总效应结果", call. = FALSE)
log_msg("总效应来源：", ifelse(any(grepl("Inverse variance|IVW", mr$results$method, ignore.case = TRUE)),
                          "IVW", "单 SNP Wald ratio（无 IVW 时的回退）"))

# ------------------------------------------------- 1) 读 IOP GWAS ----------
iop_path <- proj_path(dcfg$iop_gwas %||% "data/raw/ukb_iop_gwas.tsv.gz")
require_file(iop_path, paste0(
  "需要眼压（IOP）GWAS 汇总统计。可从 UK Biobank（Neale lab / SAIGE）、",
  "FinnGen（配镜/眼压相关端点）或已发表 IOP GWAS 补充材料获取，",
  "放到该路径并在 config$drylab$mediation$cols 里映射列名"))

cmap <- dcfg$cols %||% list()
# 完整 IOP GWAS 压缩后 ~700 MB、解压数 GB；只读需要的列以免内存爆掉
need_cols <- unique(c(cmap$snp %||% "SNP", cmap$beta %||% "Effect",
                      cmap$se %||% "StdErr", cmap$pval %||% "P-value",
                      cmap$ea %||% "EffectAllele", cmap$nea %||% "NonEffectAllele",
                      "POS"))
iop <- data.table::fread(iop_path, showProgress = FALSE, select = need_cols)
pick <- function(key, default) {
  nm <- cmap[[key]] %||% default
  if (!nm %in% names(iop)) {
    stop("IOP GWAS 中找不到列 '", nm, "'；请在 config$drylab$mediation$cols 里映射",
         call. = FALSE)
  }
  iop[[nm]]
}
iop_df <- data.frame(
  SNP = as.character(pick("snp", "SNP")),
  beta = as.numeric(pick("beta", "BETA")),
  se = as.numeric(pick("se", "SE")),
  p = as.numeric(pick("pval", "P")),
  ea = as.character(pick("ea", "A1")),
  nea = as.character(pick("nea", "A2")),
  eaf = tryCatch(as.numeric(pick("eaf", "FRQ")), error = function(e) NA_real_),
  pos = tryCatch(as.numeric(pick("pos", "POS")), error = function(e) NA_real_),
  stringsAsFactors = FALSE)
iop_df <- iop_df[is.finite(iop_df$beta) & is.finite(iop_df$se) &
                 iop_df$se > 0 & !is.na(iop_df$ea) & !is.na(iop_df$nea), ]
log_msg("IOP GWAS：", nrow(iop_df), " 个 SNP")
# 只保留后续真正用到的 SNP（显著位点 + 表达工具变量），避免千万行进 format_data
p_thr0 <- dcfg$p_threshold %||% 5e-8
iop_keep <- iop_df[iop_df$p < p_thr0 | iop_df$SNP %in% inst$SNP, ]
log_msg("IOP 精简：", nrow(iop_keep), " 个 SNP（显著 + 工具变量）")
rm(iop_df); gc()

outcome_all <- readRDS(proj_path("data", "processed", "gwas_clean.rds"))
poag_df <- data.frame(SNP = outcome_all$SNP, beta = outcome_all$beta, se = outcome_all$se,
                      p = outcome_all$p, ea = outcome_all$A1, nea = outcome_all$A2,
                      eaf = outcome_all$EAF, stringsAsFactors = FALSE)

# gwas_clean 是分析子集（显著位点窗口 + P<1e-4），IOP 工具变量未必在里面；
# 缺的从全量 sumstats 回捞，否则第二步（IOP→POAG）的 harmonise 会静默少 SNP。
miss <- setdiff(iop_keep$SNP, poag_df$SNP)
if (length(miss)) {
  log_msg("检索 IOP 工具变量的 POAG 结局：分析子集缺 ", length(miss), " 个，回捞中")
  extra <- subset_by_snp(proj_path("data", "processed", "gwas_ldsc.sumstats.gz"), miss,
                         cols = c("SNP", "A1", "A2", "BETA", "SE", "P"),
                         types = c("character", "character", "character",
                                   "numeric", "numeric", "numeric"))
  if (nrow(extra)) {
    poag_df <- rbind(poag_df, data.frame(SNP = extra$SNP, beta = extra$BETA, se = extra$SE,
                                         p = extra$P, ea = extra$A1, nea = extra$A2,
                                         eaf = NA_real_, stringsAsFactors = FALSE))
    log_msg("  回捞到 ", nrow(extra), " 个")
  }
}

mk <- function(df, type) {
  # TwoSampleMR::format_data 拒绝 data.table 类；统一转成纯 data.frame
  df <- as.data.frame(df)
  a <- list(df, type = type, snp_col = "SNP", beta_col = "beta", se_col = "se",
            effect_allele_col = "ea", other_allele_col = "nea", pval_col = "p")
  if ("eaf" %in% names(df) && any(is.finite(df$eaf))) a$eaf_col <- "eaf"
  do.call(TwoSampleMR::format_data, a)
}
# 大表格只格式化一次，避免在 101 个基因的循环里反复格式化千万行 IOP GWAS
iop_out_fmt  <- mk(iop_keep, "outcome")
poag_out_fmt <- mk(poag_df, "outcome")
harm <- function(exp_df, out_fmt, label = "") {
  e <- mk(exp_df, "exposure")
  # action=2 会因缺 EAF 丢弃回文 SNP（本 IOP GWAS 的 MAF 与效应等位基因未对齐，
  # config 有意不映射 eaf）；这里等位基因已明确对齐，用 action=1 保留这些 SNP。
  d <- TwoSampleMR::harmonise_data(e, out_fmt, action = 1)
  d <- d[d$mr_keep, ]
  if (nrow(d)) d$pair <- label
  d
}

wald <- function(b1, se1) {
  list(b = b1, se = se1, p = 2 * stats::pnorm(-abs(b1 / se1)))
}

# ------------------------------------------------- 2) 第二段：IOP → POAG ---
p_thr <- dcfg$p_threshold %||% 5e-8
iop_sig <- iop_keep[iop_keep$p < p_thr, ]
if (!nrow(iop_sig)) {
  log_msg("IOP GWAS 在 P<", p_thr, " 无 SNP，放宽到 5e-6")
  iop_sig <- iop_keep[iop_keep$p < 5e-6, ]
}
if (!nrow(iop_sig)) stop("IOP GWAS 没有可用工具变量", call. = FALSE)

# 有位置列则按真实距离做 1 Mb 粗剪；无位置则取最显著 30 个
if (any(is.finite(iop_sig$pos))) {
  iop_sig <- iop_sig[order(iop_sig$p), ]
  keep <- integer(0)
  for (i in seq_len(nrow(iop_sig))) {
    if (!length(keep)) { keep <- i; next }
    if (all(abs(iop_sig$pos[i] - iop_sig$pos[keep]) > 1e6)) keep <- c(keep, i)
  }
  iop_sig <- iop_sig[keep, ]
  log_msg("IOP 工具变量按 1 Mb 距离粗剪后：", nrow(iop_sig), " 个（建议用 plink 重做 LD 剪枝）")
} else {
  iop_sig <- iop_sig[order(iop_sig$p), ][seq_len(min(30, nrow(iop_sig))), ]
  log_msg("IOP GWAS 无位置列，仅取最显著的 ", nrow(iop_sig), " 个 SNP（弱于 LD 剪枝，需在文中说明）")
}

d2 <- harm(iop_sig[, c("SNP", "beta", "se", "p", "ea", "nea", "eaf")],
           poag_out_fmt, "iop->poag")
if (!nrow(d2)) stop("IOP 工具变量与 POAG 无可用交集", call. = FALSE)
stage2 <- TwoSampleMR::mr(d2, method_list = c("mr_ivw", "mr_egger_regression",
                                              "mr_weighted_median"))
stage2$segment <- "IOP -> POAG"
write_tsv(stage2, "mr_mediation_iop_stage2", cfg)
ivw2 <- stage2[stage2$method == "Inverse variance weighted", ]
if (!nrow(ivw2)) stop("第二段 IVW 失败", call. = FALSE)
b2 <- ivw2$b[1]; se2 <- ivw2$se[1]
log_msg("第二段 IOP→POAG：b=", round(b2, 3), "，SE=", round(se2, 3),
        "，P=", signif(ivw2$pval[1], 3), "（", nrow(d2), " 个 SNP）")

# ------------------------------------------------- 3) 第一段：表达 → IOP ---
res <- list()
for (gn in unique(inst$gene)) {
  e <- inst[inst$gene == gn, ]
  if (!nrow(e)) next
  d1 <- harm(e[, c("SNP", "beta", "se", "p", "ea", "nea", "eaf")], iop_out_fmt, "expr->iop")
  if (!nrow(d1)) { log_msg(gn, "：与 IOP GWAS 无交集，跳过"); next }
  # 单 SNP 时 mr_ivw 返回空：直接用 Wald ratio（b_out/b_exp）
  if (nrow(d1) >= 2) {
    m1 <- TwoSampleMR::mr(d1, method_list = "mr_ivw")
    b1 <- if (nrow(m1)) m1$b[1] else NA_real_
    se1 <- if (nrow(m1)) m1$se[1] else NA_real_
    p1 <- if (nrow(m1)) m1$pval[1] else NA_real_
  } else {
    wr <- wald(d1$beta.outcome[1] / d1$beta.exposure[1],
               d1$se.outcome[1] / abs(d1$beta.exposure[1]))
    b1 <- wr$b; se1 <- wr$se; p1 <- wr$p
  }

  # 总效应
  tt <- tot[tot$gene == gn, ]
  if (!nrow(tt)) next

  # 间接效应（乘积法 + delta 法 SE）
  ind <- b1 * b2
  se_ind <- sqrt(b1^2 * se2^2 + b2^2 * se1^2)
  p_ind <- 2 * stats::pnorm(-abs(ind / se_ind))
  prop <- if (is.finite(tt$beta_total[1]) && tt$beta_total[1] != 0) ind / tt$beta_total[1] else NA_real_

  # MVMR：同时放表达与 IOP，估计调整后的直接效应
  d_iop <- harm(e[, c("SNP", "beta", "se", "p", "ea", "nea", "eaf")], iop_out_fmt, "mv_iop")
  d_out <- harm(e[, c("SNP", "beta", "se", "p", "ea", "nea", "eaf")], poag_out_fmt, "mv_out")
  direct <- NA_real_; se_direct <- NA_real_; p_direct <- NA_real_; n_mv <- 0
  if (nrow(d_iop) && nrow(d_out)) {
    mv <- merge(
      d_out[, c("SNP", "beta.exposure", "beta.outcome", "se.outcome")],
      d_iop[, c("SNP", "beta.outcome")],
      by = "SNP", suffixes = c("", ".iop"))
    names(mv)[names(mv) == "beta.outcome.iop"] <- "beta.iop"
    mv <- mv[is.finite(mv$beta.exposure) & is.finite(mv$beta.iop) &
             is.finite(mv$beta.outcome) & mv$se.outcome > 0, ]
    n_mv <- nrow(mv)
    if (n_mv >= 3) {
      fit <- tryCatch(stats::lm(beta.outcome ~ beta.exposure + beta.iop,
                                data = mv, weights = 1 / se.outcome^2),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        cf <- summary(fit)$coefficients
        if ("beta.exposure" %in% rownames(cf)) {
          direct <- cf["beta.exposure", "Estimate"]
          se_direct <- cf["beta.exposure", "Std. Error"]
          p_direct <- cf["beta.exposure", "Pr(>|t|)"]
        }
      }
    }
  }

  res[[gn]] <- data.frame(
    gene = gn, tissue = e$tissue[1], n_snp_expr = nrow(d1), n_snp_mvmr = n_mv,
    beta_total = tt$beta_total[1], se_total = tt$se_total[1], p_total = tt$p_total[1],
    beta_expr_iop = b1, se_expr_iop = se1, p_expr_iop = p1,
    beta_iop_poag = b2, se_iop_poag = se2,
    beta_indirect = ind, se_indirect = se_ind, p_indirect = p_ind,
    prop_mediated = prop,
    beta_direct_mvmr = direct, se_direct_mvmr = se_direct, p_direct_mvmr = p_direct,
    stringsAsFactors = FALSE)
  log_msg(gn, "：expr→IOP b=", round(b1, 4), "（P=", signif(p1, 3), "）；中介比例 ",
          ifelse(is.finite(prop), round(prop, 2), "NA"))
}

if (!length(res)) stop("没有任何基因完成中介分析", call. = FALSE)
med <- do.call(rbind, res)
med <- med[order(med$p_indirect), ]
med$fdr_indirect <- stats::p.adjust(med$p_indirect, method = "BH")
write_tsv(med, "mr_mediation", cfg)

# ------------------------------------------------- 4) 图 -------------------
plot_d <- med[is.finite(med$prop_mediated) & med$p_total < 0.05, ]
if (nrow(plot_d) >= 2) {
  plot_d$gene <- factor(plot_d$gene, levels = rev(plot_d$gene[order(plot_d$prop_mediated)]))
  g <- ggplot2::ggplot(plot_d, ggplot2::aes(prop_mediated, gene)) +
    ggplot2::geom_vline(xintercept = c(0, 1), linetype = "dashed",
                        colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = p_indirect < 0.05), size = 2.4) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "grey55"),
                                 guide = "none") +
    ggplot2::labs(x = "proportion mediated by IOP", y = NULL,
                  title = "Two-stage MR: expression -> IOP -> POAG") +
    theme_pub()
  save_fig(g, "mr_mediation_forest", cfg, 6, max(2.6, 0.3 * nrow(plot_d) + 1.4))
}

log_msg("判读：中介比例 0.2–0.8 表示部分经 IOP 通路；>1 或 <0 说明两段方向不一致，不要硬解释；")
log_msg("      p_direct_mvmr > 0.05 且 p_indirect < 0.05 → 效应基本由 IOP 介导。")
log_msg("12 完成")

# =============================================================================
# 07_gwas_qc.R —— 课题二（1）：GWAS 汇总统计清洗、质控、LDSC 格式化
#
# 输入：config$genetics$gwas$path（FinnGen / UKB / MVP 的 sumstats，.gz 可直接读）
# 输出：data/processed/gwas_clean.rds     清洗后的数据（SNP/A1/A2/beta/se/p/N/eaf）
#       results/tables/gwas_qc_summary.tsv 各步过滤计数（写方法学时要抄这张表）
#       results/figures/gwas_qq.pdf / gwas_manhattan.pdf
#       data/processed/gwas_ldsc.sumstats.gz（可直接喂 ldsc.py --h2）
#
# 质控顺序不能乱：先统一列名与等位基因，再过滤，再算 lambda。
# 过滤掉一半 SNP 后算出来的 lambda 会偏低，别拿它当"质控良好"的证据。
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("07_gwas_qc", cfg)
load_pkgs(c("ggplot2", "data.table"))
set.seed(cfg_get(cfg, "project.seed", 1))

gcfg <- cfg_get(cfg, "genetics.gwas", list())
f_in <- proj_path(gcfg$path %||% "data/raw/gwas.txt.gz")
require_file(f_in, "请把 GWAS 汇总统计放到该路径，或在 config.yml 改 genetics.gwas.path")

# ------------------------------------------------- 1) 读取 -------------------
log_msg("读取 GWAS：", f_in)
raw <- data.table::fread(f_in, sep = gcfg$sep %||% "\t", showProgress = FALSE)
log_msg("原始 SNP 数：", nrow(raw), "，列名：", paste(names(raw), collapse = ", "))

# 列名映射（不同数据源命名差异极大，统一在这里改配置而不是改代码）
cmap <- gcfg$cols %||% list()
need <- c("snp", "beta", "se", "pval", "ea", "nea")
missing <- need[!vapply(need, function(k) !is.null(cmap[[k]]) && cmap[[k]] %in% names(raw), logical(1))]
if (length(missing)) {
  stop("以下字段未在数据中找到：", paste(missing, collapse = ", "),
       "\n  可用列：", paste(names(raw), collapse = ", "),
       "\n  请修改 config.yml 的 genetics.gwas.cols", call. = FALSE)
}

d <- data.frame(
  SNP  = as.character(raw[[cmap$snp]]),
  A1   = toupper(as.character(raw[[cmap$ea]])),   # effect allele
  A2   = toupper(as.character(raw[[cmap$nea]])),  # other allele
  beta = as.numeric(raw[[cmap$beta]]),
  se   = as.numeric(raw[[cmap$se]]),
  p    = as.numeric(raw[[cmap$pval]]),
  stringsAsFactors = FALSE
)
d$EAF <- if (!is.null(cmap$eaf) && cmap$eaf %in% names(raw)) as.numeric(raw[[cmap$eaf]]) else NA_real_
d$N   <- if (!is.null(cmap$n) && cmap$n %in% names(raw)) as.numeric(raw[[cmap$n]]) else NA_real_
if (!is.null(cmap$info) && cmap$info %in% names(raw)) d$INFO <- as.numeric(raw[[cmap$info]]) else d$INFO <- NA_real_

qc <- data.frame(step = "raw", n = nrow(d))

# 位置信息：后续 coloc / LD clumping / TWAS 都要用，这一层不带上后面全得重跑
chr_col <- grep("^(chr|CHR|chromosome|#chrom)$", names(raw), value = TRUE)[1]
pos_col <- grep("^(pos|POS|position|base_pair_location|bp)$", names(raw), value = TRUE)[1]
if (!is.na(chr_col) && !is.na(pos_col)) {
  pos_df <- data.frame(SNP = as.character(raw[[cmap$snp]]),
                       CHR = suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(raw[[chr_col]])))),
                       POS = suppressWarnings(as.numeric(raw[[pos_col]])),
                       stringsAsFactors = FALSE)
  pos_df <- pos_df[!duplicated(pos_df$SNP), ]
  d$CHR <- pos_df$CHR[match(d$SNP, pos_df$SNP)]
  d$POS <- pos_df$POS[match(d$SNP, pos_df$SNP)]
  d$pos_key <- paste0(d$CHR, "_", d$POS)
  log_msg("已挂载染色体位置：", sum(!is.na(d$CHR)), " 个 SNP")
} else {
  d$CHR <- NA_integer_; d$POS <- NA_real_; d$pos_key <- NA_character_
  log_msg("GWAS 数据无染色体/位置列：coloc 将退化为 rsID 匹配，建议先补位置")
}

# ------------------------------------------------- 2) 逐步过滤 --------------
keep <- !is.na(d$SNP) & !is.na(d$beta) & !is.na(d$se) & !is.na(d$p)
d <- d[keep, ]; qc <- rbind(qc, data.frame(step = "drop_NA_core", n = nrow(d)))

d <- d[d$se > 0 & is.finite(d$beta) & is.finite(d$se), ]
qc <- rbind(qc, data.frame(step = "drop_invalid_se", n = nrow(d)))

# 等位基因合法性：只保留 A/C/G/T 单碱基（indel 与多等位在 coloc/MR 中会破坏对齐）
d <- d[nchar(d$A1) == 1 & nchar(d$A2) == 1 &
         d$A1 %in% c("A", "C", "G", "T") & d$A2 %in% c("A", "C", "G", "T") & d$A1 != d$A2, ]
qc <- rbind(qc, data.frame(step = "keep_biallelic_SNV", n = nrow(d)))

# 回文 SNP（A/T、C/G）：MAF 接近 0.5 时链方向无法判断，MR 与 coloc 都会受污染
pal <- (d$A1 == "A" & d$A2 == "T") | (d$A1 == "T" & d$A2 == "A") |
       (d$A1 == "C" & d$A2 == "G") | (d$A1 == "G" & d$A2 == "C")
amb_ok <- !pal | is.na(d$EAF) | (d$EAF > 0.05 & d$EAF < 0.95)
d <- d[amb_ok, ]
qc <- rbind(qc, data.frame(step = "drop_ambiguous_palindromic", n = nrow(d)))

if (any(!is.na(d$INFO))) {
  d <- d[is.na(d$INFO) | d$INFO >= (gcfg$info_min %||% 0.8), ]
  qc <- rbind(qc, data.frame(step = "INFO_filter", n = nrow(d)))
}

if (any(!is.na(d$EAF))) {
  maf <- pmin(d$EAF, 1 - d$EAF)
  d <- d[is.na(maf) | maf >= (gcfg$maf_min %||% 0.01), ]
  qc <- rbind(qc, data.frame(step = "MAF_filter", n = nrow(d)))
}

if (any(!is.na(d$N))) {
  d <- d[is.na(d$N) | d$N >= (gcfg$n_min %||% 1000), ]
  qc <- rbind(qc, data.frame(step = "N_filter", n = nrow(d)))
}

# 重复位点：保留 P 最小者
d <- d[order(d$p), ]
d <- d[!duplicated(d$SNP), ]
qc <- rbind(qc, data.frame(step = "dedup_by_p", n = nrow(d)))

# 极端 beta 值（常见于低频位点或数据错误）
d <- d[abs(d$beta) < 10, ]
qc <- rbind(qc, data.frame(step = "beta_range", n = nrow(d)))

qc$retained_pct <- round(100 * qc$n / qc$n[1], 2)
write_tsv(qc, "gwas_qc_summary", cfg)

# ------------------------------------------------- 3) 基因组控制因子 --------
lambda_gc <- function(p) {
  chisq <- stats::qchisq(1 - p, df = 1)
  stats::median(chisq, na.rm = TRUE) / stats::qchisq(0.5, df = 1)
}
lam <- lambda_gc(d$p)
log_msg("基因组控制因子 lambda = ", round(lam, 4),
        if (lam > 1.1) "（>1.1，提示存在人群分层或强多效性，建议先做 LDSC 拦截）" else "")

# QQ 图
qq <- data.frame(expected = -log10(stats::ppoints(nrow(d))),
                 observed = -log10(sort(d$p)))
g_qq <- ggplot2::ggplot(qq, ggplot2::aes(expected, observed)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.3) +
  ggplot2::geom_point(size = 0.4, colour = "#1F6FB4", alpha = 0.6) +
  ggplot2::labs(x = expression("Expected " * -log[10] * "P"),
                y = expression("Observed " * -log[10] * "P"),
                title = sprintf("QQ plot (lambda = %.3f)", lam)) +
  theme_pub()
save_fig(g_qq, "gwas_qq", cfg, 4.6, 4.6)

# Manhattan 图（需要染色体与位置；没有就跳过，别用序号硬画）
if (!is.na(chr_col) && !is.na(pos_col) && "CHR" %in% names(d)) {
  dm <- data.frame(SNP = d$SNP, p = d$p, CHR = d$CHR, POS = d$POS)
  dm <- dm[!is.na(dm$CHR) & !is.na(dm$POS) & dm$CHR %in% 1:22, ]
  if (nrow(dm) > 1000) {
    dm <- dm[order(dm$CHR, dm$POS), ]
    chr_len <- tapply(dm$POS, dm$CHR, max)
    offset <- c(0, cumsum(as.numeric(chr_len))[-length(chr_len)])
    dm$x <- dm$POS + offset[match(dm$CHR, names(chr_len))]
    axis_df <- data.frame(CHR = as.integer(names(chr_len)),
                          x = as.numeric(offset) + as.numeric(chr_len) / 2)
    g_man <- ggplot2::ggplot(dm, ggplot2::aes(x, -log10(p), colour = factor(CHR %% 2))) +
      ggplot2::geom_point(size = 0.35, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c("#1F6FB4", "#8FB4D9")) +
      ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
                          colour = "#C0392B", linewidth = 0.3) +
      ggplot2::scale_x_continuous(breaks = axis_df$x, labels = axis_df$CHR) +
      ggplot2::labs(x = "chromosome", y = expression(-log[10] * "P"), title = "Manhattan") +
      theme_pub() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = 7))
    save_fig(g_man, "gwas_manhattan", cfg, 9, 3.5)
  }
} else {
  log_msg("数据中无染色体/位置列，跳过 Manhattan 图")
}

# ------------------------------------------------- 4) LDSC 格式化 -----------
# munge_sumstats 需要：SNP A1 A2 N Z（或 N_eff）。这里直接生成，省一次重跑。
z <- d$beta / d$se
n_fallback <- suppressWarnings(max(d$N, na.rm = TRUE))
if (!is.finite(n_fallback)) n_fallback <- gcfg$n_default %||% 300000
munge <- data.frame(SNP = d$SNP, A1 = d$A1, A2 = d$A2,
                    N = ifelse(is.na(d$N), n_fallback, d$N),
                    Z = z, P = d$p)
out_f <- proj_path("data", "processed", "gwas_ldsc.sumstats.gz")
data.table::fwrite(munge, out_f, sep = "\t", compress = "gzip")
log_msg("LDSC 格式已写出：", out_f)

saveRDS(d, proj_path("data", "processed", "gwas_clean.rds"))
saveRDS(list(lambda = lam, qc = qc), proj_path("data", "processed", "gwas_qc_meta.rds"))
log_msg("07 完成：保留 ", nrow(d), " 个 SNP")

# =============================================================================
# 07_gwas_qc.R —— 课题二（1）：GWAS 汇总统计清洗、质控、LDSC 格式化
#
# 输入：config$genetics$gwas$path（FinnGen / UKB / MVP 的 sumstats，.gz 可直接读）
# 输出：data/processed/gwas_clean.rds      分析子集（显著位点 ±窗口 + P<1e-4）
#       data/processed/gwas_ldsc.sumstats.gz 全量 QC 通过位点（可直接喂 ldsc.py --h2）
#       results/tables/gwas_qc_summary.tsv  各步过滤计数（写方法学时要抄这张表）
#       results/figures/gwas_qq.pdf / gwas_manhattan.pdf
#
# 质控顺序不能乱：先统一列名与等位基因，再过滤，再算 lambda。
# 过滤掉一半 SNP 后算出来的 lambda 会偏低，别拿它当"质控良好"的证据。
#
# 【为什么要流式读】
# FinnGen R12 这份 sumstats 有 2130 万行。整表 fread + as.data.frame 在本机
# （15GB 物理内存、可用不到 5GB）直接 "cannot allocate vector of size 162.7 Mb"。
# 现在用 stream_table() 一批 50 万行过两遍：
#   第 1 遍：只攒计数、全部 P 值、显著位点、绘图抽样点 —— 峰值几十 MB；
#   第 2 遍：写全量 LDSC 文件 + 抽显著位点窗口进 gwas_clean.rds。
# gwas_clean.rds 不是全量：coloc / FOCUS 只用位点附近；全量 2130 万行存成 rds，
# 下游每次 readRDS 都要 1GB 以上，既没必要也不更快。
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("07_gwas_qc", cfg)
load_pkgs(c("ggplot2", "data.table"))
set.seed(cfg_get(cfg, "project.seed", 1))

gcfg <- cfg_get(cfg, "genetics.gwas", list())
f_in <- proj_path(gcfg$path %||% "data/raw/gwas.txt.gz")
require_file(f_in, "请把 GWAS 汇总统计放到该路径，或在 config.yml 改 genetics.gwas.path")

# ------------------------------------------------- 1) 列映射 ----------------
cmap <- gcfg$cols %||% list()
hdr  <- peek_header(f_in)
need <- c("snp", "beta", "se", "pval", "ea", "nea")
missing <- need[!vapply(need, function(k) !is.null(cmap[[k]]) && cmap[[k]] %in% hdr, logical(1))]
if (length(missing)) {
  stop("以下字段未在数据中找到：", paste(missing, collapse = ", "),
       "\n  可用列：", paste(hdr, collapse = ", "),
       "\n  请修改 config.yml 的 genetics.gwas.cols", call. = FALSE)
}

chr_col   <- grep("^(chr|CHR|chrom|CHROM|chromosome)$", hdr, value = TRUE)[1]
pos_col   <- grep("^(pos|POS|position|base_pair_location|bp)$", hdr, value = TRUE)[1]
mlogp_col <- grep("^(mlogp|MLOGP|neg_log10_p|neg_log10_pval)$", hdr, value = TRUE)[1]

# 内部列名 -> 源列名；内部名同时作为 stream_table 的回传列名
src <- c(snp = cmap$snp, ea = cmap$ea, nea = cmap$nea, beta = cmap$beta,
         se = cmap$se, p = cmap$pval,
         eaf = cmap$eaf %||% NA, INFO = cmap$info %||% NA, N = cmap$n %||% NA,
         chr = chr_col %||% NA, pos = pos_col %||% NA, mlogp = mlogp_col %||% NA)
# config 里写了但文件里没有的列（如 FinnGen 的 n）在这里丢掉，不要带进读取
src <- src[!is.na(src) & nzchar(src) & src %in% hdr]
typ <- c(snp = "character", ea = "character", nea = "character",
         beta = "numeric", se = "numeric", p = "numeric", eaf = "numeric",
         INFO = "numeric", N = "numeric", chr = "character", pos = "numeric",
         mlogp = "numeric")

log_msg("读取 GWAS：", f_in, "（", round(file.size(f_in) / 1024^2), " MB）")
log_msg("列映射：", paste(names(src), src, sep = "<-", collapse = ", "))

# 样本量：FinnGen 这类 sumstats 不带 N 列，用 config 里的病例/对照数还原。
# 随手写个 30 万会让 LDSC 的 h2 与 coloc 的先验一起偏掉。
n_cases <- gcfg$n_cases; n_controls <- gcfg$n_controls
n_const <- if (!is.null(n_cases) && !is.null(n_controls)) n_cases + n_controls else (gcfg$n_default %||% 300000)
if (is.null(src["N"])) {
  log_msg("数据无 N 列，统一按 N = ", format(n_const, big.mark = ","),
          if (!is.null(n_cases)) sprintf("（%s 例 / %s 对照，来自 FinnGen R12 表型元数据）",
                                         format(n_cases, big.mark = ","),
                                         format(n_controls, big.mark = ",")) else "（config n_default）")
}
case_prop <- gcfg$case_prop %||% 0.5
log_msg("病例比例 = ", case_prop)

# ------------------------------------------------- 2) 过滤函数 --------------
# 第 1、2 遍共用同一份过滤逻辑，保证两遍口径一致
maf_min  <- gcfg$maf_min %||% 0.01
info_min <- gcfg$info_min %||% 0.8
n_min    <- gcfg$n_min %||% 1000
step_order <- c("raw", "drop_NA_core", "drop_invalid_se", "keep_biallelic_SNV",
                "drop_ambiguous_palindromic", "INFO_filter", "MAF_filter",
                "N_filter", "beta_range")

qc_filter <- function(dt) {
  cnt <- setNames(as.list(rep(0, length(step_order))), step_order)
  cnt$raw <- nrow(dt)

  # FinnGen 极显著位点的 pval 会下溢成 0，mlogp 才是真值
  if ("mlogp" %in% names(dt)) {
    fix <- !is.na(dt$p) & dt$p <= 0 & !is.na(dt$mlogp)
    if (any(fix)) dt[fix, p := 10^(-mlogp)]
  }

  keep <- !is.na(dt$snp) & !is.na(dt$beta) & !is.na(dt$se) & !is.na(dt$p) & dt$p > 0
  d <- dt[keep]; cnt$drop_NA_core <- nrow(d)

  keep <- d$se > 0 & is.finite(d$beta) & is.finite(d$se)
  d <- d[keep]; cnt$drop_invalid_se <- nrow(d)

  # 等位基因合法性：只留 A/C/G/T 单碱基（indel 与多等位在 coloc/MR 中会破坏对齐）
  a1 <- toupper(d$ea); a2 <- toupper(d$nea)
  keep <- nchar(a1) == 1 & nchar(a2) == 1 & a1 %in% c("A", "C", "G", "T") &
    a2 %in% c("A", "C", "G", "T") & a1 != a2
  d <- d[keep]; cnt$keep_biallelic_SNV <- nrow(d)

  # 回文 SNP（A/T、C/G）：MAF 接近 0.5 时链方向无法判断，MR 与 coloc 都会受污染
  a1 <- toupper(d$ea); a2 <- toupper(d$nea)
  pal <- (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
  keep <- !pal | is.na(d$eaf) | (d$eaf > 0.05 & d$eaf < 0.95)
  d <- d[keep]; cnt$drop_ambiguous_palindromic <- nrow(d)

  if ("INFO" %in% names(d) && any(!is.na(d$INFO))) {
    d <- d[is.na(d$INFO) | d$INFO >= info_min]
    cnt$INFO_filter <- nrow(d)
  }
  if ("eaf" %in% names(d) && any(!is.na(d$eaf))) {
    maf <- pmin(d$eaf, 1 - d$eaf)
    d <- d[is.na(maf) | maf >= maf_min]
    cnt$MAF_filter <- nrow(d)
  }
  if ("N" %in% names(d) && any(!is.na(d$N))) {
    d <- d[is.na(d$N) | d$N >= n_min]
    cnt$N_filter <- nrow(d)
  }
  d <- d[abs(d$beta) < 10]
  cnt$beta_range <- nrow(d)
  list(d = d, cnt = cnt)
}

# ------------------------------------------------- 3) 第 1 遍：计数与抽样 ---
cnt <- setNames(as.list(rep(0, length(step_order))), step_order)
p_chunks <- list(); lead_chunks <- list(); plot_chunks <- list()
dup_adj <- 0L

FUN1 <- local({
  k <- 0L; prev_key <- character(0)
  function(dt) {
    k <<- k + 1L
    res <- qc_filter(dt)
    d <- res$d
    for (s in step_order) cnt[[s]] <<- cnt[[s]] + res$cnt[[s]]
    if (nrow(d)) {
      key <- paste(d$chr, d$pos, toupper(d$ea), toupper(d$nea))
      dup <- duplicated(key)                       # 批内相邻/重复
      if (length(prev_key) && length(key)) dup[1] <- dup[1] || key[1] == prev_key
      if (any(dup)) { dup_adj <<- dup_adj + sum(dup); d <- d[!dup] }
      prev_key <<- key[length(key)]
    }
    p_chunks[[length(p_chunks) + 1L]] <<- d$p
    hi <- d$p < 5e-8
    if (any(hi)) lead_chunks[[length(lead_chunks) + 1L]] <<- d[hi, .(snp, chr, pos, p)]
    # Manhattan 用：显著点全留，背景按 2% 抽样（图注必须写明）
    idx <- which(!hi)
    if (length(idx)) idx <- idx[stats::runif(length(idx)) < 0.02]
    plot_chunks[[length(plot_chunks) + 1L]] <<- d[c(which(hi), idx), .(chr, pos, p)]
    NULL
  }
})

invisible(stream_table(f_in, cols = unname(src), rename = names(src),
                       types = unname(typ[names(src)]), FUN = FUN1, chunk = 5e5, verbose = TRUE))

n_raw <- cnt$raw
log_msg("原始 SNP 数：", format(n_raw, big.mark = ","))
if (dup_adj > 0) log_msg("相邻重复位点已剔除：", format(dup_adj, big.mark = ","))

qc <- data.frame(step = step_order,
                 n = vapply(step_order, function(s) as.numeric(cnt[[s]]), numeric(1)),
                 stringsAsFactors = FALSE)
qc <- qc[qc$n > 0, ]
qc$retained_pct <- round(100 * qc$n / qc$n[1], 2)
if (dup_adj > 0) {
  qc <- rbind(qc, data.frame(step = "adjacent_duplicates_removed", n = dup_adj,
                             retained_pct = round(100 * dup_adj / n_raw, 4)))
}
write_tsv(qc, "gwas_qc_summary", cfg)
n_clean <- qc$n[qc$step == "beta_range"]
log_msg("QC 后保留 ", format(n_clean, big.mark = ","), " 个 SNP（", qc$retained_pct[qc$step == "beta_range"], "%）")

# ------------------------------------------------- 4) lambda、QQ、Manhattan --
p_all <- unlist(p_chunks, use.names = FALSE)
p_chunks <- NULL
lambda_gc <- function(p) {
  chisq <- stats::qchisq(1 - p, df = 1)
  stats::median(chisq, na.rm = TRUE) / stats::qchisq(0.5, df = 1)
}
lam <- lambda_gc(p_all)
log_msg("基因组控制因子 lambda = ", round(lam, 4),
        if (lam > 1.1) "（>1.1，提示存在人群分层或强多效性，建议先做 LDSC 拦截）" else "")

# QQ：点太多画不动也没意义，等距抽 20 万个
ps <- sort(p_all)
if (length(ps) > 2e5) ps <- ps[round(seq(1, length(ps), length.out = 2e5))]
qq <- data.frame(expected = -log10(stats::ppoints(length(ps))), observed = -log10(ps))
g_qq <- ggplot2::ggplot(qq, ggplot2::aes(expected, observed)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.3) +
  ggplot2::geom_point(size = 0.4, colour = "#1F6FB4", alpha = 0.6) +
  ggplot2::labs(x = expression("Expected " * -log[10] * "P"),
                y = expression("Observed " * -log[10] * "P"),
                title = sprintf("QQ plot (lambda = %.3f)", lam)) +
  theme_pub()
save_fig(g_qq, "gwas_qq", cfg, 4.6, 4.6)
rm(ps, p_all)

dm <- data.table::rbindlist(plot_chunks)
plot_chunks <- NULL
dm[, CHR := suppressWarnings(as.integer(gsub("[^0-9]", "", chr)))]
dm <- dm[!is.na(CHR) & !is.na(pos) & CHR %in% 1:22]
if (nrow(dm) > 5000) {
  data.table::setorder(dm, CHR, pos)
  chr_len <- dm[, .(len = max(pos)), by = CHR][order(CHR)]
  chr_len[, offset := c(0, cumsum(as.numeric(len))[-.N])]
  dm <- merge(dm, chr_len[, .(CHR, offset)], by = "CHR", sort = FALSE)
  dm[, x := pos + offset]
  g_man <- ggplot2::ggplot(dm, ggplot2::aes(x, -log10(p), colour = factor(CHR %% 2))) +
    ggplot2::geom_point(size = 0.3, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c("#1F6FB4", "#8FB4D9")) +
    ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
                        colour = "#C0392B", linewidth = 0.3) +
    ggplot2::scale_x_continuous(breaks = chr_len$offset + chr_len$len / 2, labels = chr_len$CHR) +
    ggplot2::labs(x = "chromosome", y = expression(-log[10] * "P"),
                  title = "Manhattan plot (all P<1e-3; other variants 2% sampled)") +
    theme_pub() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 7))
  save_fig(g_man, "gwas_manhattan", cfg, 9, 3.5)
} else {
  log_msg("可用于 Manhattan 的位点不足，跳过")
}
rm(dm)

# ------------------------------------------------- 5) 显著位点与窗口 --------
leads <- data.table::rbindlist(lead_chunks)
lead_chunks <- NULL
leads <- leads[!is.na(chr) & !is.na(pos)]
window_bp <- (cfg_get(cfg, "genetics.gwas.window_kb", 1000) %||% 1000) * 1000
log_msg("P<5e-8 位点 ", format(nrow(leads), big.mark = ","), " 个")

if (nrow(leads)) {
  data.table::setorder(leads, p)
  picked <- list()
  for (i in seq_len(nrow(leads))) {
    hit <- FALSE
    for (pk in picked) if (pk$chr == leads$chr[i] && abs(pk$pos - leads$pos[i]) < 1e6) { hit <- TRUE; break }
    if (!hit) picked[[length(picked) + 1]] <- list(chr = leads$chr[i], pos = leads$pos[i])
  }
  w <- data.table::data.table(
    chr = vapply(picked, function(x) as.character(x$chr), character(1)),
    start = vapply(picked, function(x) as.numeric(x$pos) - window_bp, numeric(1)),
    end = vapply(picked, function(x) as.numeric(x$pos) + window_bp, numeric(1)))
  data.table::setorder(w, chr, start)
  log_msg("聚成 ", nrow(w), " 个独立位点，窗口 ±", window_bp / 1000, " kb")
} else {
  w <- data.table::data.table(chr = character(0), start = numeric(0), end = numeric(0))
}

# 窗口命中判定：按染色体建有序起点 + 累积最大终点，findInterval 二分，O(n log n)
in_window <- function(chr, pos, w) {
  out <- logical(length(pos))
  for (cc in unique(chr)) {
    ww <- w[chr == cc]
    if (!nrow(ww)) next
    j <- which(chr == cc)
    i <- findInterval(pos[j], ww$start)
    cmax <- cummax(ww$end)
    out[j] <- i > 0 & pos[j] <= cmax[pmax(i, 1)]
  }
  out
}

# ------------------------------------------------- 6) 第 2 遍：全量 LDSC + 窗口 ---
out_ldsc <- proj_path("data", "processed", "gwas_ldsc.sumstats.gz")
tmp_plain <- proj_path("data", "processed", "_gwas_ldsc_tmp.tsv")
if (file.exists(tmp_plain)) unlink(tmp_plain)
if (file.exists(out_ldsc)) unlink(out_ldsc)
clean_chunks <- list(); written <- 0L

FUN2 <- local({
  first <- TRUE; kk <- 0L; prev_key <- character(0)
  function(dt) {
    res <- qc_filter(dt)
    d <- res$d
    if (!nrow(d)) return(NULL)
    # 与第 1 遍同一套去重口径，否则两遍的位点数对不上
    key <- paste(d$chr, d$pos, toupper(d$ea), toupper(d$nea))
    dup <- duplicated(key)
    if (length(prev_key) && length(key)) dup[1] <- dup[1] || key[1] == prev_key
    if (any(dup)) d <- d[!dup]
    prev_key <<- key[length(key)]

    n_col <- rep(n_const, nrow(d))
    if ("N" %in% names(d)) n_col <- ifelse(is.na(d$N), n_const, d$N)
    munge <- data.table::data.table(
      SNP = d$snp, A1 = toupper(d$ea), A2 = toupper(d$nea),
      N = n_col, Z = d$beta / d$se, P = d$p, BETA = d$beta, SE = d$se)
    data.table::fwrite(munge, tmp_plain, sep = "\t", append = !first,
                       col.names = first, showProgress = FALSE)
    first <<- FALSE
    written <<- written + nrow(munge)

    keep <- d$p < 1e-4
    if (nrow(w) && any(!keep)) {
      rest <- which(!keep)
      keep[rest] <- in_window(d$chr[rest], d$pos[rest], w)
    }
    if (any(keep)) {
      kk <<- kk + 1L
      clean_chunks[[kk]] <<- data.table::data.table(
        SNP = d$snp[keep], A1 = toupper(d$ea[keep]), A2 = toupper(d$nea[keep]),
        beta = d$beta[keep], se = d$se[keep], p = d$p[keep],
        EAF = if ("eaf" %in% names(d)) d$eaf[keep] else NA_real_,
        N = n_col[keep],
        INFO = if ("INFO" %in% names(d)) d$INFO[keep] else NA_real_,
        CHR = d$chr[keep], POS = d$pos[keep])
    }
    NULL
  }
})

invisible(stream_table(f_in, cols = unname(src), rename = names(src),
                       types = unname(typ[names(src)]), FUN = FUN2, chunk = 5e5, verbose = TRUE))

log_msg("LDSC 格式写出 ", format(written, big.mark = ","), " 个 SNP，压缩中……")
R.utils::gzip(tmp_plain, destname = out_ldsc, overwrite = TRUE, remove = TRUE)
log_msg("LDSC 格式已写出：", out_ldsc, "（", round(file.size(out_ldsc) / 1024^2, 1), " MB）")

d <- data.table::rbindlist(clean_chunks)
clean_chunks <- NULL
d[, pos_key := paste0(CHR, "_", POS)]
d[, CHR := suppressWarnings(as.integer(gsub("[^0-9]", "", CHR)))]
data.table::setcolorder(d, c("SNP", "A1", "A2", "beta", "se", "p", "EAF", "N", "INFO", "CHR", "POS", "pos_key"))
log_msg("gwas_clean 分析子集：", format(nrow(d), big.mark = ","), " 个 SNP，", nrow(w), " 个位点窗口")

saveRDS(as.data.frame(d), proj_path("data", "processed", "gwas_clean.rds"))
saveRDS(list(lambda = lam, qc = qc, n_total = n_const, case_prop = case_prop,
             n_lead = nrow(leads), n_loci = nrow(w), window_bp = window_bp,
             source_file = basename(f_in)),
        proj_path("data", "processed", "gwas_qc_meta.rds"))

log_msg("07 完成：全量 ", format(written, big.mark = ","), " 个 SNP 写入 LDSC 文件；",
        "分析子集 ", format(nrow(d), big.mark = ","), " 个；lambda = ", round(lam, 3))

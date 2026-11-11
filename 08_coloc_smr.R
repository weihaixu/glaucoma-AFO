# =============================================================================
# 08_coloc_smr.R —— 课题二（2）：共定位（coloc）+ SMR/HEIDI 近似
#
# 回答的问题：GWAS 位点上的信号，和「哪个基因在哪个组织的表达」是同一个因果变异吗？
#
# 输入：data/processed/gwas_clean.rds（07 产出）
#       eQTL：GTEx v8 allpairs（tabix）或本地 tsv（见 config$genetics$eqtl）
# 输出：results/tables/coloc_summary.tsv      每个 gene × tissue 的 PP.H0–H4
#       results/tables/coloc_sensitivity.tsv  PP.H4 对先验的敏感性
#       results/tables/smr_heidi.tsv          SMR 效应量与 HEIDI 异质性检验
#       results/figures/coloc_locus_<gene>.pdf
#
# 判读标准：
#   PP.H4 > 0.8 才算"共享因果变异"；PP.H3 > 0.8 说明是两个独立信号，
#   很多论文只报 H4 不报 H3，这是最常见的选择性报告。
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("08_coloc_smr", cfg)
load_pkgs(c("ggplot2", "coloc", "data.table", "Rsamtools"))
set.seed(cfg_get(cfg, "project.seed", 1))

gwas <- readRDS(proj_path("data", "processed", "gwas_clean.rds"))
ecfg <- cfg_get(cfg, "genetics.eqtl", list())
ccfg <- cfg_get(cfg, "genetics.coloc", list())
win  <- (ecfg$window_kb %||% 500) * 1000
source_type <- ecfg$source %||% "local_tsv"

# ------------------------------------------------- 1) 位点定义 --------------
# 优先用配置里显式给出的位点；否则从 GWAS 显著位点自动聚类
loci <- cfg_get(cfg, "genetics.loci", NULL)
if (is.null(loci)) {
  sig <- gwas[gwas$p < 5e-8, ]
  if (!nrow(sig)) stop("GWAS 中无 P<5e-8 位点，请在 config 中显式给出 genetics.loci")
  sig <- sig[!is.na(sig$CHR) & !is.na(sig$POS), ]
  if (!nrow(sig)) stop("显著位点缺少染色体位置，无法自动定位；请在 config 给出 genetics.loci")
  sig <- sig[order(sig$p), ]
  chr_col <- if ("CHR" %in% names(sig)) "CHR" else NA
  if (is.na(chr_col)) stop("gwas_clean 缺少 CHR，无法自动定位；请在 config 中显式给出 genetics.loci（chr:start-end）")
  picked <- list()
  for (i in seq_len(nrow(sig))) {
    hit <- FALSE
    for (pk in picked) if (pk$chr == sig$CHR[i] && abs(pk$pos - sig$POS[i]) < 1e6) { hit <- TRUE; break }
    if (!hit) picked[[length(picked) + 1]] <- list(chr = sig$CHR[i], pos = sig$POS[i], snp = sig$SNP[i])
  }
  loci <- lapply(picked, function(x) list(chr = x$chr, start = x$pos - win, end = x$pos + win, lead = x$snp))
  log_msg("自动识别 ", length(loci), " 个独立位点")
}
# 每个位点 × 每个组织要发一次远程区间请求，85 个位点 × 6 个组织 ≈ 500 次，
# 一次 2–3MB，跑一晚上也跑不完。默认只做信号最强的若干个（按 lead P 排序）。
max_loci <- cfg_get(cfg, "genetics.coloc.max_loci", 20) %||% 20
if (length(loci) > max_loci) {
  log_msg("位点 ", length(loci), " 个，按 coloc.max_loci 只做前 ", max_loci, " 个")
  loci <- loci[seq_len(max_loci)]
}
log_msg("待分析位点：", length(loci), " 个")

# ------------------------------------------------- 2) eQTL 区间读取 ---------
read_eqtl_region <- function(chr, start, end, tissue, cfg) {
  st  <- cfg_get(cfg, "genetics.eqtl.source", "local_tsv")
  reg <- paste0("chr", chr, ":", start, "-", end)

  if (identical(st, "gtex")) {
    d <- cfg_get(cfg, "genetics.eqtl.gtex_dir", "data/raw/GTEx_v8")
    f <- proj_path(d, paste0(tissue, ".v8.EUR.allpairs.txt.gz"))
    if (!file.exists(f)) { log_msg("缺少 GTEx 文件：", f); return(NULL) }
    if (!file.exists(paste0(f, ".tbi"))) { log_msg("缺少 tabix 索引：", f, ".tbi"); return(NULL) }
    chunks <- tryCatch(Rsamtools::scanTabix(Rsamtools::TabixFile(f), param = reg)[[1]],
                       error = function(e) character(0))
    if (!length(chunks)) return(NULL)
    df <- data.table::fread(text = chunks, sep = "\t", header = FALSE, showProgress = FALSE)
    if (ncol(df) < 9) return(NULL)
    names(df)[1:9] <- c("variant_id", "gene_id", "tss_distance", "maf", "ma_samples",
                        "ma_count", "pval_nominal", "slope", "slope_se")
    # variant_id = chr1_12345_A_G_b38：位置键用于对齐（GTEx 不给 rsID）
    parts <- strsplit(df$variant_id, "_", fixed = TRUE)
    pv <- function(i) vapply(parts, function(x) if (length(x) >= i) x[i] else NA_character_, character(1))
    out <- data.frame(
      pos_key = paste0(sub("^chr", "", pv(1)), "_", pv(2)),
      SNP = NA_character_,
      A1 = pv(4),                        # GTEx: ALT 为效应等位基因（slope 对应 ALT）
      A2 = pv(3),                        # REF
      gene = sub("\\..*$", "", df$gene_id),
      beta = df$slope, se = df$slope_se, p = df$pval_nominal,
      maf = df$maf, tissue = tissue, stringsAsFactors = FALSE
    )
    return(out)
  }

  if (identical(st, "gtex_remote")) {
    # GTEx v8 allpairs 的公开镜像在 EBI eQTL Catalogue 的 imported/GTEx_V8/ge/，
    # 每个组织 3.5GB，整包下不动；本机 Rsamtools 的 Windows 构建不支持 http 远端
    # tabix，pysam 也装不上（PyPI 被代理拦），于是用 tools/fetch_eqtl_tabix.py
    # 自己解析 .tbi + HTTP Range 取块：一个 1Mb 窗口通常只传 1–3MB。
    base <- cfg_get(cfg, "genetics.eqtl.gtex_tabix_base",
                     "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge")
    cache <- proj_path(cfg_get(cfg, "genetics.eqtl.cache_dir", "data/raw/eqtl_cache"))
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    key <- sprintf("%s_%s_%d_%d.tsv", tissue, chr, as.integer(start), as.integer(end))
    out_f <- file.path(cache, key)
    if (!file.exists(out_f) || file.size(out_f) < 2) {
      url <- paste0(base, "/", tissue, ".tsv.gz")
      py <- Sys.which("python")
      if (!nzchar(py)) { log_msg("找不到 python，无法远程取 eQTL"); return(NULL) }
      log_msg("  远程取 eQTL：", tissue, " ", chr, ":", round(start / 1e6, 2), "-", round(end / 1e6, 2), " Mb")
      rc <- system2(py, c(proj_path("tools", "fetch_eqtl_tabix.py"),
                          "--url", url, "--index", paste0(url, ".tbi"),
                          "--chrom", as.character(chr),
                          "--start", as.character(as.integer(start)),
                          "--end", as.character(as.integer(end)),
                          "--out", out_f),
                    stdout = FALSE, stderr = FALSE)
      if (rc != 0 || !file.exists(out_f)) { log_msg("  取数失败：", tissue, " ", chr); return(NULL) }
    }
    df <- tryCatch(data.table::fread(out_f, showProgress = FALSE), error = function(e) NULL)
    if (is.null(df) || !nrow(df) || !"beta" %in% names(df)) return(NULL)
    out <- data.frame(
      SNP = as.character(df$rsid),
      pos_key = paste0(df$chromosome, "_", df$position),
      gene = sub("\\..*$", "", as.character(df$gene_id)),
      A1 = toupper(as.character(df$alt)), A2 = toupper(as.character(df$ref)),
      beta = as.numeric(df$beta), se = as.numeric(df$se),
      p = as.numeric(df$pvalue), maf = as.numeric(df$maf),
      tissue = tissue, stringsAsFactors = FALSE)
    # GTEx 每个变异给了等位基因数 an，样本量 = an/2；比在 config 里写死一个数可靠
    if ("an" %in% names(df)) attr(out, "n") <- stats::median(df$an, na.rm = TRUE) / 2
    return(out)
  }

  if (identical(st, "eqtl_catalogue")) {
    d <- cfg_get(cfg, "genetics.eqtl.catalogue_dir", "data/raw/eQTL_Catalogue")
    f <- proj_path(d, paste0(tissue, ".all.tsv.gz"))
    if (!file.exists(f)) { log_msg("缺少 eQTL Catalogue 文件：", f); return(NULL) }
    chunks <- tryCatch(Rsamtools::scanTabix(Rsamtools::TabixFile(f), param = reg)[[1]],
                       error = function(e) character(0))
    if (!length(chunks)) return(NULL)
    df <- data.table::fread(text = chunks, sep = "\t", header = TRUE, showProgress = FALSE)
    pick <- function(cands) { h <- intersect(cands, names(df)); if (length(h)) as.character(df[[h[1]]]) else rep(NA_character_, nrow(df)) }
    out <- data.frame(
      SNP = pick(c("rsid", "rsids")),
      pos_key = if (all(c("chr", "pos") %in% names(df))) paste0(df$chr, "_", df$pos) else NA_character_,
      gene = sub("\\..*$", "", pick(c("gene_id", "molecular_trait_id"))),
      A1 = pick(c("alt", "ea")), A2 = pick(c("ref", "nea")),
      beta = as.numeric(pick(c("beta", "slope"))),
      se   = as.numeric(pick(c("se", "slope_se"))),
      p    = as.numeric(pick(c("pvalue", "pval_nominal"))),
      maf  = as.numeric(pick(c("maf", "eaf"))),
      tissue = tissue, stringsAsFactors = FALSE
    )
    return(out)
  }

  # local_tsv：整表读入后按区间过滤（文件很大时请按染色体预先切分）
  f <- proj_path(cfg_get(cfg, "genetics.eqtl.local_tsv", "data/raw/eqtl_allpairs.tsv.gz"))
  if (!file.exists(f)) { log_msg("缺少本地 eQTL 文件：", f); return(NULL) }
  df <- data.table::fread(f, showProgress = FALSE)
  need <- c("SNP", "gene", "A1", "A2", "beta", "se", "p")
  if (!all(need %in% names(df))) stop("本地 eQTL 需包含列：", paste(need, collapse = ", "), call. = FALSE)
  if (!all(c("chr", "pos") %in% names(df))) stop("本地 eQTL 需包含 chr / pos 列以支持区间过滤", call. = FALSE)
  df <- df[df$chr == chr & df$pos >= start & df$pos <= end, ]
  if (!nrow(df)) return(NULL)
  df <- as.data.frame(df)
  df$pos_key <- paste0(df$chr, "_", df$pos)
  df$tissue <- tissue
  df
}

# ------------------------------------------------- 3) 等位基因对齐 ----------
complement <- function(x) chartr("ACGT", "TGCA", x)

# GWAS 与 eQTL 对齐：优先 rsID，其次 chr_pos 键（GTEx 只有位置没有 rsID）
harmonize <- function(g, q, allow_complement = TRUE) {
  q_snp_ok <- "SNP" %in% names(q) && mean(!is.na(q$SNP)) > 0.5
  if (q_snp_ok && mean(g$SNP %in% q$SNP) > 0.3) {
    g$key <- g$SNP; q$key <- q$SNP
  } else if ("pos_key" %in% names(g) && "pos_key" %in% names(q)) {
    g$key <- g$pos_key; q$key <- q$pos_key
  } else {
    return(NULL)
  }
  g <- g[!is.na(g$key), ]; q <- q[!is.na(q$key), ]
  q <- q[!duplicated(q$key), ]
  m <- merge(g, q, by = "key", suffixes = c(".g", ".q"))
  if (!nrow(m)) return(NULL)

  same <- m$A1.g == m$A1.q & m$A2.g == m$A2.q
  flip <- m$A1.g == m$A2.q & m$A2.g == m$A1.q
  if (allow_complement) {
    same <- same | (m$A1.g == complement(m$A1.q) & m$A2.g == complement(m$A2.q))
    flip <- flip | (m$A1.g == complement(m$A2.q) & m$A2.g == complement(m$A1.q))
  }
  m$beta.q[flip] <- -m$beta.q[flip]
  if ("maf" %in% names(m)) m$maf[flip] <- 1 - m$maf[flip]
  keep <- same | flip
  keep[is.na(keep)] <- FALSE          # 等位基因缺失的行直接丢弃，不能让 NA 混进 coloc
  m <- m[keep, ]
  m
}

# ------------------------------------------------- 4) coloc 主循环 ----------
tissues <- ecfg$tissues %||% c("Artery")
res_all <- list(); sens_all <- list(); smr_all <- list()

for (li in seq_along(loci)) {
  lo <- loci[[li]]
  gsub_df <- gwas[!is.na(gwas$CHR) & gwas$CHR == lo$chr &
                    gwas$POS >= lo$start & gwas$POS <= lo$end, ]
  if (nrow(gsub_df) < 50) { log_msg("位点 ", li, " GWAS SNP 过少，跳过"); next }

  for (ts in tissues) {
    q <- read_eqtl_region(lo$chr, lo$start, lo$end, ts, cfg)
    if (is.null(q) || !nrow(q)) next
    # 该组织里连名义信号都没有的基因不做 coloc：没有 eQTL 就不可能有共定位，
    # 逐个跑只是把 200 个基因 × 6 个组织白算一遍（coloc.abf 本身很便宜，
    # 贵的是每个基因一次 merge）。
    gene_minp <- tapply(q$p, q$gene, min, na.rm = TRUE)
    gset <- names(gene_minp)[is.finite(gene_minp) & gene_minp < 1e-5]
    if (!length(gset)) { log_msg("  ", ts, " @ 位点", li, "：窗口内无 P<1e-5 的 eQTL，跳过"); next }
    log_msg("  位点", li, " / ", ts, "：", length(unique(q$gene)), " 个基因中 ",
            length(gset), " 个有 eQTL 信号（P<1e-5），进入 coloc")
    for (gn in gset) {
      qq <- q[q$gene == gn, ]
      h <- harmonize(gsub_df, qq)
      if (is.null(h) || nrow(h) < 20) {
        log_msg("  [跳过] ", gn, "：harmonize 后仅 ",
                if (is.null(h)) 0 else nrow(h), " 个 SNP（<20）")
        next
      }

      # coloc 需要 SNP、beta、varbeta、MAF、N；eQTL 用 quant，GWAS 用 cc（有病例比例时）
      gtype <- cfg_get(cfg, "genetics.gwas.type", "cc")
      n_gwas <- suppressWarnings(max(h$N, na.rm = TRUE))
      if (!is.finite(n_gwas)) n_gwas <- cfg_get(cfg, "genetics.gwas.n_default", 300000)
      d1 <- list(snp = h$key, beta = h$beta.g, varbeta = h$se.g^2,
                 type = gtype, N = n_gwas)
      if (identical(gtype, "cc")) d1$s <- cfg_get(cfg, "genetics.gwas.case_prop", 0.3)
      # quant 型 GWAS（数量性状，如 IOP）必须给 MAF（或 sdY），否则 coloc 报错：
      # "must give sdY for type quant, or, if sdY unknown, MAF and N"
      if (identical(gtype, "quant")) {
        eaf_g <- if ("EAF" %in% names(h)) as.numeric(h$EAF) else rep(NA_real_, nrow(h))
        d1$MAF <- pmin(eaf_g, 1 - eaf_g)
        if (anyNA(d1$MAF)) {
          log_msg("  [跳过] ", gn, "：quant GWAS 的 EAF 含 ", sum(is.na(d1$MAF)),
                  " 个缺失，无法估计 sdY")
          next
        }
      }
      n_eqtl <- attr(q, "n")
      if (is.null(n_eqtl) || !is.finite(n_eqtl)) n_eqtl <- ecfg$n_eqtl %||% 300
      d2 <- list(snp = h$key, beta = h$beta.q, varbeta = h$se.q^2,
                 type = "quant", N = n_eqtl, MAF = h$maf)
      if (any(is.na(d1$varbeta)) || any(is.na(d2$varbeta))) next

      if (nzchar(Sys.getenv("COLOC_DEBUG"))) {
        log_msg("  [调试] ", gn, " N1=", format(d1$N), "(", class(d1$N), ")",
                " M1=", class(d1$MAF), " N2=", format(d2$N), "(", class(d2$N), ")",
                " M2=", class(d2$MAF), " b1=", class(d1$beta),
                " v1 NA=", sum(is.na(d1$varbeta)), " M1 NA=", sum(is.na(d1$MAF)),
                " M2 NA=", sum(is.na(d2$MAF)))
      }

      fit <- tryCatch(coloc::coloc.abf(d1, d2,
                        p1 = ccfg$p1 %||% 1e-4, p2 = ccfg$p2 %||% 1e-4,
                        p12 = ccfg$p12 %||% 1e-5),
                      error = function(e) {
                        log_msg("  [coloc 报错] ", gn, "/", ts, "：", conditionMessage(e), " | call: ", paste(deparse(conditionCall(e)), collapse=" "))
                        tb <- paste(capture.output(traceback(sys.frame(sys.nframe() - 1))), collapse = " | ")
                        log_msg("  [调用栈] ", substr(tb, 1, 500))
                        NULL
                      })
      if (is.null(fit)) next
      s <- fit$summary
      res_all[[length(res_all) + 1]] <- data.frame(
        locus = li, chr = lo$chr, start = lo$start, end = lo$end,
        gene = gn, tissue = ts, n_snp = s[["nsnps"]],
        PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]], PP.H2 = s[["PP.H2.abf"]],
        PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]],
        # coloc 5.x 的 fit$results 是 data.frame，行名只是 "1","2","3"…（行号），
        # 必须取它的 snp 列；旧版靠 rownames 才会得到 rsID。取不到再退回 d1$snp。
        top_snp = {
          i <- which.max(fit$results$SNP.PP.H4)
          v <- if ("snp" %in% names(fit$results)) as.character(fit$results$snp[i]) else NA_character_
          if (is.na(v) || !nzchar(v) || grepl("^[0-9]+$", v)) v <- as.character(d1$snp[i])
          v
        },
        stringsAsFactors = FALSE
      )

      # 敏感性：先验改变时 H4 是否稳定（不稳定就别写进结论）
      if (s[["PP.H4.abf"]] > 0.5) {
        for (p12 in c(1e-6, 1e-5, 1e-4)) {
          f2 <- tryCatch(coloc::coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12),
                         error = function(e) NULL)
          if (!is.null(f2)) sens_all[[length(sens_all) + 1]] <- data.frame(
            gene = gn, tissue = ts, p12 = p12, PP.H4 = f2$summary[["PP.H4.abf"]])
        }
        # 位点放大图
        # coloc 5.x 的 fit$results 只有 lABF / SNP.PP.H4，不再返回 pvalues1/2；
        # 直接由 beta/se 还原 p 值，避免依赖 coloc 的内部列名
        plot_df <- data.frame(
          snp = d1$snp,
          p1 = 2 * stats::pnorm(-abs(d1$beta / sqrt(d1$varbeta))),
          p2 = 2 * stats::pnorm(-abs(d2$beta / sqrt(d2$varbeta))),
          pp = fit$results$SNP.PP.H4[match(d1$snp, rownames(fit$results))]
        )
        plot_df <- plot_df[is.finite(plot_df$p1) & is.finite(plot_df$p2) &
                             plot_df$p1 > 0 & plot_df$p2 > 0, ]
        g1 <- ggplot2::ggplot(plot_df, ggplot2::aes(-log10(p1), -log10(p2))) +
          ggplot2::geom_point(size = 0.7, alpha = 0.6, colour = "#1F6FB4") +
          ggplot2::labs(x = "GWAS -log10(P)", y = "eQTL -log10(P)",
                        title = paste0(gn, " (", ts, ")  PP.H4=", sprintf("%.2f", s[["PP.H4.abf"]])),
                        subtitle = paste0("top SNP: ", rownames(fit$results)[which.max(fit$results$SNP.PP.H4)])) +
          theme_pub()
        save_fig(g1, paste0("coloc_locus_", gn, "_", ts), cfg, 4.6, 4.4)
      }

      # ---- SMR（单工具变量）------------------------------------------------
      # 官方 SMR 以该基因最强的 eQTL SNP 作工具：b_SMR = b_GWAS / b_eQTL，
      # se 用 delta 法。HEIDI 要靠 LD 矩阵比较多个 SNP 的异质性，没有 1000G
      # 参考面板做不了，这里直接留 NA 并标注——不能用"区域内所有 SNP 逆方差
      # 加权"冒充 HEIDI：这些 SNP 高度连锁，那样算出来的 p 会成片过度显著。
      i_top <- which.min(h$p.q)
      if (length(i_top) == 1 && is.finite(h$p.q[i_top]) &&
          h$p.q[i_top] < 5e-8 && abs(h$beta.q[i_top]) > 1e-8) {
        bq <- h$beta.q[i_top]; bg <- h$beta.g[i_top]
        b_smr <- bg / bq
        se_smr <- sqrt(h$se.g[i_top]^2 / bq^2 + bg^2 * h$se.q[i_top]^2 / bq^4)
        smr_all[[length(smr_all) + 1]] <- data.frame(
          gene = gn, tissue = ts, n_snp_region = nrow(h),
          top_snp = as.character(h$key[i_top]),
          p_eqtl = h$p.q[i_top],
          F_eqtl = (bq / h$se.q[i_top])^2,
          b_SMR = b_smr, se_SMR = se_smr,
          p_SMR = 2 * stats::pnorm(-abs(b_smr / se_smr)),
          Q_HEIDI = NA_real_, df = NA_real_, p_HEIDI = NA_real_,
          heidi_pass = NA,
          note = "HEIDI 需 1000G LD 面板，未计算",
          stringsAsFactors = FALSE)
      }
    }
  }
  log_msg("位点 ", li, " 完成")
}

# ------------------------------------------------- 5) 汇总 ------------------
# ENSG 写进表里没人看得懂；有 org.Hs.eg.db 就补一列 SYMBOL（不装也能跑）
add_symbol <- function(df) {
  if (is.null(df) || !nrow(df) || !"gene" %in% names(df)) return(df)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) return(df)
  sym <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys = unique(df$gene), column = "SYMBOL",
    keytype = "ENSEMBL", multiVals = "first"))
  df$symbol <- unname(sym[df$gene])
  df[, c("symbol", setdiff(names(df), "symbol")), drop = FALSE]
}

if (length(res_all)) {
  coloc_df <- do.call(rbind, res_all)
  coloc_df <- add_symbol(coloc_df)
  coloc_df <- coloc_df[order(-coloc_df$PP.H4), ]
  write_tsv(coloc_df, "coloc_summary", cfg)
  strong <- coloc_df[coloc_df$PP.H4 > (ccfg$pp4_threshold %||% 0.8), ]
  log_msg("PP.H4 > ", ccfg$pp4_threshold %||% 0.8, " 的 gene-tissue 对：", nrow(strong))
  if (nrow(strong)) print(strong[, c("gene", "tissue", "PP.H3", "PP.H4", "top_snp")])
} else {
  log_msg("无可分析的 gene-tissue 对：检查 eQTL 数据路径与位点坐标")
}

if (length(sens_all)) write_tsv(do.call(rbind, sens_all), "coloc_sensitivity", cfg)
if (length(smr_all)) {
  smr_df <- add_symbol(do.call(rbind, smr_all))
  smr_df <- smr_df[order(smr_df$p_SMR), ]
  smr_df$fdr_SMR <- stats::p.adjust(smr_df$p_SMR, "BH")
  write_tsv(smr_df, "smr_heidi", cfg)
  log_msg("单工具变量 SMR：", nrow(smr_df), " 个 gene-tissue 对，FDR<0.05 的 ",
          sum(smr_df$fdr_SMR < 0.05), " 个（HEIDI 未做，需 LD 面板）")

}

# ------------------------------------------------- 6) 外部 SMR 二进制 -------
# 官方 SMR/HEIDI 实现需要 BESD 格式，脚本不重复造轮子：
#   1) 用 --eqtl-flist 生成 .besd
#   2) smr --bfile <ref> --gwas-summary <file> --beqtl-summary <besd> --out <prefix> --peqtl-smr 5e-8
# 若本机有 SMR 可执行文件与 1000G 面板，跑完可核对该单工具变量估计
smr_bin <- proj_path("tools", "smr")
if (file.exists(smr_bin) && length(loci)) {
  log_msg("检测到 SMR 二进制，可按 README 中的命令行参数调用；本脚本只给单工具变量估计")
}

log_msg("08 完成")

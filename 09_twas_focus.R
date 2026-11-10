# =============================================================================
# 09_twas_focus.R —— 课题二（3）：TWAS（FUSION）+ 精细定位（FOCUS）
#
# 逻辑：eQTL 权重（SNP→表达） × GWAS z 值 → 基因表达与疾病的关联
#       FOCUS 再在共调控基因之间做精细定位，回答"这一位点到底是哪个基因"
#
# 输入：data/processed/gwas_ldsc.sumstats.gz（07 产出，FUSION 要求的 SNP A1 A2 N Z）
#       data/raw/FUSION_weights/<tissue>/<tissue>.ENSG*.wgt.RDat
#       tools/FUSION/FUSION.assoc_test.R
# 输出：results/tables/twas_<tissue>.tsv       每个基因的 TWAS.Z / TWAS.P
#       results/tables/twas_significant.tsv     Bonferroni 显著基因
#       results/tables/focus_<tissue>.tsv       后验纳入概率（PIP）
#       results/figures/twas_manhattan.pdf
#
# 注意：TWAS 是"遗传预测的表达"，不是实测表达；写文章时必须说清这一点，
#       否则会被质疑"用 GWAS 反推表达等于循环论证"。
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("09_twas_focus", cfg)
load_pkgs(c("ggplot2", "data.table"))
set.seed(cfg_get(cfg, "project.seed", 1))

tcfg <- cfg_get(cfg, "genetics.twas", list())
fusion_dir <- proj_path(tcfg$fusion_dir %||% "tools/FUSION")
sumstats <- proj_path("data", "processed", "gwas_ldsc.sumstats.gz")
require_file(sumstats, "先跑 07_gwas_qc.R 生成 LDSC 格式的 sumstats")
require_file(file.path(fusion_dir, "FUSION.assoc_test.R"),
             paste0("FUSION 下载：https://github.com/gusevlab/fusion_twas ，解压到 ", fusion_dir))

ld_ref <- proj_path(cfg_get(cfg, "genetics.ld.panel_plink", "data/raw/1000G_EUR_Phase3_plink/1000G.EUR.QC"))
if (!file.exists(paste0(ld_ref, ".bed"))) {
  log_msg("警告：未找到 LD 参考面板 ", ld_ref, "，FUSION 无法运行；请先下载 1000G plink 文件")
}

tissues <- names(tcfg$weights %||% list())
if (!length(tissues)) {
  wdir <- proj_path(tcfg$weights_dir %||% "data/raw/FUSION_weights")
  if (dir.exists(wdir)) tissues <- list.dirs(wdir, recursive = FALSE, full.names = FALSE)
}
if (!length(tissues)) stop("未找到任何 FUSION 权重目录：", tcfg$weights_dir %||% "data/raw/FUSION_weights")

# ------------------------------------------------- 1) 逐组织跑 FUSION -------
run_fusion_chr <- function(tissue, chr, cfg) {
  wdir <- proj_path(tcfg$weights_dir %||% "data/raw/FUSION_weights", tissue)
  wgts <- list.files(wdir, pattern = "\\.wgt\\.RDat$", full.names = TRUE)
  if (!length(wgts)) { log_msg(tissue, " chr", chr, "：无权重文件"); return(NULL) }
  out <- proj_path(cfg_get(cfg, "project.outdir", "results"), "tables", "fusion_raw",
                   sprintf("%s.chr%d.dat", tissue, chr))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(out)) {
    log_msg("已存在，跳过：", out)
    return(out)
  }
  args <- c(file.path(fusion_dir, "FUSION.assoc_test.R"),
            "--sumstats", sumstats,
            "--weights", paste(wgts, collapse = ","),
            "--weights_dir", wdir,
            "--ref_ld_chr", paste0(ld_ref, "."),
            "--chr", chr,
            "--out", out)
  log_msg("FUSION: ", tissue, " chr", chr, "（", length(wgts), " 个基因权重）")
  st <- system2("Rscript", args, stdout = TRUE, stderr = TRUE)
  if (!file.exists(out)) {
    log_msg("FUSION 失败，末尾输出：\n", paste(tail(st, 15), collapse = "\n"))
    return(NULL)
  }
  out
}

chrs <- tcfg$chr %||% 1:22
for (ts in tissues) {
  for (cc in chrs) run_fusion_chr(ts, cc, cfg)
}

# ------------------------------------------------- 2) 汇总 ------------------
dat_files <- list.files(proj_path(cfg_get(cfg, "project.outdir", "results"), "tables", "fusion_raw"),
                        pattern = "\\.dat$", full.names = TRUE)
if (!length(dat_files)) {
  log_msg("没有 FUSION 输出，后续步骤跳过。请检查权重文件与 LD 面板。")
} else {
  twas <- do.call(rbind, lapply(dat_files, function(f) {
    d <- tryCatch(utils::read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(NULL)
    d$tissue <- sub("\\.chr[0-9]+\\.dat$", "", basename(f))
    d
  }))
  if (!is.null(twas) && nrow(twas)) {
    need <- c("ID", "TWAS.Z", "TWAS.P")
    if (all(need %in% names(twas))) {
      twas <- twas[!is.na(twas$TWAS.P), ]
      twas <- twas[order(twas$TWAS.P), ]
      write_tsv(twas, "twas_all", cfg)

      per_ts <- split(twas, twas$tissue)
      sig_list <- lapply(names(per_ts), function(ts) {
        d <- per_ts[[ts]]
        thr <- 0.05 / nrow(d)
        s <- d[d$TWAS.P < thr, ]
        if (!nrow(s)) return(NULL)
        s$bonferroni_thr <- thr
        s
      })
      sig <- do.call(rbind, sig_list)
      if (!is.null(sig)) {
        write_tsv(sig, "twas_significant", cfg)
        log_msg("Bonferroni 显著基因：", nrow(sig))
      }

      # 曼哈顿式散点：按染色体位置（若有 P0 列）着色
      if (all(c("CHR", "P0") %in% names(twas))) {
        d <- twas[!is.na(twas$CHR) & !is.na(twas$P0), ]
        d <- d[order(d$CHR, d$P0), ]
        d$rank <- seq_len(nrow(d))
        g <- ggplot2::ggplot(d, ggplot2::aes(rank, -log10(TWAS.P), colour = factor(CHR %% 2))) +
          ggplot2::geom_point(size = 0.5, show.legend = FALSE) +
          ggplot2::scale_colour_manual(values = c("#1F6FB4", "#8FB4D9")) +
          ggplot2::geom_hline(yintercept = -log10(0.05 / nrow(d)), linetype = "dashed",
                              colour = "#C0392B", linewidth = 0.3) +
          ggplot2::labs(x = "gene ordered by genomic position",
                        y = expression(-log[10] * "TWAS P"), title = "TWAS") +
          theme_pub()
        save_fig(g, "twas_manhattan", cfg, 9, 3.6)
      }
    }
  }
}

# ------------------------------------------------- 3) FOCUS 精细定位 -------
if (requireNamespace("FOCUS", quietly = TRUE) && length(dat_files)) {
  for (ts in unique(sub("\\.chr[0-9]+\\.dat$", "", basename(dat_files)))) {
    f <- proj_path(cfg_get(cfg, "project.outdir", "results"), "tables", "fusion_raw", paste0(ts, ".chr1.dat"))
    if (!file.exists(f)) next
    wdir <- proj_path(tcfg$weights_dir %||% "data/raw/FUSION_weights", ts)
    ok <- tryCatch({
      res <- FOCUS::focus(fusion.res = read.delim(f),
                          wgt.paths = list.files(wdir, pattern = "\\.wgt\\.RDat$", full.names = TRUE),
                          ref_ld = paste0(ld_ref, "."),
                          chr = 1)
      write_tsv(as.data.frame(res), paste0("focus_", ts), cfg)
      TRUE
    }, error = function(e) {
      log_msg("FOCUS 运行失败（不同版本参数名可能不同，请以 ?FOCUS::focus 为准）：",
              conditionMessage(e))
      FALSE
    })
    if (ok) log_msg("FOCUS 完成：", ts)
  }
} else if (!requireNamespace("FOCUS", quietly = TRUE)) {
  log_msg("未安装 FOCUS，跳过精细定位。安装：remotes::install_github('kyle-gong/focus')")
}

log_msg("09 完成")

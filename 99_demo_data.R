# =============================================================================
# 99_demo_data.R —— 生成合成数据，用于在真实数据到位前跑通全链路
#
# 用途：
#   1. 验证环境、包版本、目录权限是否正常
#   2. 给你一个"正常结果长什么样"的参照（例如 coloc 的 PP.H4 应该接近 1）
#   3. 调试图表样式与表格字段
#
# 生成后运行：Rscript run_all.R --demo
# 注意：合成数据不能用于任何科学结论，只用于跑通流程。
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("99_demo_data", cfg)
load_pkgs("data.table")
set.seed(20261101)

out_dir <- ensure_dir("data", "processed")
raw_dir <- ensure_dir("data", "raw")

# ------------------------------------------------- 1) 基因集 -----------------
# 真实通路里的关键基因 + 填充基因，保证富集分析能出结果
key_genes <- c("LOX", "LOXL1", "LOXL2", "MMP2", "MMP9", "MMP14", "TIMP1", "TIMP2",
               "TGFB2", "TGFBR1", "TGFBR2", "YAP1", "WWTR1", "CTGF", "CYR61",
               "COL1A1", "COL3A1", "COL4A1", "FN1", "LAMA1", "LAMB1",
               "MYOC", "OPTN", "WDR36", "CAV1", "CAV2", "CDKN2B", "SIX6", "AFAP1",
               "ACTB", "GAPDH", "VIM", "CDH1", "RHO", "PECAM1", "ACTA2", "TAGLN",
               "IL6", "CXCL8", "CCL2", "NFKB1", "RELA", "STAT3", "JUN", "FOS")
fill <- sprintf("GENE%03d", seq_len(760))
genes <- c(key_genes, fill)

make_dataset <- function(n_case, n_ctrl, genes, effect_genes, n_batch = 2) {
  n <- n_case + n_ctrl
  base <- matrix(stats::rnorm(length(genes) * n, mean = 8, sd = 1.6),
                 nrow = length(genes), dimnames = list(genes, paste0("S", seq_len(n))))
  # 基因固有表达差异
  base <- base + stats::rnorm(length(genes), 0, 0.8)
  # 病例效应：只在 effect_genes 上，统一向上调
  # 注意必须转成无维度的向量：matrix(eff, ncol=1) 是 12x1 矩阵，
  # 加到 12 x n_case 的子矩阵上会报 non-conformable（R 只对向量做循环补齐）
  # 效应量 2.5：805 个基因同时检验时，1.2 的效应过不了 FDR<0.05（实测最小 0.09），
  # demo 就会在 02/05 连锁崩掉；合成数据要的是"有清晰信号"，不是真实效应大小
  eff <- as.numeric(stats::setNames(stats::rnorm(length(effect_genes), 2.5, 0.5), effect_genes))
  base[effect_genes, seq_len(n_case)] <- base[effect_genes, seq_len(n_case)] + eff
  # 批次效应（模拟不同扫描批次）
  batch <- rep(paste0("B", seq_len(n_batch)), length.out = n)
  if (n_batch > 1) {
    for (b in seq_len(n_batch)) {
      idx <- which(batch == paste0("B", b))
      base[, idx] <- base[, idx] + stats::rnorm(length(genes), 0, 0.6)
    }
  }
  list(expr = base, batch = batch, group = c(rep("case", n_case), rep("control", n_ctrl)))
}

# ------------------------------------------------- 2) 训练集 ---------------
eff_genes <- c("LOX", "LOXL1", "MMP2", "TGFB2", "YAP1", "CTGF", "COL1A1", "FN1",
               "TIMP1", "WWTR1", "COL3A1", "LAMA1")
d1 <- make_dataset(20, 20, genes, eff_genes, n_batch = 2)
pheno1 <- data.frame(sample = colnames(d1$expr), group = factor(d1$group, levels = c("control", "case")),
                     batch = d1$batch,
                     IOP = round(stats::rnorm(40, 17, 3) + ifelse(d1$group == "case", 5, 0), 1),
                     MD = round(stats::rnorm(40, -3, 2) - ifelse(d1$group == "case", 4, 0), 2),
                     CCT = round(stats::rnorm(40, 540, 30), 1),
                     age = round(stats::rnorm(40, 62, 9)),
                     stringsAsFactors = FALSE)
rownames(pheno1) <- pheno1$sample
saveRDS(d1$expr, file.path(out_dir, "DEMO_expr.rds"))
saveRDS(pheno1, file.path(out_dir, "DEMO_pheno.rds"))

# ------------------------------------------------- 3) 验证集 ---------------
# 用同样方向、较小效应，且只保留部分基因（模拟不同平台覆盖度）
d2 <- make_dataset(12, 12, genes, eff_genes, n_batch = 1)
pheno2 <- data.frame(sample = colnames(d2$expr), group = factor(d2$group, levels = c("control", "case")),
                     batch = NA_character_, IOP = stats::rnorm(24, 17, 3),
                     MD = stats::rnorm(24, -4, 2), CCT = stats::rnorm(24, 540, 30),
                     age = round(stats::rnorm(24, 60, 10)), stringsAsFactors = FALSE)
rownames(pheno2) <- pheno2$sample
saveRDS(d2$expr, file.path(out_dir, "DEMO_V_expr.rds"))
saveRDS(pheno2, file.path(out_dir, "DEMO_V_pheno.rds"))

saveRDS(data.frame(dataset = c("DEMO", "DEMO_V"), role = c("train", "valid"),
                   platform = "SYNTHETIC", n_sample = c(40, 24),
                   n_case = c(20, 12), n_control = c(20, 12),
                   n_gene = nrow(d1$expr), n_batch = c(2, 1), n_outlier = 0),
        file.path(out_dir, "demo_qc.rds"))

# ------------------------------------------------- 4) GWAS + eQTL ----------
# 设计一个 chr1 上的共享因果区域：同一批 SNP 同时影响 LOX 表达与 POAG 风险
n_snp <- 3000
chr <- 1L
pos <- sort(sample(500000:2000000, n_snp))
snp <- sprintf("rs%07d", seq_len(n_snp))
ref <- sample(c("A", "C", "G", "T"), n_snp, TRUE)
alt <- vapply(ref, function(r) sample(setdiff(c("A", "C", "G", "T"), r), 1), character(1))
eaf <- pmin(pmax(stats::rbeta(n_snp, 2, 5), 0.01), 0.99)

causal <- pos > 900000 & pos < 1000000          # 共享因果区间
# null 位点的 beta 必须以自身 se 为尺度（z ~ N(0,1)），否则 p 值全挤在 1 附近，
# λ_GC 会算成 0.03 这种一看就假的数。
# 注意：这里刻意保持「rnorm(n_snp) → rnorm(n_causal) → runif(n_snp)」的抽样顺序不变，
# 否则后面的表达矩阵会拿到不同的随机数，02–05 的既有结果就都得重跑。
gwas_z <- stats::rnorm(n_snp, 0, 1)
gwas_eff_causal <- stats::rnorm(sum(causal), 0.35, 0.05)
gwas_se <- stats::runif(n_snp, 0.055, 0.07)
gwas_beta <- gwas_z * gwas_se
gwas_beta[causal] <- gwas_eff_causal
gwas_p <- 2 * stats::pnorm(-abs(gwas_beta / gwas_se))

gwas_demo <- data.frame(`#chrom` = chr, pos = pos, ref = ref, alt = alt, rsids = snp,
                        af_alt = eaf, beta = gwas_beta, sebeta = gwas_se,
                        pval = gwas_p, N = 300000, check.names = FALSE)
data.table::fwrite(gwas_demo, file.path(raw_dir, "demo_gwas.tsv.gz"), sep = "\t", compress = "gzip")

# eQTL：LOX 与 GWAS 共享同一套 z（同一位点同一个因果变异 → 演示 PP.H4 高）
#       LOXL1 只有 1/4 效应（弱共定位）；MYOC 无信号（演示 H4 低、H1 高）
#
# 第一版这里给 LOX 的 beta 加了独立噪声 rnorm(0, 0.03)，结果两侧最显著 SNP
# 不是同一个（rs0000950 vs rs0000974），coloc 判 H3=0.98/H4=0.02。
# 统计上那是【正确】结论——合成 SNP 之间没有 LD，两个相邻独立信号就是两个
# 因果变异。但 fixture 的用途是验证"H4 高"这条代码路径，所以让强信号基因
# 与 GWAS 共用同一套 beta/se（z 完全一致），弱信号与无信号基因保持原设计。
eqtl_one <- function(gene, scale_factor, identical_z = FALSE) {
  b <- stats::rnorm(n_snp, 0, 0.05)
  se <- stats::runif(n_snp, 0.04, 0.06)
  b[causal] <- gwas_beta[causal] * scale_factor + stats::rnorm(sum(causal), 0, 0.03)
  if (identical_z) {                    # 共享因果变异：直接复制 GWAS 的 beta 与 se
    b[causal] <- gwas_beta[causal]
    se[causal] <- gwas_se[causal]
  }
  data.frame(SNP = snp, gene = gene, A1 = alt, A2 = ref, beta = b, se = se,
             p = 2 * stats::pnorm(-abs(b / se)), maf = pmin(eaf, 1 - eaf),
             chr = chr, pos = pos, stringsAsFactors = FALSE)
}
eqtl_demo <- rbind(eqtl_one("LOX", 1.0, TRUE), eqtl_one("LOXL1", 0.25), eqtl_one("MYOC", 0.0))
data.table::fwrite(eqtl_demo, file.path(raw_dir, "demo_eqtl.tsv.gz"), sep = "\t", compress = "gzip")

# MR 暴露：取 LOX 的最强 SNP 作为工具变量
# 剔除回文位点（A/T、C/G）：TwoSampleMR 的 harmonise 会以"链方向不明"为由
# 直接删除这类 SNP。第一版没剔除，最强的工具变量恰好是 A/T，全部被删掉，
# 10_mr.R 报"harmonise 后无可用 SNP"。合成数据不必要地复制这个现实麻烦。
top <- eqtl_demo[eqtl_demo$gene == "LOX", ]
palindromic <- (top$A1 == "A" & top$A2 == "T") | (top$A1 == "T" & top$A2 == "A") |
  (top$A1 == "C" & top$A2 == "G") | (top$A1 == "G" & top$A2 == "C")
top <- top[!palindromic, ]
top <- top[order(top$p), ][1:8, ]
top$p <- pmax(top$p, 1e-300)
exp_demo <- data.frame(SNP = top$SNP, gene = top$gene, tissue = "Artery",
                       beta = top$beta, se = top$se, ea = top$A1, nea = top$A2,
                       eaf = top$maf, p = top$p, pos = top$pos, stringsAsFactors = FALSE)
data.table::fwrite(exp_demo, file.path(raw_dir, "demo_mr_exposure.tsv.gz"), sep = "\t", compress = "gzip")

log_msg("合成数据已生成：")
log_msg("  data/processed/DEMO_expr.rds（40 例 × ", nrow(d1$expr), " 基因，2 个批次）")
log_msg("  data/processed/DEMO_V_expr.rds（24 例）")
log_msg("  data/raw/demo_gwas.tsv.gz（", n_snp, " 个 SNP，chr1 共享因果区间 900–1000 kb）")
log_msg("  data/raw/demo_eqtl.tsv.gz（LOX 强信号 / LOXL1 弱信号 / MYOC 无信号）")
log_msg("  data/raw/demo_mr_exposure.tsv.gz")
log_msg("下一步：Rscript run_all.R --demo")

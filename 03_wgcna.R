# =============================================================================
# 03_wgcna.R —— 课题一（3）：共表达网络与「表型连续变量」关联
#
# 输入：data/processed/<trainGSE>_expr.rds / _pheno.rds
# 输出：results/tables/wgcna_gene_modules.tsv     基因-模块归属 + kME
#       results/tables/wgcna_module_trait.tsv     模块-表型相关系数与 P
#       results/tables/wgcna_hub_genes.tsv        表型相关模块的 hub 基因
#       results/figures/wgcna_softthreshold.pdf
#       results/figures/wgcna_dendro.pdf
#       results/figures/wgcna_module_trait.pdf
#       data/processed/wgcna_<GSE>.rds
#
# 设计要点：
#   模块不只跟「病/非病」关联，还要跟 IOP、视野 MD、CCT、年龄这些连续变量关联。
#   二分类只能告诉你"和疾病有关"，连续变量才能告诉你"和哪个临床维度有关"，
#   这是审稿人判断你有没有真的理解疾病的地方。
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("03_wgcna", cfg)
load_pkgs(c("WGCNA", "ggplot2", "reshape2"))
set.seed(cfg_get(cfg, "project.seed", 1))

options(stringsAsFactors = FALSE)
WGCNA::allowWGCNAThreads(nThreads = max(1, parallel::detectCores() - 1))

gse <- if (isTRUE(cfg_get(cfg, "run.demo", FALSE))) "DEMO" else cfg_get(cfg, "bulk.train.gse", "GSE27276")
expr  <- readRDS(proj_path("data", "processed", paste0(gse, "_expr.rds")))
pheno <- readRDS(proj_path("data", "processed", paste0(gse, "_pheno.rds")))

# 转置成 WGCNA 要求的 样本 × 基因
datExpr <- t(expr)
storage.mode(datExpr) <- "double"

# ------------------------------------------------- 1) 样本与基因层面 QC ----
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  log_msg("剔除不合格基因 ", sum(!gsg$goodGenes), " 个、样本 ", sum(!gsg$goodSamples), " 个")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
if (nrow(datExpr) < 15) stop("样本量不足 15，WGCNA 结果不可靠；请改用 04 的通路打分策略。")

# ------------------------------------------------- 2) 软阈值 ----------------
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers,
                                networkType = cfg_get(cfg, "bulk.wgcna.network_type", "signed"),
                                verbose = 0)
fit_df <- data.frame(power = powers,
                     fit = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
                     mean_k = sft$fitIndices[, 5],
                     median_k = sft$fitIndices[, 6])

g_soft <- ggplot2::ggplot(fit_df, ggplot2::aes(power, fit)) +
  ggplot2::geom_text(ggplot2::aes(label = power), size = 2.8, colour = "grey25") +
  ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", colour = "#C0392B", linewidth = 0.4) +
  ggplot2::scale_y_continuous(breaks = seq(-0.4, 1, 0.2)) +
  ggplot2::labs(x = "soft threshold (power)",
                y = expression("scale free topology fit, signed "*R^2),
                title = paste0(gse, ": pickSoftThreshold")) +
  theme_pub()
save_fig(g_soft, paste0("wgcna_softthreshold_", gse), cfg, 5, 4)

power <- cfg_get(cfg, "bulk.wgcna.soft_power")
if (is.null(power)) {
  power <- sft$powerEstimate
  if (is.na(power)) {
    power <- fit_df$power[which.min(abs(fit_df$fit - 0.8))]
    log_msg("无 R^2>0.8 的幂次，退而选取最接近 0.8 的 power = ", power)
  }
}
log_msg("采用 soft power = ", power)

# ------------------------------------------------- 3) 建网与模块 -------------
net <- WGCNA::blockwiseModules(
  datExpr, power = power,
  # 20000 会让 13k 基因全部落进一个块：TOM 是 13k x 13k 的稠密矩阵，
  # 内存上到 3 GB 以上、实际耗时十几分钟且没有进度输出（实测在 12 分钟处仍卡住）。
  # WGCNA 官方对大矩阵的建议是 5000 分块，块内基因数小一个量级，总耗时反而更短。
  maxBlockSize = min(cfg_get(cfg, "bulk.wgcna.block_size", 5000), ncol(datExpr)),
  TOMType = cfg_get(cfg, "bulk.wgcna.tom_type", "signed"),
  networkType = cfg_get(cfg, "bulk.wgcna.network_type", "signed"),
  minModuleSize = cfg_get(cfg, "bulk.wgcna.min_module_size", 30),
  mergeCutHeight = cfg_get(cfg, "bulk.wgcna.merge_cut_height", 0.25),
  deepSplit = cfg_get(cfg, "bulk.wgcna.deep_split", 2),
  pamRespectsDendro = FALSE,
  numericLabels = FALSE,
  saveTOMs = FALSE,
  verbose = 0
)
module_color <- net$colors
names(module_color) <- colnames(datExpr)
me_raw <- net$MEs          # 保留 "ME<颜色>" 命名，signedKME 要自己从第 3 个字符开始取
me <- net$MEs
colnames(me) <- sub("^ME", "", colnames(me))
log_msg("识别模块 ", length(unique(module_color)), " 个（含 grey）")

# 树状图 + 模块颜色
f_dendro <- proj_path(cfg_get(cfg, "project.outdir", "results"), "figures",
                      paste0("wgcna_dendro_", gse, ".pdf"))
grDevices::pdf(f_dendro, width = 8, height = 4.5)
WGCNA::plotDendroAndColors(net$dendrograms[[1]],
                           module_color[net$blockGenes[[1]]],
                           "Module colors", dendroLabels = FALSE, hang = 0.03,
                           addGuide = TRUE, guideHang = 0.05, main = gse)
grDevices::dev.off()

# ------------------------------------------------- 4) 模块-表型关联 ----------
# 连续变量：配置里列出且实际存在的；二分类：group 编码 0/1
trait_list <- list()
for (nm in cfg_get(cfg, "bulk.wgcna.traits_numeric", list())) {
  if (nm %in% names(pheno) && sum(!is.na(pheno[[nm]])) >= 0.6 * nrow(pheno)) {
    trait_list[[nm]] <- as.numeric(pheno[[nm]])
  } else {
    log_msg("表型 ", nm, " 缺失或缺失率过高，跳过")
  }
}
for (nm in cfg_get(cfg, "bulk.wgcna.traits_binary", list())) {
  if (nm %in% names(pheno)) trait_list[[nm]] <- as.numeric(pheno[[nm]] == "case")
}
if (!length(trait_list)) stop("没有任何可用表型，无法做模块-表型关联")
traits <- as.data.frame(trait_list)
rownames(traits) <- pheno$sample
traits <- traits[rownames(datExpr), , drop = FALSE]

# WGCNA::cor 的 method 只接受 pearson/kendall/spearman；bicor 是独立函数
# （写成 cor(..., method="bicor") 会掉进 stats::cor 的 match.arg 报错）
modT <- WGCNA::bicor(me, traits, use = "pairwise.complete.obs")
modP <- WGCNA::corPvalueStudent(modT, nrow(datExpr))

# 长表输出：bicor 返回 行=模块(eigengene)、列=性状；
# 原来按 colnames(modT)（性状）去索引行，方向反了，必然 subscript out of bounds
long <- do.call(rbind, lapply(rownames(modT), function(mm) {
  data.frame(module = mm, trait = colnames(traits),
             r = as.numeric(modT[mm, ]), p = as.numeric(modP[mm, ]),
             n = nrow(datExpr))
}))
long <- long[long$module != "grey", ]
long$fdr <- stats::p.adjust(long$p, method = "BH")
long <- long[order(long$p), ]
write_tsv(long, "wgcna_module_trait", cfg)

# 热图：格子内标 r，括号内标 P
star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
df_hm <- data.frame(module = rep(rownames(modT), times = ncol(modT)),
                    trait = rep(colnames(modT), each = nrow(modT)),
                    r = as.vector(modT),
                    lab = paste0(sprintf("%.2f", as.vector(modT)), star(as.vector(modP))))
df_hm$module <- factor(df_hm$module, levels = rownames(modT))
df_hm$trait  <- factor(df_hm$trait, levels = colnames(modT))
g_hm <- ggplot2::ggplot(df_hm, ggplot2::aes(trait, module, fill = r)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
  ggplot2::geom_text(ggplot2::aes(label = lab), size = 2.6) +
  ggplot2::scale_fill_gradient2(low = "#1F6FB4", mid = "white", high = "#C0392B",
                                midpoint = 0, limits = c(-1, 1), name = "bicor") +
  ggplot2::labs(x = NULL, y = "module", title = paste0(gse, ": module-trait")) +
  theme_pub() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 panel.grid = ggplot2::element_blank())
save_fig(g_hm, paste0("wgcna_module_trait_", gse), cfg, 6, 0.35 * nrow(modT) + 2)

# ------------------------------------------------- 5) hub 基因 ---------------
# 必须传 me_raw（列名 "MEred"）：signedKME 内部写的是
# paste0("kME", substring(colnames(datME), 3))，传进来的列名若已是 "red"，
# 取第 3 个字符往下就只剩 "d"，于是所有列变成 kMEd / kMEown / kMEan，
# 下游按 kME<颜色> 取列全部落空，hub 表恒为空（模块相关性再强也救不回来）。
kME <- WGCNA::signedKME(datExpr, me_raw, corFnc = "bicor")
kME_df <- data.frame(gene = rownames(kME), kME, module = module_color[rownames(kME)])
write_tsv(kME_df, "wgcna_gene_modules", cfg)

# 选与任一表型 |r|>=0.4 且 P<0.05 的模块作为重点模块
focus <- unique(long$module[abs(long$r) >= 0.4 & long$p < 0.05])
log_msg("重点模块：", paste(focus, collapse = ", "))
hub <- do.call(rbind, lapply(focus, function(m) {
  col <- paste0("kME", m)
  if (!col %in% names(kME_df)) return(NULL)
  sub <- kME_df[kME_df$module == m, ]
  sub <- sub[order(-sub[[col]]), ]
  head(sub, 50)
}))
if (!is.null(hub)) {
  write_tsv(hub, "wgcna_hub_genes", cfg)
} else {
  log_msg("无满足条件的重点模块，hub 表为空")
}

saveRDS(list(datExpr = datExpr, module = module_color, ME = me, power = power,
             traits = traits, modT = modT, modP = modP, kME = kME),
        proj_path("data", "processed", paste0("wgcna_", gse, ".rds")))
log_msg("03 完成")

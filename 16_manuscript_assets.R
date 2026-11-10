# =============================================================================
# 16_manuscript_assets.R —— 把 results/ 里的碎片整理成投稿素材
#
# 这一步不是"再分析一次"，而是把 01–14 的产物按期刊要求重新组织：
#   - Table 1–6：xlsx（三线表样式）+ tsv（备份，便于再排版）
#   - figure_manifest.tsv：每张图的每个 panel 对应哪个脚本、哪个文件、是否已存在
#   - Figures/Fig2..Fig7.png：用 magick 把已有的单图拼成多 panel 图（可选）
#
# 为什么值得单独一步：
#   投稿被退回补材料时，80% 是"表 1 缺基线特征""图 3 的 panel 标号对不上"
#   这类问题。做成清单后，缺什么一眼能看到，不用回头翻 14 个脚本。
#
# 用法：
#   Rscript 16_manuscript_assets.R              # 出表 + 清单
#   Rscript 16_manuscript_assets.R --figures    # 额外拼图（需要 magick）
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("16_manuscript_assets", cfg)

args <- commandArgs(trailingOnly = TRUE)
do_figures <- "--figures" %in% args

out_dir <- ensure_dir(cfg_get(cfg, "project.outdir", "results"), "submission")
fig_dir <- ensure_dir(cfg_get(cfg, "project.outdir", "results"), "submission", "Figures")
load_pkgs(c("openxlsx", "data.table"), install = FALSE)

# ------------------------------------------------------------------ 工具 ----
pick <- function(df, cands, default = NA) {
  if (is.null(df) || !nrow(df)) return(rep(default, 0))
  hit <- intersect(cands, names(df))
  if (!length(hit)) return(rep(default, nrow(df)))
  df[[hit[1]]]
}

have_tbl <- function(name) file.exists(res_path(cfg, "tables", paste0(name, ".tsv")))

read_tbl <- function(name) {
  f <- res_path(cfg, "tables", paste0(name, ".tsv"))
  if (!file.exists(f)) { log_msg("缺 ", name, ".tsv，跳过"); return(NULL) }
  utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
}

num <- function(x) suppressWarnings(as.numeric(x))
msd <- function(x, d = 1) {              # mean ± sd，缺失返回 "—"
  x <- num(x); x <- x[is.finite(x)]
  if (!length(x)) return("—")
  paste0(formatC(mean(x), format = "f", digits = d), " ± ",
         formatC(stats::sd(x), format = "f", digits = d))
}
pct <- function(x) {
  x <- num(x); x <- x[!is.na(x)]
  if (!length(x)) return("—")
  paste0(formatC(100 * mean(x == 1 | x == "male" | x == "M"), format = "f", digits = 1), "%")
}

# --------------------------------------------------------------- Table 1 ----
# 队列基线。GEO 里年龄/性别/IOP/MD/CCT 多数缺失，缺就写"—"，不许编。
cohort_rows <- list()
add_cohort <- function(gse, role) {
  ph <- tryCatch(readRDS(proj_path("data", "processed", paste0(gse, "_pheno.rds"))),
                 error = function(e) NULL)
  if (is.null(ph)) { log_msg("无 pheno：", gse, "（先跑 01）"); return(invisible(NULL)) }
  g <- if ("group" %in% names(ph)) ph$group else rep(NA, nrow(ph))
  cohort_rows[[length(cohort_rows) + 1]] <<- data.frame(
    Cohort = gse, Role = role, N = nrow(ph),
    POAG = sum(g == "case", na.rm = TRUE),
    Control = sum(g == "control", na.rm = TRUE),
    Age = msd(pick(ph, c("age", "Age", "age_years"))),
    Male = pct(pick(ph, c("sex", "gender", "Sex"))),
    IOP = msd(pick(ph, c("IOP", "iop", "iop_mmHg"))),
    MD = msd(pick(ph, c("MD", "md", "visual_field_md"))),
    CCT = msd(pick(ph, c("CCT", "cct"))),
    stringsAsFactors = FALSE)
  invisible(NULL)
}
add_cohort(cfg_get(cfg, "bulk.train.gse", "GSE27276"), "Discovery")
for (v in cfg_get(cfg, "bulk.valid", list())) add_cohort(v$gse, "Validation")
if (!length(cohort_rows)) {
  log_msg("没有任何 pheno 文件，Table 1 只能出空表——先跑 01_geo_prepare.R")
  t1 <- data.frame(Cohort = character(0), Role = character(0))
} else t1 <- do.call(rbind, cohort_rows)

# --------------------------------------------------------------- Table 2 ----
# Top 20 差异基因（训练集）。排序规则与正文一致：先 padj，再 |logFC|。
gse <- cfg_get(cfg, "bulk.train.gse", "GSE27276")
deg <- read_tbl(paste0("deg_", gse))
t2 <- NULL
if (!is.null(deg)) {
  padj <- num(pick(deg, c("adj.P.Val", "padj", "fdr", "P.Value")))
  lfc  <- num(pick(deg, c("logFC", "log2FC")))
  t2 <- data.frame(
    Gene = pick(deg, c("gene", "Gene", "symbol", "probe")),
    logFC = round(lfc, 3),
    P = signif(num(pick(deg, c("P.Value", "pvalue", "p"))), 3),
    FDR = signif(padj, 3),
    AveExpr = round(num(pick(deg, c("AveExpr", "baseMean", "ave_expr"))), 3),
    stringsAsFactors = FALSE)
  t2 <- t2[order(t2$FDR, -abs(t2$logFC)), ][seq_len(min(20, nrow(t2))), ]
  t2$Direction <- ifelse(t2$logFC > 0, "Up", "Down")
  t2 <- t2[, c("Gene", "logFC", "AveExpr", "P", "FDR", "Direction")]
}

# --------------------------------------------------------------- Table 3 ----
# 模块—性状关联：只保留 |r| 最大的 15 行，正文按这个顺序写
mt <- read_tbl("wgcna_module_trait")
t3 <- NULL
if (!is.null(mt)) {
  r <- num(pick(mt, c("r", "cor")))
  t3 <- data.frame(Module = pick(mt, c("module", "Module")),
                   Trait = pick(mt, c("trait", "Trait")),
                   r = round(r, 3),
                   P = signif(num(pick(mt, c("p", "pvalue"))), 3),
                   FDR = signif(num(pick(mt, c("fdr", "FDR"))), 3),
                   nGenes = num(pick(mt, c("n_genes", "nGenes", "size"))),
                   stringsAsFactors = FALSE)
  t3 <- t3[order(-abs(t3$r)), ][seq_len(min(15, nrow(t3))), ]
}

# --------------------------------------------------------------- Table 4 ----
# 模型性能：训练 / 内部交叉验证 / 各外部队列。缺哪行就说明哪步没跑。
mp <- read_tbl("model_performance")
t4 <- NULL
if (!is.null(mp)) {
  # model_performance 里没有 model 列（05 只写 dataset/AUC/sens/spec），
  # 早先这里直接 pick() 得到 NULL，Table 4 的 Model 列整列是 NA——表里留一列空值很难看，
  # 所以缺列时统一填 logistic
  mdl <- pick(mp, c("model", "method"))
  if (is.null(mdl)) mdl <- rep("logistic", nrow(mp))
  t4 <- data.frame(
    Set = pick(mp, c("set", "cohort", "dataset", "Set")),
    Model = mdl,
    AUC = as.numeric(pick(mp, c("auc", "AUC"))),
    CI = ifelse(is.na(num(pick(mp, c("ci_low", "lower")))), "—",
                paste0("(", formatC(num(pick(mp, c("ci_low", "lower"))), 3), "–",
                       formatC(num(pick(mp, c("ci_high", "upper"))), 3), ")")),
    Sens = as.numeric(pick(mp, c("sensitivity", "sens"))),
    Spec = as.numeric(pick(mp, c("specificity", "spec"))),
    stringsAsFactors = FALSE)
  for (cc in c("AUC", "Sens", "Spec")) t4[[cc]] <- round(t4[[cc]], 3)
}

# bootstrap optimism 校正后的 AUC 单独占一行：审稿人看表就能判断过拟合程度，
# 不必回正文找那句话。B 用实际可用次数（n_ok），不是设定值。
mo <- read_tbl("model_optimism")
if (!is.null(mo) && !is.null(t4)) {
  ac <- num(pick(mo, c("auc_corrected", "auc_opt")))[1]
  nb <- num(pick(mo, c("n_ok", "n_boot")))[1]
  if (is.finite(ac)) {
    t4 <- rbind(t4, data.frame(
      Set = paste0("Discovery, optimism-corrected",
                   if (is.finite(nb)) paste0(" (B = ", nb, ")") else ""),
      Model = "logistic (fixed signature)",
      AUC = round(ac, 3), CI = "—", Sens = NA_real_, Spec = NA_real_,
      stringsAsFactors = FALSE))
  }
}

# --------------------------------------------------------------- Table 5 ----
# 共定位：PP.H4 由高到低；PP.H3 一起给（这是审稿人必看的）
cs <- read_tbl("coloc_summary")
t5 <- NULL
if (!is.null(cs)) {
  pp4 <- num(pick(cs, c("PP.H4", "PP.H4.abf", "pp4")))
  t5 <- data.frame(
    Gene = pick(cs, c("gene", "Gene")),
    Tissue = pick(cs, c("tissue", "Tissue")),
    nSNP = num(pick(cs, c("nSNP", "nsnp", "n_snp"))),
    PP.H0 = round(num(pick(cs, c("PP.H0", "PP.H0.abf"))), 3),
    PP.H1 = round(num(pick(cs, c("PP.H1", "PP.H1.abf"))), 3),
    PP.H2 = round(num(pick(cs, c("PP.H2", "PP.H2.abf"))), 3),
    PP.H3 = round(num(pick(cs, c("PP.H3", "PP.H3.abf"))), 3),
    PP.H4 = round(pp4, 3),
    LeadSNP = pick(cs, c("hit1", "top_snp", "SNP")),
    stringsAsFactors = FALSE)
  t5 <- t5[order(-t5$PP.H4), ]
  thr <- cfg_get(cfg, "genetics.coloc.pp4_threshold", 0.8)
  t5$Support <- ifelse(t5$PP.H4 >= thr & t5$PP.H3 < 0.5, "Shared causal variant",
                ifelse(t5$PP.H3 >= 0.5, "Distinct causal variants", "Inconclusive"))
}

# --------------------------------------------------------------- Table 6 ----
# MR 主表：只放 IVW（主分析），其余方法进补充表
mr <- read_tbl("mr_results")
t6 <- NULL
if (!is.null(mr)) {
  meth <- as.character(pick(mr, c("method", "Method")))
  keep <- grepl("Inverse variance|IVW", meth, ignore.case = TRUE)
  t6 <- data.frame(
    Exposure = pick(mr, c("exposure", "gene", "Gene")),
    Outcome = pick(mr, c("outcome", "Outcome")),
    Method = meth,
    nSNP = num(pick(mr, c("nsnp", "nSNP", "n_snp"))),
    b = round(num(pick(mr, c("b", "beta"))), 3),
    SE = round(num(pick(mr, c("se", "SE"))), 3),
    P = signif(num(pick(mr, c("pval", "p", "P"))), 3),
    OR = round(exp(num(pick(mr, c("b", "beta")))), 3),
    FDR = signif(num(pick(mr, c("fdr", "FDR"))), 3),
    stringsAsFactors = FALSE)[keep, ]
  if (nrow(t6)) t6 <- t6[order(t6$P), ][seq_len(min(30, nrow(t6))), ]
}

# ------------------------------------------------------- 写 xlsx / tsv ------
tabs <- list("Table1_cohort" = t1, "Table2_topDEG" = t2, "Table3_moduleTrait" = t3,
             "Table4_model" = t4, "Table5_coloc" = t5, "Table6_MR" = t6)

wb <- openxlsx::createWorkbook()
hdr <- openxlsx::createStyle(textDecoration = "bold", halign = "center",
                             border = "TopBottom", borderStyle = "medium")
body_style <- openxlsx::createStyle(border = "bottom", borderStyle = "thin")
for (nm in names(tabs)) {
  d <- tabs[[nm]]
  if (is.null(d) || !nrow(d)) {
    d <- data.frame(note = "该步骤尚未产出结果（对应脚本未运行或输入缺失）")
  }
  openxlsx::addWorksheet(wb, substr(nm, 1, 31))
  openxlsx::writeData(wb, substr(nm, 1, 31), d, headerStyle = hdr)
  openxlsx::addStyle(wb, substr(nm, 1, 31), body_style, rows = 2:(nrow(d) + 1),
                     cols = seq_len(ncol(d)), gridExpand = TRUE, stack = TRUE)
  openxlsx::setColWidths(wb, substr(nm, 1, 31), cols = seq_len(ncol(d)), widths = "auto")
  openxlsx::freezePane(wb, substr(nm, 1, 31), firstRow = TRUE)
  utils::write.table(d, file.path(out_dir, paste0(nm, ".tsv")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
}
openxlsx::saveWorkbook(wb, file.path(out_dir, "Tables1-6.xlsx"), overwrite = TRUE)
log_msg("表已写出：", file.path(out_dir, "Tables1-6.xlsx"))

# ---------------------------------------------------- figure manifest -------
# File 列写的是**各自脚本真正 save_fig 出去的文件名**：用 * 通配（<GSE> 这种随数据集变的
# 部分不再写死），一个 panel 有多个候选写法时用 | 分隔。
# 早先这里写的是"设想的文件名"（volcano_deg_<GSE>.png、wgcna_softpower.png…），
# 与脚本实际输出（volcano_DEMO.png、wgcna_softthreshold_DEMO.png…）对不上，
# 于是全都报"缺图"，这张清单就废了。
# Status 分四种，缺图时能区分"步骤没跑"和"根本没写绘图代码"：
#   ok / missing / no-code / hand-drawn
fig_dir_abs <- proj_path(cfg_get(cfg, "project.outdir", "results"), "figures")
resolve_fig <- function(pat) {
  if (identical(pat, "-")) return(NA_character_)
  for (p in strsplit(pat, "|", fixed = TRUE)[[1]]) {
    hit <- Sys.glob(file.path(fig_dir_abs, trimws(p)))
    if (length(hit)) return(basename(sort(hit)[1]))
  }
  NA_character_
}
fig_map <- data.frame(
  Figure = c(rep("Fig1", 3), rep("Fig2", 4), rep("Fig3", 3), rep("Fig4", 4),
             rep("Fig5", 4), rep("Fig6", 4), rep("Fig7", 3), rep("Fig8", 4),
             rep("Fig9", 3)),
  Panel  = c("A", "B", "C",  "A", "B", "C", "D",  "A", "B", "C",  "A", "B", "C", "D",
             "A", "B", "C", "D",  "A", "B", "C", "D",  "A", "B", "C",
             "A", "B", "C", "D",  "A", "B", "C"),
  Content = c(
    "Study design flow", "Data sources and sample sizes", "Analysis modules and decision rules",
    "PCA before/after ComBat", "Volcano (discovery)", "Top-50 heatmap", "Consensus gene overlap",
    "Soft-threshold fit", "Module dendrogram", "Module-trait heatmap (IOP highlighted)",
    "GO/KEGG dotplot", "GSVA Hallmark heatmap", "PROGENy pathway activity", "Cell-fraction boxplot",
    "LASSO CV curve", "SVM-RFE accuracy curve", "ROC (train + external)", "Calibration + DCA",
    "UMAP by cell type", "Marker dotplot", "CellChat differential interactions", "Pseudotime (TM/JCT)",
    "Cross-species concordance", "Human vs species logFC scatter", "Module score by condition",
    "Coloc PP4 by gene-tissue", "TWAS Manhattan", "MR forest (IVW)", "Mediation decomposition",
    "In silico perturbation score", "Permutation null distribution", "Multiverse AUC range"),
  Script = c(rep("hand-drawn", 3), rep("01", 4), rep("03", 3), rep("04", 4),
             rep("05", 4), rep("06", 4), rep("14", 3), rep("08/09/10/12", 4),
             rep("11/13", 3)),
  File = c(
    "-", "-", "-",
    "pca_*_raw.png|pca_*_clean.png", "volcano_*.png", "heatmap_top*.png", "venn_consensus.png",
    "wgcna_softthreshold_*", "wgcna_dendro_*", "wgcna_module_trait_*.png",
    "enrich_dotplot_*.png", "gsva_pathway_bar.png", "progeny_pathway.png", "deconv_box.png",
    "lasso_cv.png", "svmrfe_curve.png", "roc_train.png|roc_valid_*.png", "calibration.png",
    "sc_umap_celltype.png", "sc_dotplot_marker.png", "cellchat_compare*", "sc_pseudotime_TM.png",
    "xspecies_concordance_heatmap.png", "xspecies_logfc_scatter.png", "xspecies_modulescore.png",
    "coloc_locus_*.png", "twas_manhattan.png", "mr_forest.png", "mr_mediation_forest.png",
    "perturb_*.png", "robustness_null_auc.png", "robustness_multiverse.png"),
  stringsAsFactors = FALSE)

# 这几个 panel 目前没有绘图代码（只写了表，或分析本身没实现）：
# 补实现之前，别在正文写"Figure X? shows…"，否则图注和图对不上。
no_code <- c("Fig2D", "Fig4C", "Fig5B", "Fig9C")
fig_map$Found <- vapply(fig_map$File, resolve_fig, character(1))
fig_map$Exist <- !is.na(fig_map$Found)
fig_map$Path <- ifelse(fig_map$Exist,
                       file.path(cfg_get(cfg, "project.outdir", "results"), "figures",
                                 fig_map$Found), "-")
fig_map$Status <- ifelse(fig_map$File == "-", "hand-drawn",
                  ifelse(fig_map$Exist, "ok",
                  ifelse(paste0(fig_map$Figure, fig_map$Panel) %in% no_code, "no-code", "missing")))
utils::write.table(fig_map, file.path(out_dir, "figure_manifest.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
n_ok <- sum(fig_map$Status == "ok"); n_ms <- sum(fig_map$Status == "missing")
n_nc <- sum(fig_map$Status == "no-code"); n_hd <- sum(fig_map$Status == "hand-drawn")
log_msg("图清单：ok ", n_ok, " / ", nrow(fig_map), "；missing（步骤未跑）", n_ms,
        "；no-code（无绘图代码）", n_nc, "；手绘 ", n_hd)

# ------------------------------------------------------- 拼多 panel 图 ------
# 只有全部 panel 都在才拼，缺一个就不拼（避免出现半张图混进投稿包）。
assemble <- function(fig, cols = 2) {
  sub <- fig_map[fig_map$Figure == fig & fig_map$Exist, ]
  if (!nrow(sub) || !all(sub$Exist)) return(invisible(FALSE))
  if (!requireNamespace("magick", quietly = TRUE)) {
    log_msg("未安装 magick，跳过拼图：install.packages('magick')"); return(invisible(FALSE))
  }
  imgs <- lapply(seq_len(nrow(sub)), function(i) {
    im <- magick::image_read(proj_path(sub$Path[i]))
    magick::image_annotate(im, sub$Panel[i], size = 48, weight = 700,
                           color = "black", gravity = "northwest",
                           location = "+20+10")
  })
  rows <- split(imgs, ceiling(seq_along(imgs) / cols))
  strip <- lapply(rows, function(r) magick::image_append(do.call(c, r), stack = FALSE))
  img <- magick::image_append(do.call(c, strip), stack = TRUE)
  f <- file.path(fig_dir, paste0(fig, ".png"))
  magick::image_write(img, f, format = "png")
  log_msg("拼图完成：", f)
  invisible(TRUE)
}
if (do_figures) for (f in c("Fig2", "Fig3", "Fig4", "Fig5", "Fig6", "Fig7", "Fig8", "Fig9")) assemble(f)

# ---------------------------------------------------- 手写交付说明 --------
notes <- c(
  "# 投稿素材清单", "",
  paste0("生成时间：", format(Sys.time())),
  paste0("输出目录：", out_dir), "",
  "## 表", "",
  paste0("- Table1_cohort：", nrow(t1), " 行（队列基线）"),
  paste0("- Table2_topDEG：", if (is.null(t2)) 0 else nrow(t2), " 行（Top20 DEG）"),
  paste0("- Table3_moduleTrait：", if (is.null(t3)) 0 else nrow(t3), " 行（模块-性状）"),
  paste0("- Table4_model：", if (is.null(t4)) 0 else nrow(t4), " 行（模型性能）"),
  paste0("- Table5_coloc：", if (is.null(t5)) 0 else nrow(t5), " 行（共定位）"),
  paste0("- Table6_MR：", if (is.null(t6)) 0 else nrow(t6), " 行（MR IVW）"), "",
  "## 图", "",
  paste0("已有 panel：", sum(fig_map$Status == "ok"), " / ", nrow(fig_map),
         "（另有手绘 ", sum(fig_map$Status == "hand-drawn"), " 个）"),
  "", "## 缺什么", "",
  "- missing（绘图代码在、文件不在：多为该步没跑，或该步里这个 panel 被条件跳过）：",
  if (!any(fig_map$Status == "missing")) "  无" else
    paste0("  ", paste(paste0(fig_map$Figure, fig_map$Panel, "(", fig_map$Content, ")")[
      fig_map$Status == "missing"], collapse = "、")),
  "- no-code（根本没有绘图代码，需补实现或从图里删掉这个 panel）：",
  if (!any(fig_map$Status == "no-code")) "  无" else
    paste0("  ", paste(paste0(fig_map$Figure, fig_map$Panel, "(", fig_map$Content, ")")[
      fig_map$Status == "no-code"], collapse = "、")),
  "",
  "> 表格单元格里的 — 表示原始数据中该变量缺失（GEO 常见），",
  "> 不要在正文里补数字，宁可写 not available。"
)
writeLines(notes, file.path(out_dir, "README_assets.md"))
log_msg("完成。素材在 ", out_dir)

# =============================================================================
# 15_results_ledger.R —— 结果账本：把散落在 14 个步骤里的数字收成一张可审计的表
#
# 为什么需要这一步：
#   写论文时最容易出的错不是分析错，而是「摘要写 AUC 0.89、结果写 0.91」这类
#   数字漂移；投稿后被要求提供原始输出时又对不上。本脚本把所有进正文的数字
#   抽成 token（{{AUC_TRAIN}} 这种形式），正文只写 token，最后一次性替换。
#   好处有三：数字唯一来源、换数据集后正文自动跟着变、审稿人要原始值时能
#   直接给 source 列指到具体文件与行。
#
# 输出：
#   results/tables/results_ledger.tsv     token / display / raw / source / key / note
#   results/tables/results_ledger.md      分组可读版（贴给导师看）
#   results/tables/manuscript_token_map.tsv  两列，供 Word 查找替换
#   results/tables/missing_steps.tsv      哪些步骤没有产物（对应段落还写不了）
#   results/session_info.txt              R 与包版本（方法学可重复性段落要用）
#
# 用法：
#   Rscript 15_results_ledger.R
#   Rscript 15_results_ledger.R --check ../output/manuscript/draft.html
#       扫描稿件里所有 {{TOKEN}}，列出「账本里没有」和「账本里有但没用到」的
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("15_results_ledger", cfg)

args <- commandArgs(trailingOnly = TRUE)
check_html <- NULL
i <- which(args == "--check")
if (length(i) && i < length(args)) check_html <- args[i + 1]

# ------------------------------------------------------------------ 工具 ----

# 宽容取列：不同版本脚本的列名可能略有差别，按候选名找第一个存在的
pick <- function(df, cands) {
  if (is.null(df) || !nrow(df)) return(NULL)
  hit <- intersect(cands, names(df))
  if (!length(hit)) return(NULL)
  df[[hit[1]]]
}

# 读表；不存在则记入 missing，不报错（缺步骤是常态，不是故障）
MISSING <- character(0)
tbl <- function(name) {
  f <- res_path(cfg, "tables", paste0(name, ".tsv"))
  if (!file.exists(f)) { MISSING <<- c(MISSING, name); return(NULL) }
  d <- tryCatch(utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d)) { MISSING <<- c(MISSING, name); return(NULL) }
  attr(d, "file") <- basename(f)
  d
}

# 账本累积器
LEDGER <- list()
TK <- function(token, display, raw, source, key = "", note = "") {
  LEDGER[[length(LEDGER) + 1]] <<- data.frame(
    token = token, display = as.character(display), raw = as.character(raw),
    source = source, key = key, note = note, stringsAsFactors = FALSE)
}

fmt <- function(x, d = 2) {
  if (length(x) != 1 || !is.finite(x)) return(NA_character_)
  formatC(x, format = "f", digits = d)
}
fmt_p <- function(p) {
  if (length(p) != 1 || !is.finite(p)) return(NA_character_)
  if (p < 0.001) "< 0.001" else paste0("= ", formatC(p, format = "f", digits = 3))
}
src_of <- function(d) attr(d, "file") %||% "?"

# ============================================================ 1) 队列与样本 ==
train_gse <- cfg_get(cfg, "bulk.train.gse", "GSE27276")
ph <- tryCatch(readRDS(proj_path("data", "processed", paste0(train_gse, "_pheno.rds"))),
               error = function(e) NULL)
if (!is.null(ph) && "group" %in% names(ph)) {
  g <- ph$group
  n_all <- nrow(ph); n_case <- sum(g == "case"); n_ctrl <- sum(g == "control")
  TK("N_TRAIN", n_all, n_all, paste0(train_gse, "_pheno.rds"), train_gse, "训练集样本数")
  TK("N_TRAIN_CASE", n_case, n_case, paste0(train_gse, "_pheno.rds"), train_gse, "")
  TK("N_TRAIN_CTRL", n_ctrl, n_ctrl, paste0(train_gse, "_pheno.rds"), train_gse, "")
  TK("TRAIN_GSE", train_gse, train_gse, "config", train_gse, "")
  TK("TRAIN_RATIO", paste0(n_case, "∶", n_ctrl), paste0(n_case, "/", n_ctrl),
     paste0(train_gse, "_pheno.rds"), train_gse, "")
}

# 验证集：逐个读 pheno，拼成一句话
valid_cfg <- cfg_get(cfg, "bulk.valid", list())
valid_txt <- character(0); n_valid_total <- 0
for (v in valid_cfg) {
  gse <- v$gse
  pv <- tryCatch(readRDS(proj_path("data", "processed", paste0(gse, "_pheno.rds"))),
                 error = function(e) NULL)
  if (is.null(pv) || !"group" %in% names(pv)) next
  nc <- sum(pv$group == "case"); nk <- sum(pv$group == "control")
  n_valid_total <- n_valid_total + nrow(pv)
  valid_txt <- c(valid_txt, paste0(gse, " (", nc, " POAG / ", nk, " controls, n = ", nrow(pv), ")"))
}
if (length(valid_txt)) {
  TK("VALID_DESC", paste(valid_txt, collapse = "; "), paste(valid_txt, collapse = ";"),
     "data/processed/*_pheno.rds", "", "验证集描述，直接写进方法与结果")
  TK("N_VALID_TOTAL", n_valid_total, n_valid_total, "data/processed/*_pheno.rds", "", "")
}

# ============================================================ 2) 差异表达 ==
deg_tr <- tbl(paste0("deg_", train_gse))
cons   <- tbl("deg_consensus")

if (!is.null(deg_tr)) {
  sig_col <- pick(deg_tr, c("sig", "Sig", "significance"))
  lfc <- pick(deg_tr, c("logFC", "log2FC")); pv <- pick(deg_tr, c("adj.P.Val", "padj", "fdr"))
  if (!is.null(sig_col)) {
    is_sig <- sig_col != "NS" & !is.na(sig_col)
  } else {
    is_sig <- !is.na(pv) & pv < 0.05 & abs(lfc) > 1
  }
  n_up <- sum(is_sig & lfc > 0); n_dn <- sum(is_sig & lfc < 0)
  TK("N_DEG_TRAIN", sum(is_sig), sum(is_sig), src_of(deg_tr), train_gse, "训练集差异基因数")
  TK("N_DEG_UP", n_up, n_up, src_of(deg_tr), train_gse, "")
  TK("N_DEG_DOWN", n_dn, n_dn, src_of(deg_tr), train_gse, "")
  TK("N_GENE_TESTED", nrow(deg_tr), nrow(deg_tr), src_of(deg_tr), train_gse, "纳入检验的基因数")
}

if (!is.null(cons)) {
  lfc_c <- as.numeric(pick(cons, c("logFC_train", "logFC")))
  gene_c <- pick(cons, c("gene", "Gene", "symbol"))
  TK("N_CONSENSUS", nrow(cons), nrow(cons), src_of(cons), "", "方向一致的共识基因数")
  ord <- order(-abs(lfc_c))
  for (k in 1:3) {
    if (length(ord) >= k && !is.na(ord[k])) {
      TK(paste0("CONS_GENE_", k), gene_c[ord[k]], gene_c[ord[k]], src_of(cons), gene_c[ord[k]],
         paste0("共识基因按 |logFC| 排序第 ", k))
      TK(paste0("CONS_LOGFC_", k), fmt(lfc_c[ord[k]], 2), lfc_c[ord[k]], src_of(cons), gene_c[ord[k]])
    }
  }
  # 假说相关基因是否在列：这一条决定讨论能不能往 ECM/YAP 上引
  key_genes <- c("LOX", "LOXL1", "TGFB2", "COL1A1", "COL1A2", "MMP2", "MMP9", "TIMP1",
                 "YAP1", "WWTR1", "FN1", "SPP1", "MYOC")
  hit <- intersect(toupper(gene_c), key_genes)
  TK("KEY_GENES_HIT", if (length(hit)) paste(hit, collapse = ", ") else "none",
     paste(hit, collapse = ","), src_of(cons), "", "预设假说相关基因在共识表中的命中情况")
  TK("N_KEY_GENES_HIT", length(hit), length(hit), src_of(cons), "", "")
  # 逐基因 logFC，供讨论里点名时引用
  for (gn in hit) {
    j <- match(gn, toupper(gene_c))
    TK(paste0("LOGFC_", gn), fmt(lfc_c[j], 2), lfc_c[j], src_of(cons), gene_c[j], "")
  }
}

# 验证集方向一致比例（跨队列稳健性的核心数字）
cons_files <- list.files(res_path(cfg, "tables"), pattern = "^consensus_vs_", full.names = FALSE)
if (length(cons_files) == 0 && !is.null(cons)) {
  vcols <- grep("^logFC_", names(cons), value = TRUE)
  if (length(vcols)) {
    rates <- sapply(vcols, function(cc) {
      v <- as.numeric(cons[[cc]]); ok <- !is.na(v)
      if (!sum(ok)) return(NA_real_)
      mean(sign(v[ok]) == sign(as.numeric(pick(cons, c("logFC_train"))[ok])))
    })
    for (nm in names(rates)) {
      tk <- paste0("CONCORD_", sub("^logFC_", "", nm))
      TK(tk, paste0(fmt(100 * rates[[nm]], 1), "%"), rates[[nm]], src_of(cons), nm,
         "该验证集中与训练集方向一致的共识基因比例")
    }
  }
}

# ============================================================ 3) WGCNA ======
mt <- tbl("wgcna_module_trait")
if (!is.null(mt)) {
  r <- as.numeric(pick(mt, c("r", "cor"))); p <- as.numeric(pick(mt, c("p", "pvalue")))
  mod <- pick(mt, c("module", "Module")); tr <- pick(mt, c("trait", "Trait"))
  fdr <- pick(mt, c("fdr", "FDR"))
  n_mod <- length(unique(mod))
  TK("N_MODULES", n_mod, n_mod, src_of(mt), "", "非 grey 模块数")
  j <- which.max(abs(r))
  if (length(j)) {
    TK("TOP_MODULE", mod[j], mod[j], src_of(mt), paste(mod[j], tr[j]), "与性状相关最强的模块")
    TK("TOP_MODULE_TRAIT", tr[j], tr[j], src_of(mt), paste(mod[j], tr[j]))
    TK("TOP_MODULE_R", fmt(r[j], 2), r[j], src_of(mt), paste(mod[j], tr[j]))
    TK("TOP_MODULE_P", fmt_p(p[j]), p[j], src_of(mt), paste(mod[j], tr[j]))
    if (!is.null(fdr)) TK("TOP_MODULE_FDR", fmt_p(fdr[j]), fdr[j], src_of(mt), paste(mod[j], tr[j]))
  }
  # 与 IOP 相关最强的模块（这是本方案区别于常规 WGCNA 的关键：连续变量而非二分类）
  jj <- which(grepl("IOP|iop", tr))
  if (length(jj)) {
    k <- jj[which.max(abs(r[jj]))]
    TK("IOP_MODULE", mod[k], mod[k], src_of(mt), paste(mod[k], tr[k]))
    TK("IOP_MODULE_R", fmt(r[k], 2), r[k], src_of(mt), paste(mod[k], tr[k]))
    TK("IOP_MODULE_P", fmt_p(p[k]), p[k], src_of(mt), paste(mod[k], tr[k]))
  }
}

hub <- tbl("wgcna_hub_genes")
if (!is.null(hub)) TK("N_HUB", nrow(hub), nrow(hub), src_of(hub), "", "hub 基因总数")

# ============================================================ 4) 建模 ======
cand <- tbl("ml_candidate_genes")
if (!is.null(cand)) TK("N_CANDIDATE", nrow(cand), nrow(cand), src_of(cand), "",
                       "DEG ∩ WGCNA 候选基因数")

sig_tab <- tbl("ml_signature")
if (!is.null(sig_tab)) {
  gn_all <- pick(sig_tab, c("gene", "Gene", "symbol"))
  # summary(glm)$coefficients 第一行是 (Intercept)，它不是基因；
  # 若旧结果没有 is_intercept 列就按名字兜底过滤
  is_int <- if ("is_intercept" %in% names(sig_tab)) {
    as.logical(sig_tab$is_intercept)
  } else {
    grepl("^\\(?Intercept\\)?$", gn_all, ignore.case = TRUE)
  }
  is_int[is.na(is_int)] <- FALSE
  gn <- gsub("`", "", gn_all[!is_int & nzchar(trimws(gn_all))])
  TK("N_SIGNATURE", length(gn), length(gn), src_of(sig_tab), "", "最终 signature 基因数")
  TK("SIG_GENES", paste(gn, collapse = ", "), paste(gn, collapse = ","), src_of(sig_tab), "",
     "模型里实际进入的基因（不含截距）")
}

perf <- tbl("model_performance")
if (!is.null(perf)) {
  ds <- pick(perf, c("dataset", "cohort", "set"))
  auc <- as.numeric(pick(perf, c("AUC", "auc")))
  lo <- as.numeric(pick(perf, c("AUC_low", "low"))); hi <- as.numeric(pick(perf, c("AUC_high", "high")))
  se <- as.numeric(pick(perf, c("sensitivity", "sens"))); sp <- as.numeric(pick(perf, c("specificity", "spec")))
  # 第一行按脚本约定是训练集
  TK("AUC_TRAIN", fmt(auc[1], 3), auc[1], src_of(perf), ds[1])
  # 正文写的是 "AUC of {{AUC_TRAIN}} (95% CI {{AUC_TRAIN_CI}})"，
  # 所以这里只给区间本身，不要再带 "1.000（95%CI …）"，否则会出现 "(95% CI 1.000（95%CI 1.000–1.000）)"
  TK("AUC_TRAIN_CI", paste0(fmt(lo[1], 3), "–", fmt(hi[1], 3)),
     paste(auc[1], lo[1], hi[1], sep = "/"), src_of(perf), ds[1], "训练集 AUC 的 95% CI 区间")
  TK("SENS_TRAIN", fmt(100 * se[1], 1), se[1], src_of(perf), ds[1], "最佳阈值下敏感度（%）")
  TK("SPEC_TRAIN", fmt(100 * sp[1], 1), sp[1], src_of(perf), ds[1], "最佳阈值下特异度（%）")
  if (nrow(perf) > 1) {
    for (k in 2:nrow(perf)) {
      tk <- paste0("AUC_VALID_", ds[k])
      TK(tk, fmt(auc[k], 3), auc[k], src_of(perf), ds[k], "外部验证集 AUC")
      TK(paste0(tk, "_CI"), paste0(fmt(lo[k], 3), "–", fmt(hi[k], 3)),
         paste(auc[k], lo[k], hi[k], sep = "/"), src_of(perf), ds[k], "该验证集 AUC 的 95% CI 区间")
    }
    TK("AUC_VALID_RANGE", paste0(fmt(min(auc[-1], na.rm = TRUE), 3), "–",
                                 fmt(max(auc[-1], na.rm = TRUE), 3)),
       paste(min(auc[-1], na.rm = TRUE), max(auc[-1], na.rm = TRUE), sep = "/"),
       src_of(perf), "", "外部验证 AUC 区间（审稿人爱问上下限）")
    TK("N_VALID_SETS", nrow(perf) - 1, nrow(perf) - 1, src_of(perf), "")
  }
}

rob <- tbl("robustness_permutation")
if (!is.null(rob)) {
  met <- pick(rob, c("metric", "指标"))
  pe <- as.numeric(pick(rob, c("p_empirical", "p")))
  nm <- as.numeric(pick(rob, c("null_mean", "null")))
  q95 <- as.numeric(pick(rob, c("null_q95", "q95")))
  for (k in seq_len(nrow(rob))) {
    lab <- if (!is.null(met)) as.character(met[k]) else as.character(k)
    tag <- if (grepl("auc", lab, ignore.case = TRUE)) "AUC" else "DEG"
    TK(paste0("PERM_P_", tag), fmt_p(pe[k]), pe[k], src_of(rob), lab, "标签置换经验 P 值")
    TK(paste0("NULL_", tag, "_MEAN"), fmt(nm[k], 3), nm[k], src_of(rob), lab, "零分布均值")
    TK(paste0("NULL_", tag, "_Q95"), fmt(q95[k], 3), q95[k], src_of(rob), lab, "零分布 95 分位")
  }
}

mv <- tbl("robustness_multiverse")
if (!is.null(mv)) {
  a <- as.numeric(pick(mv, c("auc_mean", "auc", "AUC")))
  a <- a[is.finite(a)]
  TK("N_MULTIVERSE", nrow(mv), nrow(mv), src_of(mv), "", "分析管线组合数")
  if (length(a)) {
    TK("AUC_MV_RANGE", paste0(fmt(min(a), 3), "–", fmt(max(a), 3)),
       paste(min(a), max(a), sep = "/"), src_of(mv), "",
       "不同合理管线下的 AUC 区间")
  }
}

# ============================================================ 5) 遗传学 ====
gq <- tbl("gwas_qc_summary")
if (!is.null(gq)) {
  step <- pick(gq, c("step")); n <- as.numeric(pick(gq, c("n")))
  if (!is.null(step) && !is.null(n)) {
    TK("N_GWAS_RAW", format(n[1], big.mark = ","), n[1], src_of(gq), step[1], "GWAS 原始 SNP 数")
    TK("N_GWAS_CLEAN", format(n[nrow(gq)], big.mark = ","), n[nrow(gq)], src_of(gq), step[nrow(gq)],
       "质控后 SNP 数")
  }
}
lam <- tbl("gwas_lambda")
if (is.null(lam)) {
  # 07 把 λ 存在 data/processed/gwas_qc_meta.rds 里（没有单独的表），从这儿兜底取
  lam_rds <- proj_path("data", "processed", "gwas_qc_meta.rds")
  if (file.exists(lam_rds)) {
    meta <- tryCatch(readRDS(lam_rds), error = function(e) NULL)
    if (!is.null(meta) && is.finite(meta$lambda %||% NA_real_)) {
      lam <- data.frame(lambda = as.numeric(meta$lambda))
    }
  }
}
if (!is.null(lam)) TK("GWAS_LAMBDA", fmt(as.numeric(lam[[1]])[1], 3), as.numeric(lam[[1]])[1],
                      "gwas_qc_meta.rds", "", "基因组控制系数 λ（QC 后，MAF>=1%）")

col <- tbl("coloc_summary")
if (!is.null(col)) {
  pp4 <- as.numeric(pick(col, c("PP.H4", "PP4"))); pp3 <- as.numeric(pick(col, c("PP.H3", "PP3")))
  gn <- pick(col, c("symbol", "gene"))   # 有符号用符号，没有才退回 ENSG
  ts <- pick(col, c("tissue")); snp <- pick(col, c("top_snp", "lead_snp"))
  thr <- cfg_get(cfg, "genetics.coloc.pp4_threshold", 0.8)
  n_strong <- sum(pp4 > thr, na.rm = TRUE)
  TK("N_COLOC_PAIRS", nrow(col), nrow(col), src_of(col), "", "检验的 gene×tissue 对数")
  TK("N_COLOC_STRONG", n_strong, n_strong, src_of(col), "", paste0("PP.H4 > ", thr, " 的对数"))
  j <- which.max(pp4)
  if (length(j)) {
    TK("COLOC_GENE", gn[j], gn[j], src_of(col), paste(gn[j], ts[j]), "PP.H4 最高的基因")
    TK("COLOC_TISSUE", ts[j], ts[j], src_of(col), paste(gn[j], ts[j]))
    TK("COLOC_PP4", fmt(pp4[j], 3), pp4[j], src_of(col), paste(gn[j], ts[j]))
    TK("COLOC_PP3", fmt(pp3[j], 3), pp3[j], src_of(col), paste(gn[j], ts[j]),
       "同一对的 PP.H3：低才排除「两个独立信号」")
    if (!is.null(snp)) TK("COLOC_TOPSNP", snp[j], snp[j], src_of(col), paste(gn[j], ts[j]))
  }
  # 列出所有通过阈值的基因，供讨论点名
  if (n_strong > 0) {
    strong <- col[which(pp4 > thr), ]
    sgn <- if (!is.null(gn)) gn[which(pp4 > thr)] else strong$gene
    # 有 15/150 个 ENSG 在 org.Hs.eg.db 里没有 SYMBOL，别把 NA 写进正文
    bad <- is.na(sgn) | !nzchar(sgn)
    if (any(bad)) sgn[bad] <- strong$gene[bad]
    TK("COLOC_STRONG_LIST", paste(unique(sgn), collapse = ", "),
       paste(unique(sgn), collapse = ","), src_of(col), "", "通过阈值的基因清单")
  }
}

twas <- tbl("twas_significant")
if (!is.null(twas)) {
  gn <- pick(twas, c("ID", "gene", "GENE"))
  TK("N_TWAS_SIG", nrow(twas), nrow(twas), src_of(twas), "", "Bonferroni 后显著基因数")
  if (!is.null(gn)) TK("TWAS_GENES", paste(head(unique(gn), 10), collapse = ", "),
                       paste(unique(gn), collapse = ","), src_of(twas), "")
}
twas_all <- tbl("twas_all")
if (!is.null(twas_all)) {
  TK("N_TWAS_TESTED", nrow(twas_all), nrow(twas_all), src_of(twas_all), "", "TWAS 检验基因数")
  # 与 coloc 取交集：双支持的基因，是最硬的证据
  if (!is.null(col) && !is.null(twas)) {
    g1 <- unique(as.character(pick(col, c("gene"))[pp4 > cfg_get(cfg, "genetics.coloc.pp4_threshold", 0.8)]))
    g2 <- unique(as.character(pick(twas, c("ID", "gene", "GENE"))))
    both <- intersect(g1, g2)
    TK("N_DOUBLE_SUPPORT", length(both), length(both), "coloc ∩ twas",
       paste(both, collapse = ","), "共定位与 TWAS 双支持的基因数（G4 门禁）")
    if (length(both)) TK("DOUBLE_GENES", paste(both, collapse = ", "), paste(both, collapse = ","),
                         "coloc ∩ twas", "")
    else TK("DOUBLE_GENES", "none", "none", "coloc ∩ twas", "",
            "本次运行中两组证据没有共同基因；正文括注写 none 而不是留空")
  }
}

mr <- tbl("mr_results")
if (!is.null(mr)) {
  gn <- pick(mr, c("gene", "exposure")); meth <- pick(mr, c("method"))
  b <- as.numeric(pick(mr, c("b", "beta"))); p <- as.numeric(pick(mr, c("pval", "p")))
  ns <- as.numeric(pick(mr, c("n_snp", "nsnp")))
  is_ivw <- grepl("Inverse variance", meth, ignore.case = TRUE)
  sub <- if (any(is_ivw, na.rm = TRUE)) which(is_ivw) else seq_len(nrow(mr))
  j <- sub[which.min(p[sub])]
  TK("N_MR_GENES", length(unique(gn)), length(unique(gn)), src_of(mr), "", "纳入 MR 的基因数")
  TK("MR_TOP_GENE", gn[j], gn[j], src_of(mr), paste(gn[j], meth[j]))
  TK("MR_TOP_METHOD", meth[j], meth[j], src_of(mr), paste(gn[j], meth[j]))
  TK("MR_TOP_B", fmt(b[j], 3), b[j], src_of(mr), paste(gn[j], meth[j]), "效应量（log OR 或 SD 单位）")
  TK("MR_TOP_P", fmt_p(p[j]), p[j], src_of(mr), paste(gn[j], meth[j]))
  if (!is.null(ns)) TK("MR_TOP_NSNP", ns[j], ns[j], src_of(mr), paste(gn[j], meth[j]))
  # 10_mr.R 的 fdr 列在全是 Wald ratio 时会留空；用同一张表的 pval 现场做 BH，
  # 保证 MR_TOP_FDR 是可复现的真实数字而不是 NA。
  fdr <- tryCatch(as.numeric(pick(mr, c("fdr", "FDR"))), error = function(e) NULL)
  if (is.null(fdr) || !length(fdr) || all(is.na(fdr))) fdr <- stats::p.adjust(p, method = "BH")
  TK("MR_TOP_FDR", fmt_p(fdr[j]), fdr[j], src_of(mr), paste(gn[j], meth[j]))
  TK("N_MR_SIG", sum(as.numeric(p) < 0.05, na.rm = TRUE), sum(as.numeric(p) < 0.05, na.rm = TRUE),
     src_of(mr), "", "全部方法中 P<0.05 的行数（仅描述用，正式结论看 IVW+FDR）")
}

inst <- tbl("mr_instruments")
if (!is.null(inst)) {
  f <- as.numeric(pick(inst, c("F", "F_stat")))
  if (length(f) && any(is.finite(f))) TK("MR_F_MIN", fmt(min(f, na.rm = TRUE), 1), min(f, na.rm = TRUE),
                                         src_of(inst), "", "最小工具变量 F 统计量（须 > 10）")
}

med <- tbl("mr_mediation")
if (!is.null(med)) {
  gn <- pick(med, c("gene")); prop <- as.numeric(pick(med, c("prop_mediated", "prop")))
  ind <- as.numeric(pick(med, c("beta_indirect", "indirect", "ind_eff")))
  dir <- as.numeric(pick(med, c("beta_direct_mvmr", "direct", "dir_eff")))
  pp <- as.numeric(pick(med, c("p_indirect", "p")))
  j <- which.min(pp)
  if (length(j) && is.finite(prop[j])) {
    TK("MED_GENE", gn[j], gn[j], src_of(med), gn[j])
    TK("MED_PROP", paste0(fmt(100 * prop[j], 1), "%"), prop[j], src_of(med), gn[j], "中介比例（%）")
    TK("MED_INDIRECT", fmt(ind[j], 3), ind[j], src_of(med), gn[j])
    TK("MED_DIRECT", if (is.finite(dir[j])) fmt(dir[j], 3) else "not estimable (single instrument)",
       dir[j], src_of(med), gn[j], "调整 IOP 后的直接效应；单工具变量时无法用 MVMR 估计")
    TK("MED_P", fmt_p(pp[j]), pp[j], src_of(med), gn[j])
    flag <- if (prop[j] > 1 || prop[j] < 0) "direction-inconsistent; should not be read as mediation" else "consistent with partial mediation"
    if (!is.finite(dir[j]))
      flag <- paste0(flag, "; direct effect not estimable with a single instrument, so an IOP-independent route cannot be asserted")
    TK("MED_INTERP", flag, flag, src_of(med), gn[j], "措辞边界")
  }
}

# ============================================================ 6) 单细胞 ====
prop <- tbl("sc_celltype_prop")
if (!is.null(prop)) {
  # 06 写出的是宽表（细胞类型作列名、sample 作一列），不是 cell_type 长表；
  # 两种都要认，否则 N_CELLTYPES 会算成 0。
  ct <- tryCatch(pick(prop, c("cell_type", "celltype")), error = function(e) NULL)
  if (is.null(ct) || !length(ct)) ct <- setdiff(names(prop), "sample")
  ct <- unique(as.character(ct))
  TK("N_CELLTYPES", length(ct), length(ct), src_of(prop), "", "注释到的细胞类型数")
  TK("CELLTYPE_LIST", paste(ct, collapse = ", "), paste(ct, collapse = ","), src_of(prop), "")
}
# 细胞总数是注释对象的事实，不是比例表能推出来的：
# tools/write_sc_counts.R 从 scRNA_ref.rds 的 meta.data 生成 sc_celltype_counts.tsv。
cnt <- tbl("sc_celltype_counts")
if (!is.null(cnt)) {
  ct2 <- unique(as.character(pick(cnt, "cell_type")))
  if (length(ct2)) {
    TK("N_CELLTYPES", length(ct2), length(ct2), src_of(cnt), "", "注释到的细胞类型数")
    TK("CELLTYPE_LIST", paste(ct2, collapse = ", "), paste(ct2, collapse = ","), src_of(cnt), "")
  }
  ncell <- as.numeric(pick(cnt, "n"))
  nc <- sum(ncell, na.rm = TRUE)
  if (nc > 0) TK("N_CELLS", format(nc, big.mark = ","), nc, src_of(cnt), "", "单细胞总数")
}

pert <- tbl("perturb_summary")
if (!is.null(pert)) {
  gn <- pick(pert, c("gene")); ct <- pick(pert, c("celltype", "cell_type"))
  sc <- as.numeric(pick(pert, c("score"))); pp <- as.numeric(pick(pert, c("p_perm", "p")))
  j <- which.min(pp)
  TK("N_PERT", nrow(pert), nrow(pert), src_of(pert), "", "完成的 基因×细胞类型 扰动组合数")
  TK("PERT_TOP_GENE", gn[j], gn[j], src_of(pert), paste(gn[j], ct[j]))
  TK("PERT_TOP_CT", ct[j], ct[j], src_of(pert), paste(gn[j], ct[j]))
  TK("PERT_TOP_SCORE", fmt(sc[j], 3), sc[j], src_of(pert), paste(gn[j], ct[j]),
     "负值表示敲低该基因可逆转疾病表达程序")
  TK("PERT_TOP_P", fmt_p(pp[j]), pp[j], src_of(pert), paste(gn[j], ct[j]))
  # 11 的表没有 direction 列；方向的唯一依据就是 score 的符号
  # （负值 = 敲低后可逆转疾病表达程序），据此生成正文需要的动词。
  dirn <- tryCatch(pick(pert, c("direction", "dir")), error = function(e) NULL)
  if (is.null(dirn) || !length(dirn))
    dirn <- ifelse(sc < 0, "reverse", "reinforce")
  TK("PERT_TOP_DIR", dirn[j], dirn[j], src_of(pert), paste(gn[j], ct[j]))
  n_sig <- sum(pp < 0.05, na.rm = TRUE)
  TK("N_PERT_SIG", n_sig, n_sig, src_of(pert), "", "置换 P<0.05 的组合数")
}

xs <- tbl("xspecies_concordance")
if (!is.null(xs)) {
  cc <- as.numeric(pick(xs, c("concordance"))); sp <- pick(xs, c("species"))
  ct <- pick(xs, c("celltype", "cell_type")); bp <- as.numeric(pick(xs, c("binom_p")))
  j <- which.max(cc)
  TK("XS_MAX_CONC", paste0(fmt(100 * cc[j], 1), "%"), cc[j], src_of(xs), paste(sp[j], ct[j]),
     "跨物种方向一致比例最高的一组")
  TK("XS_MAX_LABEL", paste0(sp[j], " ", ct[j]), paste(sp[j], ct[j], sep = "/"), src_of(xs), "")
  TK("XS_MAX_P", fmt_p(bp[j]), bp[j], src_of(xs), paste(sp[j], ct[j]))
  TK("XS_MIN_CONC", paste0(fmt(100 * min(cc, na.rm = TRUE), 1), "%"), min(cc, na.rm = TRUE),
     src_of(xs), "", "最低一组，用于诚实交代异质性")
  TK("XS_SPECIES", paste(unique(sp), collapse = ", "), paste(unique(sp), collapse = ","), src_of(xs), "")
  ok <- which(cc >= 0.7 & bp < 0.05)
  TK("N_XS_PASS", length(ok), length(ok), src_of(xs), "", "达到保守判据（≥70% 且 P<0.05）的组合数")
}

# ================================= 6.5) 正文预留数值的补 token =============
# 正文里有几处数字当初没做成 token（占位写成 [XX]，靠人手填），这里统一补上，
# 让 17_fill_tokens.R 一次替换干净。任何一项取不到就跳过，正文仍保留 [XX]，
# 不会写进错误的值——宁可留着显眼，也不要留个看起来合理的数字。
# 这几条 token 的来源表都已在 01–14 里写出，不需要额外计算。

# --- 建模过程：三个方法各自保留的基因数 --------------------------------
las <- tbl("ml_lasso_coef")
if (!is.null(las)) TK("N_LASSO_FEAT", nrow(las), nrow(las), src_of(las), "",
                      "LASSO 在选定 lambda 下保留的非零系数个数")

rfe <- tbl("ml_svmrfe_curve")
if (!is.null(rfe)) {
  sz <- as.numeric(pick(rfe, c("size", "n_features")))
  ac <- as.numeric(pick(rfe, c("accuracy", "acc")))
  if (length(sz) && any(is.finite(ac))) {
    k <- which.max(ac)
    TK("SVM_MIN_FEAT", sz[k], sz[k], src_of(rfe), paste0("size=", sz[k]),
       "SVM-RFE 交叉验证准确率最高时对应的特征数")
  }
}

imp <- tbl("ml_rf_importance")
if (!is.null(imp)) {
  dec <- pick(imp, c("decision", "Decision"))
  n_b <- if (!is.null(dec)) sum(dec == "Confirmed", na.rm = TRUE) else nrow(imp)
  TK("N_BORUTA", n_b, n_b, src_of(imp), if (is.null(dec)) "ranger 前 20" else "Boruta Confirmed",
     if (is.null(dec)) "未跑 Boruta，取 ranger 重要性前 20" else "Boruta 确认的基因数")
}

# --- 校准与决策曲线：Hosmer-Lemeshow、十分位偏差、净收益阈值区间 ---------
# 这三组数由 05 的 model_calibration.tsv 提供（本文件只读，不重算）
cal <- tbl("model_calibration")
if (!is.null(cal)) {
  hl  <- as.numeric(pick(cal, c("hoslem_p", "hl_p")))[1]
  gap <- as.numeric(pick(cal, c("max_abs_gap", "max_gap")))[1]
  lo  <- as.numeric(pick(cal, c("dca_thr_lo", "thr_lo")))[1]
  hi  <- as.numeric(pick(cal, c("dca_thr_hi", "thr_hi")))[1]
  dg  <- as.numeric(pick(cal, c("hoslem_df")))[1]
  if (is.finite(hl))  TK("CALIB_HL_P", fmt_p(hl), hl, src_of(cal), "",
                         "Hosmer-Lemeshow 拟合优度 P 值；大于 0.05 才能说校准可接受")
  if (is.finite(dg))  TK("CALIB_HL_DF", dg, dg, src_of(cal), "", "H-L 卡方自由度")
  if (is.finite(gap)) TK("CALIB_MAX_GAP", fmt(gap, 3), gap, src_of(cal), "",
                         "十分位内预测概率与实测比例的绝对偏差最大值")
  if (is.finite(lo))  TK("DCA_THR_LO", fmt(lo, 2), lo, src_of(cal), "",
                         "净收益优于全治/不治的阈值下界")
  if (is.finite(hi))  TK("DCA_THR_HI", fmt(hi, 2), hi, src_of(cal), "",
                         "净收益优于全治/不治的阈值上界")
}

# --- Bootstrap 校正 AUC（05 的 model_optimism.tsv）-------------------------
opt <- tbl("model_optimism")
if (!is.null(opt)) {
  a_c <- as.numeric(pick(opt, c("auc_corrected", "auc_opt")))[1]
  opm <- as.numeric(pick(opt, c("optimism_mean", "optimism")))[1]
  nb  <- as.numeric(pick(opt, c("n_ok", "n_boot")))[1]
  if (is.finite(a_c)) {
    TK("AUC_CORRECTED", fmt(a_c, 3), a_c, src_of(opt), "",
       "bootstrap optimism 校正后的训练集 AUC（系数层校正，未含基因筛选）")
  }
  if (is.finite(opm)) TK("AUC_OPTIMISM", fmt(opm, 3), opm, src_of(opt), "",
                         "表观 AUC 与校正 AUC 之差（optimism）")
  if (is.finite(nb))  TK("AUC_BOOT_N", nb, nb, src_of(opt), "",
                         "可用 bootstrap 重抽次数")
}

# --- 单细胞 TM/JCT 对 bulk 方向的复现数 --------------------------------
scd <- tbl("sc_deg_TM")
if (!is.null(scd) && !is.null(cons)) {
  gs  <- toupper(as.character(pick(scd, c("gene"))))
  ls_ <- as.numeric(pick(scd, c("avg_log2FC", "avg_logFC", "logFC")))
  gc_ <- toupper(as.character(pick(cons, c("gene"))))
  lc_ <- as.numeric(pick(cons, c("logFC_train", "logFC")))
  j <- match(gc_, gs)
  ok <- !is.na(j) & is.finite(ls_[j]) & is.finite(lc_)
  if (sum(ok) > 0) {
    n_same <- sum(sign(ls_[j][ok]) == sign(lc_[ok]))
    TK("TM_CONCORD_N", n_same, n_same, paste(src_of(scd), "+", src_of(cons)), "",
       "单细胞 TM/JCT 中与 bulk 同向变化的共识基因数")
    TK("TM_CONCORD_PCT", paste0(fmt(100 * n_same / sum(ok), 1), "%"), n_same / sum(ok),
       paste(src_of(scd), "+", src_of(cons)), "",
       "分母是两平台都检测到的共识基因数，不是 N_CONSENSUS")
    TK("TM_CONCORD_DENOM", sum(ok), sum(ok), paste(src_of(scd), "+", src_of(cons)), "")
  }
}

# --- 跨物种：Spearman rho 与模块打分 AUC -------------------------------
if (!is.null(xs)) {
  cc3 <- as.numeric(pick(xs, c("concordance")))
  rho <- as.numeric(pick(xs, c("spearman_rho", "rho")))
  rp  <- as.numeric(pick(xs, c("spearman_p", "rho_p")))
  if (!is.null(rho) && length(cc3) && !is.null(sp) && !is.null(ct)) {
    j3 <- which.max(cc3)
    TK("XS_RHO", fmt(rho[j3], 2), rho[j3], src_of(xs), paste(sp[j3], ct[j3]),
       "该组 logFC 的 Spearman 相关系数")
    if (!is.null(rp)) TK("XS_RHO_P", fmt_p(rp[j3]), rp[j3], src_of(xs), paste(sp[j3], ct[j3]))
  }
}
msc <- tbl("xspecies_modulescore")
if (!is.null(msc)) {
  au <- as.numeric(pick(msc, c("auc", "AUC")))
  gs2 <- pick(msc, c("gene_set", "set")); ctm <- pick(msc, c("celltype", "cell_type"))
  if (length(au) && any(is.finite(au))) {
    j4 <- which.max(au)
    TK("XS_AUC", fmt(au[j4], 3), au[j4], src_of(msc), paste(gs2[j4], ctm[j4]),
       "人源模块打分在该物种细胞上的判别 AUC")
    TK("XS_AUC_LABEL", paste(gs2[j4], ctm[j4]), paste(gs2[j4], ctm[j4], sep = "/"), src_of(msc), "")
  }
}

# --- FOCUS 精细定位后的缩小集 -----------------------------------------
focus_files <- list.files(res_path(cfg, "tables"), pattern = "^focus_.*\\.tsv$", full.names = TRUE)
if (length(focus_files)) {
  tot <- 0L; n_hi <- 0L
  for (ff in focus_files) {
    d <- tryCatch(utils::read.delim(ff, check.names = FALSE, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    if (is.null(d) || !nrow(d)) next
    pip <- pick(d, c("PIP", "pip"))
    if (is.null(pip)) next
    tot <- tot + nrow(d); n_hi <- n_hi + sum(as.numeric(pip) > cfg_get(cfg, "genetics.twas.focus_min_pp", 0.5),
                                             na.rm = TRUE)
  }
  if (tot > 0) {
    TK("N_FOCUS_TESTED", tot, tot, "focus_*.tsv", "", "FOCUS 精细定位的基因×组织数")
    TK("N_FOCUS_PIP", n_hi, n_hi, "focus_*.tsv", "",
       paste0("PIP > ", cfg_get(cfg, "genetics.twas.focus_min_pp", 0.5), " 的基因×组织数"))
  }
}

# --- 扰动排名与遗传学证据的重合度（正文 3.8 那句 "top five" 用）--------
if (!is.null(pert)) {
  gn3 <- as.character(pick(pert, c("gene")))
  sc3 <- as.numeric(pick(pert, c("score")))
  top5 <- unique(gn3[order(sc3)])[seq_len(min(5, length(unique(gn3))))]
  gcol <- character(0); gmr <- character(0)
  if (!is.null(col)) {
    p4 <- as.numeric(pick(col, c("PP.H4", "PP4"))); g4 <- as.character(pick(col, c("gene")))
    gcol <- unique(g4[which(p4 > cfg_get(cfg, "genetics.coloc.pp4_threshold", 0.8))])
  }
  if (!is.null(mr)) {
    gm <- as.character(pick(mr, c("gene", "exposure"))); pm <- as.numeric(pick(mr, c("pval", "p")))
    gmr <- unique(gm[which(pm < 0.05)])
  }
  n_ov <- sum(top5 %in% union(gcol, gmr))
  TK("N_PERT_GENETIC", n_ov, n_ov, "perturb_summary ∩ (coloc ∪ mr)",
     paste(top5, collapse = ","),
     "排前五的扰动基因中同时有共定位或 MR 支持的个数；注意两者不独立")
  TK("PERT_TOP5", paste(top5, collapse = ", "), paste(top5, collapse = ","), src_of(pert), "")
}

# --- MR-Egger 截距：有水平多效性嫌疑的基因数 ---------------------------
# 列在 mr_sensitivity.tsv（10_mr.R 写的 sensitivity 汇总），不是 mr_results.tsv；
# mr_results 是按 method 逐行的方法学结果，没有基因级的截距列。
# 正文那句 "all but [XX] genes" 取这里的数。
mrs <- tbl("mr_sensitivity")
ep_src <- if (!is.null(mrs) && any(c("egger_intercept_p", "intercept_p") %in% names(mrs))) {
  mrs
} else {
  mr
}
if (!is.null(ep_src)) {
  ep <- as.numeric(pick(ep_src, c("egger_intercept_p", "intercept_p")))
  if (length(ep) && any(is.finite(ep))) {
    n_ep <- sum(ep < 0.05, na.rm = TRUE)
    TK("N_EGGER_OUT", n_ep, n_ep, src_of(ep_src), "",
       "Egger 截距 P<0.05 的基因数（>0 时讨论里必须交代多效性）")
    TK("N_EGGER_TESTED", sum(is.finite(ep)), sum(is.finite(ep)), src_of(ep_src), "")
    TK("EGGER_P_MIN", fmt_p(min(ep[is.finite(ep)])), min(ep[is.finite(ep)]), src_of(ep_src), "")
  }
}

# --- 组间显著改变的配体-受体对 -----------------------------------------
# 06 按条件各写一张 cellchat_<cond>.tsv。以 source|target|interaction_name 为键
# 比较两侧：一侧显著（P<0.05）且 |Δprob|>0.05 记为一对"改变"。
# 这是个可复现的操作性定义，正文的写法必须与之一致，不要写成"CellChat 差异分析 N 对"。
ccf <- list.files(res_path(cfg, "tables"), pattern = "^cellchat_.*\\.tsv$", full.names = TRUE)
ccf <- ccf[!grepl("cellchat_compare", ccf)]
if (length(ccf) == 2) {
  rd <- lapply(ccf, function(f) tryCatch(
    utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL))
  need <- c("source", "target", "interaction_name", "prob", "pval")
  if (all(vapply(rd, function(d) !is.null(d) && all(need %in% names(d)), logical(1)))) {
    key <- function(d) paste(d$source, d$target, d$interaction_name, sep = "|")
    k1 <- key(rd[[1]]); k2 <- key(rd[[2]])
    m1 <- rd[[1]]; m2 <- rd[[2]]
    idx <- match(intersect(k1, k2), k1); idx2 <- match(intersect(k1, k2), k2)
    dp <- m2$prob[idx2] - m1$prob[idx]
    sig <- (m1$pval[idx] < 0.05) | (m2$pval[idx2] < 0.05)
    n_alt <- sum(abs(dp) > 0.05 & sig, na.rm = TRUE)
    TK("N_LR_PAIRS", length(unique(c(k1, k2))), length(unique(c(k1, k2))), "cellchat_*.tsv", "",
       "两个条件下检出的配体-受体对并集")
    TK("N_LR_ALTERED", n_alt, n_alt, "cellchat_*.tsv", "",
       "|Δprob|>0.05 且至少一侧 P<0.05 的对数（操作性定义，正文需照此描述）")
    # 正文 Figure 6C 那句 "the largest decreases involved [XX]" 取这里：
    # 在满足上述操作性定义的配对里挑 Δprob 最负的一个，连同细胞对一起给出，
    # 避免作者凭图凭印象手挑一个对不上的名字
    if (n_alt > 0) {
      cand <- which(abs(dp) > 0.05 & sig)
      j <- cand[which.min(dp[cand])]
      lab <- paste0(m2$interaction_name[idx2][j], " signalling (",
                    m2$source[idx2][j], " to ", m2$target[idx2][j],
                    ", delta prob ", fmt(dp[j], 2), ")")
      TK("LR_TOP_DECREASE", lab, dp[j], "cellchat_*.tsv", "",
         "下降幅度最大的配体-受体对（delta prob 最负），供正文 Figure 6C 句使用")
    }
  }
}

# --- CMap 反向富集命中数 ----------------------------------------------
cm <- tbl("cmap_connectivity")
if (!is.null(cm)) {
  cs2 <- as.numeric(pick(cm, c("cs", "score")))
  n_sig_cm <- if (length(cs2) && any(is.finite(cs2))) sum(cs2 < 0, na.rm = TRUE) else nrow(cm)
  TK("N_CMAP", format(nrow(cm), big.mark = ",", trim = TRUE, scientific = FALSE), nrow(cm), src_of(cm), "", "药物连通性分析的化合物总数")
  TK("N_CMAP_NEG", format(n_sig_cm, big.mark = ",", trim = TRUE, scientific = FALSE), n_sig_cm, src_of(cm), "", "连通性得分为负（可能逆转疾病特征）的化合物数")
}

# --- 去卷积组成差异（04 的 deconv_group_diff）--------------------------
dd <- tbl("deconv_group_diff")
if (!is.null(dd) && nrow(dd) && "fdr" %in% names(dd)) {
  fdrv <- as.numeric(dd$fdr)
  ddv  <- as.numeric(dd$diff)
  ctv  <- as.character(dd$cell_type)
  fin  <- is.finite(fdrv) & is.finite(ddv)
  up <- ctv[fin & ddv > 0 & fdrv < 0.05]
  dn <- ctv[fin & ddv < 0 & fdrv < 0.05]
  # 类型代码 -> 正文可读标签（未收录的原样输出）
  lab <- c("TM_Fibro" = "trabecular-meshwork fibroblasts",
           "Schlemm_Endo" = "Schlemm canal endothelium",
           "Macrophage" = "macrophages", "Lymphocyte" = "lymphocytes",
           "Mast" = "mast cells", "Pericyte" = "pericytes",
           "Schwann" = "Schwann cells", "Ciliary_Muscle" = "ciliary muscle",
           "Uveal_Fibro" = "uveal fibroblasts", "Sclera_Fibro" = "scleral fibroblasts",
           "K_Fibro" = "keratocyte-like fibroblasts", "FibroX" = "fibroblast-like cells",
           "K_Endo" = "corneal endothelium", "VascEndo" = "vascular endothelium",
           "Lymphatic_Endo" = "lymphatic endothelium",
           "CC_VenEndo" = "collector-channel venous endothelium",
           "Uveal_Melanocyte" = "uveal melanocytes",
           "Conj_Melanocyte" = "conjunctival melanocytes",
           "Goblet" = "goblet cells",
           "Conj_Epi-Superficial" = "superficial conjunctival epithelium",
           "Conj_Epi-Wing" = "wing conjunctival epithelium",
           "Conj_Epi-Basal" = "basal conjunctival epithelium",
           "K_Epi-Superficial" = "superficial corneal epithelium",
           "K_Epi-Wing" = "wing corneal epithelium",
           "K_Epi-Basal" = "basal corneal epithelium",
           "Limbal_Epi-Superficial" = "superficial limbal epithelium",
           "Limbal_Epi-Wing" = "wing limbal epithelium",
           "Limbal_Epi-Basal" = "basal limbal epithelium")
  lab_fn <- function(v) { l <- unname(lab[v]); ifelse(is.na(l), v, l) }
  fmt_v <- function(v) if (length(v)) paste(lab_fn(v), collapse = ", ") else "none"
  TK("DECONV_HIGHER", fmt_v(up), paste(up, collapse = ","), src_of(dd), "",
     "POAG 中估计占比更高的细胞类型（FDR<0.05）")
  TK("DECONV_LOWER", fmt_v(dn), paste(dn, collapse = ","), src_of(dd), "",
     "POAG 中估计占比更低的细胞类型（FDR<0.05）")
}

# --- SCENIC 差异 regulon 数 ------------------------------------------------
# 06 不做 SCENIC（网络推断在 python 的 pySCENIC 里），但正文写了差异 regulon 数。
# tools/scenic_regulon_count.R 把 pySCENIC 的 AUC 矩阵做成 scenic_regulon_diff.tsv，
# 这里只读结果、不重算；文件不存在就让正文留着 {{N_SCENIC_REGULON}}。
scf <- res_path(cfg, "tables", "scenic_regulon_diff.tsv")
if (file.exists(scf)) {
  sc <- utils::read.delim(scf, check.names = FALSE, stringsAsFactors = FALSE)
  pa <- as.numeric(pick(sc, c("padj", "p_adj", "FDR")))
  if (length(pa) && any(is.finite(pa))) {
    n_sc <- sum(pa < 0.05, na.rm = TRUE)
    TK("N_SCENIC_REGULON", n_sc, n_sc, "scenic_regulon_diff.tsv", "",
       "BH 校正后 padj<0.05 的差异活动 regulon 数")
    TK("N_SCENIC_TESTED", sum(is.finite(pa)), sum(is.finite(pa)),
       "scenic_regulon_diff.tsv", "", "参与检验的 regulon 总数")
  }
}

# ============================================================ 7) 版本信息 ==
si <- utils::capture.output(utils::sessionInfo())
dir.create(res_path(cfg), recursive = TRUE, showWarnings = FALSE)
writeLines(si, res_path(cfg, "session_info.txt"))
rv <- paste(R.version$major, R.version$minor, sep = ".")
TK("R_VERSION", rv, rv, "session_info.txt", "", "R 版本，写进方法学")
# 关键包版本：方法段落要求写明
for (p in c("Seurat", "limma", "WGCNA", "TwoSampleMR", "coloc", "clusterProfiler", "CellChat")) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  if (!is.na(v)) TK(paste0("VER_", toupper(p)), v, v, "session_info.txt", p, "包版本")
}

# ============================================================ 8) 写出 =====
if (!length(LEDGER)) stop("账本为空：没有任何结果表可读。先跑 01–14，或确认 project.outdir 指对了",
                          call. = FALSE)
led <- do.call(rbind, LEDGER)
led <- led[!duplicated(led$token), ]
write_tsv(led, "results_ledger", cfg)

# 供 Word 查找替换：{{TOKEN}} -> 值
map <- data.frame(find = paste0("{{", led$token, "}}"), replace = led$display,
                  stringsAsFactors = FALSE)
utils::write.table(map, res_path(cfg, "tables", "manuscript_token_map.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("已写出查找替换表：manuscript_token_map.tsv（", nrow(map), " 条）")

# 缺步骤审计
# key_table 取「该步只要跑通就必然产出」的主表，不要取可选分支的产物：
# 04 的反卷积依赖单细胞参考（或 xCell），缺参考时富集照样完成，不能因此把 04 判成没跑。
step_files <- c(
  "01" = "qc_summary", "02" = "deg_consensus", "03" = "wgcna_module_trait",
  "04" = "enrich_go_Consensus", "05" = "model_performance", "06" = "sc_celltype_prop",
  "07" = "gwas_qc_summary", "08" = "coloc_summary", "09" = "twas_all", "10" = "mr_results",
  "11" = "perturb_summary", "12" = "mr_mediation", "13" = "robustness_permutation",
  "14" = "xspecies_concordance")
# 缺的时候到底缺什么，写清楚，省得每次都去翻脚本
step_notes <- c(
  "01" = "", "02" = "", "03" = "",
  "04" = "富集已完成；deconv_fraction 需 data/processed/scRNA_ref.rds 或 xCell 包",
  "05" = "", "06" = "需单细胞对象 data/processed/scRNA_ref.rds（含 cell_type 列）",
  "07" = "", "08" = "", "09" = "需 FUSION 权重或 FOCUS 输入，未配置时跳过",
  "10" = "", "11" = "需单细胞对象 + 扰动基因集（scRNA_ref.rds）",
  "12" = "需中介暴露 GWAS（如 IOP）与暴露-结局配对",
  "13" = "", "14" = "需跨物种（小鼠/猕猴）表达矩阵")
have <- vapply(step_files, function(f) file.exists(res_path(cfg, "tables", paste0(f, ".tsv"))), logical(1))
ms <- data.frame(step = names(step_files), key_table = unname(step_files),
                 has_output = unname(have),
                 consequence = ifelse(have, "", "对应结果段落暂时无法撰写，正文保留 {{TOKEN}} 占位"),
                 note = ifelse(have, "", unname(step_notes[names(step_files)])),
                 stringsAsFactors = FALSE)
utils::write.table(ms, res_path(cfg, "tables", "missing_steps.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
if (any(!have)) log_msg("缺产物的步骤：", paste(names(step_files)[!have], collapse = ", "))

# 可读版 md：按前缀分组
grp <- sub("_.*$", "", led$token)
md <- c("# 结果账本（自动生成，勿手改）", "",
        paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M"), "，共 ", nrow(led), " 个 token"), "")
for (g in unique(grp)) {
  md <- c(md, paste0("## ", g), "", "| token | 值 | 来源 | 备注 |", "|---|---|---|---|")
  sub <- led[grp == g, ]
  md <- c(md, sprintf("| {{%s}} | %s | %s | %s |", sub$token, sub$display, sub$source,
                      ifelse(nzchar(sub$note), sub$note, "")))
  md <- c(md, "")
}
writeLines(md, res_path(cfg, "results_ledger.md"))
log_msg("账本：", nrow(led), " 个 token，md 版见 results/results_ledger.md")

# 可选：校验稿件里的 token
if (!is.null(check_html)) {
  if (!file.exists(check_html)) stop("稿件不存在：", check_html, call. = FALSE)
  txt <- paste(readLines(check_html, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  used <- unique(unlist(regmatches(txt, gregexpr("\\{\\{[A-Z0-9_]+\\}\\}", txt))))
  used <- gsub("[{}]", "", used)
  unknown <- setdiff(used, led$token)
  unused <- setdiff(led$token, used)
  cat("\n==== 稿件 token 校验 ====\n")
  cat("稿件引用：", length(used), " 个；账本共 ", nrow(led), " 个\n", sep = "")
  cat("账本中不存在（拼写错误或该步骤未跑）：",
      if (length(unknown)) paste(unknown, collapse = ", ") else "无", "\n")
  cat("账本有但稿件未用（可选补充）：",
      if (length(unused)) paste(head(unused, 20), collapse = ", ") else "无", "\n")
  utils::write.table(data.frame(token = unknown), res_path(cfg, "tables", "token_unknown.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
}

log_msg("完成")

# =============================================================================
# 04_enrich_deconv.R —— 课题一（4）：通路富集、通路活性、免疫/细胞反卷积
#
# 输入：data/processed/<trainGSE>_expr.rds / _pheno.rds
#       results/tables/deg_<trainGSE>.tsv
#       data/processed/wgcna_<trainGSE>.rds（可选，有则做模块-细胞比例关联）
#       data/processed/scRNA_ref.rds（可选，Seurat 对象；有则自建 signature）
# 输出：results/tables/enrich_go_<set>.tsv / enrich_kegg_<set>.tsv / enrich_reactome_<set>.tsv
#       results/tables/gsva_pathway_diff.tsv / gsea_pathway.tsv
#         （集合来源见 set_src：有网用 MSigDB Hallmark，离线兜底用本地 GO BP 词条）
#       results/tables/progeny_diff.tsv
#       results/tables/deconv_fraction.tsv
#       results/tables/deconv_group_diff.tsv
#       results/figures/enrich_dotplot_<set>.pdf, deconv_box.pdf
#
# 注意：
#   - CIBERSORTx 的 LM22 是血液肿瘤芯片训练的矩阵，直接套眼组织会被审稿人质疑；
#     这里默认用单细胞参考自建 signature + NNLS，属于可解释的替代方案
#   - GSVA/GSEA 的基因集版本要写进方法学（msigdbr 版本号会影响结果）
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("04_enrich_deconv", cfg)
load_pkgs(c("ggplot2", "clusterProfiler", "org.Hs.eg.db", "ReactomePA",
            "GSVA", "limma", "reshape2"))
# msigdbr 是可选依赖：只影响 GSVA/GSEA（Hallmark），缺了不应该阻断 ORA 与反卷积
has_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
if (!has_msigdbr) log_msg("[WARN] 未安装 msigdbr，将跳过 Hallmark 相关的 GSVA/GSEA")
set.seed(cfg_get(cfg, "project.seed", 1))

gse <- if (isTRUE(cfg_get(cfg, "run.demo", FALSE))) "DEMO" else cfg_get(cfg, "bulk.train.gse", "GSE27276")
expr  <- readRDS(proj_path("data", "processed", paste0(gse, "_expr.rds")))
pheno <- readRDS(proj_path("data", "processed", paste0(gse, "_pheno.rds")))
tt    <- utils::read.delim(res_path(cfg, "tables", paste0("deg_", gse, ".tsv")),
                           stringsAsFactors = FALSE)

p_cut <- cfg_get(cfg, "bulk.enrich.p_cut", 0.05)
q_cut <- cfg_get(cfg, "bulk.enrich.q_cut", 0.2)

# --------------------------------------------------------- 1) ORA 富集 ------

gene_sets <- list(
  Up   = tt$gene[tt$sig == "Up"],
  Down = tt$gene[tt$sig == "Down"]
)
if (file.exists(res_path(cfg, "tables", "deg_consensus.tsv"))) {
  cons <- utils::read.delim(res_path(cfg, "tables", "deg_consensus.tsv"), stringsAsFactors = FALSE)
  gene_sets[["Consensus"]] <- cons$gene
}

ens <- function(g) {
  ids <- suppressWarnings(clusterProfiler::bitr(g, fromType = "SYMBOL", toType = "ENTREZID",
                                                OrgDb = org.Hs.eg.db))
  unique(ids$ENTREZID)
}

ora_all <- list()
for (nm in names(gene_sets)) {
  g <- unique(gene_sets[[nm]])
  if (length(g) < 5) { log_msg(nm, " 基因数不足 5，跳过富集"); next }
  eid <- ens(g)
  if (length(eid) < 5) { log_msg(nm, " 映射到 ENTREZ 不足 5，跳过"); next }

  go <- clusterProfiler::enrichGO(eid, OrgDb = org.Hs.eg.db, ont = "BP",
                                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                                  readable = TRUE)

  # KEGG 每次都要联网读 rest.kegg.jp（无本地库），离线/受限网络下会卡到超时。
  # 已有结果表就直接复用；没有才联网，且失败只记 WARN，不中断整条流程。
  kegg_f <- res_path(cfg, "tables", paste0("enrich_kegg_", nm, ".tsv"))
  if (file.exists(kegg_f)) {
    kdf <- utils::read.delim(kegg_f, stringsAsFactors = FALSE)
    log_msg("复用已有 KEGG 结果：", basename(kegg_f), "（", nrow(kdf), " 行）")
  } else {
    kegg <- tryCatch(
      clusterProfiler::enrichKEGG(eid, organism = "hsa",
                                  pAdjustMethod = "BH", pvalueCutoff = 1),
      error = function(e) {
        log_msg("[WARN] KEGG 富集失败（", conditionMessage(e), "）；跳过 ", nm, " 的 KEGG")
        NULL
      })
    kdf <- if (!is.null(kegg) && nrow(as.data.frame(kegg))) {
      d <- as.data.frame(kegg)
      d$Description <- sub(" - Homo sapiens \\(human\\)$", "", d$Description)
      d
    } else NULL
  }

  reac <- ReactomePA::enrichPathway(eid, organism = "human",
                                    pAdjustMethod = "BH", pvalueCutoff = 1)

  if (!is.null(go) && nrow(as.data.frame(go))) {
    write_tsv(as.data.frame(go), paste0("enrich_go_", nm), cfg)
    ora_all[[paste0("GO_", nm)]] <- as.data.frame(go)
  }
  if (!is.null(kdf) && nrow(kdf)) {
    write_tsv(kdf, paste0("enrich_kegg_", nm), cfg)
    ora_all[[paste0("KEGG_", nm)]] <- kdf
  }
  if (!is.null(reac) && nrow(as.data.frame(reac))) {
    write_tsv(as.data.frame(reac), paste0("enrich_reactome_", nm), cfg)
    ora_all[[paste0("Reactome_", nm)]] <- as.data.frame(reac)
  }

  # dotplot（取 GO BP 前 15 条）
  if (!is.null(go) && nrow(as.data.frame(go)) >= 3) {
    gd <- head(as.data.frame(go), 15)
    gd$Description <- factor(gd$Description, levels = rev(gd$Description))
    gd$ratio <- gd$Count / as.numeric(sub("/.*", "", gd$BgRatio))
    g_dot <- ggplot2::ggplot(gd, ggplot2::aes(ratio, Description, size = Count, colour = p.adjust)) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_gradient(low = "#C0392B", high = "#1F6FB4", name = "FDR") +
      ggplot2::labs(x = "gene ratio", y = NULL, title = paste0(nm, ": GO BP (", gse, ")")) +
      theme_pub()
    save_fig(g_dot, paste0("enrich_dotplot_", nm), cfg, 7, 5)
  }
}

# --------------------------------------------------------- 2) msigdb 集合 ----
# msigdbr 10.x 把 category 改名为 collection，这里做双向兼容
get_msig <- function(cat = "H") {
  if (!has_msigdbr) return(NULL)
  if (!isTRUE(cfg_get(cfg, "bulk.enrich.msigdb", TRUE))) {
    log_msg("bulk.enrich.msigdb=false，跳过 Hallmark（GSVA/GSEA 不产出）")
    return(NULL)
  }
  # msigdbr 首次要联网下载基因集；不让它无限重试（默认 timeout 可能是 300s 且带重试）
  old_to <- getOption("timeout")
  options(timeout = 60)
  on.exit(options(timeout = old_to), add = TRUE)
  df <- tryCatch({
    if ("collection" %in% names(formals(msigdbr::msigdbr))) {
      msigdbr::msigdbr(species = "Homo sapiens", collection = cat)
    } else {
      msigdbr::msigdbr(species = "Homo sapiens", category = cat)
    }
  }, error = function(e) {
    # msigdbr 首次调用会去 zenodo 下载基因集，离线/受限网络下必然失败。
    # 这不是致命错误：ORA 与反卷积不依赖 Hallmark，不能因为它整段流程停掉。
    log_msg("[WARN] msigdbr 调用失败（", conditionMessage(e),
            "）；跳过 GSVA/GSEA，ORA 与反卷积继续")
    NULL
  })
  if (is.null(df) || !nrow(df)) {
    log_msg("[WARN] msigdbr 未返回任何基因集：", cat, "；跳过 Hallmark 相关分析")
    return(NULL)
  }
  split(df$gene_symbol, df$gs_name)
}

hallmark <- get_msig(cfg_get(cfg, "bulk.enrich.msigdb_category", "H"))
set_src <- "Hallmark (MSigDB)"

# 离线兜底：msigdbr 不可用（未装 / 未联网）时，改用本地 GO BP 词条做 GSVA。
# 不用 ORA 的显著词条来建集合——那是用「富集显著的通路」算通路活性，循环论证。
# 这里直接从 org.Hs.eg.db + GO.db 取全部 BP 词条（去 IEA 电子注释），在当前
# 数据集里 10–500 个基因的才进 GSVA。集合来源写进图标题，投稿时如实描述。
if (!length(hallmark) &&
    requireNamespace("GO.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE) &&
    isTRUE(cfg_get(cfg, "bulk.enrich.gsva_fallback", TRUE))) {
  log_msg("[兜底] msigdbr 不可用，改从本地 org.Hs.eg.db 取 GO BP 词条")
  go_map <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL"),
    columns = "GO", keytype = "SYMBOL")
  go_map <- go_map[!is.na(go_map$GO) & go_map$EVIDENCE != "IEA", c("SYMBOL", "GO")]
  go_meta <- AnnotationDbi::select(GO.db::GO.db, keys = unique(go_map$GO),
                                   columns = c("ONTOLOGY", "TERM"))
  go_meta <- go_meta[go_meta$ONTOLOGY == "BP", c("GOID", "TERM")]
  go_pair <- unique(merge(go_map, go_meta, by.x = "GO", by.y = "GOID"))
  gs <- split(go_pair$SYMBOL, go_pair$TERM)
  gs <- lapply(gs, function(g) intersect(unique(g), rownames(expr)))
  gs <- gs[lengths(gs) >= 10 & lengths(gs) <= 500]
  # 兜底是为真实数据准备的（2 万基因宇宙下有几千个 BP 词条）；demo 只有 500 个
  # 基因，达标词条会寥寥无几——少于 5 个时不如不做，否则图上孤零零一条
  if (length(gs) < 5) {
    log_msg("[兜底] GO BP 达标词条仅 ", length(gs), " 个（基因宇宙太小），跳过 GSVA")
    gs <- list()
  }
  if (length(gs)) {
    hallmark <- gs
    set_src <- "GO BP terms (local org.Hs.eg.db)"
    log_msg("[兜底] GO BP 可用词条：", length(gs), " 个")
  }
}

if (length(hallmark)) {
  log_msg("GSVA 基因集（", set_src, "）：", length(hallmark), " 个")

  # ------------------------------------------------------- 3) GSVA ----------
  # GSVA 2.x 改成了 param 对象式 API：旧签名 gsva(mat, sets, ...) 会报
  # "unable to find an inherited method for function 'gsva' for signature
  # 'param = \"matrix\"'"。按版本分派，1.x 走旧签名。
  gsva_mat <- if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
    GSVA::gsva(GSVA::gsvaParam(as.matrix(expr), hallmark, kcdf = "Gaussian"))
  } else {
    GSVA::gsva(as.matrix(expr), hallmark, method = "gsva",
               kcdf = "Gaussian", verbose = FALSE)
  }
  design <- stats::model.matrix(~ 0 + group, data = pheno)
  colnames(design) <- sub("^group", "", colnames(design))
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(gsva_mat, design),
                                            limma::makeContrasts(case - control, levels = design)))
  gtt <- limma::topTable(fit, number = Inf, sort.by = "P")
  gtt$pathway <- rownames(gtt)
  write_tsv(gtt, "gsva_pathway_diff", cfg)

  g_bar <- head(gtt[order(gtt$P.Value), ], 20)
  # 词条名可能很长（GO BP 的 TERM 动辄 80+ 字符），超长折行会让条形图变表格
  g_bar$pathway <- ifelse(nchar(g_bar$pathway) > 46,
                          paste0(substr(g_bar$pathway, 1, 44), "…"), g_bar$pathway)
  # GO BP 词条名前 44 字符可能相同（"positive regulation of…"一族），截断后必须去重，
  # 否则 factor 会把两条合并成一条，20 条变 19 条还不报错
  g_bar$pathway <- make.unique(g_bar$pathway)
  g_bar$pathway <- factor(g_bar$pathway, levels = rev(g_bar$pathway))
  g_bar$dir <- ifelse(g_bar$logFC > 0, "Up in POAG", "Down in POAG")
  g_gsva <- ggplot2::ggplot(g_bar, ggplot2::aes(logFC, pathway, fill = dir)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = c("Up in POAG" = "#C0392B", "Down in POAG" = "#1F6FB4"), name = NULL) +
    ggplot2::labs(x = "GSVA score difference (case - control)", y = NULL,
                  title = paste0("Pathway activity (", set_src, ", ", gse, ")")) +
    theme_pub()
  save_fig(g_gsva, "gsva_pathway_bar", cfg, 7.5, 5.5)

  # ------------------------------------------------------- 4) GSEA ----------
  ranks <- tt$t
  names(ranks) <- tt$gene
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
  gsea <- clusterProfiler::GSEA(ranks, TERM2GENE = data.frame(
    term = rep(names(hallmark), lengths(hallmark)),
    gene = unlist(hallmark, use.names = FALSE)), pvalueCutoff = 1, seed = TRUE)
  if (!is.null(gsea) && nrow(as.data.frame(gsea))) {
    write_tsv(as.data.frame(gsea), "gsea_pathway", cfg)
  }
} else {
  log_msg("GSVA/GSEA 基因集不可用，跳过；后续 PROGENy 与反卷积不受影响")
}

# --------------------------------------------------------- 5) PROGENy -------
prog_ok <- requireNamespace("progeny", quietly = TRUE)
if (prog_ok) {
  ps <- progeny::progeny(expr, scale = TRUE, organism = "Human",
                         top = 100, perm = 1)
  # progeny 不同版本返回值方向不一致：先统一成「样本 × 通路」，再转置给 limma
  if (nrow(ps) != ncol(expr)) ps <- t(ps)
  pdesign <- stats::model.matrix(~ 0 + group, data = pheno)
  colnames(pdesign) <- sub("^group", "", colnames(pdesign))
  pfit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(t(ps), pdesign),
                                             limma::makeContrasts(case - control, levels = pdesign)))
  ptt <- limma::topTable(pfit, number = Inf, sort.by = "P")
  ptt$pathway <- rownames(ptt)
  write_tsv(ptt, "progeny_diff", cfg)

  # Fig4C：PROGENy 通路活性差异。14 条通路全画——footprint 方法的好处就是
  # 通路集合固定，跨数据集可以直接对比；* 为 BH 校正后 <0.05
  # NA 的通路先剔除：demo 这类小基因宇宙里，footprint top-100 与矩阵重叠过少
  # 会导致得分零方差、limma 系数为 NA；画出来是空槽位，图例还会冒出 "NA" 一档
  n_na <- sum(is.na(ptt$logFC))
  if (n_na) log_msg("[WARN] PROGENy 有 ", n_na, "/", nrow(ptt),
                    " 条通路 logFC 为 NA（footprint 基因在矩阵中覆盖不足），图中不画")
  pbar <- ptt[!is.na(ptt$logFC), ]
  if (nrow(pbar)) {
  pbar$pathway <- factor(pbar$pathway, levels = pbar$pathway)
  pbar$dir <- ifelse(pbar$logFC > 0, "Up in POAG", "Down in POAG")
  pbar$star <- ifelse(pbar$adj.P.Val < 0.05, "*", "")
  g_prog <- ggplot2::ggplot(pbar, ggplot2::aes(logFC, pathway, fill = dir)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_text(ggplot2::aes(label = star, hjust = ifelse(logFC > 0, -0.25, 1.25)),
                       size = 4.2, colour = "grey25") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c("Up in POAG" = "#C0392B", "Down in POAG" = "#1F6FB4"), name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
    ggplot2::labs(x = "PROGENy activity difference (case - control)", y = NULL,
                  title = "Pathway activity (PROGENy footprint)") +
    theme_pub() +
    ggplot2::theme(legend.position = "top")
  save_fig(g_prog, "progeny_pathway", cfg, 6.8, 4.6)
  }
} else {
  log_msg("未安装 progeny，跳过通路活性推断")
}

# ------------------------------------------------ 6) 细胞比例反卷积 ---------
# 6a. 自建 signature：单细胞参考 → 每个细胞类型高特异性高表达基因
build_signature <- function(sc, cell_types, max_genes = 300, min_expr = 0.5) {
  stopifnot(inherits(sc, "Seurat"))
  mat <- tryCatch(Seurat::GetAssayData(sc, assay = "RNA", layer = "data"),
                  error = function(e) Seurat::GetAssayData(sc, assay = "RNA", slot = "data"))
  if (!"cell_type" %in% colnames(sc@meta.data)) stop("单细胞对象缺少 cell_type 列")
  ct <- as.character(sc$cell_type)
  means <- sapply(cell_types, function(k) Matrix::rowMeans(mat[, ct == k, drop = FALSE]))
  means <- as.matrix(means)
  others <- apply(means, 1, function(r) sort(r, decreasing = TRUE)[2])
  sig <- list()
  for (k in cell_types) {
    spec <- means[, k] - others
    cand <- names(sort(spec[means[, k] > min_expr], decreasing = TRUE))
    sig[[k]] <- head(cand, max_genes)
  }
  genes <- unique(unlist(sig))
  S <- means[genes, , drop = FALSE]
  S[S < 0] <- 0
  attr(S, "markers") <- sig
  S
}

deconv_nnls <- function(expr, S) {
  if (!requireNamespace("nnls", quietly = TRUE)) {
    stop("需要 nnls 包：install.packages('nnls')", call. = FALSE)
  }
  common <- intersect(rownames(S), rownames(expr))
  if (length(common) < 50) stop("signature 与 bulk 交集基因过少（", length(common), "），无法反卷积")
  S <- S[common, , drop = FALSE]
  S <- S[rowSums(S) > 0, , drop = FALSE]
  common <- rownames(S)
  X <- expr[common, , drop = FALSE]
  frac <- matrix(NA_real_, ncol(X), ncol(S), dimnames = list(colnames(X), colnames(S)))
  for (j in seq_len(ncol(X))) {
    y <- X[, j]
    fit <- nnls::nnls(S, y)$x
    if (sum(fit) > 0) frac[j, ] <- fit / sum(fit)
  }
  as.data.frame(frac)
}

sc_path <- proj_path("data", "processed", "scRNA_ref.rds")
frac <- NULL
if (file.exists(sc_path) && requireNamespace("Seurat", quietly = TRUE)) {
  log_msg("使用单细胞参考构建 signature 并做 NNLS 反卷积")
  sc <- readRDS(sc_path)
  cts <- cfg_get(cfg, "bulk.deconv.signature_cells",
                 c("TM", "JCT", "SC", "immune", "fibroblast", "endothelial"))
  cts <- intersect(cts, unique(as.character(sc$cell_type)))
  S <- build_signature(sc, cts,
                       cfg_get(cfg, "bulk.deconv.max_genes_per_celltype", 300),
                       cfg_get(cfg, "bulk.deconv.min_expr", 0.5))
  frac <- deconv_nnls(expr, S)
} else if (requireNamespace("xCell", quietly = TRUE)) {
  log_msg("无单细胞参考，改用 xCell（注意：xCell 基于组织表达谱，对眼组织为间接推断）")
  xc <- xCell::xCellAnalysis(expr)
  keep <- grep("aDC|DC|Macrophage|Monocyte|Neutrophil|B cell|CD4|CD8|T cell|NK|Fibroblast|Endothelial|Mesangial|Pericyte|Adipocyte|Smooth muscle",
               rownames(xc), value = TRUE)
  frac <- as.data.frame(t(xc[keep, , drop = FALSE]))
} else {
  log_msg("既无单细胞参考也无 xCell，跳过反卷积")
}

if (!is.null(frac)) {
  frac$sample <- rownames(frac)
  frac <- merge(pheno[, c("sample", "group")], frac, by = "sample")
  write_tsv(frac, "deconv_fraction", cfg)

  cell_cols <- setdiff(names(frac), c("sample", "group"))
  diff_rows <- lapply(cell_cols, function(cc) {
    x <- frac[[cc]][frac$group == "case"]; y <- frac[[cc]][frac$group == "control"]
    if (length(x) < 3 || length(y) < 3) return(NULL)
    wt <- stats::wilcox.test(x, y)
    data.frame(cell_type = cc, mean_case = mean(x), mean_control = mean(y),
               diff = mean(x) - mean(y), p = wt$p.value)
  })
  ddiff <- do.call(rbind, diff_rows)
  if (!is.null(ddiff)) {
    ddiff$fdr <- stats::p.adjust(ddiff$p, "BH")
    write_tsv(ddiff, "deconv_group_diff", cfg)
  }

  long <- reshape2::melt(frac[, c("sample", "group", cell_cols)], id.vars = c("sample", "group"))
  g_box <- ggplot2::ggplot(long, ggplot2::aes(group, value, fill = group)) +
    ggplot2::geom_boxplot(outlier.size = 0.4, width = 0.6) +
    ggplot2::facet_wrap(~ variable, scales = "free_y", ncol = 4) +
    ggplot2::scale_fill_manual(values = c(control = "#1F6FB4", case = "#C0392B"), guide = "none") +
    ggplot2::labs(x = NULL, y = "estimated fraction") +
    theme_pub()
  save_fig(g_box, "deconv_box", cfg, 9, 6)

  # 模块 eigengene 与细胞比例的相关（把"哪个细胞类型驱动哪个模块"讲清楚）
  wgcna_file <- proj_path("data", "processed", paste0("wgcna_", gse, ".rds"))
  if (file.exists(wgcna_file)) {
    w <- readRDS(wgcna_file)
    me <- w$ME[frac$sample, , drop = FALSE]
    cc <- WGCNA::cor(me, frac[, cell_cols, drop = FALSE], use = "p", method = "bicor")
    cc_df <- data.frame(module = rownames(cc), cell_type = colnames(cc),
                        r = as.vector(cc),
                        p = as.vector(WGCNA::corPvalueStudent(cc, nrow(me))))
    write_tsv(cc_df[order(cc_df$p), ], "deconv_module_correlation", cfg)
  }
}

log_msg("04 完成")

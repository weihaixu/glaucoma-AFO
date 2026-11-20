# =============================================================================
# 06_scrna_cellchat.R —— 课题一（6）：单细胞定位、细胞通讯、拟时序
#
# 目的：把 bulk 得到的"分子轴"落到具体细胞类型上，回答
#       「TGF-β2 / ECM 相关信号是谁发给谁的」
#
# 输入：config 中 scrna 段指定的 Seurat rds，或 CellRanger 输出目录列表
# 输出：data/processed/scRNA_annotated.rds
#       results/tables/sc_celltype_prop.tsv        各样本细胞类型比例
#       results/tables/sc_deg_TM.tsv               TM 区室差异基因
#       results/tables/cellchat_<cond>.tsv         条件内显著 L-R 对
#       results/tables/cellchat_compare.tsv        疾病 vs 对照 通讯差异
#       results/figures/sc_umap_celltype.pdf / cellchat_circle_<cond>.pdf
#
# 依赖：Seurat、harmony、CellChat、monocle3（可选）、SCENIC（可选，走命令行）
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("06_scrna_cellchat", cfg)
load_pkgs(c("ggplot2", "Seurat", "harmony", "CellChat"))
set.seed(cfg_get(cfg, "project.seed", 1))
options(future.globals.maxSize = 8 * 1024^3)

# Seurat v5 用 layer=；v4 的 slot= 在新版已 defunct，先探测 layers 再决定
get_mat <- function(obj, what = "data") {
  a <- Seurat::GetAssay(obj, assay = "RNA")
  layers <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  if (what %in% layers) {
    Seurat::GetAssayData(obj, assay = "RNA", layer = what)
  } else if ("data" %in% layers) {
    Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
  } else if ("counts" %in% layers) {
    Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
  } else {
    stop("RNA assay 没有可用的 data/counts layer: ", paste(layers, collapse = ","))
  }
}

sc_cfg  <- cfg_get(cfg, "scrna", list())
mode    <- sc_cfg$mode %||% "rds"
group_regex <- sc_cfg$group_regex %||% "POAG|glaucoma|treated"

# ------------------------------------------------- 1) 读入 -------------------
# 断点续跑：scRNA_annotated.rds 含 scale.data 等冗余层，约 4.2GB，读入后
# 剩余内存不足以为 CellChat 分配连续矩阵。改为从 10x 导出目录重建精简对象
# （counts + 元数据，带细胞类型/分组），只做 NormalizeData，内存占用小得多。
exp_dir <- sc_cfg$export_dir %||% "data/raw/scRNA_export"
resume_f <- proj_path("data", "processed", "scRNA_annotated.rds")
sc_resumed <- FALSE
if (!isTRUE(sc_cfg$no_resume) && dir.exists(proj_path(exp_dir)) &&
    file.exists(proj_path(exp_dir, "matrix.mtx"))) {
  # 优先用 tools/stream_mtx.R 生成的精简矩阵（52MB，已按细胞类型下采样）
  small <- proj_path("data", "processed", "scRNA_cellchat_small2.rds")
  if (file.exists(small) && file.exists(paste0(small, ".meta.csv"))) {
    m <- readRDS(small)
    md <- utils::read.csv(paste0(small, ".meta.csv"), row.names = 1,
                          stringsAsFactors = FALSE)
    common <- intersect(colnames(m), rownames(md))
    m <- m[, common, drop = FALSE]; md <- md[common, , drop = FALSE]
    log_msg("读入精简矩阵 ", nrow(m), " × ", ncol(m), "（流式下采样产物）")
  } else {
    stop("缺少 scRNA_cellchat_small.rds：先运行 tools/stream_mtx.R 生成精简表达矩阵")
  }
  sc <- Seurat::CreateSeuratObject(counts = m, meta.data = md, project = "scRNA")
  rm(m)
  sc <- Seurat::NormalizeData(sc, verbose = FALSE)
  if (!"cell_type" %in% names(sc@meta.data)) stop("obs.csv 缺 cell_type 列")
  if (!"condition" %in% names(sc@meta.data)) sc$condition <- "control"
  sc_resumed <- TRUE
  log_msg("断点续跑：从 10x 重建精简对象 ", ncol(sc), " 细胞 × ", nrow(sc),
          " 基因（含 cell_type/condition），跳过完整 rds/QC/聚类")
} else if (identical(mode, "rds")) {
  f <- proj_path(sc_cfg$rds %||% "data/processed/scRNA_ref.rds")
  require_file(f, "请把单细胞 Seurat 对象放到该路径，或把 scrna.mode 改成 10x")
  sc <- readRDS(f)
  log_msg("读入 Seurat 对象：", ncol(sc), " 细胞 × ", nrow(sc), " 基因")
} else {
  dirs <- sc_cfg$tenx_dirs %||% list()
  if (!length(dirs)) stop("scrna.tenx_dirs 为空：请填入 CellRanger filtered_feature_bc_matrix 目录")
  objs <- lapply(dirs, function(d) {
    m <- Seurat::Read10X(d)
    o <- Seurat::CreateSeuratObject(m, project = basename(d), min.cells = 3, min.features = 200)
    o$sample <- basename(d)
    o
  })
  sc <- merge(objs[[1]], y = objs[-1], add.cell.ids = vapply(objs, function(o) o$sample[1], character(1)))
}

# ------------------------------------------------- 2) QC 与标准化 ------------
if (!sc_resumed) {
if (!"percent.mt" %in% names(sc@meta.data)) sc[["percent.mt"]] <- Seurat::PercentageFeatureSet(sc, pattern = "^MT-")
qc <- sc_cfg$qc %||% list()
sc <- subset(sc, subset = nFeature_RNA > (qc$min_features %||% 200) &
                            nFeature_RNA < (qc$max_features %||% 8000) &
                            percent.mt < (qc$max_mt %||% 20))
log_msg("QC 后保留 ", ncol(sc), " 个细胞")

sc <- Seurat::NormalizeData(sc, verbose = FALSE)
sc <- Seurat::FindVariableFeatures(sc, nfeatures = 2000, verbose = FALSE)
sc <- Seurat::ScaleData(sc, verbose = FALSE)
sc <- Seurat::RunPCA(sc, npcs = 30, verbose = FALSE)

# 按样本做 Harmony，去掉个体/批次效应（不加这一步，聚类会按供体分开）
if ("sample" %in% names(sc@meta.data) && length(unique(sc$sample)) > 1) {
  sc <- harmony::RunHarmony(sc, group.by.vars = "sample", verbose = FALSE)
  red <- "harmony"
} else {
  red <- "pca"
}
sc <- Seurat::FindNeighbors(sc, reduction = red, dims = 1:20, verbose = FALSE)
sc <- Seurat::FindClusters(sc, resolution = sc_cfg$resolution %||% 0.6, verbose = FALSE)
sc <- Seurat::RunUMAP(sc, reduction = red, dims = 1:20, verbose = FALSE)

# ------------------------------------------------- 3) 细胞注释 --------------
# 有 cell_type 列就直接用；否则按 marker 打分做粗注释，再人工核对
if (!"cell_type" %in% names(sc@meta.data)) {
  markers <- list(
    TM        = c("MYOC", "MGP", "TIMD4", "ANGPTL7", "SFRP1"),
    JCT       = c("MYOC", "MGP", "CLU", "SPARCL1", "PECAM1"),
    SC        = c("PECAM1", "KDR", "FLT1", "CDH5", "VWF"),
    immune    = c("PTPRC", "CD68", "LYZ", "HLA-DRA", "AIF1"),
    fibroblast= c("DCN", "LUM", "COL1A1", "PDGFRA"),
    endothelial = c("PECAM1", "VWF", "CLDN5"),
    pericyte  = c("RGS5", "ACTA2", "PDGFRB"),
    RGC       = c("RBPMS", "NEFL", "SNCG", "THY1")
  )
  markers <- lapply(markers, function(g) intersect(g, rownames(sc)))
  markers <- markers[lengths(markers) > 0]
  sc <- Seurat::AddModuleScore(sc, features = markers, name = "score", seed = 1)
  score_cols <- paste0("score", seq_along(markers))
  names(score_cols) <- names(markers)
  sm <- sc@meta.data[, score_cols, drop = FALSE]
  colnames(sm) <- names(markers)
  sc$cell_type <- colnames(sm)[max.col(as.matrix(sm), ties.method = "first")]
  log_msg("按 marker 打分完成粗注释，请用 DotPlot 人工复核后写入 cell_type 列")
}

if (!"condition" %in% names(sc@meta.data)) {
  cols <- sc@meta.data[, seq_len(min(4, ncol(sc@meta.data))), drop = FALSE]
  txt <- do.call(paste, c(as.list(cols), sep = " | "))
  sc$condition <- ifelse(grepl(group_regex, txt, ignore.case = TRUE), "case", "control")
}
saveRDS(sc, proj_path("data", "processed", "scRNA_annotated.rds"))

# UMAP
g_umap <- Seurat::DimPlot(sc, group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggplot2::labs(title = "Cell types")
save_fig(g_umap, "sc_umap_celltype", cfg, 7, 6)

# 各样本细胞类型比例
prop <- as.data.frame.matrix(table(sc$sample, sc$cell_type))
prop <- prop / rowSums(prop)
prop$sample <- rownames(prop)
write_tsv(prop, "sc_celltype_prop", cfg)
}   # <- else 分支（非断点续跑）到这里结束

# ------------------------------------------------- 4) TM 区室差异基因 ------
tm_types <- intersect(sc_cfg$tm_celltypes %||% c("TM", "JCT"), unique(sc$cell_type))
if (length(tm_types) && length(unique(sc$condition)) > 1) {
  sub <- subset(sc, subset = cell_type %in% tm_types)
  Idents(sub) <- "condition"
  deg <- Seurat::FindMarkers(sub, ident.1 = "case", ident.2 = "control",
                             logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
  deg$gene <- rownames(deg)
  write_tsv(deg, "sc_deg_TM", cfg)
  log_msg("TM 区室差异基因 ", sum(deg$p_val_adj < 0.05), " 个（FDR<0.05）")
}

# ------------------------------------------------- 5) CellChat -------------
# 大类 52k 细胞的 data 层 1.6GB，createCellChat/计算通讯概率会需要 >4GB 连续内存，
# 在 8-16GB 的机器上会失败。按细胞类型分层下采样到每类至多 N 个（保留全部类型），
# 这是 CellChat 官方对 >10k 细胞数据的建议做法。
run_cellchat <- function(obj, label) {
  max_per_type <- sc_cfg$cellchat_max_per_type %||% 300L
  set.seed(1)
  keep <- unlist(lapply(split(colnames(obj), obj$cell_type), function(x) {
    if (length(x) > max_per_type) sample(x, max_per_type) else x
  }), use.names = FALSE)
  log_msg("CellChat[", label, "] 下采样到 ", length(keep), " 细胞（每类 ≤", max_per_type,
          "，共 ", length(unique(obj$cell_type)), " 类）")
  # 必须用 GetAssayData 而不是 a@layers$data：Seurat v5 的 layer 是裸稀疏矩阵，
  # 不带 dimnames（基因名在 Assay@features、细胞名在 Assay@cells）。
  # 直接取 layer 会让 CellChatDB 一个基因都匹配不上，data.signaling 变 0 行。
  dat <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
  stopifnot(length(rownames(dat)) > 0, length(colnames(dat)) > 0)
  data.input <- dat[, keep, drop = FALSE]
  lab <- as.character(obj$cell_type)[match(keep, colnames(obj))]
  meta <- data.frame(labels = lab, row.names = keep)
  cc <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "labels")
  # 某些安装下 CellChatDB.human 惰性加载是空的（subsetData 后 data.signaling 0 行），
  # 必须 data() 显式加载再塞进 @DB
  suppressWarnings(utils::data("CellChatDB.human", package = "CellChat", envir = environment()))
  cc@DB <- CellChatDB.human
  cc <- CellChat::subsetData(cc)
  # presto（CRAN）提供快速 Wilcoxon；没装的话退回标准实现（52k 细胞约慢 20 倍）
  cc <- CellChat::identifyOverExpressedGenes(cc,
         do.fast = requireNamespace("presto", quietly = TRUE))
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, raw.use = TRUE, population.size = TRUE)
  cc <- CellChat::filterCommunication(cc, min.cells = 10)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)
  cc
}

cc_list <- list()
for (cond in unique(sc$condition)) {
  obj <- subset(sc, subset = condition == cond)
  if (length(unique(obj$cell_type)) < 2) { log_msg(cond, " 细胞类型不足 2 类，跳过 CellChat"); next }
  cc <- run_cellchat(obj, cond)
  cc_list[[cond]] <- cc

  lr <- CellChat::subsetCommunication(cc)
  if (nrow(lr)) write_tsv(lr, paste0("cellchat_", cond), cfg)

  f <- proj_path(cfg_get(cfg, "project.outdir", "results"), "figures",
                 paste0("cellchat_circle_", cond, ".pdf"))
  grDevices::pdf(f, width = 7, height = 7)
  CellChat::netVisual_circle(cc@net$count, vertex.weight = as.numeric(table(cc@idents)),
                             weight.scale = TRUE, label.edge = FALSE,
                             title.name = paste0("Interaction count — ", cond))
  grDevices::dev.off()
}

if (length(cc_list) == 2) {
  merged <- CellChat::mergeCellChat(cc_list, add.names = names(cc_list))
  cmp <- CellChat::rankNet(merged, mode = "comparison", stacked = TRUE,
                           do.stat = TRUE, return.data = TRUE)
  if (!is.null(cmp$signaling.contribution)) {
    write_tsv(as.data.frame(cmp$signaling.contribution), "cellchat_compare", cfg)
  } else {
    write_tsv(as.data.frame(cmp[[1]]), "cellchat_compare", cfg)
  }

  f2 <- proj_path(cfg_get(cfg, "project.outdir", "results"), "figures", "cellchat_compare_rank.pdf")
  grDevices::pdf(f2, width = 6, height = 6)
  print(CellChat::rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE))
  grDevices::dev.off()
  saveRDS(merged, proj_path("data", "processed", "cellchat_merged.rds"))
}

# ------------------------------------------------- 6) 拟时序（可选） -------
if (requireNamespace("monocle3", quietly = TRUE) && length(tm_types)) {
  sub <- subset(sc, subset = cell_type %in% tm_types)
  cds <- suppressWarnings(monocle3::new_cell_data_set(
    get_mat(sub, "counts"),
    cell_metadata = sub@meta.data,
    gene_metadata = data.frame(gene_short_name = rownames(sub), row.names = rownames(sub))))
  cds <- monocle3::preprocess_cds(cds, num_dim = 30)
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
  cds <- monocle3::cluster_cells(cds)
  cds <- monocle3::learn_graph(cds, use_partition = FALSE)
  if (interactive()) {
    cds <- monocle3::order_cells(cds)          # 交互式点选根节点
  } else {                                      # 批处理：以对照组的第一个细胞为根
    root <- colnames(sub)[which(sub$condition == "control")[1]]
    if (is.na(root)) root <- colnames(sub)[1]
    cds <- monocle3::order_cells(cds, root_cells = root)
  }
  g_ps <- monocle3::plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
  save_fig(g_ps, "sc_pseudotime_TM", cfg, 7, 5.5)
  saveRDS(cds, proj_path("data", "processed", "monocle3_TM.rds"))
  log_msg("提示：order_cells 默认交互式，命令行批处理请指定 root_cells")
}

# ------------------------------------------------- 7) SCENIC 导出 ----------
# R 端不跑 pySCENIC（Arboreto 依赖）。15_results_ledger.R 只读
# results/tables/scenic_regulon_diff.tsv（由 tools/scenic_regulon_count.R 生成），
# 并不读这里的表达矩阵，因此默认不导出：as.matrix() 一个 59413×8049 的稀疏矩阵
# 需要 ~2-4GB 连续内存并写出数 GB 的 TSV，是 06 的主要 OOM 来源。
# 只有确实要在命令行跑 pySCENIC 时，才把 scrna.export_scenic_input 设为 true，
# 此时按基因分块写出，避免一次性 as.matrix()。
if (isTRUE(sc_cfg$export_scenic_input)) {
  scen_dir <- res_path(cfg, "scenic")
  if (!dir.exists(scen_dir)) dir.create(scen_dir, recursive = TRUE, showWarnings = FALSE)
  cnt <- get_mat(sc, "counts")
  keep <- rowSums(cnt > 0) >= max(1, 0.01 * ncol(cnt))
  cnt <- cnt[keep, , drop = FALSE]
  log_msg("SCENIC 输入导出：过滤后 ", nrow(cnt), " 基因 × ", ncol(cnt), " 细胞（分块写）")
  ef <- file.path(scen_dir, "expr_matrix.tsv")
  con <- file(ef, open = "wt")
  writeLines(paste("", paste(colnames(cnt), collapse = "\t")), con)
  chunk <- 2000L
  for (i in seq_len(ceiling(nrow(cnt) / chunk))) {
    idx <- ((i - 1) * chunk + 1):min(i * chunk, nrow(cnt))
    blk <- as.matrix(cnt[idx, , drop = FALSE])
    write.table(blk, con, sep = "\t", quote = FALSE, col.names = FALSE,
                row.names = rownames(cnt)[idx])
  }
  close(con)
  utils::write.table(data.frame(cell = colnames(sc), cell_type = sc$cell_type),
                     file.path(scen_dir, "cell_info.tsv"), sep = "\t",
                     quote = FALSE, row.names = FALSE)
  log_msg("SCENIC 输入已导出到 ", scen_dir)
} else {
  log_msg("SCENIC 表达矩阵导出已跳过（scrna.export_scenic_input=false）；如需 pySCENIC 输入请设为 true")
}

log_msg("06 完成")

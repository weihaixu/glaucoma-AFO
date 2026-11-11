# =============================================================================
# 00_check_env.R —— 环境自检
#
# 装完 R 和依赖包后，先跑这个，再跑 --demo。它把"环境问题"和"数据问题"
# 分开报出来，省得在 01 里对着一个下载报错猜半天。
#
# 用法（在 glaucoma_R_pipeline/ 目录下）：
#   Rscript 00_check_env.R
# 输出：
#   控制台表格 + logs/env_check.log + results/tables/env_check.tsv
# =============================================================================

options(warn = 1, timeout = 600)

root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(root, "config", "config.yml"))) {
  stop("请在 glaucoma_R_pipeline/ 目录下运行（当前：", root, "）", call. = FALSE)
}

rows <- list()
add <- function(group, item, ok, detail = "") {
  rows[[length(rows) + 1]] <<- data.frame(
    group = group, item = item,
    status = if (isTRUE(ok)) "OK" else "MISSING",
    detail = as.character(detail), stringsAsFactors = FALSE)
}

# ------------------------------------------------------------ 1) 基础环境 ----
add("基础", "R 版本 >= 4.2", getRversion() >= "4.2.0",
    paste("当前", as.character(getRversion()), R.version$platform))
add("基础", "项目根目录", TRUE, root)
add("基础", "cairo PDF 设备（中文字体）", isTRUE(capabilities("cairo")),
    if (isTRUE(capabilities("cairo"))) "可用" else "缺 cairo，中文 PDF 可能出豆腐块")
add("基础", "libcurl（https 下载）", isTRUE(capabilities("libcurl")), "")
add("基础", "可用内存", TRUE,
    tryCatch({
      if (requireNamespace("ps", quietly = TRUE)) {
        paste0(round(ps::ps_system_memory()$total / 1024^3, 1), " GB")
      } else {
        "未装 ps 包，无法读取（单细胞建议 >= 32 GB）"
      }
    }, error = function(e) "未知"))

# ------------------------------------------------------------ 2) 目录与写权限 ----
dirs <- c("data/raw", "data/annotation", "data/processed", "results/tables",
          "results/figures", "results/scenic", "logs", "tools")
for (d in dirs) {
  p <- file.path(root, d)
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  add("目录", d, dir.exists(p), "")
}
wtest <- file.path(root, "logs", ".write_test")
ok_w <- tryCatch({ writeLines("x", wtest); file.remove(wtest); TRUE },
                 error = function(e) FALSE)
add("目录", "logs/ 可写", ok_w, if (!ok_w) "检查目录权限或杀毒软件锁定" else "")

# ------------------------------------------------------------ 3) R 包 ----------
pkg_spec <- list(
  "基础"   = c("yaml", "data.table", "ggplot2", "ggrepel", "reshape2", "openxlsx"),
  "统计"   = c("limma", "sva", "preprocessCore", "BiocParallel"),
  "富集"   = c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "GSVA", "msigdbr",
               "ComplexHeatmap"),
  "机器学习" = c("glmnet", "e1071", "ranger", "Boruta", "pROC", "rms"),
  "单细胞"  = c("Seurat", "SeuratObject", "harmony", "CellChat", "monocle3"),
  "遗传"   = c("coloc", "susieR", "Rsamtools", "TwoSampleMR", "MRPRESSO",
               "ieugwasr", "genetics.binaRies")
)

missing_pkgs <- character(0)
for (grp in names(pkg_spec)) {
  for (p in pkg_spec[[grp]]) {
    ok <- requireNamespace(p, quietly = TRUE)
    ver <- if (ok) tryCatch(as.character(utils::packageVersion(p)),
                            error = function(e) "?") else ""
    if (!ok) missing_pkgs <- c(missing_pkgs, p)
    add(grp, p, ok, ver)
  }
}

# ------------------------------------------------------------ 4) 外部工具 ------
tool_spec <- list(
  "外部工具" = list(
    list("Rscript", "Rscript", "run_all.R 用子进程调用，必须有"),
    list("python",  "python",  "ldsc / FOCUS 用"),
    list("plink",   "plink",   "FUSION 与 LD 剪枝用（1.9 或 2.0）"),
    list("plink2",  "plink2",  "可选，PLINK2 更快"),
    list("smr",     "smr",     "官方 SMR 二进制，复核 08 的近似结果")
  )
)
for (t in tool_spec[["外部工具"]]) {
  p <- Sys.which(t[[2]])
  add("外部工具", paste0(t[[1]], "（", t[[2]], "）"), nzchar(p),
      if (nzchar(p)) p else t[[3]])
}
# ldsc 需要 Python 2.7，单独提示（这是最常见的卡点）
ldsc_dir <- file.path(root, "tools", "ldsc")
add("外部工具", "tools/ldsc", dir.exists(ldsc_dir),
    if (dir.exists(ldsc_dir)) "存在（注意：官方 ldsc 只支持 Python 2.7）" else "未下载")

# ------------------------------------------------------------ 5) 配置文件 ------
# 注意：这里必须用 helpers 的 read_config，不能用裸的 yaml::read_yaml——
# 在 GBK locale（中文 Windows 默认）下 libyaml 读到中文注释会【静默】返回
# 空文档不报错；read_config 的第二层会把注释剥掉再解析。之前这里误报
# "config.yml 不可解析"就是这个原因。
cfg <- tryCatch({
  source(file.path(root, "R", "helpers.R"))
  read_config(file.path(root, "config", "config.yml"))
}, error = function(e) NULL)
add("配置", "config/config.yml 可解析", !is.null(cfg), "")
if (!is.null(cfg)) {
  add("配置", "run.demo = false（真实数据模式）",
      !isTRUE(cfg$run$demo),
      if (isTRUE(cfg$run$demo)) "当前是 demo 模式，跑真实数据前改成 false" else "")
  add("配置", "project.outdir", TRUE, as.character(cfg$project$outdir))
}

# ------------------------------------------------------------ 6) 数据文件 ------
if (!is.null(cfg)) {
  want <- c(
    "GWAS 汇总统计" = cfg$genetics$gwas$path,
    "MR 暴露工具变量" = cfg$genetics$mr$exposure,
    "单细胞对象" = cfg$scrna$rds
  )
  for (nm in names(want)) {
    p <- want[[nm]]
    if (is.null(p) || identical(p, "none") || !nzchar(p)) next
    full <- file.path(root, p)
    add("数据", nm, file.exists(full),
        if (file.exists(full)) paste0(round(file.size(full) / 1024^2, 1), " MB") else p)
  }
  gtex_dir <- file.path(root, cfg$genetics$eqtl$gtex_dir)
  for (ts in cfg$genetics$eqtl$tissues) {
    f <- file.path(gtex_dir, paste0(ts, ".v8.EUR.allpairs.txt.gz"))
    add("数据", paste0("GTEx ", ts), file.exists(f),
        if (file.exists(f)) "" else "缺文件（或 .tbi 索引）")
  }
}

# ------------------------------------------------------------ 7) 汇总输出 ------
rep <- do.call(rbind, rows)
n_bad <- sum(rep$status != "OK")

show <- rep[rep$status != "OK" | rep$group %in% c("基础", "配置"), ]
cat("\n================ 环境自检 ================\n")
print(show, row.names = FALSE)
cat("\n总计", nrow(rep), "项，未通过", n_bad, "项。\n")

if (n_bad) {
  cat("\n未通过项按组统计：\n")
  print(table(rep$group[rep$status != "OK"]))
}

if (length(missing_pkgs)) {
  bioc <- c("limma", "sva", "preprocessCore", "BiocParallel", "clusterProfiler",
            "org.Hs.eg.db", "ReactomePA", "GSVA", "ComplexHeatmap", "coloc",
            "susieR", "Rsamtools", "monocle3")
  cat("\n缺失的 R 包：", paste(missing_pkgs, collapse = ", "), "\n")
  cat("一键安装（Bioc 与 CRAN 分开）：\n")
  cr <- setdiff(missing_pkgs, bioc)
  bc <- intersect(missing_pkgs, bioc)
  if (length(cr)) cat('  install.packages(c(', paste0('"', cr, '"', collapse = ", "), '))\n')
  if (length(bc)) cat('  BiocManager::install(c(', paste0('"', bc, '"', collapse = ", "), '), ask = FALSE, update = FALSE)\n')
  cat("  remotes::install_github('jinworks/CellChat')\n")
  cat("  remotes::install_github('cole-trapnell-lab/monocle3')\n")
}

out_dir <- file.path(root, "results", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(rep, file.path(out_dir, "env_check.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
log_dir <- file.path(root, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
utils::capture.output(print(rep, row.names = FALSE),
                      file = file.path(log_dir, "env_check.log"))

cat("\n报告已写入：results/tables/env_check.tsv, logs/env_check.log\n")
cat("下一步：", if (n_bad == 0) "Rscript 99_demo_data.R && Rscript run_all.R --demo",
    "\n", sep = "")
if (n_bad > 0) cat("先把上面 MISSING 的项补齐，再跑 demo。\n")

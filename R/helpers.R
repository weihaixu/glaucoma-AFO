# =============================================================================
# R/helpers.R —— 公共函数库
# 所有分析脚本第一行 source("R/helpers.R")，不要重复定义函数。
# 设计原则：不依赖 tidyverse 管道，全部 base R + 少量常用包，降低版本踩坑概率。
# =============================================================================

# ---------------------------------------------------------------- 基础工具 ----

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# 从当前目录向上找含 config/config.yml 的目录作为项目根
find_root <- function(start = getwd(), max_up = 5) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_up)) {
    if (file.exists(file.path(d, "config", "config.yml"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("找不到项目根目录（含 config/config.yml）。请先 setwd() 到 glaucoma_R_pipeline/ 再运行。")
}

PROJ_ROOT <- find_root()

# 支持用环境变量 GLAUCOMA_CONFIG 指定配置文件（run_all.R --demo 用这个切换）
read_config <- function(path = NULL) {
  if (is.null(path)) {
    env_path <- Sys.getenv("GLAUCOMA_CONFIG", "")
    path <- if (nzchar(env_path)) env_path else file.path(PROJ_ROOT, "config", "config.yml")
  }
  if (!file.exists(path)) stop("配置文件不存在：", path, call. = FALSE)
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("需要 yaml 包：install.packages('yaml')")
  }
  # 两层解析。第一层直接解析（macOS/Linux 与 UTF-8 locale 下没问题）；
  # 第二层是给中文 Windows（GBK locale）的兜底：实测 libyaml 读到配置里的
  # 中文注释时会【静默】返回空文档、不报错——所有 cfg_get 落回默认值，
  # 脚本"能跑"但全是默认参数，比崩溃更危险。剥掉注释后剩余键值都是 ASCII，
  # 解析必然成功。注意：值里不要写中文（demo/正式配置的值都是 ASCII）。
  txt <- paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  cfg <- tryCatch(yaml::read_yaml(text = txt), error = function(e) NULL)
  if (is.null(cfg$project) || is.null(cfg$project$name)) {
    ln <- readLines(path, encoding = "UTF-8", warn = FALSE)
    ln <- sub("#.*$", "", ln)               # 剥掉整行与行尾注释
    txt2 <- paste(ln[nzchar(trimws(ln))], collapse = "\n")
    cfg <- tryCatch(yaml::read_yaml(text = txt2), error = function(e) NULL)
  }
  if (is.null(cfg$project) || is.null(cfg$project$name)) {
    stop("配置解析结果为空（", path, "）：请检查文件编码（应为无 BOM 的 UTF-8），",
         "或确认 YAML 语法无误（可粘贴到 yamllint 验证）。", call. = FALSE)
  }
  cfg$.root <- PROJ_ROOT
  cfg <- coerce_config_numbers(cfg)
  cfg
}

# YAML 1.1 的坑：浮点必须带小数点，`1e-4` / `5e-8` 会被 libyaml 当成【字符串】。
# 后果很隐蔽：`nsnps * p1` 报 non-numeric（coloc 直接失败），
# 而 `p < 5e-8` 这种比较会退化成字符串字典序比较——不报错，结论悄悄变。
# 所以在读入后统一把"看起来是数字的字符串"还原成数值，一次修好所有配置。
coerce_config_numbers <- function(x) {
  if (is.list(x)) return(lapply(x, coerce_config_numbers))
  if (!is.character(x) || length(x) != 1L) return(x)
  if (grepl("^[-+]?(\\d+\\.?\\d*|\\.\\d+)([eE][-+]?\\d+)?$", x)) {
    v <- suppressWarnings(as.numeric(x))
    if (!is.na(v)) return(v)
  }
  x
}

# 嵌套取值：cfg_get(cfg, "bulk.deg.fdr", 0.05)
cfg_get <- function(cfg, key, default = NULL) {
  parts <- strsplit(key, ".", fixed = TRUE)[[1]]
  v <- cfg
  for (p in parts) {
    if (!is.list(v) || is.null(v[[p]])) return(default)
    v <- v[[p]]
  }
  if (is.null(v)) default else v
}

proj_path <- function(...) file.path(PROJ_ROOT, ...)

# 结果目录路径：跟随 project.outdir（demo 用 results_demo，真实数据用 results）
res_path <- function(cfg, ...) {
  proj_path(cfg_get(cfg, "project.outdir", "results"), ...)
}

ensure_dir <- function(...) {
  d <- proj_path(...)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# 日志：同时打印与写入 logs/<script>.log
LOG_FILE <- NULL
log_init <- function(script_name, cfg = NULL) {
  dir <- if (!is.null(cfg)) proj_path(cfg_get(cfg, "project.log_dir", "logs")) else proj_path("logs")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  LOG_FILE <<- file.path(dir, paste0(script_name, ".log"))
  log_msg("===== ", script_name, " 开始 @ ", format(Sys.time()), " =====")
}

log_msg <- function(...) {
  txt <- paste0(...)
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", txt)
  message(line)
  if (!is.null(LOG_FILE)) cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
  invisible(NULL)
}

# 显式检查输入文件，缺失时给出「该放什么文件、放哪」的提示，而不是让 read 报错
require_file <- function(path, hint = "") {
  if (!file.exists(path)) {
    stop(sprintf("缺少输入文件：%s\n  %s", path, hint), call. = FALSE)
  }
  invisible(TRUE)
}

# 已存在且 overwrite=FALSE 时跳过重算
step_done <- function(outfile, cfg) {
  file.exists(outfile) && !isTRUE(cfg_get(cfg, "run.overwrite", FALSE))
}

# ---------------------------------------------------------------- 包管理 ------

load_pkgs <- function(pkgs, install = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && install) {
    bioc <- c("GEOquery", "limma", "sva", "Biobase", "clusterProfiler",
              "org.Hs.eg.db", "ReactomePA", "GSVA", "msigdbr", "preprocessCore",
              "ComplexHeatmap", "coloc", "snpStats", "BiocParallel")
    for (m in missing) {
      if (m %in% bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(m, ask = FALSE, update = FALSE)
      } else {
        install.packages(m)
      }
    }
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  }
  if (length(missing)) {
    stop("缺少 R 包（先运行 install = TRUE 或手动安装）：\n  ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  for (pk in pkgs) suppressPackageStartupMessages(library(pk, character.only = TRUE))
  invisible(TRUE)
}

# ---------------------------------------------------------------- 出图/出表 ----

theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey90"),
      axis.text = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      legend.key = ggplot2::element_blank()
    )
}

# 同时存 PDF（投稿）与 PNG（预览）
save_fig <- function(plot, name, cfg, width = NULL, height = NULL, dir = "figures") {
  d <- ensure_dir(cfg_get(cfg, "project.outdir", "results"), dir)
  w <- width %||% cfg_get(cfg, "project.fig_width", 7)
  h <- height %||% cfg_get(cfg, "project.fig_height", 5)
  dpi <- cfg_get(cfg, "project.dpi", 300)
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(file.path(d, paste0(name, ".pdf")), plot, width = w, height = h, device = pdf_dev)
  # PNG 必须走 cairo：Windows 默认的 type="windows" 设备会把字符串转成本地编码
  # （这台机器是 C locale），UTF-8 的 ± – — 等字符在图里变成 "B1" 这种残骸；
  # cairo 设备按 UTF-8 渲染，PDF/PNG 两边一致
  if (capabilities("cairo")) {
    ggplot2::ggsave(file.path(d, paste0(name, ".png")), plot, width = w, height = h,
                    dpi = dpi, type = "cairo")
  } else {
    ggplot2::ggsave(file.path(d, paste0(name, ".png")), plot, width = w, height = h, dpi = dpi)
  }
  invisible(file.path(d, paste0(name, ".pdf")))
}

write_tsv <- function(df, name, cfg, dir = "tables", row.names = FALSE) {
  d <- ensure_dir(cfg_get(cfg, "project.outdir", "results"), dir)
  f <- file.path(d, paste0(name, ".tsv"))
  utils::write.table(df, f, sep = "\t", quote = FALSE, row.names = row.names)
  log_msg("写出表格：", f, " (", nrow(df), " 行)")
  invisible(f)
}

# ---------------------------------------------------------------- GEO 相关 ----

# GSE27276 -> GSE27nnn
geo_series_dir <- function(gse) paste0(substr(gse, 1, nchar(gse) - 3), "nnn")
# GPL6883 -> GPLnnn
geo_platform_dir <- function(gpl) paste0(substr(gpl, 1, nchar(gpl) - 3), "nnn")

# --- 下载：多通道自适应 -----------------------------------------------------
# 为什么要写这么长：在 Windows 上 utils::download.file 会栽在证书吊销检查上
# （schannel 连不上 CRL 服务器），报 "SSL connect error" 或 wininet 的
# "不能连接到吊销服务器"，同一个 URL 换 `curl --ssl-no-revoke` 就通。
# 这种报错看起来像"被墙"，极易误判成网络不可用、进而放弃跑真实数据。
# 所以按顺序试几条通道，谁通用谁，并把成功的那条写进日志——换台机器再跑时，
# 日志能直接说明哪条路走得通。
dl_transports <- function() {
  out <- c("R-libcurl", "R-wininet")
  if (nzchar(Sys.which("curl"))) out <- c(out, "curl-no-revoke", "curl-direct")
  out
}

dl_once <- function(url, dest, how, timeout = 3600) {
  unlink(dest)
  if (how == "R-libcurl") {
    tryCatch(utils::download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl"),
             error = function(e) NULL)
  } else if (how == "R-wininet") {
    tryCatch(utils::download.file(url, dest, mode = "wb", quiet = TRUE, method = "wininet"),
             error = function(e) NULL)
  } else {
    args <- c("--fail", "--location", "--silent", "--show-error",
              "--retry", "3", "--retry-delay", "2",
              "--connect-timeout", "60", "--max-time", timeout,
              "--output", dest)
    # 只在不做吊销检查的那条通道上加这个开关；老版本 curl 不认识会直接报错，
    # 所以不能无条件加。
    if (how == "curl-no-revoke" && .Platform$OS.type == "windows") {
      args <- c(args, "--ssl-no-revoke")
    }
    if (how == "curl-direct") args <- c(args, "--noproxy", "*")
    args <- c(args, url)
    suppressWarnings(system2(unname(Sys.which("curl")), args, stdout = TRUE, stderr = TRUE))
  }
  file.exists(dest) && file.size(dest) > 0
}

# 远端文件字节数（HEAD）。拿不到就返回 NA —— 网络不通时不该阻塞已有缓存。
dl_remote_size <- function(url, timeout = 60) {
  cur <- Sys.which("curl")
  if (!nzchar(cur)) return(NA_real_)
  args <- c("-sS", "--noproxy", "*", "--head",
            "--connect-timeout", "30", "--max-time", timeout)
  if (.Platform$OS.type == "windows") args <- c(args, "--ssl-no-revoke")
  h <- suppressWarnings(system2(unname(cur), c(args, url), stdout = TRUE, stderr = FALSE))
  ln <- grep("^[Cc]ontent-[Ll]ength:", h, value = TRUE)
  if (!length(ln)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(sub("^[^:]*:", "", ln[1]))))
}

download_if_missing <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 0) {
    # 缓存也要核对：曾有一次 soft 下载被中断，残片以 <GPL>_annot.txt.gz 的名字
    # 留下，后续所有运行都把它当注释表读，探针只映射上 1/4，而日志里看不出来。
    want <- dl_remote_size(url)
    if (is.na(want) || file.size(dest) == want) return(invisible(dest))
    log_msg("缓存不完整（", file.size(dest), " / ", want, " 字节），重新下载：", basename(dest))
    unlink(dest)
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  log_msg("下载：", url)
  want <- dl_remote_size(url)
  if (!is.na(want)) log_msg("  远端大小 ", format(want, big.mark = ","), " 字节")

  tried <- character(0)
  for (how in dl_transports()) {
    if (dl_once(url, dest, how)) {
      if (!is.na(want) && file.size(dest) != want) {
        log_msg("  通道 ", how, " 截断（", file.size(dest), " / ", want, " 字节），丢弃")
        unlink(dest)
        tried <- c(tried, paste0(how, "(截断)"))
        next
      }
      log_msg("  通道 ", how, " 成功（", format(file.size(dest), big.mark = ","), " 字节）")
      return(invisible(dest))
    }
    tried <- c(tried, how)
    log_msg("  通道 ", how, " 失败")
  }
  unlink(dest)
  stop("下载失败：", url, "\n  已尝试的通道：", paste(tried, collapse = " / "),
       "\n  若只有 R-libcurl 与 R-wininet 报 SSL/吊销 错误，说明本机需要 ",
       "curl 命令（Windows 10 以上自带）或检查代理设置。", call. = FALSE)
}

# 手动解析 series matrix（比 GEOquery 更稳，且不依赖 Bioconductor 版本）
# 返回 list(expr = matrix, pheno = data.frame)
read_series_matrix <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)

  meta_idx <- grep("^!Sample_", lines)
  meta <- list()
  for (i in meta_idx) {
    f <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    key <- sub("^!Sample_", "", f[1])
    vals <- gsub('^"|"$', "", f[-1])
    if (key %in% names(meta)) {           # characteristics_ch1 常出现多次
      k <- 2
      while (paste0(key, ".", k) %in% names(meta)) k <- k + 1
      key <- paste0(key, ".", k)
    }
    meta[[key]] <- vals
  }
  pheno <- as.data.frame(meta, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(pheno$geo_accession)) rownames(pheno) <- pheno$geo_accession

  i1 <- grep("^!series_matrix_table_begin", lines)
  i2 <- grep("^!series_matrix_table_end", lines)
  if (length(i1) != 1 || length(i2) != 1) stop("series matrix 表区定位失败：", path)
  block <- lines[(i1 + 1):(i2 - 1)]
  tab <- utils::read.delim(text = block, header = TRUE, check.names = FALSE,
                           quote = "\"", na.strings = c("null", "NA", ""),
                           stringsAsFactors = FALSE)
  id_col <- names(tab)[1]
  expr <- as.matrix(tab[, -1, drop = FALSE])
  rownames(expr) <- as.character(tab[[id_col]])
  storage.mode(expr) <- "double"
  colnames(expr) <- sub("\\.\\d+$", "", colnames(expr))   # 重复列名去后缀
  list(expr = expr, pheno = pheno, id_col = id_col)
}

# 平台 -> Bioconductor 注释包。只在 GEO 拿不到 annot 时才用：
# GPL5175（Affy Human Exon 1.0 ST）的 soft 有 4.5 GB，annot 是 404，
# 而 huex10sttranscriptcluster.db 只有 12 MB，探针号就是转录簇 ID。
gpl_probe_db <- function(gpl) {
  reg <- c(GPL5175 = "huex10sttranscriptcluster.db",
           GPL570  = "hgu133plus2.db",
           GPL6883 = "illuminaHumanv3.db")
  if (!gpl %in% names(reg)) return(NULL)
  reg[[gpl]]
}

# 用 Bioconductor 注释包产出 data.frame(probe, symbol)
annotation_from_db <- function(gpl) {
  pkg <- gpl_probe_db(gpl)
  if (is.null(pkg)) return(NULL)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_msg("注释包 ", pkg, " 未安装，跳过")
    return(NULL)
  }
  db <- getExportedValue(pkg, pkg)
  keys <- AnnotationDbi::keys(db, keytype = "PROBEID")
  m <- suppressMessages(AnnotationDbi::select(db, keys = keys,
                                              columns = "SYMBOL", keytype = "PROBEID"))
  m <- m[!is.na(m$SYMBOL) & m$SYMBOL != "", , drop = FALSE]
  n_multi <- sum(duplicated(m$PROBEID))
  # 一号多基因（如 2315554 -> TTLL10 + MIR200B/200A/429）取注释里第一条，
  # 即平台的 primary assignment；按字母排序会挑到 MIR，那是错的。
  m <- m[!duplicated(m$PROBEID), , drop = FALSE]
  log_msg("注释来源：Bioconductor ", pkg, "（", nrow(m), " 探针；",
          n_multi, " 个探针一对多，取主注释）")
  data.frame(probe = as.character(m$PROBEID), symbol = m$SYMBOL,
             stringsAsFactors = FALSE)
}

# 读取 GPL 注释表，返回 data.frame(probe, symbol)
# 优先本地 data/annotation/<GPL>_annot.txt.gz，其次自动下载注释文件
# 顺序按可用性排：GEO 上 <GPL>.annot.gz 是绝大多数平台唯一存在的注释
# （旧代码只试 <GPL>_full.txt.gz，该文件对 GPL2507/GPL570/GPL6883 全是 404），
# 最后才回退 soft —— soft 里平台表与样本表混在一起，只有平台表那段能用。
read_gpl_annotation <- function(gpl, cfg) {
  local <- proj_path("data", "annotation", paste0(gpl, "_annot.txt.gz"))
  if (!file.exists(local)) {
    base <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/platforms/",
                   geo_platform_dir(gpl), "/", gpl)
    cand <- c(paste0(base, "/annot/", gpl, ".annot.gz"),
              paste0(base, "/annot/", gpl, "_full.txt.gz"))
    # soft 是最后的兜底，但有的平台 soft 大得离谱（GPL5175 有 4.5 GB），
    # 而 Bioconductor 有等价注释包时就完全没必要下 —— 直接跳过。
    if (is.null(gpl_probe_db(gpl))) {
      cand <- c(cand, paste0(base, "/soft/", gpl, "_family.soft.gz"))
    }
    for (u in cand) {
      ok <- tryCatch({download_if_missing(u, local); TRUE},
                     error = function(e) { log_msg("注释源不可用：", basename(u)); FALSE })
      if (ok) { log_msg("注释来源：", basename(u)); break }
    }
    # 第三条路：GEO 上没有 annot、soft 又大到不合理时，用 Bioconductor 注释包
    if (!file.exists(local)) {
      a <- annotation_from_db(gpl)
      if (!is.null(a)) return(a)
      stop("拿不到 ", gpl, " 的注释文件（GEO 无 annot，Bioconductor 也无对应包）",
           call. = FALSE)
    }
  }
  con <- gzfile(local, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  # soft 里若混进多个表块，只取第一个 !platform_table_begin..end
  b1 <- grep("^!platform_table_begin", lines)
  e1 <- grep("^!platform_table_end", lines)
  if (length(b1) && length(e1)) lines <- lines[(b1[1] + 1):(e1[1] - 1)]
  keep <- !grepl("^[#!^]", lines)
  lines <- lines[keep]
  tab <- utils::read.delim(text = lines, header = TRUE, check.names = FALSE,
                           quote = "\"", stringsAsFactors = FALSE, na.strings = c("", "NA"))
  cand <- c("Gene symbol", "GENE_SYMBOL", "Gene_Symbol", "Symbol", "GENE",
            "ILMN_Gene", "gene_assignment", "Gene Symbol")
  sym_col <- intersect(cand, names(tab))[1]
  if (is.na(sym_col)) {
    # 表里没有基因符号列 = 这文件不是注释表（多半是别的平台文件误命名落在这里）
    log_msg(basename(local), " 没有基因符号列，丢弃该文件")
    close(con); on.exit(NULL, add = FALSE)
    unlink(local)
    a <- annotation_from_db(gpl)
    if (!is.null(a)) return(a)
    stop("注释文件格式不可识别且无替代来源：", local, call. = FALSE)
  }
  sym <- as.character(tab[[sym_col]])
  if (sym_col == "gene_assignment") {           # Affymetrix: "NM_x // SYMBOL // ..."
    sym <- vapply(strsplit(sym, " // ", fixed = TRUE),
                  function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  }
  out <- data.frame(probe = as.character(tab[[1]]), symbol = trimws(sym),
                    stringsAsFactors = FALSE)
  out <- out[!is.na(out$symbol) & out$symbol != "" & out$symbol != "---", ]
  out[!duplicated(out$probe), ]
}

# 探针 -> 基因：同基因多探针取平均表达最高者（避免 merge 造成的倍数稀释）
collapse_probes <- function(expr, annot) {
  idx <- match(rownames(expr), annot$probe)
  ok <- !is.na(idx)
  expr <- expr[ok, , drop = FALSE]
  sym <- annot$symbol[idx[ok]]
  # 一个探针映射多个基因时只保留第一个
  sym <- sub(" ?///.*$", "", sym)
  sym <- sub(" ?//.*$", "", sym)
  mu <- rowMeans(expr, na.rm = TRUE)
  ord <- order(sym, -mu)
  expr <- expr[ord, , drop = FALSE]
  sym <- sym[ord]
  keep <- !duplicated(sym)
  expr <- expr[keep, , drop = FALSE]
  rownames(expr) <- sym[keep]
  log_msg("探针合并：", sum(ok), " 探针 -> ", nrow(expr), " 基因")
  expr
}

# 是否已经在 log2 尺度上：log2 后的芯片数据 99% 分位通常在 4–16，
# 原始强度则在数百到数万。原来这里写的是 `q > thr`，函数名说"是 log"、
# 判断却返回"是原始"，调用处再取反，结果两组数据全判反：
# 已是 log2 的 GSE27276 / GSE45570 被二次 log2（动态范围压成 1.2），
# 而原始强度的 GSE4316 完全没转换（limma 报出 logFC 6105）。
is_log_scale <- function(expr, thr = 50) {
  q <- stats::quantile(expr, 0.99, na.rm = TRUE)
  is.finite(q) && q <= thr
}

zscore_rows <- function(mat) t(scale(t(mat)))

# 安全除法
safe_div <- function(a, b) ifelse(b == 0 | is.na(b), NA_real_, a / b)

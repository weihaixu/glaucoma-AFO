#!/usr/bin/env Rscript
# 修 Seurat 的依赖链：一次性解析全部硬依赖并装齐缺失项，直到 library(Seurat) 成功
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"),
        install.packages.check.lock = TRUE)

strip_v <- function(x) trimws(sub("\\s*\\(.*$", "", x))

get_deps <- function(pkg) {
  d <- suppressWarnings(packageDescription(pkg))
  if (is.null(d) || is.atomic(d) || !length(d$Package) || is.na(d$Package)) return(character(0))
  paste(c(d$Depends, d$Imports, d$LinkingTo), collapse = ", ") |>
    (\(s) gsub("\\s", "", s))() |>
    (\(s) strsplit(s, ",")[[1]])() |>
    (\(v) trimws(v[nzchar(v)]))() |>
    (\(v) strip_v(v))() |>
    (\(v) v[v != "R"])()
}

rounds <- 0
repeat {
  ok <- tryCatch({
    suppressPackageStartupMessages(library(Seurat))
    TRUE
  }, error = function(e) FALSE)
  if (ok) { cat("SEURAT_OK after", rounds, "rounds\n"); break }
  all_deps <- unique(unlist(lapply(c("Seurat", "SeuratObject"), get_deps)))
  # 传递依赖再扫一层
  all_deps <- unique(c(all_deps, unlist(lapply(all_deps, get_deps))))
  missing <- all_deps[!vapply(all_deps, function(p) p %in% rownames(installed.packages()), logical(1))]
  missing <- setdiff(missing, rownames(installed.packages()))
  cat("round", rounds + 1, "missing:", paste(missing, collapse = ", "), "\n")
  if (!length(missing)) {
    # 都在但加载仍失败：逐个 requireNamespace 找加载不了的
    bad <- all_deps[!vapply(all_deps, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
    if (!length(bad)) { cat("STUCK: no missing deps but Seurat still fails\n"); break }
    cat("load-broken:", paste(bad, collapse = ", "), "\n")
    # 对加载坏的重装
    install.packages(unique(bad), type = "win.binary")
  } else {
    install.packages(missing, type = "win.binary")
  }
  rounds <- rounds + 1
  if (rounds >= 8) { cat("TOO_MANY_ROUNDS\n"); break }
}
r <- tryCatch({ suppressPackageStartupMessages(library(Seurat)); "SEURAT_OK" },
              error = function(e) paste("SEURAT_FAIL:", conditionMessage(e)))
cat(r, "\n")

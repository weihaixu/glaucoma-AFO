#!/usr/bin/env Rscript
# 安装 06/11/12/14 步缺的两个包：CellChat 与 monocle3
# 顺序：CRAN 硬依赖 -> CellChat -> monocle3；日志同步打到 stdout 和 logs/install_sc_pkgs.log
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"),
        Ncpus = max(1L, parallel::detectCores() %/% 2L))

log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# --- 1) CRAN 依赖（NMF 是 CellChat 的硬依赖，pryr 是 monocle3 的）---
log("step1: CRAN deps (NMF, pryr)")
install.packages(c("NMF", "pryr"), type = "win.binary")
for (p in c("NMF", "pryr")) {
  log("  ", p, ":", p %in% rownames(installed.packages()))
}

# --- 2) CellChat：先试 CRAN（若已回归），失败走 GitHub ---
if (!requireNamespace("CellChat", quietly = TRUE)) {
  log("step2: CellChat via CRAN")
  try(install.packages("CellChat", type = "win.binary"), silent = TRUE)
}
if (!requireNamespace("CellChat", quietly = TRUE)) {
  log("step2b: CellChat via GitHub (jinworks/CellChat)")
  try(remotes::install_github("jinworks/CellChat", upgrade = "never"), silent = TRUE)
}
log("CellChat:", requireNamespace("CellChat", quietly = TRUE))

# --- 3) monocle3：只在 GitHub ---
log("step3: monocle3 via GitHub (cole-trapnell-lab/monocle3)")
try(remotes::install_github("cole-trapnell-lab/monocle3", upgrade = "never"), silent = TRUE)
log("monocle3:", requireNamespace("monocle3", quietly = TRUE))

# --- 4) 加载终检 ---
for (p in c("CellChat", "monocle3")) {
  r <- tryCatch({ suppressPackageStartupMessages(library(p, character.only = TRUE)); "LOAD_OK" },
                error = function(e) paste("LOAD_FAIL:", conditionMessage(e)))
  log(p, "->", r)
}
log("DONE")

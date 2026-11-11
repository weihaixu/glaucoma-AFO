# =============================================================================
# 17_fill_tokens.R —— 把账本里的数字回填进稿件，并检查有没有漏/多
#
# 前提：15_results_ledger.R 已生成 results/tables/results_ledger.tsv
# 稿件里所有会变的数字都写成 {{TOKEN}}（例如 {{AUC_TRAIN}}），本脚本一次性替换。
#
# 三种状态，都要处理，不能只看"替换成功"：
#   1) 稿件有、账本有        → 替换，记 OK
#   2) 稿件有、账本没有      → 记 MISSING_TO_TOKEN（通常是拼错 token 名，或该步没跑）
#   3) 账本有、稿件没用到    → 记 UNUSED（可能是漏写了一段，也可能是确实不需要）
#
# 用法：
#   Rscript 17_fill_tokens.R
#   Rscript 17_fill_tokens.R --in ../output/manuscript/formatted-manuscript.html \
#                            --out ../output/manuscript/manuscript_filled.html
#   Rscript 17_fill_tokens.R --docx        # 额外调用 html-to-doc 转 docx（需配 converter 路径）
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("17_fill_tokens", cfg)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) && i < length(args)) args[i + 1] else default
}
in_html  <- get_arg("--in",  proj_path("..", "output", "manuscript", "formatted-manuscript.html"))
out_html <- get_arg("--out", sub("\\.html$", "_filled.html", in_html))
make_docx <- "--docx" %in% args

require_file(in_html, "稿件模板不存在；先用 --in 指定路径")

# ------------------------------------------------------------------ 账本 ----
ledger_f <- res_path(cfg, "tables", "results_ledger.tsv")
if (!file.exists(ledger_f)) {
  stop("找不到 results_ledger.tsv。请先运行：Rscript 15_results_ledger.R", call. = FALSE)
}
# 用 encoding="UTF-8"（只标记、不转码），不要用 fileEncoding：
# Windows 中文 locale 下 fileEncoding 会尝试转成本地 GBK，遇到 "λ" "–" 这类字符直接报
# invalid input；而若不指定，字符串的编码标记是 unknown，和 UTF-8 的 html 混在一起做
# gsub 会按字节错位，替换值里的中文顿号会把紧随其后的 HTML 一起吞掉（"(TIMP1c)" 那种残句）
led <- utils::read.delim(ledger_f, stringsAsFactors = FALSE, check.names = FALSE,
                         encoding = "UTF-8")
need <- c("token", "display")
if (!all(need %in% names(led))) stop("账本缺少 token / display 列，请检查 15 的输出", call. = FALSE)
led <- led[!is.na(led$token) & nzchar(led$token), ]
led <- led[!duplicated(led$token), ]          # 同名 token 只保留第一次（先写先得）
led$display <- enc2utf8(as.character(led$display))

# ------------------------------------------------------------- 执行替换 ----
html <- paste(readLines(in_html, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
html <- enc2utf8(html)

# 稿件里的 {{...}} 有两种：自动 token（{{AUC_TRAIN}}）和作者手填占位（{{[bootstrap AUC]}}）。
# 后者不匹配 [A-Z0-9_] 形态，如果只用严格正则会被静默漏掉，最后以为"全填完了"
ph_all <- function(x) {
  v <- unique(regmatches(x, gregexpr("\\{\\{[^{}]+\\}\\}", x))[[1]])
  gsub("^\\{\\{|\\}\\}$", "", v)
}
ph <- ph_all(html)
is_tok <- grepl("^[A-Z0-9_]+$", ph)
manual_ph <- ph[!is_tok]
found <- ph[is_tok]
log_msg("稿件中共 ", length(found), " 个自动 token，另有 ", length(manual_ph), " 处人工占位")

rep_n <- 0L
for (i in seq_len(nrow(led))) {
  tk <- led$token[i]
  val <- as.character(led$display[i])
  if (is.na(val) || !nzchar(val)) next
  pat <- paste0("{{", tk, "}}")
  n <- grepl(pat, html, fixed = TRUE)
  if (!n) next
  cnt <- lengths(regmatches(html, gregexpr(pat, html, fixed = TRUE)))
  html <- gsub(pat, val, html, fixed = TRUE)
  rep_n <- rep_n + cnt
}
log_msg("替换 ", rep_n, " 处")

# 检查残留
left <- ph_all(html)
left <- left[grepl("^[A-Z0-9_]+$", left)]
status <- data.frame(
  token = c(found[found %in% led$token], left, setdiff(led$token, found), manual_ph),
  status = c(rep("OK", sum(found %in% led$token)),
             rep("MISSING_IN_LEDGER", length(left)),
             rep("UNUSED_IN_TEXT", length(setdiff(led$token, found))),
             rep("MANUAL_PLACEHOLDER", length(manual_ph))),
  stringsAsFactors = FALSE)
status <- status[!duplicated(paste(status$token, status$status)), ]

# 给缺失 token 一个可执行的排查方向：
# 「该步没跑」和「15 没实现这个 token」是两件完全不同的事，混在一起看会白花时间
ms_f <- res_path(cfg, "tables", "missing_steps.tsv")
# 必须按字符读：read.delim 的 type.convert 会把 step 列的 "06" 变成整数 6，
# 后面按 "06" 匹配就全部落空，于是所有缺失都被误判成「已产出但账本没这个 token」
# 注意不要加 encoding = "UTF-8"：那会给 note 打上 UTF-8 标记，
# 再和源码里的中文字面量 paste0 时，R 会按本地编码做一次多余转换，
# 结果写出 "<U+9700>" 这种转义（Excel 里全是问号）。保持 unknown 让字节直通。
ms <- if (file.exists(ms_f)) {
  utils::read.delim(ms_f, stringsAsFactors = FALSE, colClasses = "character")
} else NULL
step_of <- function(tk) {
  rules <- c(
    "^CONCORD_" = "02",
    "^N_GENE_TESTED$|^N_DEG|^N_CONSENSUS$|^KEY_GENES|^N_KEY_GENES" = "02",
    "^N_HUB$" = "03",
    "^N_CELLS$|^N_CELLTYPES$|^CELLTYPE_LIST$|^TM_CONCORD|^N_LR|^LR_" = "06",
    "^N_SCENIC|^SCENIC_" = "06",
    "^TWAS_|^N_TWAS|^N_FOCUS_PIP$|^N_DOUBLE_SUPPORT$|^DOUBLE_GENES$" = "09",
    "^PERT_|^N_PERT|^N_CMAP" = "11",
    "^MED_" = "12",
    "^XS_|^N_XS_" = "14")
  for (p in names(rules)) if (grepl(p, tk)) return(unname(rules[p]))
  ""
}
hint_of <- function(tk) {
  # 两类「步骤跑过了但数字仍空」的固定原因，先说清楚，省得去翻代码
  if (grepl("^CONCORD_", tk)) {
    return("验证集名称与本次运行不一致（demo 下是 {{CONCORD_DEMO_V}}）；换真实数据集后会自然对上")
  }
  if (identical(tk, "N_HUB")) {
    return("03 已跑，但没有模块达到 |r|>=0.4 且 P<0.05，hub 定义为空；真实数据通常非空")
  }
  if (grepl("^N_SCENIC|^SCENIC_", tk)) {
    return("SCENIC 不在 06 内：先跑 tools/scenic_regulon_count.R 生成 scenic_regulon_diff.tsv，再重跑 15")
  }
  if (grepl("^N_LR|^LR_", tk)) {
    return("需要 06 产出的 cellchat_<条件>.tsv 两张表；单细胞未跑时无法给出")
  }
  st <- step_of(tk)
  if (!nzchar(st)) return("账本无此 token：确认 15 是否已实现，再查稿件里的拼写")
  if (!is.null(ms) && st %in% ms$step) {
    i <- match(st, ms$step)
    if (!isTRUE(as.logical(ms$has_output[i]))) {
      note <- if ("note" %in% names(ms)) ms$note[i] else ""
      return(paste0("步骤 ", st, " 未产出", if (nzchar(note)) paste0("：", note) else ""))
    }
  }
  paste0("步骤 ", st, " 已产出但账本无此 token：查 15 的实现，或稿件里 token 名写错")
}
status$hint <- ""
mi <- status$status == "MISSING_IN_LEDGER"
if (any(mi)) status$hint[mi] <- vapply(status$token[mi], hint_of, character(1))
mp <- status$status == "MANUAL_PLACEHOLDER"
if (any(mp)) status$hint[mp] <- "作者手填占位（非自动 token）：需人工核对后替换"

if (length(left)) {
  log_msg("!! 以下 token 在稿件里但没有对应数值（正文这几处还是空的）：")
  for (t in left) log_msg("   - ", t)
}
unused <- setdiff(led$token, found)
if (length(unused)) log_msg("提示：", length(unused), " 个 token 未被正文使用（不一定错，核对即可）")

# 用二进制写，避免 Windows 下按本地编码回写导致中文再次损坏
writeBin(charToRaw(enc2utf8(html)), out_html)
log_msg("已写出：", out_html)
# 报告里有中文 hint（"步骤 14 未产出：…"）。这台机器 R 的本地编码是 C，
# 一旦给 fileEncoding 就会多转一层，写出 <U+9700> 这种转义；
# 不加反而字节直通（15 写的 missing_steps.tsv 就是这么来的），所以这里刻意不写 fileEncoding。
rep_f <- res_path(cfg, "submission", "token_fill_report.tsv")
# 首次在干净目录里跑到这一步时 results/submission/ 还不存在，write.table 会直接
# 报 "cannot open the connection" 并 halt —— 而此时 HTML 已经写出去了，
# 看起来像"回填成功但报告没生成"，实际是整个脚本失败。先建目录。
dir.create(dirname(rep_f), recursive = TRUE, showWarnings = FALSE)
utils::write.table(status, rep_f, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("写出表格：", rep_f, " (", nrow(status), " 行)")

# ------------------------------------------------------------- 正文字数 ----
# 字数必须在 token 全部替换之后才算：一个 {{AUC_TRAIN}} 算一个词，换成
# "0.87 (0.81–0.92)" 就是三个词。手填的数字跑一次真实数据就作废，所以这里
# 交给 tools/fill_wordcount.py 重算，并同一条命令里同步 HTML 与 docx。
# 口径：1. Introduction 到 6. Declarations 之前（不含摘要、声明、参考文献、表格、图注）。
py <- Sys.which("python")
wc_script <- proj_path("tools", "fill_wordcount.py")
if (nzchar(py) && file.exists(wc_script)) {
  cmd <- sprintf("\"%s\" \"%s\" --html \"%s\"", py, wc_script, out_html)
  docx_same <- sub("\\.html$", ".docx", out_html)
  if (file.exists(docx_same)) cmd <- paste0(cmd, sprintf(" --docx \"%s\"", docx_same))
  log_msg("统计正文字数：", cmd)
  st_wc <- system(cmd)
  if (st_wc != 0) log_msg("字数脚本返回非 0（HTML 与 docx 字数不一致或含旧值），请查看上方输出")
} else {
  log_msg("跳过字数回填：缺 python 或 tools/fill_wordcount.py")
}

# --------------------------------------------------------------- 转 docx ----
if (make_docx) {
  conv <- cfg_get(cfg, "manuscript.converter", "")
  if (!nzchar(conv) || !file.exists(conv)) {
    log_msg("未配置 manuscript.converter，跳过 docx 转换。")
    log_msg("可在 config.yml 里加：manuscript: { converter: \"<插件>/skills/html-to-doc/scripts/convert_doc.py\" }")
  } else {
    out_docx <- sub("\\.html$", ".docx", out_html)
    cmd <- sprintf("python \"%s\" --html \"%s\" --output \"%s\"", conv, out_html, out_docx)
    log_msg("执行：", cmd)
    st <- system(cmd)
    if (st == 0) log_msg("docx 已生成：", out_docx) else log_msg("转换失败，见上方输出")
  }
}

cat("\n---- 汇总 ----\n")
print(table(status$status))
cat("替换后仍残留占位符：", length(left), " 个\n")
if (length(left)) cat(paste0("  ", left, collapse = "\n"), "\n")

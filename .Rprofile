# =============================================================================
# .Rprofile —— 在解析任何脚本之前把 LC_CTYPE 切到 UTF-8
#
# 为什么必须放在这里：中文 Windows 的默认 locale 是 C，R 会把源码里的
# 多字节字符（— ± μ ×）当成多个单字节字符，字符串的 Encoding() 是
# "unknown"。cairo 设备救不了这种情况——它拿到的就是一个没被标记成 UTF-8
# 的字节串，画到图里变成 "GSE27276 b" 这种残骸。
# locale 必须在【解析脚本之前】设好，所以只能放 .Rprofile：它比脚本加载更早，
# 而 Rscript 是边解析边执行，改完即时生效。
# =============================================================================
invisible(tryCatch({
  cur <- Sys.getlocale("LC_CTYPE")
  if (!grepl("UTF-8|utf8", cur, ignore.case = TRUE)) {
    for (loc in c("en_US.UTF-8", "C.UTF-8", "zh_CN.UTF-8",
                  "English_United States.65001")) {
      r <- suppressWarnings(Sys.setlocale("LC_CTYPE", loc))
      if (grepl("UTF-8|utf8|65001", r, ignore.case = TRUE)) break
    }
  }
}, error = function(e) NULL))

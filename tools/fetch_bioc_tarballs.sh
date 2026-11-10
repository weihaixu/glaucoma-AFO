#!/bin/bash
# fetch_bioc_tarballs.sh —— 断点续传拉 Bioc tarball，装完本地安装
# 用法：bash tools/fetch_bioc_tarballs.sh   （建议后台跑，国际带宽慢，可能要几小时）
cd "$(dirname "$0")/../data/raw/_tarballs" || exit 1
B="https://bioconductor.org/packages/3.22"
declare -A FILES=(
  [clusterProfiler_4.18.4.tar.gz]="$B/bioc/src/contrib"
  [ReactomePA_1.54.0.tar.gz]="$B/bioc/src/contrib"
  [org.Hs.eg.db_3.22.0.tar.gz]="$B/data/annotation/src/contrib"
  [GO.db_3.22.0.tar.gz]="$B/data/annotation/src/contrib"
  [reactome.db_1.95.0.tar.gz]="$B/data/annotation/src/contrib"
  [org.Mm.eg.db_3.22.0.tar.gz]="$B/data/annotation/src/contrib"
)
for round in 1 2 3 4 5 6 7 8; do
  all_ok=1
  for f in "${!FILES[@]}"; do
    # 校验：tar -tzf 能列目录才算完整
    if [ -f "$f" ] && tar -tzf "$f" >/dev/null 2>&1; then continue; fi
    rm -f "$f"
    echo "[round $round] 下载 $f"
    curl -sL -C - -m 7200 --retry 3 -o "$f" "${FILES[$f]}/$f"
    if tar -tzf "$f" >/dev/null 2>&1; then echo "  OK $(stat -c%s "$f")B"; else all_ok=0; rm -f "$f"; fi
  done
  [ "$all_ok" = 1 ] && break
  sleep 120
done
echo "=== 下载阶段结束，开始本地安装 ==="
"/g/R/R-4.5.2/bin/Rscript.exe" -e 'for (f in list.files("data/raw/_tarballs", pattern="tar.gz$", full.names=TRUE)) { cat("install:", f, "\n"); try(install.packages(f, repos=NULL, type="source")) }; for (p in c("clusterProfiler","ReactomePA","org.Hs.eg.db","GO.db","reactome.db","org.Mm.eg.db")) cat(p, requireNamespace(p, quietly=TRUE), "\n")' > ../../logs/bioc_local_install.log 2>&1
tail -8 ../../logs/bioc_local_install.log

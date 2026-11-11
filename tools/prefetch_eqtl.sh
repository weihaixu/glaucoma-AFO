#!/usr/bin/env bash
# prefetch_eqtl.sh —— 并发预取 GTEx 区间（08_coloc_smr.R 会直接读缓存）
#
# 为什么单独写：coloc 要对 位点 × 组织 逐个取数，串行跑一次几十秒，
# 120 组就是两小时。这里用 xargs -P 起 6 个并发，十几分钟拿完，
# 之后跑 08 全程只读本地缓存。
#
# 用法：bash tools/prefetch_eqtl.sh [n_loci] [parallel]
set -u
cd "$(dirname "$0")/.." || exit 1

N_LOCI=${1:-20}
PAR=${2:-6}
BASE="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge"
CACHE="data/raw/eqtl_cache"
mkdir -p "$CACHE"

TISSUES=("Artery_Tibial" "Nerve_Tibial" "Muscle_Skeletal" "Adipose_Subcutaneous" \
         "Skin_Sun_Exposed_Lower_leg" "Brain_Frontal_Cortex_BA9")

TASKS=$(mktemp)
tail -n +2 data/processed/coloc_loci.tsv | head -n "$N_LOCI" | while read -r chr start end lead_p lead_snp; do
  for t in "${TISSUES[@]}"; do
    printf '%s %s %s %s\n' "$t" "$chr" "$start" "$end"
  done
done > "$TASKS"

echo "任务数：$(wc -l < "$TASKS")，并发：$PAR"

fetch_one() {
  t=$1; c=$2; s=$3; e=$4
  out="$CACHE/${t}_${c}_${s}_${e}.tsv"
  if [ -s "$out" ]; then return 0; fi
  for attempt in 1 2 3; do
    python tools/fetch_eqtl_tabix.py --url "$BASE/$t.tsv.gz" --index "$BASE/$t.tsv.gz.tbi" \
      --chrom "$c" --start "$s" --end "$e" --out "$out" 2>/dev/null && return 0
    sleep 3
  done
  echo "FAILED $t $c $s $e" >&2
  return 1
}
export -f fetch_one
export BASE CACHE

xargs -P "$PAR" -L 1 bash -c 'fetch_one "$@"' _ < "$TASKS"
echo "缓存文件数：$(ls -1 "$CACHE" | wc -l)"
du -sh "$CACHE"
rm -f "$TASKS"

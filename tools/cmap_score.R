# tools/cmap_score.R -- CMap 阶段 2：逐桶聚合 + 带符号 KS 打分
# 读取 cmap_work/part_*.tsv（dr/gi 为全局 int 编码），聚合后逐 drug 打分。
suppressMessages(library(data.table))
setDTthreads(4)

src      <- "data/raw/lincs/Cmap_QCeddata"
work     <- "data/raw/lincs/cmap_work"
sig_file <- "results/tables/deg_GSE27276.tsv"
out_file <- "results/tables/cmap_connectivity.tsv"
NBUCKET  <- 12L

panel <- c(
  "lung_A549_ph1_QCed.txt.gz", "lung_A549_ph2_QCed.txt.gz",
  "skin_A375_ph1_QCed.txt.gz", "skin_A375_ph2_QCed.txt.gz",
  "liver_HEPG2_ph1_QCed.txt.gz", "liver_HEPG2_ph2_QCed.txt.gz",
  "lung_HCC515_ph1_QCed.txt.gz", "lung_HCC515_ph2_QCed.txt.gz",
  "digestive_HT29_ph1_QCed.txt.gz", "digestive_HT29_ph2_QCed.txt.gz",
  "breast_MCF7_ph1_QCed.txt.gz", "breast_MCF7_ph2_QCed.txt.gz",
  "prostate_PC3_ph1_QCed.txt.gz", "prostate_PC3_ph2_QCed.txt.gz",
  "prostate_VCAP_ph1_QCed.txt.gz", "skin_SKMEL28_ph1_QCed.txt.gz"
)
bfiles <- file.path(work, sprintf("part_%02d.tsv", seq_len(NBUCKET) - 1L))

# ---- 全局字典（与 fill 阶段一致）----
drug_levels <- character(0)
for (g in panel) {
  h <- fread(file.path(src, g), nrows = 0, encoding = "Latin-1", showProgress = FALSE)
  cs <- names(h)[names(h) != "Gene"]
  pf <- iconv(cs, from = "Latin1", to = "UTF-8")
  drug_levels <- union(drug_levels, sub("^[^_]+_([^_]+)_.*$", "\\1", pf))
}
G <- character(0)
for (g in panel) {
  d <- fread(file.path(src, g), select = "Gene", encoding = "Latin-1", showProgress = FALSE)
  G <- union(G, d$Gene)
  rm(d); invisible(gc())
}
cat("global drugs:", length(drug_levels), " genes:", length(G), "\n")

# ---- 疾病签名 ----
cons <- fread(sig_file)
padj_col <- if ("adj.P.Val" %in% names(cons)) "adj.P.Val" else grep("padj|FDR", names(cons), value = TRUE)[1]
lfc_col <- grep("^logFC", names(cons), value = TRUE)[1]
sig <- cons[is.finite(get(padj_col)) & get(padj_col) < 0.05]
cons_lfc <- setNames(as.numeric(sig[[lfc_col]]), toupper(sig$gene))
UP_T <- toupper(names(cons_lfc)[cons_lfc > 0])
DN_T <- toupper(names(cons_lfc)[cons_lfc < 0])
cat("signature:", length(UP_T), "up /", length(DN_T), "down\n")

ks_cs <- function(gene, zscore) {
  o <- order(-zscore)
  nr <- length(o)
  g <- toupper(gene[o])
  es <- function(hits) {
    if (length(hits) < 5) return(NA_real_)
    r <- seq_along(hits)
    dev <- r / length(hits) - hits / nr
    dev[which.max(abs(dev))]
  }
  e_up <- es(which(g %in% UP_T))
  e_dn <- es(which(g %in% DN_T))
  if (is.na(e_up) || is.na(e_dn)) return(NA_real_)
  e_up - e_dn
}

# ---- 逐桶聚合 + 打分 ----
t0 <- Sys.time()
res_list <- list()
for (bb in seq_len(NBUCKET)) {
  bf <- bfiles[bb]
  if (!file.exists(bf) || file.size(bf) < 10) next
  t <- fread(bf, showProgress = FALSE)
  t <- t[, .(zsum = sum(zsum), n = sum(n)), by = .(dr, gi)]
  sc <- t[, list(cs = ks_cs(G[gi], zsum / n)), by = dr]
  res_list[[length(res_list) + 1L]] <- sc
  cat(sprintf("bucket %2d/%d: %d drugs, union=%s rows, elapsed=%s\n",
              bb, NBUCKET, nrow(sc), format(nrow(t), big.mark = ","),
              format(Sys.time() - t0)), sep = "")
  rm(t, sc); invisible(gc(full = TRUE))
}

allres <- rbindlist(res_list)
allres[, drug := drug_levels[dr]]
allres <- allres[is.finite(cs)][order(cs)]
fwrite(allres[, .(drug, cs)], out_file, sep = "\t", quote = FALSE)
cat("compounds scored:", nrow(allres), "\n")
cat("cs<0 (inverse enrichment):", sum(allres$cs < 0), "\n")
cat("most negative 5:", paste(head(allres$drug, 5), collapse = ", "), "\n")
cat("written:", out_file, "\n")

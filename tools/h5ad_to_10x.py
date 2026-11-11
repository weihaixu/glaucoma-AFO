#!/usr/bin/env python
"""GSE263230 h5ad -> 每样本 10x 风格目录（matrix.mtx/barcodes/genes）+ obs.csv

06_scrna_cellchat.R 期望 Seurat rds 或 CellRanger 输出；anndata rds 直读要
reticulate + SeuratDisk，在 Windows 上都不稳。这里导出成最朴素的格式，
再在 R 里 Read10X() 组装，注释列从 obs.csv 拿回来。
"""
import os
import sys
import argparse

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", required=True, help="样本名（如 GTM_1672_P3）")
    a = ap.parse_args()

    A = ad.read_h5ad(a.h5ad)
    X = A.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = sp.csr_matrix(X)  # cells x genes

    os.makedirs(a.outdir, exist_ok=True)
    sio.mmwrite(os.path.join(a.outdir, "matrix.mtx"), X.T)  # 10x 是 genes x cells

    # CELLxGENE 的 var_names 是 ENSG；CellChat/11 的目标基因全是 symbol，
    # 必须换 feature_name（同名 symbol 用 make.unique 去重）
    if "feature_name" in A.var.columns:
        gene_names = [str(x) for x in A.var["feature_name"]]
    else:
        gene_names = list(np.asarray(A.var_names).astype(str))
    # 去重（"/" 等非法字符 read.delim 不挑，直接 make.unique 风格加后缀）
    seen = {}
    uniq = []
    for g in gene_names:
        if g in seen:
            seen[g] += 1
            uniq.append("%s.%d" % (g, seen[g]))
        else:
            seen[g] = 0
            uniq.append(g if g else "blank")
    genes = pd.DataFrame({
        "gene": uniq,
        "foo": np.arange(X.shape[1]) + 1,
    })
    genes.to_csv(os.path.join(a.outdir, "genes.tsv"), sep="\t",
                 index=False, header=False)

    with open(os.path.join(a.outdir, "barcodes.tsv"), "w") as f:
        for b in A.obs_names:
            f.write(str(b) + "\n")

    obs = A.obs.copy()
    obs.insert(0, "barcode", [str(x) for x in A.obs_names])
    obs.insert(1, "sample", a.sample)
    obs.to_csv(os.path.join(a.outdir, "obs.csv"), index=False)

    print("sample=%s cells=%d genes=%d obs_cols=%s" %
          (a.sample, X.shape[0], X.shape[1], ",".join(map(str, A.obs.columns))))


if __name__ == "__main__":
    main()

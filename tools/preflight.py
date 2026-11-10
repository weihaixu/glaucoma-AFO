#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tools/preflight.py —— 不依赖 R 的「开跑前体检」

为什么需要这个脚本：
    15_results_ledger.R 的第 8 节会写 missing_steps.tsv，但它的前提是
    results/tables/ 下已经有东西；而且跑它需要 R。实际最常见的情形是
    「R 还没装 / 数据还没下 / 什么都没跑」，这时 15 会直接 stop()：
        stop("账本为空：没有任何结果表可读。先跑 01–14 ……")
    也就是说，缺步清单在最需要它的时候恰恰拿不到。

本脚本做三件不依赖 R 的事：
    1. 逐步检查 01–14 的关键结果表是否存在（与 15 的 step_files 完全同一张映射表）
    2. 检查 config.yml 里声明的外部数据是否在磁盘上（这些是要人工下载的大件）
    3. 输出「下一步该做什么」——缺数据 vs 缺运行，处置方式完全不同

输出（默认 results/preflight/）：
    missing_steps.tsv      step / key_table / has_output / run_status / consequence
    data_readiness.tsv     config_key / expected_path / exists / used_by / how_to_get
    summary.txt            给人看的一段结论

用法：
    python tools/preflight.py
    python tools/preflight.py --root . --outdir results/preflight

注意：本脚本只读，不写 results/tables/，不会和 15 的输出混淆。
      真正投稿用的 missing_steps.tsv 仍以 15 生成的为准。
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from datetime import datetime

# ---- 与 15_results_ledger.R 第 8 节保持同一张映射表（改一处要同步另一处）----
STEP_TABLES = [
    ("01", "qc_summary", "数据准备与 QC"),
    ("02", "deg_consensus", "差异表达 + 共识基因"),
    ("03", "wgcna_module_trait", "WGCNA 模块-性状"),
    # 04 的主产物是富集表；反卷积只是可选分支（要 scRNA 参考或 xCell），
    # 拿 deconv_fraction 当判据会把跑成功的 04 误判成没跑
    ("04", "enrich_go_Consensus", "富集 / PROGENy（反卷积见 6a，可选）"),
    ("05", "model_performance", "机器学习建模"),
    ("06", "sc_celltype_prop", "单细胞注释 / CellChat"),
    ("07", "gwas_qc_summary", "GWAS 质控"),
    ("08", "coloc_summary", "共定位"),
    ("09", "twas_all", "TWAS / FOCUS"),
    ("10", "mr_results", "孟德尔随机化"),
    ("11", "perturb_summary", "计算扰动（数字敲除）"),
    ("12", "mr_mediation", "两阶段 / 中介 MR"),
    ("13", "robustness_permutation", "置换零分布 + multiverse + CMap"),
    ("14", "xspecies_concordance", "跨物种保守"),
]

# ---- config.yml 里声明的外部数据（路径 -> 由谁使用 / 怎么拿）----
DATA_INPUTS = [
    ("genetics.gwas.path", "data/raw/finngen_R12_H7_GLAUCOMA.gz", "07/08/09/10",
     "FinnGen R12 官网按 endpoint 下载（H7_GLAUCOMA）"),
    ("genetics.ld.panel_plink", "data/raw/1000G_EUR_Phase3_plink/1000G.EUR.QC.bed", "07/08",
     "1000G EUR Phase3 PLINK 面板（LDSC 仓库提供）"),
    ("genetics.eqtl.gtex_dir", "data/raw/GTEx_v8/Artery_Tibial.v8.EUR.allpairs.txt.gz", "08/09/10",
     "GTEx v8 eQTL allpairs（按组织逐个下）；眼组织需 EyeGEx"),
    ("genetics.pqtl.path", "data/raw/UKBPPP_pqtl_selected.tsv.gz", "08/10",
     "UKB-PPP pQTL（Synapse 申请）"),
    ("genetics.mr.exposure", "data/raw/eqtl_instruments.tsv.gz", "10",
     "由 08 的 eQTL 结果整理，或直接用 eQTL Catalogue"),
    ("genetics.twas.weights_dir", "data/raw/FUSION_weights", "09",
     "FUSION 预计算权重（GTEx v8 / 眼组织）"),
    ("scrna.rds", "data/processed/scRNA_ref.rds", "06/11/14",
     "人 TM/JCT 单细胞（GSE313736）读入后存 Seurat 对象；也可改 mode=10x"),
    ("drylab.mediation.iop_gwas", "data/raw/iop_gwas.tsv.gz", "12",
     "IOP GWAS（UKB / FinnGen / 已发表 meta）"),
    ("drylab.robustness.lincs_path", "", "13", "LINCS/CMap 签名（可选，不填则跳过 CMap）"),
]

RUN_STEPS_DEFAULT = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10"]
DRYLAB_STEPS = ["11", "12", "13", "14"]


def read_cfg_text(root: str, rel: str = "config/config.yml") -> str:
    p = rel if os.path.isabs(rel) else os.path.join(root, rel)
    if not os.path.exists(p):
        return ""
    with open(p, encoding="utf-8") as f:
        return f.read()


def cfg_get(text: str, key: str, default: str = ""):
    """极简取值：只处理 `key: "value"` / `key: value` 这种单行标量。"""
    m = re.search(r"^\s*" + re.escape(key.split(".")[-1]) + r":\s*(\"[^\"]*\"|'[^']*'|[^#\n]+)",
                  text, re.M)
    if not m:
        return default
    v = m.group(1).strip().strip("\"'").strip()
    return v or default


def has_r() -> tuple[bool, str]:
    """找 Rscript：先 PATH，再按常见安装位置扫目录。

    R 装在非 PATH 位置（如 G:\\R\\R-4.5.2）是常态，只查 PATH 会把"已装好"误报成
    "未安装——硬阻塞"，用户第一次跑就会被带偏。
    """
    exe = "Rscript.exe" if os.name == "nt" else "Rscript"
    for d in os.environ.get("PATH", "").split(os.pathsep):
        cand = os.path.join(d, exe)
        if os.path.exists(cand):
            return True, cand

    # R_HOME 直接指到安装目录的情况
    rh = os.environ.get("R_HOME")
    if rh:
        for sub in ("bin", os.path.join("bin", "x64")):
            cand = os.path.join(rh, sub, exe)
            if os.path.exists(cand):
                return True, cand

    if os.name == "nt":
        bases = [
            r"C:\Program Files\R",
            r"C:\Program Files (x86)\R",
            os.path.expandvars(r"%LOCALAPPDATA%\Programs\R"),
            os.path.expandvars(r"%PROGRAMFILES%\R"),
        ]
        # 非系统盘的常见手装位置：G:\R\R-4.5.2、D:\R\R-4.x …
        for drive in "CDEFGH":
            bases.append(drive + r":\R")
        found = []
        for base in bases:
            if not os.path.isdir(base):
                continue
            for name in os.listdir(base):
                for sub in ("bin", os.path.join("bin", "x64")):
                    cand = os.path.join(base, name, sub, exe)
                    if os.path.exists(cand):
                        found.append(cand)
        if found:
            # 有多版本时取版本号最大的
            def ver(p):
                m = re.search(r"R-(\d+)\.(\d+)\.(\d+)", p)
                return tuple(int(x) for x in m.groups()) if m else (0, 0, 0)
            return True, sorted(found, key=ver)[-1]
    return False, ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".", help="项目根目录（含 config/ 的那一层）")
    ap.add_argument("--outdir", default=None, help="默认 <root>/results/preflight")
    ap.add_argument("--config", default="config/config.yml",
                    help="用哪份配置；跑 demo 时传 config/config.demo.yml，"
                         "否则会去查 results/ 而不是 results_demo/")
    a = ap.parse_args()
    root = os.path.abspath(a.root)

    cfg_text = read_cfg_text(root, a.config)
    # outdir 跟着配置里的 project.outdir 走，别硬编码 results
    proj_out = cfg_get(cfg_text, "outdir", "results") or "results"
    tbl_dir = os.path.join(root, proj_out, "tables")
    outdir = os.path.abspath(a.outdir or os.path.join(root, proj_out, "preflight"))
    os.makedirs(outdir, exist_ok=True)

    r_ok, r_where = has_r()

    # ---------------- 1) 缺步清单 ----------------
    rows = []
    for step, key, desc in STEP_TABLES:
        f = os.path.join(tbl_dir, key + ".tsv")
        exists = os.path.exists(f)
        if exists:
            status, conseq = "已产出", ""
        elif not r_ok:
            status = "阻塞：本机无 R"
            conseq = "安装 R 后执行 Rscript run_all.R --only %s" % step
        elif step in DRYLAB_STEPS:
            status = "阻塞：需单独触发"
            conseq = "不在 run.steps 内，需 Rscript run_all.R --only %s" % step
        else:
            status = "未运行"
            conseq = "Rscript run_all.R --only %s" % step
        rows.append((step, key, desc, "yes" if exists else "no", status, conseq))

    with open(os.path.join(outdir, "missing_steps.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("step\tkey_table\tstep_desc\thas_output\trun_status\tconsequence\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    # ---------------- 2) 数据就绪度 ----------------
    drows = []
    for key, rel, used, how in DATA_INPUTS:
        if not rel:
            drows.append((key, "(未配置)", "n/a", used, how))
            continue
        p = os.path.join(root, rel.replace("/", os.sep))
        # 目录型输入（GTEx / FUSION_weights）判目录是否存在
        if "." not in os.path.basename(rel) or os.path.basename(rel).startswith("1000G"):
            ok = os.path.exists(p)
        else:
            ok = os.path.exists(p)
        drows.append((key, rel, "yes" if ok else "no", used, how))

    with open(os.path.join(outdir, "data_readiness.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("config_key\texpected_path\texists\tused_by\thow_to_get\n")
        for r in drows:
            f.write("\t".join(r) + "\n")

    # ---------------- 3) 结论 ----------------
    n_missing = sum(1 for r in rows if r[3] == "no")
    n_data_missing = sum(1 for r in drows if r[2] == "no")
    lines = [
        "# 开跑前体检（preflight）",
        "",
        "生成时间：%s" % datetime.now().strftime("%Y-%m-%d %H:%M"),
        "项目根目录：%s" % root,
        "",
        "## 环境",
        "R / Rscript：%s" % ("已找到（%s）" % r_where if r_ok else "未安装 —— 这是当前唯一的硬阻塞"),
        "",
        "## 结果表",
        "已产出 %d / %d 步；缺 %d 步" % (len(rows) - n_missing, len(rows), n_missing),
        "%s/tables/ 目录：%s" % (proj_out, "存在" if os.path.isdir(tbl_dir) else "不存在或为空"),
        "",
        "## 外部数据",
        "%s 声明的 %d 项输入，缺失 %d 项" % (a.config, len(drows), n_data_missing),
        "",
        "## 判读",
    ]
    if n_missing == len(rows) and not r_ok:
        lines += [
            "当前不是「缺 11/12/14」的状态，而是**一步都没跑**：14 步全部没有产物，",
            "因为本机没有 R，results/ 也还没有任何表。",
            "顺序必须是：装 R → 下载 07–14 依赖的外部数据（见 data_readiness.tsv）",
            "→ 01–10 跑通 → 单细胞对象就位后跑 06/11/14 → 12 需要 IOP GWAS",
            "→ 最后 15 才会写出真正的 missing_steps.tsv 与账本。",
        ]
    elif n_missing:
        miss = [r[0] for r in rows if r[3] == "no"]
        lines.append("缺产物步骤：%s（处置见 missing_steps.tsv 的 consequence 列）" % ", ".join(miss))
    else:
        lines.append("01–14 全部有产物，可以跑 15 → 16 → 17 出稿。")

    with open(os.path.join(outdir, "summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print("\n".join(lines))
    print("\n已写出：")
    for fn in ("missing_steps.tsv", "data_readiness.tsv", "summary.txt"):
        print("  " + os.path.join(outdir, fn))
    return 0


if __name__ == "__main__":
    sys.exit(main())

# 青光眼房水外流阻力机制 —— 课题一 / 课题二 R 分析骨架

配套开题方案《小梁网机械微环境重塑驱动房水外流阻力增加的分子机制研究》。

- **课题一**（bulk 转录组 + 单细胞）：`01`–`06`
- **课题二**（GWAS × eQTL × pQTL → 共定位 / TWAS / MR）：`07`–`10`

脚本是**骨架**：流程、参数、判读标准、输出字段都是按可发表的标准写的，但**不含任何真实结果**。换数据、改配置即可跑，不需要动代码逻辑。

> **先读这四份文档，再动脚本**
> - `docs/操作手册-从环境到成图.md`：按天执行的操作手册（装环境 → demo → 拿数据 → 01–10 每步怎么判 → 成图写作 → 周计划 → 故障速查）
> - `docs/实验SOP-课题三.md`：细胞与动物的详细实验方案（分组、原代 HTM、牵张、水凝胶、qPCR/WB/IF、磁流模型、统计）
> - `docs/纯干实验版-方案与操作流程.md`：**不做湿实验**时的替代设计与 16 周流程（`11`–`14` 服务于这一版）
> - `docs/成稿手册-从结果到投稿.md`：从分析产物到投稿稿件（token 回填、图表规划、逐节写法、期刊匹配、STROBE-MR 清单、投稿前 30 项、投稿信模板、审稿应答预案）
>
> 环境装完先跑 `Rscript 00_check_env.R`，通过后再跑 demo。

---

## 一、三步上手

```bash
# 0. 装依赖（首次，约 10–30 分钟，Bioconductor 包较多）
Rscript -e 'source("R/helpers.R"); load_pkgs(c("yaml","ggplot2","limma","sva","WGCNA","clusterProfiler","org.Hs.eg.db","ReactomePA","GSVA","msigdbr","glmnet","e1071","ranger","pROC","data.table","Seurat","harmony","CellChat","coloc","Rsamtools","TwoSampleMR"), install = TRUE)'

# 1. 合成数据跑通全链路（不需要下载任何真实数据）
Rscript 99_demo_data.R
Rscript run_all.R --demo

# 2. 换成真实数据：改 config/config.yml 里的路径，然后
Rscript run_all.R
```

第一次务必先跑 `--demo`：它会在 5–10 分钟内把 01→10 全过一遍，把环境问题（包版本、字体、内存、写权限）一次性暴露出来。合成数据里 chr1 的 900–1000 kb 区间被设计成「GWAS 与 LOX 表达共享因果变异」，所以 `coloc` 的 PP.H4 应该接近 1、`MR` 的 IVW 应该显著——如果 demo 跑出来不是这样，说明代码或环境有问题，此时不要急着换真实数据。

---

## 二、目录结构

```
glaucoma_R_pipeline/
├── run_all.R                    总调度（--demo / --only / --from / --config）
├── 99_demo_data.R               合成数据生成
├── R/helpers.R                  公共函数（配置、日志、GEO 解析、出图）
├── config/
│   ├── config.yml               真实数据配置（改这里，不改脚本）
│   └── config.demo.yml          合成数据配置
├── 01_geo_prepare.R             下载 / 探针重注释 / 标准化 / 批次校正 / QC
├── 02_deg_limma.R               差异表达 + 火山图 + 热图 + 共识基因
├── 03_wgcna.R                   共表达模块 + 与 IOP/MD/CCT 连续变量关联
├── 04_enrich_deconv.R           GO/KEGG/Reactome + GSVA + GSEA + PROGENy + 反卷积
├── 05_ml_signature.R            LASSO + SVM-RFE + RF/Boruta → 模型 → 外部验证
├── 06_scrna_cellchat.R          Seurat → 注释 → CellChat → 拟时序 → SCENIC 导出
├── 07_gwas_qc.R                 GWAS 清洗 / lambda / QQ / Manhattan / LDSC 格式
├── 08_coloc_smr.R               共定位（coloc）+ SMR/HEIDI 近似
├── 09_twas_focus.R              FUSION TWAS + FOCUS 精细定位
├── 10_mr.R                      工具变量 / IVW / Egger / 加权中位数 / PRESSO
├── 11_insilico_perturbation.R   单细胞网络上的计算扰动（数字敲除）
├── 12_mediation_mr.R            两阶段 MR：表达 → IOP → POAG + MVMR
├── 13_robustness_cmap.R         标签置换零分布 + 多管线 multiverse + CMap
├── 14_cross_species.R           人/猴/鼠跨物种一致性（替代动物模型验证）
├── 15_results_ledger.R          结果账本：正文数字全部 token 化、可溯源
├── 16_manuscript_assets.R       Table 1–6 成表 + 图 panel 清单（+ 可选拼图）
├── 17_fill_tokens.R             把账本数字回填进稿件，并查漏/查多余（回填后自动重算字数）
├── tools/                       装包辅助 + 稿件收尾脚本
│   ├── fill_wordcount.py        正文字数：HTML 与 docx 分别计数并同步表头，不一致就报警
│   └── fill_placeholders.py     ORCID / 仓库 URL / DOI：格式校验后同时写进 HTML 与 docx
├── data/{raw,annotation,processed}/
├── results/{tables,figures,scenic}/
└── logs/                        每步独立日志
```

---

## 三、数据放置清单

| 内容 | 放哪里 | 从哪来 |
|---|---|---|
| GEO series matrix | `data/raw/GSE27276_series_matrix.txt.gz` | 脚本自动下载；也可手动从 GEO FTP 拉 |
| GPL 注释 | `data/annotation/GPLxxxx_annot.txt.gz` | 脚本自动下载；公司给的注释文件也可直接放 |
| GWAS 汇总统计 | `data/raw/finngen_R12_H7_GLAUCOMA.gz` | FinnGen R12 端点 `H7_GLAUCOMA` |
| GTEx eQTL | `data/raw/GTEx_v8/Artery.v8.EUR.allpairs.txt.gz` + `.tbi` | GTEx Portal（需含 tabix 索引） |
| eQTL Catalogue | `data/raw/eQTL_Catalogue/<dataset>.all.tsv.gz` | eQTL Catalogue FTP |
| LD 参考面板 | `data/raw/1000G_EUR_Phase3_plink/1000G.EUR.QC.{bed,bim,fam}` | 1000G Phase3 |
| FUSION | `tools/FUSION/FUSION.assoc_test.R` | github.com/gusevlab/fusion_twas |
| FUSION 权重 | `data/raw/FUSION_weights/<tissue>/*.wgt.RDat` | FUSION 官网预计算权重 |
| 单细胞参考 | `data/processed/scRNA_ref.rds` | Seurat 对象，需含 `cell_type` 列 |
| 单细胞原始 | `scrna.tenx_dirs` 指向 CellRanger 输出目录 | `filtered_feature_bc_matrix/` |

GWAS 的列名各家差异极大，统一在 `config.yml` 的 `genetics.gwas.cols` 里映射（`snp/beta/se/pval/ea/nea/eaf/n`），脚本不写死列名。

---

## 四、每步在做什么、产出什么

### 课题一

| 步 | 核心动作 | 主要产出 | 最容易出错的地方 |
|---|---|---|---|
| 01 | 探针→基因（同基因取平均表达最高探针）、log2 判定、quantile 标准化、ComBat、PCA 离群标记 | `<GSE>_expr.rds` / `_pheno.rds`、`qc_summary.tsv` | 直接拿 GEO 矩阵不做重注释；`auto_remove_outliers` 默认 **false**，先看 PCA 再决定 |
| 02 | limma + BH，双阈值；验证集只做方向一致性 | `deg_<GSE>.tsv`、`deg_consensus.tsv` | 用训练集阈值去卡验证集会过度筛选 |
| 03 | `pickSoftThreshold` → `blockwiseModules` → 模块与 **IOP/MD/CCT/年龄** 关联（bicor） | `wgcna_module_trait.tsv`、`wgcna_hub_genes.tsv` | 只跟「病/非病」关联，信息量不足；表型列名要能在 pheno 里匹配到 |
| 04 | ORA（GO/KEGG/Reactome）、GSVA（Hallmark）、GSEA、PROGENy、NNLS 反卷积 | `enrich_*.tsv`、`gsva_hallmark_diff.tsv`、`deconv_fraction.tsv` | 用 LM22 套眼组织；msigdbr 版本要写进方法学 |
| 05 | LASSO（CV）/ SVM-RFE（Guyon 算法手写）/ RF-Boruta → 交集 → logistic → ROC/校准/DCA/列线图 → 外部队列验证 | `ml_signature.tsv`、`model_performance.tsv`、ROC/校准/DCA 图 | 特征筛选与评估用同一批样本（数据泄漏）；只报 AUC 不报 CI |
| 06 | Harmony 去批次 → 聚类 → marker 注释 → CellChat 条件对比 → Monocle3 | `cellchat_compare.tsv`、`sc_deg_TM.tsv`、`sc_pseudotime_TM.pdf` | 不做 Harmony 会按供体聚类；注释必须 DotPlot 人工复核 |

### 课题二

| 步 | 核心动作 | 主要产出 | 最容易出错的地方 |
|---|---|---|---|
| 07 | 列名统一 → 回文 SNP / MAF / INFO / N 过滤 → dedup → lambda → QQ/Manhattan → LDSC 格式 | `gwas_clean.rds`、`gwas_qc_summary.tsv`、`gwas_ldsc.sumstats.gz` | 先过滤再算 lambda；indel 未剔除导致等位基因对不齐 |
| 08 | 位点 ±500 kb 取 eQTL → 等位基因对齐（rsID 优先，GTEx 用 chr_pos）→ `coloc.abf` → 敏感性 → SMR/HEIDI 近似 | `coloc_summary.tsv`、`smr_heidi.tsv` | 只报 PP.H4 不报 H3；先验 p12 不做敏感性分析 |
| 09 | FUSION 逐染色体 assoc test → 汇总 → Bonferroni → FOCUS | `twas_all.tsv`、`twas_significant.tsv`、`focus_<tissue>.tsv` | 把 TWAS 结果说成「实测表达」；权重与 GWAS 的 build/等位基因不一致 |
| 10 | F 统计量 → LD 剪枝 → harmonise → 多方法 MR → 多效性/异质性/留一法 → FDR | `mr_results.tsv`、`mr_sensitivity.tsv`、森林图 | 弱工具变量（F<10）、样本重叠未交代、多重检验未校正 |

### 课题三（纯干实验版，用 `--only 11,12,13` 触发）

| 步 | 核心动作 | 主要产出 | 最容易出错的地方 |
|---|---|---|---|
| 11 | 单细胞建基因网络 → 目标基因置零传播 → 与疾病特征比方向 + 度数匹配置换 | `perturb_summary.tsv`、`perturb_<gene>_<celltype>.tsv`、散点图 | 说成"因果"；细胞类型标签与 config 不一致 |
| 12 | 两阶段 MR（表达→IOP→POAG）+ delta 法算间接效应 + MVMR 算直接效应 | `mr_mediation.tsv`、`mr_mediation_forest.pdf` | 中介比例 >1 硬解释；IOP GWAS 无位置列导致剪枝退化 |
| 13 | 标签置换零分布 → 8 条管线 multiverse → CMap 药物连通性（可选） | `robustness_permutation.tsv`、`robustness_multiverse.tsv`、`cmap_connectivity.tsv` | 只报最好看的 AUC；把 CMap 预测写成疗效 |
| 14 | 物种同源映射 → TM/JCT 伪批量差异 → 方向一致比例 + Spearman → 人源基因集模块打分 | `xspecies_concordance.tsv`、`xspecies_de.tsv`、`xspecies_modulescore.tsv`、三张图 | 写"跨物种验证"而实际只做 toupper 近似映射；concordance<0.7 仍称保守 |

---

## 五、结果与论文图表的对应

| 论文图 | 来源 | 建议呈现 |
|---|---|---|
| 图 1 研究设计 | 手绘 | 数据流 + 假说示意 |
| 图 2 QC 与 PCA | `pca_*_raw/clean.pdf` | 校正前后并排，标注离群样本 |
| 图 3 差异谱 | `volcano_*`、`heatmap_top50_*` | 火山图 + 共识基因韦恩图 |
| 图 4 模块-表型 | `wgcna_module_trait_*` | 格子内 `r (P)`，星号标显著性 |
| 图 5 通路 | `enrich_dotplot_*`、`gsva_hallmark_bar` | GO BP + Hallmark 各一 |
| 图 6 细胞比例 | `deconv_box`、`deconv_module_correlation` | 箱线图 + 模块相关性热图 |
| 图 7 模型 | `roc_train`、`roc_valid_*`、`calibration`、`dca`、`nomogram` | 训练 + 至少一个外部验证 |
| 图 8 单细胞 | `sc_umap_celltype`、`cellchat_circle_*`、`cellchat_compare_rank` | UMAP + 通讯圈图 + 差异通路 |
| 图 9 遗传定位 | `coloc_locus_*`、`twas_manhattan`、`mr_forest` | coloc 散点 + TWAS + 森林图 |
| 表 S1 | `qc_summary.tsv` | 数据集、平台、样本量、批次 |
| 表 S2 | `gwas_qc_summary.tsv` | 每步保留 SNP 数（方法学必备） |

---

## 六、常见报错与处理

| 报错 | 原因 | 处理 |
|---|---|---|
| `找不到项目根目录` | 不在项目目录下运行 | `setwd()` 到 `glaucoma_R_pipeline/` |
| `缺少输入文件：...` | 数据没放对位置 | 按第三节的清单放；路径写进 `config.yml` |
| `series matrix 表区定位失败` | 下载不完整 | 删掉 `data/raw/*.gz` 重跑 01 |
| `需要 yaml 包` / 其他包缺失 | 依赖未装 | 用第一节第 0 步的 `install = TRUE` |
| WGCNA `样本量不足 15` | 数据集太小 | 该数据集只做 DEG + GSVA，不做 WGCNA |
| `coloc` 无输出 | eQTL 区间内 SNP < 20 | 加大 `window_kb`，或换组织 |
| FUSION 失败 | 缺 LD 面板/权重 | 检查 `ref_ld_chr` 前缀与 `.bed` 是否存在 |
| CellChat 报 `净计数矩阵` 相关错误 | 细胞类型只有 1 类 | 提高聚类分辨率或放宽 QC |

---

## 七、必须自己补的三件事

1. **真实临床变量**：`IOP / MD / CCT / age` 在 GEO 里多数缺失，需要从原文补充材料或课题组病例库补。没有连续变量，03 的价值会掉一半。
2. **单细胞注释的人工复核**：`06` 里的 marker 打分只是初稿，必须 `DotPlot` + 文献 marker 表逐群核对，尤其是 TM / JEC / SC 的边界。
3. **实验验证的设计**：脚本只到「候选」层面。qPCR / WB / 免疫荧光 / 原代 HTM 牵张实验的方案在开题方案「课题三」，不在本目录。

---

## 八、已知限制（写文章时要在 Limitations 里交代）

- 反卷积默认用单细胞参考自建 signature；没有参考时退回 xCell，属于间接推断。
- SMR/HEIDI 是 R 内近似实现（Wald 比 + 加权 Q 检验），发表时建议用官方 SMR 二进制复核（`08` 末尾有命令行提示）。
- FUSION 权重依赖第三方预计算文件，缺失时 `09` 会跳过而不是报错。
- TWAS 的「遗传预测表达」不等于实测表达，结论表述要留余地。

---

## 九、成稿阶段（15 → 16 → 17）

分析跑完后，正文数字不要手工抄。整条链路是「账本 → 素材 → 回填」：

```bash
Rscript 15_results_ledger.R          # 把所有进正文的数字抽成 token，并记录来源文件与字段
Rscript 16_manuscript_assets.R       # Table 1–6 成 xlsx；figure_manifest 标注每张 panel 是否已生成
Rscript 17_fill_tokens.R             # 把 token 替换成数字；报告 MISSING / UNUSED
```

产出：

| 文件 | 用途 |
|---|---|
| `results/tables/results_ledger.tsv` | 每个 token 的显示值、原始值、来源文件、字段、备注 |
| `results/tables/missing_steps.tsv` | 哪些步骤还没有产物（对应段落还写不了） |
| `results/tables/results_ledger.md` | 分组可读版，可直接发给导师看 |
| `results/submission/Tables1-6.xlsx` | 三线表样式的成表 |
| `results/submission/figure_manifest.tsv` | 9 张图、32 个 panel 的文件状态 |
| `results/submission/token_fill_report.tsv` | OK / MISSING_IN_LEDGER / UNUSED_IN_TEXT |
| `results/session_info.txt` | R 与包版本（方法学段落与补充材料 S8 要用） |

稿件模板在 `../output/manuscript/formatted-manuscript.html`，转换好的 Word 在 `../output/manuscript/POAG_TM_manuscript_draft.docx`。模板里所有 `{{TOKEN}}` 由 `17` 回填；所有 `[XX]` 和方括号需人工填写，清单见 `docs/成稿手册-从结果到投稿.md` 第十二节。

`17` 的用法：

```bash
Rscript 17_fill_tokens.R                                   # 默认读 config 里的 manuscript.template
Rscript 17_fill_tokens.R --in <模板.html> --out <回填.html>
Rscript 17_fill_tokens.R --docx                            # 回填后直接转 docx（需在 config 填 manuscript.converter）
```

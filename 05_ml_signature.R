# =============================================================================
# 05_ml_signature.R —— 课题一（5）：三算法交叉筛选 + 诊断模型 + 外部验证
#
# 输入：results/tables/deg_consensus.tsv（可缺）
#       results/tables/wgcna_gene_modules.tsv（可缺）
#       data/processed/<trainGSE>_expr.rds / _pheno.rds
#       data/processed/<validGSE>_expr.rds / _pheno.rds
# 输出：results/tables/ml_candidate_genes.tsv
#       results/tables/ml_lasso_coef.tsv / ml_svmrfe_curve.tsv / ml_rf_importance.tsv
#       results/tables/ml_signature.tsv           最终基因组合及其系数
#       results/tables/model_performance.tsv      训练/验证集 AUC、敏感度、特异度
#       results/tables/model_calibration.tsv      H-L 检验、十分位偏差、DCA 阈值区间
#       results/tables/model_optimism.tsv         bootstrap optimism 与校正后 AUC
#       results/figures/roc_train.pdf / roc_valid_<GSE>.pdf
#       results/figures/calibration.pdf / dca.pdf / nomogram.pdf
#
# 方法学要点（避免被质疑"数据泄漏"）：
#   1. 特征筛选只在训练集做；验证集从头到尾不参与选基因
#   2. 三种算法取交集，单用 LASSO 会被认为过拟合
#   3. LASSO 的 lambda 用交叉验证选，但交叉验证必须包在训练集内部
#   4. 报告 AUC 的 95% CI，而不是只报一个点估计
# =============================================================================

source("R/helpers.R")
cfg <- read_config()
log_init("05_ml_signature", cfg)
load_pkgs(c("ggplot2", "glmnet", "e1071", "ranger", "pROC"))
set.seed(cfg_get(cfg, "project.seed", 1))

demo <- isTRUE(cfg_get(cfg, "run.demo", FALSE))
train_gse <- if (demo) "DEMO" else cfg_get(cfg, "bulk.train.gse", "GSE27276")

expr  <- readRDS(proj_path("data", "processed", paste0(train_gse, "_expr.rds")))
pheno <- readRDS(proj_path("data", "processed", paste0(train_gse, "_pheno.rds")))
y <- factor(pheno$group, levels = c("control", "case"))
if (any(table(y) < 5)) stop("分组样本过少（每组需 >=5），不适合建模")

# ------------------------------------------------- 1) 候选基因 --------------
cand <- list()
f_deg <- res_path(cfg, "tables", "deg_consensus.tsv")
if (file.exists(f_deg)) {
  cand$deg <- utils::read.delim(f_deg, stringsAsFactors = FALSE)$gene
} else {
  tt <- utils::read.delim(res_path(cfg, "tables", paste0("deg_", train_gse, ".tsv")),
                          stringsAsFactors = FALSE)
  cand$deg <- tt$gene[tt$sig != "NS"]
}
f_w <- res_path(cfg, "tables", "wgcna_hub_genes.tsv")
if (file.exists(f_w)) {
  hub <- utils::read.delim(f_w, stringsAsFactors = FALSE)
  cand$wgcna <- unique(hub$gene)
}
cand$all <- Reduce(intersect, cand)
if (length(cand$all) < 10) {
  log_msg("交集候选基因仅 ", length(cand$all), " 个，回退为 DEG ∩ 表达矩阵")
  cand$all <- intersect(cand$deg, rownames(expr))
}
cand$all <- intersect(cand$all, rownames(expr))
if (!length(cand$all)) {
  stop("候选基因为 0：常见原因是 DEG 阈值过严（bulk.deg.fdr/logfc）或数据无信号。",
       "先看 results/tables/deg_", train_gse, ".tsv 里 sig 列非 NS 的行数。", call. = FALSE)
}
log_msg("候选基因 ", length(cand$all), " 个")
write_tsv(data.frame(gene = cand$all, source = "DEG∩WGCNA"), "ml_candidate_genes", cfg)

X <- t(expr[cand$all, , drop = FALSE])
Xz <- scale(X)   # 建模统一用 z-score，便于跨平台迁移

# ------------------------------------------------- 2) LASSO -----------------
cvfit <- glmnet::cv.glmnet(Xz, y, family = "binomial", alpha = 1,
                           nfolds = min(cfg_get(cfg, "bulk.ml.lasso_folds", 10), sum(y == "case")),
                           type.measure = "deviance", standardize = FALSE)
lam <- if (cfg_get(cfg, "bulk.ml.lasso_lambda", "lambda.min") == "lambda.1se") cvfit$lambda.1se else cvfit$lambda.min
coef_mat <- as.matrix(stats::coef(cvfit, s = lam))
lasso_df <- data.frame(gene = rownames(coef_mat), coef = as.numeric(coef_mat[, 1]))
lasso_df <- lasso_df[lasso_df$coef != 0 & lasso_df$gene != "(Intercept)", ]
lasso_df <- lasso_df[order(-abs(lasso_df$coef)), ]
write_tsv(lasso_df, "ml_lasso_coef", cfg)
log_msg("LASSO 选中 ", nrow(lasso_df), " 个基因（lambda = ", signif(lam, 4), "）")

g_lasso <- ggplot2::ggplot(data.frame(lambda = cvfit$lambda, cvm = cvfit$cvm, cvsd = cvfit$cvsd),
                           ggplot2::aes(log(lambda), cvm)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = cvm - cvsd, ymax = cvm + cvsd), width = 0.02, colour = "grey60") +
  ggplot2::geom_point(size = 0.8, colour = "#C0392B") +
  ggplot2::geom_vline(xintercept = log(lam), linetype = "dashed", linewidth = 0.4) +
  ggplot2::labs(x = expression(log(lambda)), y = "binomial deviance", title = "LASSO 10-fold CV") +
  theme_pub()
save_fig(g_lasso, "lasso_cv", cfg, 5, 4)

# ------------------------------------------------- 3) SVM-RFE --------------
# 经典 Guyon 算法：线性 SVM 权重平方排序，分批剔除最差特征，逐规模记录 CV 准确率
svm_rfe <- function(X, y, sizes, v = 5, drop_frac = 0.2) {
  keep <- colnames(X)
  sizes <- sort(unique(sizes[sizes <= length(keep)]), decreasing = TRUE)
  if (!length(sizes)) sizes <- min(5, length(keep))

  cv_acc <- function(cols) {
    set.seed(1)
    folds <- sample(rep(seq_len(v), length.out = nrow(X)))
    preds <- rep(NA_character_, nrow(X))
    for (k in seq_len(v)) {
      tr <- folds != k
      m <- e1071::svm(X[tr, cols, drop = FALSE], y[tr], kernel = "linear", cost = 1, scale = FALSE)
      preds[!tr] <- as.character(stats::predict(m, X[!tr, cols, drop = FALSE]))
    }
    mean(preds == as.character(y), na.rm = TRUE)
  }

  acc <- data.frame(size = integer(0), accuracy = numeric(0))
  gene_sets <- list()
  for (target in sizes) {
    while (length(keep) > target) {
      m <- e1071::svm(X[, keep, drop = FALSE], y, kernel = "linear", cost = 1, scale = FALSE)
      r <- drop(crossprod(m$coefs, m$SV))^2          # 线性核权重平方
      n_drop <- min(length(keep) - target, max(1, floor(drop_frac * length(keep))))
      keep <- keep[-order(r)[seq_len(n_drop)]]       # 权重小的先剔除
    }
    a <- cv_acc(keep)
    acc <- rbind(acc, data.frame(size = length(keep), accuracy = a))
    gene_sets[[as.character(length(keep))]] <- keep
    log_msg("  SVM-RFE size=", length(keep), " CV acc=", round(a, 3))
  }
  best <- acc$size[which.max(acc$accuracy)]
  list(accuracy = acc, genes = gene_sets[[as.character(best)]], best_size = best)
}
sizes <- cfg_get(cfg, "bulk.ml.svm_rfe_subsets", c(2, 3, 4, 5, 6, 8, 10, 12, 15, 20))
rfe <- svm_rfe(Xz, y, sizes)
write_tsv(rfe$accuracy, "ml_svmrfe_curve", cfg)
best_size <- rfe$best_size
svm_genes <- rfe$genes
log_msg("SVM-RFE 最优规模 ", best_size, "（CV 准确率 ",
        round(max(rfe$accuracy$accuracy), 3), "）")

# Fig5B：SVM-RFE 准确率曲线。没有这张图，「为什么停在这个基因数」在正文里
# 只是一句结论；有了它，审稿人能自己看到曲线在高规模处是平的还是掉下来的。
if (nrow(rfe$accuracy) >= 2) {
  acc_rng <- range(rfe$accuracy$accuracy, na.rm = TRUE)
  g_rfe <- ggplot2::ggplot(rfe$accuracy, ggplot2::aes(size, accuracy)) +
    ggplot2::geom_line(colour = "#1F6FB4", linewidth = 0.6) +
    ggplot2::geom_point(colour = "#1F6FB4", size = 1.9) +
    ggplot2::geom_vline(xintercept = best_size, linetype = "dashed",
                        colour = "#C0392B", linewidth = 0.4) +
    # 标注放左上角而不是贴着竖线：贴线时若最优规模是最大 size，文本会被面板裁掉
    ggplot2::annotate("text", x = min(rfe$accuracy$size), y = acc_rng[2],
                      label = paste0("selected size: n = ", best_size),
                      hjust = 0, vjust = 1.3, size = 3.1, colour = "#C0392B") +
    ggplot2::scale_x_continuous(breaks = rfe$accuracy$size) +
    ggplot2::scale_y_continuous(limits = c(max(0, acc_rng[1] - 0.04), min(1, acc_rng[2] + 0.03))) +
    ggplot2::labs(x = "Number of genes retained", y = "Cross-validated accuracy",
                  title = "SVM-RFE feature selection") +
    theme_pub()
  save_fig(g_rfe, "svmrfe_curve", cfg, 5.2, 4)
}

# ------------------------------------------------- 4) 随机森林 / Boruta -----
rf_genes <- character(0)
if (isTRUE(cfg_get(cfg, "bulk.ml.use_boruta", TRUE)) && requireNamespace("Boruta", quietly = TRUE)) {
  set.seed(cfg_get(cfg, "project.seed", 1))
  bor <- Boruta::Boruta(x = as.data.frame(Xz), y = y, ntree = 500, doTrace = 0)
  rf_genes <- Boruta::getSelectedAttributes(bor, withTentative = FALSE)
  imp <- as.data.frame(Boruta::attStats(bor))
  imp$gene <- rownames(imp)
  write_tsv(imp[order(-imp$medianImp), ], "ml_rf_importance", cfg)
  log_msg("Boruta 确认 ", length(rf_genes), " 个基因")
} else {
  set.seed(cfg_get(cfg, "project.seed", 1))
  rf <- ranger::ranger(x = as.data.frame(Xz), y = y, num.trees = cfg_get(cfg, "bulk.ml.rf_ntree", 1000),
                       importance = "permutation", classification = TRUE)
  imp <- data.frame(gene = names(rf$variable.importance),
                    importance = as.numeric(rf$variable.importance))
  imp <- imp[order(-imp$importance), ]
  write_tsv(imp, "ml_rf_importance", cfg)
  rf_genes <- head(imp$gene, 20)
  log_msg("ranger 取重要性前 20 个基因")
}

# ------------------------------------------------- 5) 交集与建模 ------------
sets <- list(LASSO = lasso_df$gene, SVMRFE = svm_genes, RF = rf_genes)
sig <- Reduce(intersect, sets)
if (length(sig) < 2) {
  log_msg("三法交集不足 2 个，退化为「至少被两种方法选中」")
  tb <- table(unlist(sets))
  sig <- names(tb)[tb >= 2]
  if (length(sig) < 2) {
    log_msg("仍不足，改用 LASSO 前若干基因")
    sig <- head(lasso_df$gene, 5)
  }
}
log_msg("最终 signature：", paste(sig, collapse = ", "))

dat <- data.frame(Xz[, sig, drop = FALSE], group = y)
mod <- stats::glm(group ~ ., data = dat, family = stats::binomial())
coefs <- summary(mod)$coefficients
sig_tab <- data.frame(gene = rownames(coefs),
                      is_intercept = rownames(coefs) == "(Intercept)",
                      estimate = coefs[, 1],
                      se = coefs[, 2], z = coefs[, 3], p = coefs[, 4], row.names = NULL)
# is_intercept 是给下游（15_results_ledger）用的：summary(glm) 第一行是截距，
# 不是基因，不能进 signature 基因数，也不能写进正文的基因列表
write_tsv(sig_tab, "ml_signature", cfg)

prob_train <- as.numeric(stats::predict(mod, type = "response"))
roc_train <- pROC::roc(y, prob_train, quiet = TRUE)
ci_train <- as.numeric(pROC::ci.auc(roc_train))

# ------------------------------------------------- 6) 评价曲线 -------------
plot_roc <- function(roc, title) {
  df <- data.frame(fpr = 1 - roc$specificities, tpr = roc$sensitivities)
  df <- df[order(df$fpr, df$tpr), ]
  ggplot2::ggplot(df, ggplot2::aes(fpr, tpr)) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_step(colour = "#C0392B", linewidth = 0.8) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "1 - specificity", y = "sensitivity",
                  title = title,
                  subtitle = sprintf("AUC = %.3f (95%%CI %.3f-%.3f)", ci_train[2], ci_train[1], ci_train[3])) +
    theme_pub()
}
save_fig(plot_roc(roc_train, paste0("Training set: ", train_gse)), "roc_train", cfg, 4.6, 4.6)

# 校准曲线（bootstrap）
calib <- function(y, p, g = 5) {
  br <- unique(stats::quantile(p, probs = seq(0, 1, length.out = g + 1)))
  grp <- cut(p, br, include.lowest = TRUE, labels = FALSE)
  obs <- tapply(y == "case", grp, mean)
  prd <- tapply(p, grp, mean)
  data.frame(pred = as.numeric(prd), obs = as.numeric(obs))
}
cd <- calib(y, prob_train)
g_cal <- ggplot2::ggplot(cd, ggplot2::aes(pred, obs)) +
  ggplot2::geom_abline(linetype = "dashed", colour = "grey55") +
  ggplot2::geom_line(colour = "#1F6FB4", linewidth = 0.8) +
  ggplot2::geom_point(size = 2, colour = "#1F6FB4") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(x = "predicted probability", y = "observed proportion", title = "Calibration") +
  theme_pub()
save_fig(g_cal, "calibration", cfg, 4.4, 4.2)

# 决策曲线（net benefit，手写避免额外依赖）
net_benefit <- function(y, p, thresholds = seq(0.01, 0.99, 0.01)) {
  prev <- mean(y == "case")
  out <- data.frame(threshold = thresholds, model = NA_real_, all = NA_real_, none = 0)
  for (i in seq_along(thresholds)) {
    t <- thresholds[i]
    tp <- sum(p >= t & y == "case"); fp <- sum(p >= t & y == "control")
    n <- length(y)
    out$model[i] <- tp / n - fp / n * t / (1 - t)
    out$all[i] <- prev - (1 - prev) * t / (1 - t)
  }
  out
}
nb <- net_benefit(y, prob_train)
nb_long <- rbind(
  data.frame(threshold = nb$threshold, nb = nb$model, strategy = "Model"),
  data.frame(threshold = nb$threshold, nb = nb$all, strategy = "Treat all"),
  data.frame(threshold = nb$threshold, nb = nb$none, strategy = "Treat none")
)
g_dca <- ggplot2::ggplot(nb_long, ggplot2::aes(threshold, nb, colour = strategy)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::scale_colour_manual(values = c(Model = "#C0392B", `Treat all` = "grey40", `Treat none` = "grey75"), name = NULL) +
  ggplot2::coord_cartesian(ylim = c(-0.05, max(nb$model, na.rm = TRUE) * 1.1)) +
  ggplot2::labs(x = "threshold probability", y = "net benefit", title = "Decision curve") +
  theme_pub()
save_fig(g_dca, "dca", cfg, 5, 4)

# ------------------------------------------------- 6b) 可引用的校准数字 ----
# 图能看，但正文里要写"P = 0.41""偏差 0.03""阈值 0.12 到 0.68"，
# 这些必须由脚本给出而不是手填。Hosmer-Lemeshow 用十分位手写实现，
# 不引入 ResourceSelection 依赖（公式与 hoslem.test 的 g=10 一致）。
cal_num <- local({
  g <- 10
  br <- unique(stats::quantile(prob_train, probs = seq(0, 1, length.out = g + 1)))
  grp <- cut(prob_train, br, include.lowest = TRUE, labels = FALSE)
  obs_n <- tapply(as.integer(y == "case"), grp, sum)
  exp_n <- tapply(prob_train, grp, sum)
  n_g <- tapply(as.integer(y == "case"), grp, length)
  pbar <- exp_n / n_g
  keep <- is.finite(pbar) & pbar > 0 & pbar < 1 & n_g > 0
  chi <- sum((obs_n[keep] - exp_n[keep])^2 / (n_g[keep] * pbar[keep] * (1 - pbar[keep])))
  df <- sum(keep) - 2
  list(chi2 = chi, df = df,
       p = stats::pchisq(chi, df, lower.tail = FALSE),
       max_gap = max(abs(obs_n[keep] / n_g[keep] - pbar[keep])),
       n_groups = sum(keep))
})

# 决策曲线上"优于全治与不治"的连续区间：取最长的一段，
# 而不是取所有满足点的极值——否则跨度会被两端的零散点撑大
better <- nb$model > pmax(nb$all, nb$none) + 1e-9
rl <- rle(better)
thr_lo <- thr_hi <- NA_real_
if (any(rl$values)) {
  k <- which(rl$values)[which.max(rl$lengths[rl$values])]
  ends <- cumsum(rl$lengths)
  i1 <- if (k == 1) 1 else ends[k - 1] + 1
  i2 <- ends[k]
  thr_lo <- nb$threshold[i1]; thr_hi <- nb$threshold[i2]
}

cal_tab <- data.frame(
  hoslem_chi2 = round(cal_num$chi2, 3), hoslem_df = cal_num$df,
  hoslem_p = signif(cal_num$p, 4), max_abs_gap = round(cal_num$max_gap, 4),
  dca_thr_lo = round(thr_lo, 2), dca_thr_hi = round(thr_hi, 2),
  n_calib_groups = cal_num$n_groups, n = length(y))
write_tsv(cal_tab, "model_calibration", cfg)
log_msg("校准：H-L P = ", signif(cal_num$p, 3), "，十分位最大偏差 ",
        round(cal_num$max_gap, 3), "；DCA 优势区间 ", round(thr_lo, 2), "–", round(thr_hi, 2))

# ------------------------------------------------- 6c) Bootstrap 校正 AUC ---
# 正文那句 "a bootstrap-corrected estimate of X" 必须由脚本给出，不能手填。
# 口径（Harrell optimism）：每次有放回重抽，在重抽样本上**重拟合系数**，
# 得到该样本的表观 AUC，再把这个模型拿回原样本算 AUC，两者之差即该次 optimism；
# 取均值后从表观 AUC 里扣掉。
# 只校正系数层的过拟合，不校正 LASSO/Boruta 的基因筛选——筛选那层靠外部验证队列兜底，
# 所以 note 列把这句口径写进表里，正文和补充材料按它表述，别让人误以为全流程都校正过。
# 配置里历史上叫 boot_n（真实配置设的 1000），统一按 bootstrap 读，保留旧键兜底
n_boot <- as.integer(cfg_get(cfg, "bulk.ml.bootstrap",
                             cfg_get(cfg, "bulk.ml.boot_n", 500)))
set.seed(cfg_get(cfg, "project.seed", 1))
auc_app <- as.numeric(pROC::auc(roc_train))
auc_in <- auc_out <- opt <- rep(NA_real_, n_boot)
for (b in seq_len(n_boot)) {
  idx <- sample.int(length(y), replace = TRUE)
  yb <- y[idx]
  if (length(unique(yb)) < 2) next          # 极端重抽只含一类，AUC 无定义
  db <- dat[idx, , drop = FALSE]
  fit_b <- try(suppressWarnings(stats::glm(group ~ ., data = db, family = stats::binomial())),
               silent = TRUE)
  if (inherits(fit_b, "try-error")) next
  p_in <- try(as.numeric(stats::predict(fit_b, newdata = db, type = "response")), silent = TRUE)
  p_out <- try(as.numeric(stats::predict(fit_b, newdata = dat, type = "response")), silent = TRUE)
  if (inherits(p_in, "try-error") || inherits(p_out, "try-error")) next
  a_in <- try(as.numeric(pROC::auc(pROC::roc(yb, p_in, quiet = TRUE))), silent = TRUE)
  a_out <- try(as.numeric(pROC::auc(pROC::roc(y, p_out, quiet = TRUE))), silent = TRUE)
  if (inherits(a_in, "try-error") || inherits(a_out, "try-error")) next
  if (!is.finite(a_in) || !is.finite(a_out)) next
  auc_in[b] <- a_in; auc_out[b] <- a_out; opt[b] <- a_in - a_out
}
ok <- is.finite(opt)
n_ok <- sum(ok)
opt_tab <- data.frame(
  n_boot = n_boot, n_ok = n_ok,
  auc_apparent = round(auc_app, 3),
  auc_boot_mean = if (n_ok) round(mean(auc_in[ok]), 3) else NA_real_,
  optimism_mean = if (n_ok) round(mean(opt[ok]), 4) else NA_real_,
  optimism_sd = if (n_ok > 1) round(stats::sd(opt[ok]), 4) else NA_real_,
  auc_corrected = if (n_ok) round(auc_app - mean(opt[ok]), 3) else NA_real_,
  n_params = length(sig) + 1L,
  note = "fixed signature; coefficient-level optimism only, gene selection not corrected")
write_tsv(opt_tab, "model_optimism", cfg)
if (n_ok >= 20) {
  log_msg("Bootstrap（", n_boot, " 次，成功 ", n_ok, " 次）：表观 AUC ", round(auc_app, 3),
          "，optimism ", round(mean(opt[ok]), 4), "，校正后 ", round(auc_app - mean(opt[ok]), 3))
} else {
  log_msg("[警告] bootstrap 仅 ", n_ok, " 次可用（<20），校正 AUC 无法给出；",
          "常见原因是样本量太小导致重抽后完全分离。正文该处会保留 [XX]")
}


# 列线图（rms，可选）
if (requireNamespace("rms", quietly = TRUE)) {
  d_nomo <- data.frame(Xz[, sig, drop = FALSE], group = as.integer(y == "case"))
  names(d_nomo) <- make.names(names(d_nomo))
  dd <- rms::datadist(d_nomo); options(datadist = "dd")
  fml <- stats::as.formula(paste("group ~", paste(setdiff(names(d_nomo), "group"), collapse = " + ")))
  fit_rms <- rms::lrm(fml, data = d_nomo)
  nom <- rms::nomogram(fit_rms, fun = stats::plogis)
  f_nomo <- proj_path(cfg_get(cfg, "project.outdir", "results"), "figures", "nomogram.pdf")
  grDevices::pdf(f_nomo, width = 8, height = 6)
  plot(nom, xfrac = 0.35)
  grDevices::dev.off()
  options(datadist = NULL)
}

# ------------------------------------------------- 7) 外部验证 -------------
perf <- data.frame(dataset = train_gse, n = length(y),
                   n_case = sum(y == "case"), n_control = sum(y == "control"),
                   AUC = ci_train[2], AUC_low = ci_train[1], AUC_high = ci_train[3],
                   cutoff = as.numeric(pROC::coords(roc_train, "best", ret = "threshold")[1]),
                   sensitivity = as.numeric(pROC::coords(roc_train, "best", ret = "sensitivity")[1]),
                   specificity = as.numeric(pROC::coords(roc_train, "best", ret = "specificity")[1]))

valid_gse <- if (demo) c("DEMO_V") else vapply(cfg_get(cfg, "bulk.valid", list()), function(x) x$gse, character(1))
for (v in valid_gse) {
  fv <- proj_path("data", "processed", paste0(v, "_expr.rds"))
  if (!file.exists(fv)) next
  ev <- readRDS(fv); pv <- readRDS(proj_path("data", "processed", paste0(v, "_pheno.rds")))
  gv <- intersect(sig, rownames(ev))
  if (length(gv) < length(sig)) log_msg(v, " 缺失 ", length(sig) - length(gv), " 个 signature 基因，按可用基因预测")
  if (length(gv) < 2) { log_msg(v, " 可用基因不足，跳过"); next }
  Xv <- scale(t(ev[gv, , drop = FALSE]))
  colnames(Xv) <- gv
  newd <- as.data.frame(Xv)
  # 缺失基因用 0（即训练集均值）代替，等价于「该基因取人群平均表达」
  for (g in setdiff(sig, gv)) newd[[g]] <- 0
  pv_prob <- as.numeric(stats::predict(mod, newdata = newd[, sig, drop = FALSE], type = "response"))
  yv <- factor(pv$group, levels = c("control", "case"))
  rocv <- pROC::roc(yv, pv_prob, quiet = TRUE)
  civ <- as.numeric(pROC::ci.auc(rocv))
  perf <- rbind(perf, data.frame(dataset = v, n = length(yv),
                                 n_case = sum(yv == "case"), n_control = sum(yv == "control"),
                                 AUC = civ[2], AUC_low = civ[1], AUC_high = civ[3],
                                 cutoff = perf$cutoff[1],
                                 sensitivity = mean(pv_prob[yv == "case"] >= perf$cutoff[1]),
                                 specificity = mean(pv_prob[yv == "control"] < perf$cutoff[1])))
  # 验证集 ROC 单独出图（副标题用验证集自己的 CI）
  dfv <- data.frame(fpr = 1 - rocv$specificities, tpr = rocv$sensitivities)
  dfv <- dfv[order(dfv$fpr, dfv$tpr), ]
  g_v <- ggplot2::ggplot(dfv, ggplot2::aes(fpr, tpr)) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_step(colour = "#1F6FB4", linewidth = 0.8) + ggplot2::coord_equal() +
    ggplot2::labs(x = "1 - specificity", y = "sensitivity", title = paste0("Validation: ", v),
                  subtitle = sprintf("AUC = %.3f (95%%CI %.3f-%.3f)", civ[2], civ[1], civ[3])) +
    theme_pub()
  save_fig(g_v, paste0("roc_valid_", v), cfg, 4.6, 4.6)
}
write_tsv(perf, "model_performance", cfg)
print(perf)

saveRDS(list(model = mod, genes = sig, coefs = sig_tab, lambda = lam),
        proj_path("data", "processed", "diagnostic_model.rds"))
log_msg("05 完成")

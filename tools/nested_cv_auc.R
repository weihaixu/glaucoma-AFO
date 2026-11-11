# 嵌套交叉验证：每折内部重新做基因筛选 + 拟合，折外预测，衡量真实泛化能力。
# 现有的 bootstrap optimism 只在固定基因集上重复拟合，选择偏倚没进循环，
# 所以 36 例、89 个候选基因、AUC=1.000 也能给出 optimism≈0。
suppressMessages({library(glmnet); library(pROC)})
source("R/helpers.R")
cfg <- read_config()
expr <- readRDS("data/processed/GSE27276_expr.rds")
pheno <- readRDS("data/processed/GSE27276_pheno.rds")
cand <- read.delim("results/tables/ml_candidate_genes.tsv")$gene
y <- factor(pheno$group, levels = c("control", "case"))
X <- scale(t(expr[intersect(cand, rownames(expr)), , drop = FALSE]))

set.seed(20261111)
K <- 5; R <- 20
aucs <- numeric(0)
for (r in seq_len(R)) {
  fold <- sample(rep(seq_len(K), length.out = nrow(X)))
  for (k in seq_len(K)) {
    tr <- fold != k; te <- fold == k
    if (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2) next
    cv <- try(glmnet::cv.glmnet(X[tr, ], y[tr], family = "binomial", alpha = 1,
                                nfolds = 5, standardize = FALSE), silent = TRUE)
    if (inherits(cv, "try-error")) next
    pr <- try(as.numeric(predict(cv, X[te, , drop = FALSE], s = cv$lambda.min, type = "response")), silent = TRUE)
    if (inherits(pr, "try-error")) next
    a <- try(as.numeric(pROC::auc(pROC::roc(y[te], pr, quiet = TRUE, direction = "<"))), silent = TRUE)
    if (!inherits(a, "try-error")) aucs <- c(aucs, a)
  }
}
cat(sprintf("嵌套 CV：%d 折 x %d 次，有效折 %d\n", K, R, length(aucs)))
cat(sprintf("折外 AUC 中位数 %.3f（IQR %.3f–%.3f），均值 %.3f\n",
            median(aucs), quantile(aucs, .25), quantile(aucs, .75), mean(aucs)))
write.table(data.frame(rep = seq_along(aucs), auc = aucs),
            "results/tables/model_nested_cv.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

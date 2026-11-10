# =============================================================================
# run_all.R —— 总调度：按顺序调用各步骤脚本，失败即停
#
# 用法：
#   Rscript run_all.R                  # 跑 config.yml 里 run.steps 指定的步骤
#   Rscript run_all.R --demo           # 用合成数据跑通全链路（推荐第一次用这个）
#   Rscript run_all.R --only 02,03     # 只跑指定步骤
#   Rscript run_all.R --from 04        # 从第 04 步开始往后跑
#   Rscript run_all.R --config config/my.yml
#
# 为什么用子进程而不是 source()：每步独立内存、独立报错，某一步崩了不会污染
# 后面的全局环境；也方便你单独调试某一步。
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i) || i == length(args)) return(default)
  args[i + 1]
}

demo <- "--demo" %in% args
if (demo) Sys.setenv(GLAUCOMA_CONFIG = "config/config.demo.yml")
user_cfg <- get_arg("--config")
if (!is.null(user_cfg)) Sys.setenv(GLAUCOMA_CONFIG = user_cfg)

source("R/helpers.R")
cfg <- read_config()

scripts <- c(
  "01" = "01_geo_prepare.R",
  "02" = "02_deg_limma.R",
  "03" = "03_wgcna.R",
  "04" = "04_enrich_deconv.R",
  "05" = "05_ml_signature.R",
  "06" = "06_scrna_cellchat.R",
  "07" = "07_gwas_qc.R",
  "08" = "08_coloc_smr.R",
  "09" = "09_twas_focus.R",
  "10" = "10_mr.R",
  # 课题三（干实验版）：不在默认 run.steps 里，用 --only 11,12,13 触发
  "11" = "11_insilico_perturbation.R",
  "12" = "12_mediation_mr.R",
  "13" = "13_robustness_cmap.R",
  "14" = "14_cross_species.R",
  # 成稿阶段：15 结果账本 → 16 图表素材 → 17 数字回填
  # 这三步不在默认 run.steps 里，用 --only 15,16 或 --only 17 单独触发
  "15" = "15_results_ledger.R",
  "16" = "16_manuscript_assets.R",
  "17" = "17_fill_tokens.R"
)

steps <- cfg_get(cfg, "run.steps", names(scripts))
only <- get_arg("--only")
if (!is.null(only)) steps <- trimws(strsplit(only, ",")[[1]])
from <- get_arg("--from")
if (!is.null(from)) steps <- steps[as.integer(steps) >= as.integer(from)]
steps <- intersect(steps, names(scripts))

log_init("run_all", cfg)
log_msg("配置：", Sys.getenv("GLAUCOMA_CONFIG", "config/config.yml"))
log_msg("执行步骤：", paste(steps, collapse = ", "))

timing <- data.frame(step = character(0), script = character(0),
                     status = character(0), seconds = numeric(0))

# R 不一定在 PATH 里（例如装在 G:\R\R-4.5.2），必须用当前 R 自己的 Rscript，
# 否则用户双击/在别的终端里跑 run_all 会全部 FAIL(127)
rscript_bin <- file.path(R.home("bin"),
                         if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!file.exists(rscript_bin)) rscript_bin <- "Rscript"
log_msg("使用 Rscript：", rscript_bin)

for (s in steps) {
  f <- scripts[[s]]
  if (!file.exists(f)) { log_msg("脚本缺失，跳过：", f); next }
  log_msg("───── [", s, "] ", f, " ─────")

  # 17 默认把 filled HTML 写在模板旁边；投稿交付要落到 <outdir>/submission 下
  extra <- character(0)
  if (s == "17") {
    src_html <- proj_path("..", "output", "manuscript", "formatted-manuscript.html")
    if (file.exists(src_html)) {
      dst_dir <- ensure_dir(cfg_get(cfg, "project.outdir", "results"), "submission")
      extra <- c("--in", src_html,
                 "--out", file.path(dst_dir, "formatted-manuscript_filled.html"))
    } else {
      log_msg("[17] 未找到稿件模板，将使用脚本默认路径：", src_html)
    }
  }

  t0 <- Sys.time()
  status <- tryCatch({
    code <- system2(rscript_bin, c("--vanilla", shQuote(c(f, extra))),
                    stdout = "", stderr = "")
    if (code == 0) "OK" else paste0("FAIL(", code, ")")
  }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timing <- rbind(timing, data.frame(step = s, script = f, status = status, seconds = round(dt, 1)))
  log_msg("[", s, "] ", status, "，耗时 ", round(dt, 1), " 秒")
  if (!identical(status, "OK")) {
    log_msg("步骤 ", s, " 失败，中止。单步调试：Rscript --vanilla ", f)
    print(timing)
    quit(status = 1)
  }
}

print(timing)
utils::write.table(timing, proj_path(cfg_get(cfg, "project.log_dir", "logs"), "run_all_timing.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("全部完成，总耗时 ", round(sum(timing$seconds), 1), " 秒")

source("R/helpers.R")
for (f in c("config/config.demo.yml", "config/config.yml")) {
  cat("=== ", f, " ===\n", sep = "")
  cfg <- read_config(f)
  cat("project.name : ", as.character(cfg$project$name %||% NA), "\n")
  st <- cfg$run$steps
  cat("run.steps    : ", paste(st, collapse = ","), " (class=", class(st), ", len=", length(st), ")\n", sep = "")
  cat("run.demo     : ", as.character(cfg$run$demo %||% NA), "\n")
  cat("deg.fdr      : ", as.character(cfg$bulk$deg$fdr %||% NA), " (class=", class(cfg$bulk$deg$fdr), ")\n", sep = "")
  cat("coloc.p12    : ", as.character(cfg$genetics$coloc$p12 %||% NA), "\n")
  cat("gwas.n       : ", as.character(cfg$genetics$gwas$n %||% NA), "\n")
}

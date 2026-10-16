## Runs the real 05_fit.R and 06_combine.R scripts end to end on the synthetic fixture (make fixture).
## Exits with status 1 unless every planted program lands in its expected group.
root <- normalizePath(file.path("analysis", "hcmi_persistence", "fixture", "out"), mustWork = FALSE)
unlink(root, recursive = TRUE)
Sys.setenv(GRACKLE_DATA = file.path(root, "data"), GRACKLE_SCRATCH = file.path(root, "scratch"), GRACKLE_RUN_ID = "fixture")
source("analysis/hcmi_persistence/config.R")
dir.create(file.path(cfg$processed, "gene_graphs"), recursive = TRUE, showWarnings = FALSE)
fx <- simulatePersistenceFixture()
saveRDS(fx$Y, file.path(cfg$processed, "Y.rds"))
saveRDS(fx$samples, file.path(cfg$processed, "samples_included.rds"))
saveRDS(fx$k, file.path(cfg$processed, "k.rds"))
saveRDS(fx$gene_graph, file.path(cfg$processed, "gene_graphs", "coess_2d.rds"))
saveRDS(list(BLOCK_P1 = colnames(fx$Y)[fx$blocks$P1]), file.path(cfg$processed, "hallmark.rds"))
grid <- expand.grid(alpha_2 = c(0, 0.1), alpha_1 = c(0, 0.1))[, c("alpha_1", "alpha_2")]
grid$setting_id <- sprintf("a%02d", seq_len(nrow(grid)))
t <- merge(grid, data.frame(fold = 0:5))
t <- t[order(t$setting_id, t$fold), ]
tasks <- data.frame(task_id = sprintf("coess_2d_%s_f%d", t$setting_id, t$fold), arm = "coess_2d", mode = "main",
                    setting_id = t$setting_id, alpha_1 = t$alpha_1, alpha_2 = t$alpha_2, fold = t$fold,
                    k_mult = 1, variant = NA)
tasks_file <- file.path(root, "tasks_main_coess_2d.csv")
utils::write.csv(tasks, tasks_file, row.names = FALSE, na = "")
rscript <- file.path(R.home("bin"), "Rscript")
status <- vapply(seq_len(nrow(tasks)), function(i)
  system2(rscript, c("analysis/hcmi_persistence/05_fit.R", "--tasks", tasks_file, "--index", i)), integer(1))
if (any(status != 0)) stop("05_fit.R failed for tasks: ", paste(which(status != 0), collapse = ", "))
if (system2(rscript, c("analysis/hcmi_persistence/06_combine.R", "--arm", "coess_2d", "--stage", "main",
                       "--tasks", tasks_file, "--n_perm", "999", "--n_boot", "200", "--n_perm_match", "200")) != 0) {
  stop("06_combine.R failed")
}
out <- readRDS(file.path(cfg$run, "results", "coess_2d", "combined.rds"))
calls <- fixtureCalls(out, fx)
print(rbind(expected = fx$expected, observed = calls[names(fx$expected)]))
if (!identical(unname(calls[names(fx$expected)]), unname(fx$expected))) quit(save = "no", status = 1)
cat("fixture PASS\n")

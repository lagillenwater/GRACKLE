## Writes the main task tables: 36 alpha settings x (reference + 5 folds) for each arm (216 rows each).
arms <- c("coess_2d", "coess_3d", "string_mech", "string_all")
alphas <- c(0, 0.01, 0.03, 0.1, 0.3, 1)
grid <- expand.grid(alpha_2 = alphas, alpha_1 = alphas)[, c("alpha_1", "alpha_2")]
grid$setting_id <- sprintf("a%02d", seq_len(nrow(grid)))
dir.create("analysis/hcmi_persistence/tasks", showWarnings = FALSE)
invisible(lapply(arms, function(arm) {
  t <- merge(grid, data.frame(fold = 0:5))
  t <- t[order(t$setting_id, t$fold), ]
  out <- data.frame(task_id = sprintf("%s_%s_f%d", arm, t$setting_id, t$fold), arm = arm, mode = "main",
                    setting_id = t$setting_id, alpha_1 = t$alpha_1, alpha_2 = t$alpha_2, fold = t$fold,
                    k_mult = 1, variant = NA)
  utils::write.csv(out, sprintf("analysis/hcmi_persistence/tasks/main_%s.csv", arm), row.names = FALSE, na = "")
}))

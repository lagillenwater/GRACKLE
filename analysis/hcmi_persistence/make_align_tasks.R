## Writes the pair-alignment task table: sample graph (pair-only or combined) x alpha_1 x alpha_2 on the
## coess_2d gene graph, with sex-fingerprint genes dropped, each as a reference fit plus 5 folds.
## alpha_1 = 0 has no sample graph, so it appears once ("none") for each alpha_2.
a1 <- c(0.01, 0.03, 0.1, 0.3, 1)
a2 <- c(0, 0.1, 1)
grid <- rbind(data.frame(sample_graph = "none", alpha_1 = 0, alpha_2 = a2),
              expand.grid(sample_graph = c("pair", "combined"), alpha_1 = a1, alpha_2 = a2, stringsAsFactors = FALSE))
grid$setting_id <- sprintf("%s-%g-%g", grid$sample_graph, grid$alpha_1, grid$alpha_2)
t <- merge(grid, data.frame(fold = 0:5))
t <- t[order(t$setting_id, t$fold), ]
out <- data.frame(task_id = sprintf("coess_2d_%s_f%d", t$setting_id, t$fold), arm = "coess_2d", mode = "align",
                  setting_id = t$setting_id, alpha_1 = t$alpha_1, alpha_2 = t$alpha_2, fold = t$fold, k_mult = 1,
                  variant = NA, sample_graph = t$sample_graph, drop_sex = TRUE)
utils::write.csv(out, "analysis/hcmi_persistence/tasks/align_coess_2d.csv", row.names = FALSE, na = "")

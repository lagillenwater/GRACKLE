## Writes the pair-alignment task tables: sample graph (pair-only or combined) x alpha_1 x alpha_2 for one
## gene-graph arm, with the sex-fingerprint genes dropped, each as a reference fit plus 5 folds.
## alpha_1 = 0 has no sample graph, so it appears once ("none") for each alpha_2.
alignTasks <- function(arm, a1, a2, graphs = c("pair", "combined"), with_none = TRUE, name = arm) {
  grid <- expand.grid(sample_graph = graphs, alpha_1 = a1, alpha_2 = a2, stringsAsFactors = FALSE)
  if (with_none) grid <- rbind(data.frame(sample_graph = "none", alpha_1 = 0, alpha_2 = a2), grid)
  grid$setting_id <- sprintf("%s-%g-%g", grid$sample_graph, grid$alpha_1, grid$alpha_2)
  t <- merge(grid, data.frame(fold = 0:5))
  t <- t[order(t$setting_id, t$fold), ]
  out <- data.frame(task_id = sprintf("%s_%s_f%d", arm, t$setting_id, t$fold), arm = arm, mode = "align",
                    setting_id = t$setting_id, alpha_1 = t$alpha_1, alpha_2 = t$alpha_2, fold = t$fold, k_mult = 1,
                    variant = NA, sample_graph = t$sample_graph, drop_sex = TRUE)
  utils::write.csv(out, sprintf("analysis/hcmi_persistence/tasks/align_%s.csv", name), row.names = FALSE, na = "")
}
## coess_2d (run 20261022_075f8ae): the full grid, including the plain-NMF baseline (alpha_1 = alpha_2 = 0).
alignTasks("coess_2d", c(0.01, 0.03, 0.1, 0.3, 1), c(0, 0.1, 1))
## The pair-only graph still paired better at alpha_1 = 1, the top of the coess_2d grid: extend it.
alignTasks("coess_2d", c(3, 10), 0, graphs = "pair", with_none = FALSE, name = "coess_2d_ext")
## STRING experiments graph (03c), scored by GSEA: alone and with the pair-only sample graph at the
## strengths where the positive control still held.
alignTasks("string_exp", c(0.1, 0.3), c(0.1, 0.3), graphs = "pair")

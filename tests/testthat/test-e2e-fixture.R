test_that("the fixture has the planted structure", {
  fx <- simulatePersistenceFixture()
  expect_equal(dim(fx$Y), c(75 + 60, 500))
  expect_equal(sum(fx$samples$role == "organoid"), 60)
  expect_true(all(fx$Y >= 0 & fx$Y <= 1))
  expect_equal(nrow(pairTable(fx$samples)), 60)
})

test_that("fixture programs land in their expected persistence groups", {
  skip_if_not(Sys.getenv("GRACKLE_RUN_E2E") == "true", "set GRACKLE_RUN_E2E=true (make e2e) to run")
  skip_if_no_tf()
  fx <- simulatePersistenceFixture()
  grid <- expand.grid(alpha_2 = c(0, 0.1), alpha_1 = c(0, 0.1))[, c("alpha_1", "alpha_2")]
  grid$setting_id <- sprintf("a%02d", seq_len(nrow(grid)))
  tasks <- merge(grid, data.frame(fold = 0:5))
  results <- lapply(seq_len(nrow(tasks)), function(i) {
    t <- tasks[i, ]
    list(meta = list(task_id = paste(t$setting_id, t$fold), arm = "fixture", mode = "main",
                     setting_id = t$setting_id, alpha_1 = t$alpha_1, alpha_2 = t$alpha_2, fold = t$fold),
         fit = fitTask(fx$Y, fx$samples, fx$gene_graph, fx$k, t$alpha_1, t$alpha_2, t$fold),
         error = NULL)
  })
  out <- combineArm(results, grid, fx$samples, pairTable(fx$samples),
                    gene_sets = list(BLOCK_P1 = colnames(fx$Y)[fx$blocks$P1]),
                    n_perm = 999, n_boot = 200, n_perm_match = 200)
  expect_equal(fixtureCalls(out, fx), fx$expected)
})

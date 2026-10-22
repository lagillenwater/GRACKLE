align_samples <- function() {
  data.frame(profile_id = c(paste0("t", 1:6), paste0("o", 1:3)), case_id = c(paste0("c", 1:6), paste0("c", 1:3)),
             sample_id = c(paste0("t", 1:6), paste0("o", 1:3)), parent_sample_id = c(rep(NA, 6), paste0("t", 1:3)),
             role = rep(c("tumor", "organoid"), c(6, 3)), cancer_type = "X", stringsAsFactors = FALSE)
}

test_that("pairingScore is 1 when each organoid's program profile matches its own tumor's", {
  s <- align_samples()
  W <- withr::with_seed(3, matrix(stats::rexp(6 * 4), 6, 4, dimnames = list(s$profile_id[1:6], NULL)))
  W <- rbind(W, W[1:3, ] * 1.0001)
  rownames(W) <- s$profile_id
  r <- pairingScore(W, s)
  expect_equal(nrow(r), 3)
  expect_equal(r$score, rep(1, 3))
})

test_that("pairingScore gives an organoid at its own tumor's far end of the tree a low score", {
  s <- align_samples()
  W <- rbind(t1 = c(10, 0), t2 = c(9, 1), t3 = c(8, 2), t4 = c(0, 10), t5 = c(1, 9), t6 = c(2, 8),
             o1 = c(0, 10), o2 = c(9, 1), o3 = c(8, 2))
  r <- pairingScore(W, s)
  expect_lt(r$score[r$organoid == "o1"], 0.5)
  expect_gt(r$score[r$organoid == "o2"], 0.5)
})

test_that("enrichedShare counts the programs with a significant gene set", {
  genes <- paste0("g", 1:200)
  H <- rbind(LV1 = c(rep(10, 10), rep(0, 190)), LV2 = c(rep(0, 190), rep(10, 10)))
  colnames(H) <- genes
  sets <- list(A = genes[1:20])
  expect_equal(enrichedShare(H, sets, method = "ora"), 0.5)
})

## Smoothed programs spread their loadings, which a top-5% overlap test penalizes; GSEA uses every gene's rank.
test_that("enrichedShare with GSEA finds the program whose highest loadings sit on a gene set", {
  skip_if_not_installed("fgsea")
  genes <- paste0("g", 1:200)
  rest <- withr::with_seed(5, stats::runif(180, 0, 0.5))
  H <- rbind(LV1 = c(seq(10, 5, length.out = 20), rest), LV2 = c(rest[1:20], rest[21:180], seq(10, 5, length.out = 20))[1:200])
  colnames(H) <- genes
  sets <- list(A = genes[1:20], B = genes[101:120])
  g <- gseaPrograms(H, sets)
  expect_named(g, c("program", "gene_set", "nes", "p", "fdr"))
  expect_lt(g$fdr[g$program == "LV1" & g$gene_set == "A"], 0.05)
  expect_equal(enrichedShare(H, sets, method = "gsea"), 0.5)
})

test_that("balanceSelect keeps settings that beat the baseline's pairing and picks the one nearest the ideal", {
  tab <- data.frame(setting_id = c("base", "A", "B", "C", "D"), pairing = c(0.6, 0.7, 0.65, 0.62, 0.55),
                    enrichment = c(0.5, 0.4, 0.6, 0.55, 0.9), baseline = c(TRUE, FALSE, FALSE, FALSE, FALSE),
                    eligible = TRUE)
  out <- balanceSelect(tab)
  expect_equal(out$floor_ok, c(FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(out$pareto, c(FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(out$setting_id[out$selected], "B")
})

test_that("balanceSelect selects nothing when no setting beats the baseline", {
  tab <- data.frame(setting_id = c("base", "A"), pairing = c(0.6, 0.5), enrichment = c(0.5, 0.9),
                    baseline = c(TRUE, FALSE), eligible = TRUE)
  expect_false(any(balanceSelect(tab)$selected))
})

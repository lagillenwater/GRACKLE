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

## A strong pair graph raised pairing by folding tumor-only programs into shared ones (pair-only alpha_1 >= 1:
## microenvironment programs no longer dropped out more often than tumor-intrinsic ones).
test_that("controlHolds is TRUE only when microenvironment programs are lost more often than tumor-intrinsic ones", {
  p <- data.frame(program = paste0("LV", 1:6), stable = TRUE, microenv = rep(c(TRUE, FALSE), each = 3),
                  group = c("drops_out", "drops_out", "persists", "persists", "persists", "reduced_concordant"))
  expect_true(controlHolds(p))                        # 2 of 3 against 1 of 3
  p$group <- c("drops_out", "persists", "persists", "drops_out", "reduced_concordant", "persists")
  expect_false(controlHolds(p))                       # 1 of 3 against 2 of 3
  p$stable[4:5] <- FALSE
  expect_true(controlHolds(p))                        # unstable programs do not count: 1 of 3 against 0 of 1
})

## The share of enriched programs saturates under GSEA (95-100% in every setting), so it cannot rank settings.
test_that("enrichmentStrength is the median over programs of each program's best NES", {
  g <- data.frame(program = c("LV1", "LV1", "LV2", "LV2", "LV3"), gene_set = c("A", "B", "A", "B", "A"),
                  nes = c(1, 2, 1.5, 0.5, 3), p = 0.01, fdr = 0.02)
  expect_equal(enrichmentStrength(g), 2)              # best per program 2, 1.5, 3
})

## Programs that transfer in a sample cluster are read off the cluster's tumor and organoid means, with
## each program standardized over all samples so a program every organoid loses reads low in organoids.
test_that("clusterActivity averages standardized activity by cluster and role and marks transfers and losses", {
  A <- cbind(P1 = c(9, 9, 0, 0, 9, 9, 0, 0), P2 = c(0, 0, 9, 9, 0, 0, 0, 0), P3 = 1)
  role <- rep(c("tumor", "organoid"), each = 4)
  cl <- c(1, 1, 2, 2, 1, 1, 2, 2)
  out <- clusterActivity(A, role, cl)
  expect_named(out, c("cluster", "program", "tumor_z", "organoid_z", "status"))
  Z <- scale(log1p(A))
  expect_equal(out$tumor_z[out$cluster == "2" & out$program == "P2"], mean(Z[3:4, "P2"]))
  expect_equal(out$organoid_z[out$cluster == "1" & out$program == "P1"], mean(Z[5:6, "P1"]))
  expect_equal(out$status[out$program == "P1"], c("transfers", ""))
  expect_equal(out$status[out$program == "P2"], c("", "lost"))
  expect_equal(out$tumor_z[out$program == "P3"], c(0, 0))   # no variance reads as 0, not NA
  expect_equal(unique(clusterActivity(A, role, c(1, 1, 2, 2, 1, 1, 1, 1))$cluster), "1")  # cluster 2 has no organoids
})

## The sample tree in steps 08 and 09 standardizes each program within tumors and within organoids, which
## removes the shift every organoid shares (its lost microenvironment programs).
test_that("standardizeActivity scales log1p activity over all samples or within each role, keeping row order", {
  A <- cbind(P1 = c(1, 3, 5, 10, 20, 40), P2 = c(2, 2, 2, 4, 8, 6))
  role <- c("tumor", "organoid", "tumor", "organoid", "tumor", "organoid")
  all <- standardizeActivity(A)
  expect_equal(all, scale(log1p(A)), ignore_attr = TRUE)
  w <- standardizeActivity(A, role)
  expect_equal(dimnames(w), dimnames(A))
  expect_equal(w[c(1, 3, 5), ], scale(log1p(A[c(1, 3, 5), ])), ignore_attr = TRUE)
  expect_equal(w[c(2, 4, 6), ], scale(log1p(A[c(2, 4, 6), ])), ignore_attr = TRUE)
  A[1, 1] <- NA
  expect_false(anyNA(standardizeActivity(A, role)))    # unmatched in a sample's fold reads as 0
  expect_equal(standardizeActivity(A, role)[1, "P1"], 0, ignore_attr = TRUE)
})

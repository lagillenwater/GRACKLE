#' simulatePersistenceFixture
#'
#' Synthetic tumor-organoid data with four planted programs (spec Section 8.5):
#' - P1 and P2 are patient-specific and shared by tumor and organoid (expected: persists).
#' - P3 is tumor-only microenvironment, scaled by 1 - purity (expected: drops_out).
#' - P4 is a lineage constant within cancer type, plus noise (expected: present_not_patient_specific).
#' A dense baseline profile is the fifth component. The gene graph links genes within each
#' block and adds random edges.
#' @return A list with Y (min-max scaled), samples (with folds), gene_graph, blocks, k, expected.
#' @export
simulatePersistenceFixture <- function(n_pairs_per_type = 20, n_types = 3, n_unpaired = 15, n_genes = 500,
                                       noise_sd = 0.2, seed = 20261015) {
  fx <- withr::with_seed(seed, {
    types <- c("COAD", "PAAD", "BRCA", "SKCM", "GBM")[seq_len(n_types)]
    n_paired <- n_pairs_per_type * n_types
    n_cases <- n_paired + n_unpaired
    case_id <- sprintf("HCM-TEST-%04d-C00", seq_len(n_cases))
    cancer_type <- c(rep(types, each = n_pairs_per_type), rep(types, length.out = n_unpaired))
    purity <- stats::runif(n_cases, 0.3, 0.9)
    a1 <- stats::rgamma(n_cases, 2, 1)
    a2 <- stats::rgamma(n_cases, 2, 1)
    lineage <- c(0.5, 1.5, 3, 2, 1)[match(cancer_type, types)]
    tumors <- data.frame(profile_id = paste0(case_id, "-01A"), case_id = case_id, sample_id = paste0(case_id, "-01A"),
                         parent_sample_id = NA_character_, role = "tumor", cancer_type = cancer_type,
                         tumor_descriptor = "primary", prior_treatment = stats::runif(n_cases) < 0.4,
                         age = round(stats::runif(n_cases, 30, 80)), sex = sample(c("F", "M"), n_cases, TRUE),
                         ancestry = "EUR", purity_dna = purity, stringsAsFactors = FALSE)
    p <- seq_len(n_paired)
    orgs <- tumors[p, ]
    orgs$profile_id <- paste0(case_id[p], "-85A")
    orgs$sample_id <- orgs$profile_id
    orgs$parent_sample_id <- paste0(case_id[p], "-01A")
    orgs$role <- "organoid"
    orgs$purity_dna <- NA_real_
    samples <- rbind(tumors, orgs)
    n <- nrow(samples)
    ci <- match(samples$case_id, case_id)
    is_t <- samples$role == "tumor"
    jitter <- function() stats::rnorm(n, sd = 0.1)
    act <- cbind(P1 = pmax(a1[ci] + jitter(), 0), P2 = pmax(a2[ci] + jitter(), 0),
                 P3 = ifelse(is_t, 2 * (1 - purity[ci]), 0), P4 = pmax(lineage[ci] + jitter(), 0), base = 1)
    blocks <- list(P1 = 1:50, P2 = 51:100, P3 = 101:150, P4 = 151:200)
    genes <- sprintf("G%03d", seq_len(n_genes))
    H_true <- matrix(0, 5, n_genes, dimnames = list(colnames(act), genes))
    H_true[cbind(rep(1:4, each = 50), unlist(blocks))] <- 1
    H_true["base", ] <- stats::runif(n_genes, 0.5, 1.5)
    Y <- pmax(act %*% H_true + matrix(stats::rnorm(n * n_genes, sd = noise_sd), n), 0)
    dimnames(Y) <- list(samples$profile_id, genes)
    G <- matrix(0, n_genes, n_genes, dimnames = list(genes, genes))
    G[do.call(rbind, lapply(blocks, function(b) cbind(b, sample(b))))] <- 1
    G[matrix(sample.int(n_genes, 2 * n_genes, TRUE), ncol = 2)] <- 1
    list(Y = Y, samples = samples, G = G, blocks = blocks)
  })
  Y <- min_max_scale(fx$Y, apply(fx$Y, 2, min), apply(fx$Y, 2, max))
  G <- pmax(fx$G, t(fx$G))
  diag(G) <- 0
  list(Y = Y, samples = makeFolds(fx$samples, seed = seed), gene_graph = G, blocks = fx$blocks, k = 5L,
       expected = c(P1 = "persists", P2 = "persists", P3 = "drops_out", P4 = "present_not_patient_specific"))
}

#' fixtureCalls
#'
#' Maps each planted block to the reference program with the highest mean loading on it and
#' returns that program's group (NA if the program is not stable).
#' @param combined Output of [combineArm()].
#' @param fixture Output of [simulatePersistenceFixture()].
#' @return Named character vector of groups for P1..P4.
#' @export
fixtureCalls <- function(combined, fixture) {
  H <- combined$ref_H
  prog <- vapply(fixture$blocks, function(b) rownames(H)[which.max(rowMeans(H[, b, drop = FALSE]))], character(1))
  stats::setNames(combined$programs$group[match(prog, combined$programs$program)], names(fixture$blocks))
}

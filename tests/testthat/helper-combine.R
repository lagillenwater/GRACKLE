## 60 pairs in 2 cancer types and 5 folds, plus 10 unpaired tumors. Three true programs:
## LV1 shared by tumor and organoid, LV2 absent from organoids, LV3 independent.
fake_world <- function(seed = 14, k = 3, g = 60) {
  withr::with_seed(seed, {
    n_pairs <- 60
    case <- sprintf("c%03d", seq_len(n_pairs + 10))
    tumors <- data.frame(profile_id = paste0(case, "-T"), case_id = case, sample_id = paste0(case, "-T"),
                         parent_sample_id = NA_character_, role = "tumor",
                         cancer_type = rep(c("COAD", "PAAD"), length.out = length(case)),
                         purity_dna = stats::runif(length(case), 0.3, 0.9), stringsAsFactors = FALSE)
    orgs <- data.frame(profile_id = paste0(case[1:n_pairs], "-O"), case_id = case[1:n_pairs],
                       sample_id = paste0(case[1:n_pairs], "-O"), parent_sample_id = paste0(case[1:n_pairs], "-T"),
                       role = "organoid", cancer_type = tumors$cancer_type[1:n_pairs], purity_dna = NA_real_,
                       stringsAsFactors = FALSE)
    samples <- makeFolds(rbind(tumors, orgs), seed = seed)
    a <- stats::rgamma(length(case), 2, 1)
    ci <- match(samples$case_id, case)
    is_t <- samples$role == "tumor"
    act <- cbind(a[ci] + stats::rnorm(nrow(samples), sd = 0.1),
                 ifelse(is_t, 1 + stats::rnorm(nrow(samples), sd = 0.2), 0),
                 stats::rgamma(nrow(samples), 2, 1))
    act <- pmax(act, 0)
    rownames(act) <- samples$profile_id
    H_ref <- normalizeRows(matrix(stats::rexp(k * g), k, dimnames = list(paste0("LV", 1:k), paste0("g", 1:g))))
  })
  list(samples = samples, act = act, H_ref = H_ref, pairs = pairTable(samples))
}

fake_results <- function(world, setting_id, alpha_1, noise, recon, converged = rep(TRUE, 6),
                         missing_folds = integer(0), seed = 1) {
  withr::with_seed(seed, lapply(setdiff(0:5, missing_folds), function(f) {
    perm <- if (f == 0) seq_len(nrow(world$H_ref)) else sample(nrow(world$H_ref))
    H <- normalizeRows(world$H_ref[perm, ] + matrix(abs(stats::rnorm(length(world$H_ref), sd = noise)), nrow(world$H_ref)))
    rownames(H) <- paste0("LV", seq_len(nrow(H)))
    held <- world$samples$profile_id[world$samples$fold == f]
    W_test <- if (f == 0) NULL else world$act[held, perm, drop = FALSE]
    if (!is.null(W_test)) colnames(W_test) <- rownames(H)
    list(meta = list(task_id = sprintf("x_%s_f%d", setting_id, f), arm = "x", mode = "main", setting_id = setting_id,
                     alpha_1 = alpha_1, alpha_2 = 0, fold = f, k_mult = 1, variant = NA, k = nrow(H)),
         fit = list(H = H, W_test = W_test, converged = converged[f + 1], n_iter = 10L,
                    final_loss = data.frame(total = 1), recon_test = if (f == 0) NA_real_ else recon[f]),
         error = NULL)
  }))
}

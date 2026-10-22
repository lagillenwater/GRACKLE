## Step 05: one cross-fitting task (spec Section 8.1); one Slurm array index per task.
##   Rscript analysis/hcmi_persistence/05_fit.R --tasks <tasks.csv> --index <1-based row> [--iterations <cap>]
## A task whose result exists without an error is skipped, so failed indices can be resubmitted alone.
## With --iterations, a result that stopped unconverged below that cap is refit (needsFit()).
source("analysis/hcmi_persistence/config.R")
args <- parseArgs(list(tasks = "", index = as.character(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")) + 1L),
                       iterations = ""))
cap <- if (nzchar(args$iterations)) list(iterations = as.integer(args$iterations)) else list()
tasks <- as.data.frame(fread(args$tasks))
t <- tasks[as.integer(args$index), , drop = FALSE]
if (!nrow(t) || is.na(t$task_id)) stop("no task at index ", args$index, " of ", args$tasks)
out <- file.path(cfg$run, "fits", t$arm, t$mode, paste0(t$task_id, ".rds"))
if (!needsFit(readResult(out), if (length(cap)) cap$iterations else NA)) {
  message("already done: ", out)
  quit(save = "no", status = 0)
}
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
Y <- readRDS(file.path(cfg$processed, "Y.rds"))
samples <- readRDS(file.path(cfg$processed, "samples_included.rds"))
samples <- samples[match(rownames(Y), samples$profile_id), , drop = FALSE]
k <- as.integer(round(readRDS(file.path(cfg$processed, "k.rds")) * t$k_mult))
## Optional task columns (pair-alignment tables): drop sex-fingerprint genes; pair-only sample graph.
hgnc_file <- file.path(cfg$raw, "annotation", "hgnc_complete_set_2024-10-01.txt")
drop_sex <- isTRUE(as.logical(t$drop_sex))
if (drop_sex) Y <- Y[, !colnames(Y) %in% sexFingerprintGenes(as.data.frame(fread(hgnc_file, quote = ""))), drop = FALSE]
weights <- if (identical(t$sample_graph, "pair")) c(pair = 1, annot = 0, purity = 0) else c(pair = 1, annot = 1, purity = 1)
graph_file <- file.path(cfg$processed, "gene_graphs",
                        if (t$arm == "coess_3d") sprintf("coess_3d_fold%d.rds", t$fold) else paste0(t$arm, ".rds"))
G <- readRDS(graph_file)[colnames(Y), colnames(Y)]
transfer <- t$mode == "transfer"
meta <- c(as.list(t), k = k, cap)
res <- tryCatch(
  list(fit = do.call(fitTask, c(list(Y, samples, G, k, t$alpha_1, t$alpha_2, t$fold,
                                     train_role = if (transfer) "tumor" else "all",
                                     test_role = if (transfer) "organoid" else "all",
                                     graph_weights = weights), cap)),
       error = NULL),
  error = function(e) list(fit = NULL, error = conditionMessage(e)))
saveResult(c(list(meta = meta), res), out)
writeProvenance(out, inputs = c(file.path(cfg$processed, c("Y.rds", "samples_included.rds", "k.rds")), graph_file,
                                if (drop_sex) hgnc_file),
                extra = list(task = meta))
if (!is.null(res$error)) {
  message("fit failed: ", res$error)
  quit(save = "no", status = 1)
}
message(sprintf("%s: n_iter %d, converged %s", t$task_id, res$fit$n_iter, res$fit$converged))

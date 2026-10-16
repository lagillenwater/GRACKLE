## Step 06: combine one arm (--stage main | sensitivity) or compare arms (--stage compare).
## Spec Sections 7.3, 7.4, 8.2-8.4, 9 and 10.
source("analysis/hcmi_persistence/config.R")
args <- parseArgs(list(arm = "", stage = "main", tasks = "", n_perm = "10000", n_boot = "2000",
                       n_perm_match = "1000", cores = Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
np <- as.integer(args$n_perm); nb <- as.integer(args$n_boot); nm <- as.integer(args$n_perm_match)
cores <- as.integer(args$cores)
res_dir <- file.path(cfg$run, "results")
## Results for the task ids in `ids` only, so files left over from an earlier selection are never mixed in.
readFits <- function(arm, mode, ids) {
  fs <- file.path(cfg$run, "fits", arm, mode, paste0(ids, ".rds"))
  lapply(fs[file.exists(fs)], readRDS)
}
optional <- function(path, reader) if (file.exists(path)) reader(path) else NULL
samples <- readRDS(file.path(cfg$processed, "samples_included.rds"))
pairs <- pairTable(samples)
gene_sets <- optional(file.path(cfg$raw, "annotation", "hallmark.rds"), readRDS)
if (is.null(gene_sets)) gene_sets <- optional(file.path(cfg$processed, "hallmark.rds"), readRDS)
if (is.null(gene_sets)) gene_sets <- list()
mut <- optional(file.path(cfg$processed, "mutations.rds"), readRDS)
est <- optional(file.path(cfg$processed, "estimate.csv"), function(p) as.data.frame(fread(p)))
writeArm <- function(out, d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, file.path(d, "combined.rds"))
  fwrite(out$selection, file.path(d, "selection.csv"))
  fwrite(out$programs, file.path(d, "programs.csv"))
  fwrite(out$convergence, file.path(d, "convergence.csv"))
  if (!is.null(out$maps_long)) fwrite(out$maps_long, file.path(d, "maps_all_settings.csv"))
  fwrite(out$per_type, file.path(d, "per_type.csv"))
  if (!is.null(out$enrichment)) fwrite(out$enrichment, file.path(d, "enrichment.csv"))
  if (!is.null(out$genotype)) fwrite(out$genotype, file.path(d, "genotype.csv"))
  writeProvenance(file.path(d, "combined.rds"), extra = list(selected = out$selected, n_perm = np, n_boot = nb))
}

if (args$stage == "main") {
  tasks_file <- if (nzchar(args$tasks)) args$tasks else file.path(cfg$analysis, "tasks", sprintf("main_%s.csv", args$arm))
  tasks <- as.data.frame(fread(tasks_file))
  expected <- unique(tasks[, c("setting_id", "alpha_1", "alpha_2")])
  out <- combineArm(readFits(args$arm, "main", tasks$task_id), expected, samples, pairs, gene_sets, mut, est,
                    n_perm = np, n_boot = nb, n_perm_match = nm, seed = cfg$seed, cores = cores)
  writeArm(out, file.path(res_dir, args$arm))
  sel <- out$selection[out$selection$selected, ]
  sens <- merge(data.frame(k_mult = c(0.5, 2)), data.frame(fold = 0:5))
  sens <- data.frame(task_id = sprintf("%s_%s_k%s_f%d", args$arm, sel$setting_id, sens$k_mult, sens$fold),
                     arm = args$arm, mode = "sensitivity", setting_id = paste0(sel$setting_id, "_k", sens$k_mult),
                     alpha_1 = sel$alpha_1, alpha_2 = sel$alpha_2, fold = sens$fold, k_mult = sens$k_mult, variant = NA)
  tr <- merge(data.frame(variant = c("selected", "baseline")), data.frame(fold = 1:5))
  tr <- data.frame(task_id = sprintf("%s_%s_transfer_%s_f%d", args$arm, sel$setting_id, tr$variant, tr$fold),
                   arm = args$arm, mode = "transfer", setting_id = paste0(sel$setting_id, "_", tr$variant),
                   alpha_1 = sel$alpha_1, alpha_2 = ifelse(tr$variant == "baseline", 0, sel$alpha_2),
                   fold = tr$fold, k_mult = 1, variant = tr$variant)
  td <- file.path(cfg$run, "tasks"); dir.create(td, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(sens, file.path(td, sprintf("sensitivity_%s.csv", args$arm)), row.names = FALSE, na = "")
  utils::write.csv(tr, file.path(td, sprintf("transfer_%s.csv", args$arm)), row.names = FALSE, na = "")
  message(args$arm, ": selected ", out$selected, "; stable programs ", sum(out$programs$stable))
}

if (args$stage == "sensitivity") {
  tasks <- as.data.frame(fread(file.path(cfg$run, "tasks", sprintf("sensitivity_%s.csv", args$arm))))
  results <- readFits(args$arm, "sensitivity", tasks$task_id)
  invisible(lapply(unique(tasks$k_mult), function(km) {
    sub <- Filter(function(r) r$meta$k_mult == km, results)
    expected <- unique(tasks[tasks$k_mult == km, c("setting_id", "alpha_1", "alpha_2")])
    d <- file.path(res_dir, args$arm, sprintf("sensitivity_k%s", km))
    unlink(d, recursive = TRUE)   # one outcome per run: an earlier ERROR.txt or result never survives
    out <- tryCatch(combineArm(sub, expected, samples, pairs, gene_sets, mut, est, n_perm = np, n_boot = nb,
                               n_perm_match = nm, seed = cfg$seed, cores = cores),
                    error = function(e) { dir.create(d, recursive = TRUE, showWarnings = FALSE)
                      writeLines(conditionMessage(e), file.path(d, "ERROR.txt")); NULL })
    if (!is.null(out)) writeArm(out, d)
  }))
}

if (args$stage == "compare") {
  d <- file.path(res_dir, "compare"); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  transfer <- unlist(lapply(cfg$arms, function(a) {
    tt <- as.data.frame(fread(file.path(cfg$run, "tasks", sprintf("transfer_%s.csv", a))))
    readFits(a, "transfer", tt$task_id)
  }), recursive = FALSE)
  fwrite(transferTable(transfer), file.path(d, "transfer.csv"))
  a <- readRDS(file.path(res_dir, "coess_2d", "combined.rds"))
  b <- readRDS(file.path(res_dir, "coess_3d", "combined.rds"))
  ag <- mapAgreement(a$ref_H, b$ref_H, a$programs, b$programs, tau = max(a$tau_selected, b$tau_selected))
  fwrite(ag$summary, file.path(d, "map_agreement.csv"))
  fwrite(ag$pairs, file.path(d, "map_agreement_pairs.csv"))
  structure <- data.frame(cancer_type = samples$cancer_type[match(pairs$tumor, samples$profile_id)],
                          fold = samples$fold[match(pairs$tumor, samples$profile_id)])
  fwrite(concordancePower(structure, seed = cfg$seed), file.path(d, "power_curve.csv"))
  writeProvenance(file.path(d, "map_agreement.csv"))
  print(ag$summary)
}

## Step 03: expression matrix, gene selection, scaling, ESTIMATE scores (spec Section 4.4).
source("analysis/hcmi_persistence/config.R")
samples <- readRDS(file.path(cfg$processed, "samples_included.rds"))
gdc <- as.data.frame(fread(file.path(cfg$raw, "gdc", "gdc_files.tsv")))
hgnc <- as.data.frame(fread(file.path(cfg$raw, "annotation", "hgnc_complete_set_2024-10-01.txt"), quote = ""))
files <- file.path(cfg$raw, "gdc", "star", samples$profile_id, gdc$file_name[match(samples$profile_id, gdc$profile_id)])
names(files) <- samples$profile_id
logexpr <- buildExpressionMatrix(files, hgnc)
genes <- selectGenes(logexpr, n = 5000)
Y <- logexpr[, genes, drop = FALSE]
Y <- min_max_scale(Y, apply(Y, 2, min), apply(Y, 2, max))
saveRDS(logexpr, file.path(cfg$processed, "expression_log2tpm.rds"))
saveRDS(Y, file.path(cfg$processed, "Y.rds"))
writeLines(genes, file.path(cfg$processed, "genes.txt"))
tumor_ids <- samples$profile_id[samples$role == "tumor"]
fwrite(estimateScores(t(logexpr[tumor_ids, , drop = FALSE])), file.path(cfg$processed, "estimate.csv"))
writeProvenance(file.path(cfg$processed, "Y.rds"), inputs = file.path(cfg$processed, "samples_included.rds"),
                extra = list(n_profiles = nrow(Y), n_genes = ncol(Y)))
message("Y: ", nrow(Y), " profiles x ", ncol(Y), " genes")

## Step 02: the sample table, folds, filter counts and the mutation matrix (spec Sections 4.2, 4.3, 4.6).
source("analysis/hcmi_persistence/config.R")
gdc <- as.data.frame(fread(file.path(cfg$raw, "gdc", "gdc_files.tsv")))
long <- fread(file.path(cfg$raw, "cbioportal", "clinical_sample_long.tsv"))
plong <- fread(file.path(cfg$raw, "cbioportal", "clinical_patient_long.tsv"))
first <- function(x) as.character(x)[1L]   # an attribute listed twice for a sample keeps its first value
wide_s <- dcast(long, sampleId + patientId ~ clinicalAttributeId, value.var = "value", fun.aggregate = first, fill = NA_character_)
wide_p <- dcast(plong, patientId ~ clinicalAttributeId, value.var = "value", fun.aggregate = first, fill = NA_character_)
patient_only <- setdiff(names(wide_p), names(wide_s))
wide <- merge(wide_s, wide_p[, c("patientId", patient_only), with = FALSE], by = "patientId", all.x = TRUE)
pick <- function(col) if (col %in% names(wide)) as.character(wide[[col]]) else rep(NA_character_, nrow(wide))
clinical <- data.frame(sample_id = wide$sampleId, growth_pattern = pick("GROWTH_PATTERN"),
                       culture_media = pick("CULTURE_MEDIA"), treatment_status = pick("TREATMENT_STATUS"),
                       cancer_type = pick("DIAGNOSIS_TCGA_HCMI_CODE"), age = pick("AGE"), sex = pick("SEX"),
                       race_category = pick("RACE_CATEGORY"), purity_absolute = pick("PURITY_ABSOULTE"),
                       purity_ascat = pick("PURITY_ASCAT"), purity_hatchet = pick("PURITY_HACHET"),
                       purity_purple = pick("PURITY_PURPLE"), consensus_purity = pick("CONSENSUS_PURITY"),
                       stringsAsFactors = FALSE)
built <- buildSampleTable(gdc, clinical)
s <- built$samples
inc <- makeFolds(s[s$included, , drop = FALSE], seed = cfg$seed)
s$fold <- inc$fold[match(s$profile_id, inc$profile_id)]

tum <- inc$role == "tumor"
cand <- s$sample_type %in% "Next Generation Cancer Model"
coverage <- data.frame(
  field = c("purity_dna, included tumors", "model_type, model samples", "cancer_type, included samples"),
  n = c(sum(tum), sum(cand), nrow(inc)),
  n_available = c(sum(tum & !is.na(inc$purity_dna)), sum(cand & !is.na(s$model_type)), sum(!is.na(inc$cancer_type))))
## cBioPortal's CONSENSUS_PURITY is not a purity estimate here (see spec 4.6); report how well the
## four DNA callers agree with each other instead.
callers <- clinical[match(inc$sample_id[tum], clinical$sample_id),
                    c("purity_absolute", "purity_ascat", "purity_hatchet", "purity_purple")]
cm <- suppressWarnings(stats::cor(vapply(callers, as.numeric, numeric(nrow(callers))), use = "pairwise.complete.obs"))
coverage$dna_callers_median_pairwise_r <- c(stats::median(cm[upper.tri(cm)], na.rm = TRUE), NA, NA)
print(coverage)
fwrite(coverage, file.path(cfg$processed, "field_coverage.csv"))
if (coverage$n_available[1] == 0) stop("no DNA-based tumor purity: stop and ask the user (spec Section 4.6)")

onco <- as.data.frame(fread(file.path(cfg$raw, "annotation", "oncokb_cancerGeneList.txt")))
mut <- as.data.frame(fread(file.path(cfg$raw, "cbioportal", "mutations.tsv")))
seq_ids <- readLines(file.path(cfg$raw, "cbioportal", "sequenced_samples.txt"))
tumors <- inc[tum, , drop = FALSE]
M <- mutationMatrix(mut, tumors$profile_id, tumors$sample_id, unique(onco[["Hugo Symbol"]]), seq_ids)
saveRDS(M, file.path(cfg$processed, "mutations.rds"))

keep_cols <- c("profile_id", "case_id", "case_submitter_id", "sample_id", "sample_type", "role", "cancer_type",
               "tumor_descriptor", "preservation_method", "model_id", "model_type", "medium", "passage",
               "derivation_center", "parent_sample_id", "pair_flag", "purity_dna", "consensus_purity",
               "prior_treatment", "age", "sex", "ancestry", "fold", "included", "exclusion_reason")
fwrite(s[, keep_cols], file.path(cfg$processed, "sample_table.csv"))
fwrite(unique(inc[, c("case_id", "fold")]), file.path(cfg$processed, "folds.csv"))
fwrite(built$filter_counts, file.path(cfg$processed, "filter_counts.csv"))
fwrite(pairTable(inc), file.path(cfg$processed, "pairs.csv"))
saveRDS(inc, file.path(cfg$processed, "samples_included.rds"))
writeProvenance(file.path(cfg$processed, "samples_included.rds"),
                inputs = c(file.path(cfg$raw, "gdc", "gdc_files.tsv"), file.path(cfg$raw, "cbioportal", "clinical_sample_long.tsv")))
print(built$filter_counts)
message("pairs: ", nrow(pairTable(inc)), "; flagged organoids: ", sum(!is.na(inc$pair_flag)))

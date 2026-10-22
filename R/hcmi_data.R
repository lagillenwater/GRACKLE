## HCMI inputs: GDC STAR counts and the sample table (spec Sections 4.2 and 4.4).

.tumor_sample_types <- c("Primary Tumor", "Metastatic", "Recurrent Tumor", "Additional Metastatic",
                         "Human Tumor Original Cells", "Neoplasms of Uncertain and Unknown Behavior",
                         "Slides", "FFPE Scrolls")

#' readStarCounts
#' @param path A GDC STAR - Counts file.
#' @return Data frame with gene_id, gene_name, gene_type, tpm_unstranded (N_* summary rows removed).
#' @export
readStarCounts <- function(path) {
  x <- data.table::fread(path, sep = "\t", skip = "gene_id", data.table = FALSE)
  x[!startsWith(x$gene_id, "N_"), c("gene_id", "gene_name", "gene_type", "tpm_unstranded")]
}

#' buildExpressionMatrix
#'
#' log2(TPM + 1) for protein-coding genes mapped to HGNC symbols through Ensembl ids
#' (version suffix removed; PAR_Y copies dropped). If several Ensembl ids map to one symbol,
#' the one with the highest mean TPM is kept.
#' @param files Named character vector of STAR files (names are profile ids).
#' @param hgnc Data frame with symbol and ensembl_gene_id.
#' @return Profiles x symbols matrix, columns sorted.
#' @export
buildExpressionMatrix <- function(files, hgnc) {
  first <- readStarCounts(files[[1]])
  keep <- first$gene_type == "protein_coding" & !grepl("_PAR_Y$", first$gene_id)
  sym <- hgnc$symbol[match(sub("\\..*$", "", first$gene_id[keep]), hgnc$ensembl_gene_id)]
  tpm <- vapply(files, function(f) {
    x <- readStarCounts(f)
    if (!identical(x$gene_id, first$gene_id)) stop("buildExpressionMatrix: gene order differs in ", f)
    x$tpm_unstranded[keep]
  }, numeric(sum(keep)))
  tpm <- matrix(tpm, nrow = sum(keep))
  ok <- !is.na(sym)
  tpm <- tpm[ok, , drop = FALSE]
  sym <- sym[ok]
  o <- order(sym, -rowMeans(tpm))
  sel <- o[!duplicated(sym[o])]
  out <- log2(t(tpm[sel, , drop = FALSE]) + 1)
  dimnames(out) <- list(names(files), sym[sel])
  out[, order(colnames(out)), drop = FALSE]
}

#' selectGenes
#' @param logexpr Profiles x genes matrix.
#' @param n Number of genes. Default 5000.
#' @return The n genes with the highest variance across profiles.
#' @export
selectGenes <- function(logexpr, n = 5000) {
  v <- matrixStats::colVars(logexpr)
  names(v) <- colnames(logexpr)
  names(sort(v, decreasing = TRUE))[seq_len(min(n, length(v)))]
}

#' buildSampleTable
#'
#' Applies the inclusion rules of spec Section 4.2 in order, recording the counts after each
#' step: expanded models and non-tumor samples out; 3D organoids only; the FFPE rule; one
#' profile per sample; parent-tumor pairing by tumor descriptor; organoids without a tumor
#' partner out. Cancer type is set per case (most frequent code, else the ICD-O site in the HCMI id).
#' Organoids inherit their parent tumor's annotations; tumors get DNA purity as the median of the
#' four DNA callers.
#' @param gdc GDC file metadata: profile_id, case_id, case_submitter_id, sample_submitter_id, sample_type, tumor_descriptor, preservation_method, passage_count.
#' @param clinical cBioPortal fields per sample.
#' @return A list with samples (every profile, with included and exclusion_reason) and filter_counts.
#' @export
buildSampleTable <- function(gdc, clinical) {
  s <- gdc
  s$role <- ifelse(s$sample_type %in% .tumor_sample_types, "tumor",
                   ifelse(s$sample_type %in% "Next Generation Cancer Model", "organoid", NA_character_))
  s$sample_id <- s$sample_submitter_id
  s$exclusion_reason <- NA_character_
  cl <- clinical[match(s$sample_id, clinical$sample_id, incomparables = NA), , drop = FALSE]
  counts <- list()
  record <- function(step) {
    inc <- is.na(s$exclusion_reason)
    counts[[length(counts) + 1L]] <<- data.frame(
      step = step, n_profiles = sum(inc), n_cases = length(unique(s$case_id[inc])),
      n_tumors = sum(inc & s$role %in% "tumor"), n_organoids = sum(inc & s$role %in% "organoid"))
  }
  exclude <- function(cond, reason) {
    hit <- cond & is.na(s$exclusion_reason)
    s$exclusion_reason[hit] <<- reason
  }
  record("all RNA profiles")
  exclude(s$sample_type %in% "Expanded Next Generation Cancer Model", "expanded model (duplicate of an original model)")
  exclude(is.na(s$role), "not a tumor or model sample")
  record("tumors and original models")
  s$model_type <- ifelse(s$role %in% "organoid", cl$growth_pattern, NA_character_)
  exclude(s$role %in% "organoid" & is.na(s$model_type), "model type unknown")
  exclude(s$role %in% "organoid" & !(s$model_type %in% "3D Organoid"), "model type not 3D organoid")
  record("3D organoids and tumors")
  tumor_key <- paste(s$case_id, s$tumor_descriptor)
  frozen <- unique(tumor_key[s$role %in% "tumor" & is.na(s$exclusion_reason) & !(s$preservation_method %in% "FFPE")])
  exclude(s$role %in% "tumor" & s$preservation_method %in% "FFPE" & tumor_key %in% frozen,
          "FFPE with a frozen profile of the same tumor")
  record("FFPE rule")
  inc_idx <- which(is.na(s$exclusion_reason))
  o <- inc_idx[order(s$sample_id[inc_idx], s$profile_id[inc_idx])]
  exclude(seq_len(nrow(s)) %in% o[duplicated(s$sample_id[o])], "duplicate profile of the same sample")
  record("one profile per sample")
  inc_t <- s[s$role %in% "tumor" & is.na(s$exclusion_reason), c("case_id", "tumor_descriptor", "sample_id")]
  key_t <- paste(inc_t$case_id, inc_t$tumor_descriptor)
  n_key <- table(key_t)
  key_o <- paste(s$case_id, s$tumor_descriptor)
  single <- key_o %in% names(n_key)[n_key == 1]
  s$parent_sample_id <- ifelse(s$role %in% "organoid" & single, inc_t$sample_id[match(key_o, key_t)], NA_character_)
  s$pair_flag <- ifelse(s$role %in% "organoid" & is.na(s$parent_sample_id),
                        "no single tumor with the model's descriptor: paired with all tumors of the case", NA_character_)
  pt <- pairTable(s[is.na(s$exclusion_reason), , drop = FALSE])
  exclude(s$role %in% "organoid" & !(s$profile_id %in% pt$organoid), "no included tumor partner")
  record("paired organoids and all tumors")
  ## Cancer type is a property of the case: the most frequent code among its samples, else the
  ## ICD-O site carried in the HCMI id (e.g. "site:C18"). makeFolds() needs one value per case.
  case_code <- tapply(cl$cancer_type, s$case_id, function(v) {
    v <- v[!is.na(v) & nzchar(v)]
    if (length(v)) names(sort(table(v), decreasing = TRUE))[1] else NA_character_
  })
  case_type <- unname(case_code[s$case_id])
  site <- paste0("site:", sub("^HCM-[A-Z]{4}-[0-9]{4}-([A-Z][0-9]{2}).*$", "\\1", s$case_submitter_id))
  case_type <- ifelse(is.na(case_type), site, case_type)
  ann <- data.frame(cancer_type = case_type,
                    prior_treatment = ifelse(is.na(cl$treatment_status), NA, cl$treatment_status == "Treated"),
                    age = suppressWarnings(as.numeric(cl$age)), sex = cl$sex, ancestry = cl$race_category,
                    stringsAsFactors = FALSE)
  parent_row <- match(s$parent_sample_id, s$sample_id, incomparables = NA)
  from <- ifelse(s$role %in% "organoid" & !is.na(parent_row), parent_row, seq_len(nrow(s)))
  s[c("cancer_type", "prior_treatment", "age", "sex", "ancestry")] <- ann[from, ]
  is_org <- s$role %in% "organoid"
  s$medium <- ifelse(is_org, cl$culture_media, NA_character_)
  s$passage <- ifelse(is_org, s$passage_count, NA)
  s$derivation_center <- ifelse(is_org, sub("^HCM-([A-Z]{4})-.*$", "\\1", s$case_submitter_id), NA_character_)
  s$model_id <- ifelse(is_org, s$sample_id, NA_character_)
  callers <- suppressWarnings(matrix(as.numeric(as.matrix(cl[, c("purity_absolute", "purity_ascat",
                                                                 "purity_hatchet", "purity_purple")])), nrow(s)))
  med <- matrixStats::rowMedians(callers, na.rm = TRUE)
  med[!is.finite(med)] <- NA_real_
  s$purity_dna <- ifelse(s$role %in% "tumor", med, NA_real_)
  s$consensus_purity <- ifelse(s$role %in% "tumor", suppressWarnings(as.numeric(cl$consensus_purity)), NA_real_)
  s$included <- is.na(s$exclusion_reason)
  rownames(s) <- NULL
  list(samples = s, filter_counts = do.call(rbind, counts))
}

#' sexFingerprintGenes
#'
#' Genes that mark a donor's sex rather than tumor biology: every Y-linked gene plus XIST and TSIX.
#' Pseudoautosomal genes (located on both X and Y) are kept. They would pair an organoid with its
#' tumor by sex alone.
#' @param hgnc HGNC table with `symbol` and `location` columns.
#' @return Gene symbols.
#' @export
sexFingerprintGenes <- function(hgnc) {
  hgnc$symbol[grepl("^Y", hgnc$location) | hgnc$symbol %in% c("XIST", "TSIX")]
}

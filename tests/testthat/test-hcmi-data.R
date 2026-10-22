write_star <- function(path, tpm = c(10.5, 0, 0, 40.2, 2)) {
  lines <- c("# gene-model: GENCODE v36",
             "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded",
             "N_unmapped\t\t\t10\t10\t10\t\t\t",
             "N_multimapping\t\t\t5\t5\t5\t\t\t",
             paste("ENSG00000000003.15", "TSPAN6", "protein_coding", 100, 50, 50, tpm[1], 3, 3, sep = "\t"),
             paste("ENSG00000000005.6", "TNMD", "protein_coding", 0, 0, 0, tpm[2], 0, 0, sep = "\t"),
             paste("ENSG00000002586.20_PAR_Y", "CD99", "protein_coding", 0, 0, 0, tpm[3], 0, 0, sep = "\t"),
             paste("ENSG00000002586.20", "CD99", "protein_coding", 500, 250, 250, tpm[4], 12, 13, sep = "\t"),
             paste("ENSG00000000419.13", "DPM1", "lncRNA", 20, 10, 10, tpm[5], 1, 1, sep = "\t"))
  writeLines(lines, path)
  path
}
hgnc_toy <- data.frame(symbol = c("TSPAN6", "TNMD", "CD99", "DPM1"),
                       ensembl_gene_id = c("ENSG00000000003", "ENSG00000000005", "ENSG00000002586", "ENSG00000000419"),
                       entrez_id = c(7105, 64102, 4267, 8813), stringsAsFactors = FALSE)

test_that("STAR counts parse and become a log2 TPM matrix of protein-coding HGNC genes", {
  f1 <- write_star(tempfile()); f2 <- write_star(tempfile(), tpm = c(1, 3, 0, 7, 9))
  x <- readStarCounts(f1)
  expect_false(any(startsWith(x$gene_id, "N_")))
  m <- buildExpressionMatrix(c(p1 = f1, p2 = f2), hgnc_toy)
  expect_equal(rownames(m), c("p1", "p2"))
  expect_equal(colnames(m), c("CD99", "TNMD", "TSPAN6"))
  expect_equal(unname(m["p1", "CD99"]), log2(40.2 + 1))
  expect_equal(unname(m["p2", "TSPAN6"]), log2(1 + 1))
})

test_that("selectGenes keeps the most variable genes", {
  m <- cbind(a = c(1, 1, 1), b = c(0, 5, 10), c = c(1, 2, 3))
  expect_equal(selectGenes(m, 2), c("b", "c"))
})

hcmi_inputs <- function() {
  gdc <- data.frame(
    profile_id = sprintf("f%02d", 1:14),
    case_id = c("A", "A", "A", "B", "C", "C", "D", "E", "E", "F", "G", "H", "H", "H"),
    case_submitter_id = sprintf("HCM-TEST-%04d-C18", c(1, 1, 1, 2, 3, 3, 4, 5, 5, 6, 7, 8, 8, 8)),
    sample_submitter_id = c("A-01A", "A-01Z", "A-85A", "B-01A", "C-01A", "C-85A", "D-85A", "E-01A", "E-01A",
                            "F-85A", "G-85A", "H-01A", "H-06A", "H-85A"),
    sample_type = c("Primary Tumor", "Slides", "Next Generation Cancer Model", "Primary Tumor", "Primary Tumor",
                    "Next Generation Cancer Model", "Next Generation Cancer Model", "Primary Tumor", "Primary Tumor",
                    "Next Generation Cancer Model", "Expanded Next Generation Cancer Model", "Primary Tumor",
                    "Metastatic", "Next Generation Cancer Model"),
    tumor_descriptor = c("primary", "primary", "primary", "primary", "primary", "primary", "primary", "primary",
                         "primary", "primary", "primary", "primary", "metastatic", "recurrence"),
    preservation_method = c("Frozen", "FFPE", "Frozen", "FFPE", rep("Frozen", 10)),
    passage_count = c(NA, NA, 6, NA, NA, 4, 5, NA, NA, 7, 8, NA, NA, 3), stringsAsFactors = FALSE)
  clinical <- data.frame(
    sample_id = c("A-01A", "A-85A", "B-01A", "C-01A", "C-85A", "D-85A", "E-01A", "F-85A", "H-01A", "H-06A", "H-85A"),
    growth_pattern = c("3D Organoid", "3D Organoid", "3D Organoid", "2D Cell Line", "2D Cell Line", NA, "3D Organoid",
                       "3D Organoid", "3D Organoid", "3D Organoid", "3D Organoid"),
    culture_media = c(NA, "OPAC", NA, NA, "CM", "X", NA, "Y", NA, NA, "Z"),
    treatment_status = c("Treated", "Treated", "Treatment Naive", NA, NA, NA, "Treated", NA, "Treatment Naive",
                         "Treatment Naive", "Treatment Naive"),
    cancer_type = c("COAD", "COAD", "PAAD", "BRCA", "BRCA", "SKCM", "COAD", "GBM", "PAAD", "PAAD", "PAAD"),
    age = c("60", "60", "55", "70", "70", "40", "62", "33", "50", "50", "50"),
    sex = c("F", "F", "M", "F", "F", "M", "M", "F", "M", "M", "M"), race_category = "White",
    purity_absolute = c("0.5", "1", "0.7", "0.6", "1", NA, "0.4", NA, "0.8", "0.3", "1"),
    purity_ascat = c("0.6", "1", NA, NA, "1", NA, "0.5", NA, NA, "0.35", "1"),
    purity_hatchet = NA_character_, purity_purple = NA_character_, consensus_purity = NA_character_,
    stringsAsFactors = FALSE)
  list(gdc = gdc, clinical = clinical)
}

test_that("buildSampleTable applies every inclusion rule", {
  x <- hcmi_inputs()
  b <- buildSampleTable(x$gdc, x$clinical)
  s <- stats::setNames(split(b$samples, seq_len(nrow(b$samples))), b$samples$profile_id)
  reason <- function(id) s[[id]]$exclusion_reason
  expect_match(reason("f02"), "FFPE")
  expect_true(s[["f04"]]$included)                           # FFPE kept when no frozen profile exists
  expect_match(reason("f06"), "not 3D organoid")
  expect_match(reason("f07"), "model type unknown")
  expect_match(reason("f09"), "duplicate")
  expect_match(reason("f10"), "no included tumor partner")
  expect_match(reason("f11"), "expanded")
  expect_equal(s[["f03"]]$parent_sample_id, "A-01A")
  expect_true(is.na(s[["f14"]]$parent_sample_id)); expect_false(is.na(s[["f14"]]$pair_flag))
  expect_equal(sort(b$samples$profile_id[b$samples$included & b$samples$role == "organoid"]), c("f03", "f14"))
})

test_that("organoids inherit their parent tumor's annotations; tumors get DNA purity", {
  x <- hcmi_inputs()
  s <- buildSampleTable(x$gdc, x$clinical)$samples
  row <- function(id) s[s$profile_id == id, ]
  expect_equal(row("f03")$cancer_type, "COAD"); expect_true(row("f03")$prior_treatment)
  expect_equal(row("f03")$medium, "OPAC"); expect_equal(row("f03")$passage, 6)
  expect_equal(row("f03")$derivation_center, "TEST")
  expect_true(is.na(row("f03")$purity_dna))
  expect_equal(row("f01")$purity_dna, 0.55)
  expect_equal(row("f13")$purity_dna, 0.325)
  expect_false(row("f14")$prior_treatment)
})

test_that("filter counts never increase", {
  x <- hcmi_inputs()
  fc <- buildSampleTable(x$gdc, x$clinical)$filter_counts
  expect_true(all(diff(fc$n_profiles) <= 0))
  expect_equal(utils::tail(fc$n_organoids, 1), 2)
})

test_that("an organoid's pair partner comes from the included tumors", {
  x <- hcmi_inputs()
  s <- buildSampleTable(x$gdc, x$clinical)$samples
  pt <- pairTable(s[s$included, ])
  expect_setequal(paste(pt$tumor, pt$organoid), c("f01 f03", "f12 f14", "f13 f14"))
})

test_that("cancer type is one value per case, falling back to the ICD site in the HCMI id", {   # makeFolds needs one type per case
  x <- hcmi_inputs()
  x$clinical$cancer_type[x$clinical$sample_id == "A-01A"] <- NA          # case A: only the organoid carries COAD
  x$clinical$cancer_type[x$clinical$sample_id %in% c("H-01A", "H-06A", "H-85A")] <- NA   # case H: no code anywhere
  s <- buildSampleTable(x$gdc, x$clinical)$samples
  per_case <- tapply(s$cancer_type, s$case_id, function(v) length(unique(v)))
  expect_true(all(per_case == 1))
  expect_equal(unique(s$cancer_type[s$case_id == "A"]), "COAD")
  expect_equal(unique(s$cancer_type[s$case_id == "H"]), "site:C18")
  expect_false(anyNA(s$cancer_type))
})

test_that("sexFingerprintGenes lists Y-linked genes plus XIST and TSIX, not pseudoautosomal ones", {
  hgnc <- data.frame(symbol = c("RPS4Y1", "XIST", "TSIX", "CD99", "GAPDH", "KDM5D"),
                     location = c("Yp11.2", "Xq13.2", "Xq13.2", "Xp22.33 and Yp11.2", "12p13.31", "Yq11.223"))
  expect_setequal(sexFingerprintGenes(hgnc), c("RPS4Y1", "XIST", "TSIX", "KDM5D"))
})

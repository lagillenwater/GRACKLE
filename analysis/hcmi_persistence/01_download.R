## Step 01: download every scripted input and verify the manual ones (spec Section 4.1).
## Writes raw/ and raw/manifest.csv (path, sha256, source, retrieved_utc).
source("analysis/hcmi_persistence/config.R")
suppressPackageStartupMessages(library(httr2))

manifest <- list()
note <- function(path, source) {
  manifest[[length(manifest) + 1L]] <<- data.frame(path = path, sha256 = sha256(path), source = source,
                                                    retrieved_utc = format(Sys.time(), tz = "UTC"))
  invisible(path)
}
fetch <- function(url, dest) {
  if (!file.exists(dest)) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    request(url) |> req_retry(max_tries = 5) |> req_perform(path = dest)
  }
  note(dest, url)
}

## 1. GDC STAR counts and their metadata ----------------------------------------------------
gdc_dir <- file.path(cfg$raw, "gdc"); dir.create(gdc_dir, recursive = TRUE, showWarnings = FALSE)
status <- request("https://api.gdc.cancer.gov/status") |> req_perform() |> resp_body_json()
fields <- c("file_id", "file_name", "md5sum", "cases.case_id", "cases.submitter_id",
            "cases.samples.submitter_id", "cases.samples.sample_type", "cases.samples.tumor_descriptor",
            "cases.samples.preservation_method", "cases.samples.passage_count")
filters <- list(op = "and", content = list(
  list(op = "in", content = list(field = "cases.project.project_id", value = list("HCMI-CMDC"))),
  list(op = "in", content = list(field = "data_type", value = list("Gene Expression Quantification"))),
  list(op = "in", content = list(field = "analysis.workflow_type", value = list("STAR - Counts")))))
hits <- request("https://api.gdc.cancer.gov/files") |>
  req_body_json(list(filters = filters, fields = paste(fields, collapse = ","), format = "JSON", size = 10000)) |>
  req_perform() |> resp_body_json(simplifyVector = FALSE)
hits <- hits$data$hits
stopifnot(all(vapply(hits, function(h) length(h$cases[[1]]$samples), integer(1)) == 1L))
get <- function(x) if (is.null(x)) NA else x
gdc <- rbindlist(lapply(hits, function(h) {
  cs <- h$cases[[1]]; sm <- cs$samples[[1]]
  data.table(profile_id = h$file_id, file_name = h$file_name, md5sum = h$md5sum, case_id = cs$case_id,
             case_submitter_id = cs$submitter_id, sample_submitter_id = get(sm$submitter_id),
             sample_type = get(sm$sample_type), tumor_descriptor = get(sm$tumor_descriptor),
             preservation_method = get(sm$preservation_method), passage_count = get(sm$passage_count))
}))
fwrite(gdc, file.path(gdc_dir, "gdc_files.tsv"), sep = "\t")
writeLines(status$data_release, file.path(gdc_dir, "GDC_DATA_RELEASE.txt"))
star_dir <- file.path(gdc_dir, "star")
todo <- gdc[!file.exists(file.path(star_dir, gdc$profile_id, gdc$file_name))]
batches <- split(todo$profile_id, ceiling(seq_len(nrow(todo)) / 100))
for (b in batches) {   # download batches (I/O), not a loop over data
  if (length(b) == 1) {
    dest <- file.path(star_dir, b, todo$file_name[todo$profile_id == b])
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    request(paste0("https://api.gdc.cancer.gov/data/", b)) |> req_retry(max_tries = 5) |> req_perform(path = dest)
  } else {
    tgz <- tempfile(fileext = ".tar.gz")
    request("https://api.gdc.cancer.gov/data") |> req_body_json(list(ids = as.list(b))) |>
      req_retry(max_tries = 5) |> req_perform(path = tgz)
    utils::untar(tgz, exdir = star_dir)
    unlink(tgz)
  }
}
paths <- file.path(star_dir, gdc$profile_id, gdc$file_name)
bad <- unname(tools::md5sum(paths)) != gdc$md5sum
if (any(bad)) stop("md5 mismatch for ", sum(bad), " STAR files, for example ", gdc$profile_id[bad][1])
note(file.path(gdc_dir, "gdc_files.tsv"), paste("GDC files endpoint, data release", status$data_release))
message("GDC: ", nrow(gdc), " STAR files, md5 verified, release ", status$data_release)

## 2. cBioPortal clinical data and mutations (REST API) --------------------------------------
cb <- "https://www.cbioportal.org/api"
cb_dir <- file.path(cfg$raw, "cbioportal"); dir.create(cb_dir, recursive = TRUE, showWarnings = FALSE)
clin <- function(type) {
  x <- request(sprintf("%s/studies/pancan_hcmi_2025/clinical-data", cb)) |>
    req_url_query(clinicalDataType = type) |> req_perform() |> resp_body_json(simplifyVector = TRUE)
  x <- as.data.table(x)
  x[, intersect(c("sampleId", "patientId", "clinicalAttributeId", "value"), names(x)), with = FALSE]
}
fwrite(clin("SAMPLE"), file.path(cb_dir, "clinical_sample_long.tsv"), sep = "\t")
fwrite(clin("PATIENT"), file.path(cb_dir, "clinical_patient_long.tsv"), sep = "\t")
## GET /mutations needs a gene id; the fetch endpoint takes the whole sample list.
muts <- request(sprintf("%s/molecular-profiles/pancan_hcmi_2025_mutations/mutations/fetch", cb)) |>
  req_url_query(projection = "DETAILED", pageSize = 10000000) |>
  req_body_json(list(sampleListId = "pancan_hcmi_2025_sequenced")) |>
  req_perform() |> resp_body_json(simplifyVector = TRUE)
fwrite(data.table(sample_id = muts$sampleId, Hugo_Symbol = muts$gene$hugoGeneSymbol,
                  Variant_Classification = muts$mutationType), file.path(cb_dir, "mutations.tsv"), sep = "\t")
writeLines(request(sprintf("%s/sample-lists/pancan_hcmi_2025_sequenced/sample-ids", cb)) |>
             req_perform() |> resp_body_json(simplifyVector = TRUE), file.path(cb_dir, "sequenced_samples.txt"))
invisible(lapply(file.path(cb_dir, c("clinical_sample_long.tsv", "clinical_patient_long.tsv", "mutations.tsv",
                                     "sequenced_samples.txt")), note, source = "cBioPortal API pancan_hcmi_2025"))

## 3. STRING v12 ------------------------------------------------------------------------------
string_dir <- file.path(cfg$raw, "string")
invisible(lapply(c("protein.links.full.v12.0", "protein.links.v12.0", "protein.aliases.v12.0"), function(f)
  fetch(sprintf("https://stringdb-downloads.org/download/%s/9606.%s.txt.gz", f, f),
        file.path(string_dir, sprintf("9606.%s.txt.gz", f)))))

## 4. HGNC, OncoKB, MSigDB Hallmark -----------------------------------------------------------
ann_dir <- file.path(cfg$raw, "annotation"); dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)
fetch("https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2024-10-01.txt",
      file.path(ann_dir, "hgnc_complete_set_2024-10-01.txt"))
fetch("https://www.oncokb.org/api/v1/utils/cancerGeneList.txt", file.path(ann_dir, "oncokb_cancerGeneList.txt"))
h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
saveRDS(split(h$gene_symbol, h$gs_name), file.path(ann_dir, "hallmark.rds"))
note(file.path(ann_dir, "hallmark.rds"), paste("msigdbr", packageVersion("msigdbr"), "collection H"))

## 5. DepMap. The portal blocks scripted downloads, so the 2D screens and NextGen metadata come
## from DepMap's figshare deposits (pinned file ids and md5s); only the genome-wide NextGen gene
## effects (portal file set "NextGen Model Manuscript 2026") are a manual download.
figshare <- data.frame(
  file = c("CRISPRGeneEffect.csv", "Model.csv", "AchillesNonessentialControls.csv",
           "nextgen_model_metadata.csv", "nextgen_screen_metadata.csv"),
  id = c(51064667, 51065297, 51063566, 64494438, 64494441),
  md5 = c("6edf7ade09b9b34199210b559d4745d3", "675210d17675f3517b0ce39a3c274f16", "ca6a757db1cae49b913423ce73b818cf",
          "97b26ec52b8a05e52df7abd507b8139f", "3f31162458b449dfd1dd1db0417d1e88"),
  source = c(rep("DepMap 24Q4 Public, figshare+ doi 10.25452/figshare.plus.27993248", 3),
             rep("NextGen manuscript metadata, figshare doi 10.6084/m9.figshare.29472362", 2)),
  stringsAsFactors = FALSE)
dm_dir <- file.path(cfg$raw, "depmap"); dir.create(dm_dir, recursive = TRUE, showWarnings = FALSE)
invisible(lapply(seq_len(nrow(figshare)), function(i) {
  dest <- file.path(dm_dir, figshare$file[i])
  if (!file.exists(dest) || unname(tools::md5sum(dest)) != figshare$md5[i]) {
    request(sprintf("https://ndownloader.figshare.com/files/%d", figshare$id[i])) |>
      req_retry(max_tries = 5) |> req_perform(path = dest)
  }
  if (unname(tools::md5sum(dest)) != figshare$md5[i]) stop("md5 mismatch for ", figshare$file[i])
  note(dest, figshare$source[i])
}))
ng_tab <- file.path(cfg$analysis, "depmap_nextgen_files.csv")
if (!file.exists(ng_tab)) stop("NextGen 3D gene effects are a manual download (the DepMap portal blocks scripts): see plan Task 22")
ng <- as.data.frame(fread(ng_tab))
if (!"gene_effect" %in% ng$role) stop("depmap_nextgen_files.csv needs a gene_effect row")
ng_paths <- file.path(cfg$manual, "depmap_nextgen", ng$filename)
if (!all(file.exists(ng_paths))) stop("missing NextGen files on Alpine: ", paste(ng$filename[!file.exists(ng_paths)], collapse = ", "))
invisible(lapply(ng_paths, note, source = paste("manual DepMap portal download:", unique(ng$release))))
sg_tab <- file.path(cfg$analysis, "sanger_files.csv")
if (file.exists(sg_tab)) {
  sg <- as.data.frame(fread(sg_tab))
  sg_paths <- file.path(cfg$manual, "sanger", sg$filename)
  if (!all(file.exists(sg_paths))) stop("missing Sanger files on Alpine: ", paste(sg$filename[!file.exists(sg_paths)], collapse = ", "))
  invisible(lapply(sg_paths, note, source = paste("manual Sanger download", unique(sg$release))))
} else {
  message("no sanger_files.csv: coess_3d will use DepMap 3D screens only")
}

fwrite(rbindlist(manifest), file.path(cfg$raw, "manifest.csv"))
message("01_download complete: ", length(manifest), " files in raw/manifest.csv")

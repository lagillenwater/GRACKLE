## Shared configuration for the HCMI persistence pipeline. Every script sources this file from
## the repository root. The roots come from environment variables set by scripts/alpine/env.sh
## on Alpine, or by fixture/run_fixture.R locally.
need_env <- function(v) {
  x <- Sys.getenv(v)
  if (!nzchar(x)) stop("environment variable ", v, " is not set (source scripts/alpine/env.sh)")
  x
}
cfg <- list(
  raw = file.path(need_env("GRACKLE_DATA"), "raw"),
  manual = file.path(need_env("GRACKLE_DATA"), "raw", "manual"),
  processed = file.path(need_env("GRACKLE_DATA"), "processed"),
  run = file.path(need_env("GRACKLE_SCRATCH"), Sys.getenv("GRACKLE_RUN_ID", "dev")),
  analysis = file.path("analysis", "hcmi_persistence"),
  seed = 20261015L,
  arms = c("coess_2d", "coess_3d", "string_mech", "string_all"))
invisible(lapply(cfg[c("raw", "processed", "run")], dir.create, recursive = TRUE, showWarnings = FALSE))

## --name value pairs; unknown names are an error.
parseArgs <- function(defaults) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) %% 2) stop("arguments must come in --name value pairs")
  if (length(a)) {
    keys <- sub("^--", "", a[c(TRUE, FALSE)])
    unknown <- setdiff(keys, names(defaults))
    if (length(unknown)) stop("unknown argument(s): ", paste(unknown, collapse = ", "))
    defaults[keys] <- a[c(FALSE, TRUE)]
  }
  defaults
}

sha256 <- function(paths) {
  tool <- if (nzchar(Sys.which("sha256sum"))) "sha256sum" else c("shasum", "-a", "256")
  vapply(paths, function(p) strsplit(system2(tool[1], c(tool[-1], shQuote(p)), stdout = TRUE), "\\s+")[[1]][1],
         character(1), USE.NAMES = FALSE)
}

## Sidecar <path>.provenance.json: commit, job, run, inputs' sha256, session.
writeProvenance <- function(path, inputs = character(0), extra = list()) {
  rec <- c(list(output = basename(path),
                output_sha256 = if (file.exists(path)) sha256(path) else NA,
                code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                slurm_job_id = Sys.getenv("SLURM_JOB_ID", NA),
                run_id = Sys.getenv("GRACKLE_RUN_ID", "dev"),
                created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                inputs = as.list(stats::setNames(sha256(inputs), basename(inputs))),
                session = utils::capture.output(utils::sessionInfo())), extra)
  jsonlite::write_json(rec, paste0(path, ".provenance.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
}

suppressPackageStartupMessages({
  library(data.table)
  devtools::load_all(quiet = TRUE)
})

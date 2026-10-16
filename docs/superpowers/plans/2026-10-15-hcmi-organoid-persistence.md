# HCMI Tumor-to-Organoid Persistence Map Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and run a leakage-free cross-fitting pipeline on CU Alpine. It maps which GRACKLE gene programs persist from HCMI parent tumors into matched 3D organoids, under four gene-graph priors.

**Architecture:**
- **Package functions in `R/`** hold all the logic, each with tests: graphs, tuning, matching, scoring, annotation, plus the cores of the fit and combine steps.
- **Thin numbered scripts** in `analysis/hcmi_persistence/` connect files on disk to those functions.
- **Slurm job scripts** in `scripts/alpine/` run the scripts on Alpine, submitted only through `ralpine`.
- **Local first:** every piece is checked on the laptop against a synthetic fixture before anything runs on Alpine.

**Tech Stack:**
- R 4.5 and TensorFlow (CPU) through `reticulate`, from one conda-forge environment built with micromamba on both the laptop and Alpine.
- R packages: `nnls`, `cluster`, `clue`, `irlba`, `data.table`, `matrixStats`, `withr`, `testthat` 3, `ggplot2`, `patchwork`.
- Slurm on Alpine.

**Spec:** `docs/superpowers/specs/2026-10-15-hcmi-organoid-persistence-design.md`. Read it before starting. Section references below (for example "spec 8.3") point into it.

## Global Constraints

**Git**
- Work on branch `organoid-transfer`. Never commit to `main` or `alpine`.
- Ask the user before any force-push, `git reset --hard`, branch deletion, or action touching `main` or `alpine`.
- The working tree has the user's uncommitted `tests/HTP_v2.R` and untracked `GRACKLE_data/`, `REPOSITORY_REVIEW.md` and `tests/HTP_k_sweep.R`. **Never stage them.** Always `git add` explicit paths.
- Commit messages end with: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

**User coding rules** (from `~/.claude/CLAUDE.md`)
- No emojis anywhere.
- Avoid nested loops; prefer vectorized code. A single `lapply` or `vapply` over independent units is fine.

**Model**
- The GRACKLE objective and multiplicative updates are unchanged (spec 5).

**Data and seeds**
- Seed `20261015` for folds, parallel analysis and every permutation null.
- Folds: 5, stratified by cancer type, grouped by case.
- Genes: protein-coding only, top 5,000 by variance of log2(TPM + 1), min-max scaled per gene. No gene is dropped for lacking graph edges.

**Fits and tuning**
- GRACKLE fits use `dtype = "float32"`, `iterations = 2000`, `diff_threshold = 1e-6`.
- Alpha grid: `alpha_1, alpha_2 in {0, 0.01, 0.03, 0.1, 0.3, 1}`, giving 36 settings.

**Sample graph**
- `k_nn = 10`; parts scaled to mean degree 1, then weighted pair : annot : purity = 1 : 1 : 1.
- Purity edges are tumor-tumor only, using DNA-based purity only.

**Gene graphs**
- Co-essentiality: a gene is informative if its SD exceeds the 95th percentile of the non-essential controls' SDs; GLS whitening; GLS correlation; same-arm pairs and negative values set to 0; `k_nn = 10`.
- STRING: prior 0.041.

**Matching**
- Threshold `tau` is the 99th percentile of the best-match null with gene loadings permuted, 1,000 permutations per fold.
- A program is stable if it matches in at least 4 of 5 folds.

**Persistence**
- 10,000 permutations within cancer type x fold, 2,000 bootstraps, BH FDR < 0.05.
- A program is "reduced" if its retention FDR < 0.05 and its median log2 ratio is <= -1.
- The per-type view needs at least 10 pairs.

**Annotation**
- Enrichment uses the top 5% of genes, a hypergeometric test against Hallmark, and the final gene set as background.
- Genotype uses OncoKB genes with non-silent mutations in at least 10 tumors and at least 2 wild-type tumors.

**Alpine**
- Access only through `scripts/alpine/ralpine`; never raw `ssh`.
- Partition `acpu`, QoS `cpu-normal`; memory = cores x 3,840 MB.
- Logs go to `logs/%x-%j.out` (`%x-%A_%a.out` for arrays).
- Run outputs go to `/scratch/alpine/$USER/grackle/hcmi_persistence/<run_id>/`, with `run_id = <UTC YYYYMMDD>_<short commit>`.

**Results**
- Run results stay uncommitted (gitignored) until the user has reviewed them; then they are committed with `git add -f`.

**Commands**
- Local R runs through `make` targets, which use `$(HOME)/.local/envs/grackle/bin/Rscript`.
- Run every command from the repository root.

## Review Focus

These are the inputs most likely to break the pipeline that no task's main tests cover. Each one gets a test in the task named.

1. **NNLS activities with many exact zeros** (an organoid side that is all 0, or tied values). Concordance must return `rho = 0` and `p = 1`, never `NA` or an error. Test in Task 11.
2. **A cancer type x fold stratum holding a single pair.** The permutation must leave that pair in place, and code must never call `sample()` on a length-1 integer vector. Test in Task 11.
3. **Missing or constant annotations** (NA age, NA sex, a constant ancestry column, NA purity). The sample graph must build without `NA`, and a tumor with NA purity must get no purity edges. Test in Task 5.
4. **Gene identifiers that map many-to-many or not at all** (STRING aliases, DepMap symbols). The graph must have no `NA`, duplicate gene pairs must be aggregated by max, and unmapped genes must have zero rows. Tests in Tasks 6 and 7.
5. **A fit missing or failed for one fold of a setting.** Combine must exclude that setting and list it, never score from partial folds and never error. Test in Task 14.

## File Structure

| Path | Responsibility | Task |
|---|---|---|
| `environment.yml`, `environment.lock.<platform>.txt` | Pinned conda-forge environment | 1, 19 |
| `scripts/dev/setup_env.sh` | Builds the environment with micromamba (laptop and Alpine) | 1 |
| `Makefile` | `env`, `document`, `test`, `e2e`, `fixture` targets | 1 |
| `DESCRIPTION`, `NAMESPACE`, `R/Benchmarks.R`, `R/simulation_functions.R`, `R/Similarity.R` | Package loads in a clean environment | 2 |
| `R/GRACKLE.R` | Loss/convergence fixes, `svdInit()` | 3 |
| `R/folds.R` | `makeFolds()` | 4 |
| `R/sample_graph.R` | Pair, annotation and purity parts; `buildSampleGraph()` | 5 |
| `R/gene_graph.R` | STRING graphs (Task 6), co-essentiality graphs (Task 7) | 6, 7 |
| `R/rank_selection.R` | `parallelAnalysisRank()` | 8 |
| `R/tuning.R` | `lambdaFromAlpha()`, `selectSetting()` | 9 |
| `R/program_matching.R` | Hungarian matching, null threshold, stability | 10 |
| `R/persistence.R` | Concordance, retention, nulls, groups | 11 |
| `R/annotation.R` | Enrichment, purity/ESTIMATE, genotype | 12 |
| `R/cross_fitting.R` | `fitTask()`, the core of `05_fit.R` | 13 |
| `R/combine.R` | `combineArm()`, the core of `06_combine.R` | 14 |
| `R/compare.R` | Transfer table, map agreement | 15 |
| `R/fixture.R` | Synthetic fixture and its expected calls | 16 |
| `R/hcmi_data.R`, `R/depmap_data.R` | Parsing GDC, cBioPortal and DepMap inputs | 17 |
| `scripts/alpine/ralpine`, `scripts/alpine/*.sbatch`, `scripts/alpine/env.sh`, `scripts/alpine/submit_chain.sh` | Alpine access and jobs | 18, 19, 20, 21 |
| `analysis/hcmi_persistence/*.R`, `analysis/hcmi_persistence/tasks/*.csv` | Pipeline scripts and task tables | 20, 21 |
| `tests/testthat/*` | Unit tests, fixtures, end-to-end test | every code task |

---

### Task 1: Pinned environment on the laptop

**Files:**
- Create: `environment.yml`
- Create: `scripts/dev/setup_env.sh`
- Create: `Makefile`
- Create: `logs/.gitkeep`
- Modify: `.gitignore`

**Interfaces:**
- Produces: `make env`, which builds `$HOME/.local/envs/grackle`. `make test FILTER=<name>` runs `devtools::test()` in that environment. Later tasks run every test through these targets.

- [ ] **Step 1: Write `environment.yml`**

```yaml
name: grackle
channels:
  - conda-forge
  - nodefaults
dependencies:
  - python=3.12
  - r-base=4.5
  - tensorflow=2.21
  - r-tensorflow=2.20
  - r-reticulate
  - r-nnls
  - r-igraph
  - r-cluster
  - r-clue
  - r-irlba
  - r-data.table
  - r-r.utils
  - r-matrixstats
  - r-matrix
  - r-glmnet
  - r-fnn
  - r-ppcor
  - r-rcolorbrewer
  - r-dplyr
  - r-tidyr
  - r-tibble
  - r-ggplot2
  - r-patchwork
  - r-jsonlite
  - r-httr2
  - r-msigdbr
  - r-readxl
  - r-withr
  - r-devtools
  - r-roxygen2
  - r-testthat
  - r-remotes
```

- [ ] **Step 2: Write `scripts/dev/setup_env.sh` and make it executable**

```bash
#!/usr/bin/env bash
# Builds the GRACKLE conda-forge environment with a standalone micromamba binary.
# The same script runs on the laptop and on Alpine (scripts/alpine/setup_env.sbatch).
#
#   scripts/dev/setup_env.sh <env-prefix> <micromamba-dir>
#
# Laptop:  scripts/dev/setup_env.sh "$HOME/.local/envs/grackle" "$HOME/.local/bin"
# Alpine:  scripts/dev/setup_env.sh "/projects/$USER/software/envs/grackle" "/projects/$USER/software/bin"
#
# Writes environment.lock.<platform>.txt (an explicit package list) into the repository root.
set -euo pipefail
PREFIX="${1:?env prefix required}"
MAMBA_DIR="${2:?micromamba directory required}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
case "$(uname -s)-$(uname -m)" in
  Darwin-arm64) PLATFORM=osx-arm64 ;;
  Linux-x86_64) PLATFORM=linux-64 ;;
  *) echo "unsupported platform $(uname -s)-$(uname -m)" >&2; exit 1 ;;
esac
MAMBA="$MAMBA_DIR/micromamba"
if [[ ! -x "$MAMBA" ]]; then
  mkdir -p "$MAMBA_DIR"
  curl -Ls "https://micro.mamba.pm/api/micromamba/$PLATFORM/latest" \
    | tar -xj -C "$MAMBA_DIR" --strip-components=1 bin/micromamba
fi
export MAMBA_ROOT_PREFIX="$(dirname "$PREFIX")/.mamba"
if [[ -d "$PREFIX/conda-meta" ]]; then
  "$MAMBA" install -y -p "$PREFIX" -f "$REPO/environment.yml"
else
  "$MAMBA" create -y -p "$PREFIX" -f "$REPO/environment.yml"
fi
# tidyestimate is CRAN-only; pin it.
"$PREFIX/bin/Rscript" -e 'if (!requireNamespace("tidyestimate", quietly = TRUE) || as.character(packageVersion("tidyestimate")) != "1.1.1") remotes::install_version("tidyestimate", version = "1.1.1", repos = "https://cloud.r-project.org", upgrade = "never")'
"$MAMBA" env export -p "$PREFIX" --explicit > "$REPO/environment.lock.$PLATFORM.txt"
RETICULATE_PYTHON="$PREFIX/bin/python" "$PREFIX/bin/Rscript" -e 'library(tensorflow); x <- tf$constant(1); stopifnot(as.numeric(x$numpy()) == 1); cat("tensorflow", tf$version$VERSION, "ok\n")'
echo "environment ready: $PREFIX ($PLATFORM)"
```

Run: `chmod +x scripts/dev/setup_env.sh`

- [ ] **Step 3: Write `Makefile`**

```make
# Development entry points. GRACKLE_ENV is the micromamba environment built by
# scripts/dev/setup_env.sh; on Alpine, scripts/alpine/env.sh exports it.
GRACKLE_ENV ?= $(HOME)/.local/envs/grackle
RSCRIPT ?= $(GRACKLE_ENV)/bin/Rscript
RETICULATE_PYTHON ?= $(GRACKLE_ENV)/bin/python
export RETICULATE_PYTHON
FILTER ?=

.PHONY: env document test e2e fixture
env:
	scripts/dev/setup_env.sh "$(GRACKLE_ENV)" "$(HOME)/.local/bin"
document:
	$(RSCRIPT) -e 'devtools::document()'
test:
	$(RSCRIPT) -e 'f <- "$(FILTER)"; devtools::test(filter = if (nzchar(f)) f else NULL, stop_on_failure = TRUE)'
e2e:
	GRACKLE_RUN_E2E=true $(RSCRIPT) -e 'devtools::test(filter = "e2e", stop_on_failure = TRUE)'
fixture:
	$(RSCRIPT) analysis/hcmi_persistence/fixture/run_fixture.R
```

- [ ] **Step 4: Update `.gitignore` and create `logs/.gitkeep`**

Append to `.gitignore`:

```
GRACKLE_data/
logs/*
!logs/.gitkeep
analysis/hcmi_persistence/results/
analysis/hcmi_persistence/fixture/out/
_moved_aside/
*.R~
.#*
\#*#
```

Run: `mkdir -p logs && touch logs/.gitkeep`

`logs/` has to exist in a fresh checkout. Slurm does not create the output directory, and a job whose log path is missing fails without any output.

- [ ] **Step 5: Build the environment**

Run: `make env`
Expected: the last lines are `tensorflow 2.21.<x> ok` and `environment ready: /Users/<user>/.local/envs/grackle (osx-arm64)`.

If the solve fails on `tensorflow=2.21` or `r-base=4.5`, change that line in `environment.yml` to the newest version conda-forge has for **both** osx-arm64 and linux-64. Checked 2026-10-15: tensorflow 2.19.1 and 2.21.0, r-base 4.5.x and 4.6.x. Record the change in the commit message.

- [ ] **Step 6: Commit**

```bash
git add environment.yml environment.lock.osx-arm64.txt scripts/dev/setup_env.sh Makefile .gitignore logs/.gitkeep
git commit -m "build: pinned conda-forge environment and make targets

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Package loads in a clean environment

**Files:**
- Modify: `R/Benchmarks.R` (`runNMF`, `runPLIER`)
- Modify: `R/simulation_functions.R` (`simulateExpression`)
- Modify: `R/GRACKLE.R:8` (remove `#' @import dynutils`)
- Modify: `R/Similarity.R:6` (remove `#' @useDynLib GRACKLE`)
- Delete: `R/RcppExports.R`
- Modify: `DESCRIPTION`
- Regenerate: `NAMESPACE`, `man/`
- Create: `tests/testthat.R`, `tests/testthat/test-load.R`

**Interfaces:**
- Produces: `devtools::load_all()` and `devtools::test()` work without `PLIER`, `sgnesR`, `NMF` or `dynutils` installed. Every later task depends on this.

- [ ] **Step 1: Write the failing test**

`tests/testthat.R`:

```r
library(testthat)
library(GRACKLE)
test_check("GRACKLE")
```

`tests/testthat/test-load.R`:

```r
test_that("the namespace does not import optional packages", {
  imports <- names(getNamespaceImports("GRACKLE"))
  expect_false(any(c("PLIER", "sgnesR", "NMF", "dynutils") %in% imports))
  expect_true(is.function(asNamespace("GRACKLE")$GRACKLE))
})

test_that("the compiled-code stubs are gone", {
  expect_false(exists("nmf_with_graph_regularization", envir = asNamespace("GRACKLE"), inherits = FALSE))
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `make test FILTER=load`
Expected: FAIL. `load_all` errors with "there is no package called 'PLIER'" (or `sgnesR`), because the current `NAMESPACE` imports it.

- [ ] **Step 3: Make the optional dependencies optional**

In `R/Benchmarks.R`:
- Delete the roxygen line `#' @import NMF` (line 8).
- Replace the body of `runNMF` with:

```r
runNMF <- function(expression_matrix, k = 5, method = 'lee', seed = 42) {
  if (!requireNamespace("NMF", quietly = TRUE)) stop("runNMF needs the NMF package (install it separately)")
  res <- NMF::nmf(expression_matrix, rank = k, seed = seed)
  W <- NMF::basis(res)
  H <- NMF::coef(res)
  return(list(W = W, H = H))
}
```

- Delete the roxygen line `#' @import PLIER`.
- Replace the body of `runPLIER` with:

```r
runPLIER <- function(expression_matrix, k) {
  if (!requireNamespace("PLIER", quietly = TRUE)) stop("runPLIER needs the PLIER package (GitHub: wgmao/PLIER)")
  env <- new.env()
  utils::data("canonicalPathways", package = "PLIER", envir = env)
  res <- PLIER::PLIER(data = expression_matrix, priorMat = env$canonicalPathways, k, rseed = 42)
  return(list(W = res$B, H = res$Z))
}
```

In `R/simulation_functions.R`, in the `simulateExpression` block:
- Delete the roxygen line `#' @import sgnesR ` (it has a trailing space).
- Make the first line of the function body:

```r
  if (!requireNamespace("sgnesR", quietly = TRUE)) stop("simulateExpression needs the sgnesR package")
```

- Replace `rp<-new("rsgns.param",time=0,stop_time=time,readout_interval= time)` with:

```r
    rp <- methods::new(methods::getClass("rsgns.param", where = asNamespace("sgnesR")), time = 0, stop_time = time, readout_interval = time)
```

- Replace `rsg <- new("rsgns.data",network=g, rconst=rc)` with:

```r
    rsg <- methods::new(methods::getClass("rsgns.data", where = asNamespace("sgnesR")), network = g, rconst = rc)
```

- Replace `xx <- rsgns.rn(rsg, rp, timeseries=FALSE, sample=1)` with:

```r
    xx <- sgnesR::rsgns.rn(rsg, rp, timeseries = FALSE, sample = 1)
```

In `R/GRACKLE.R`, delete the line `#' @import dynutils`. In `R/Similarity.R`, delete the line `#' @useDynLib GRACKLE`.

Then: `git rm R/RcppExports.R`

- [ ] **Step 4: Rewrite `DESCRIPTION`**

Keep the Title, Authors, Description and License lines exactly as they are; licensing is out of scope. Replace everything from `Imports:` to the end with:

```
Imports:
    clue,
    cluster,
    data.table,
    dplyr,
    FNN,
    glmnet,
    igraph,
    irlba,
    jsonlite,
    Matrix,
    matrixStats,
    methods,
    nnls,
    parallel,
    ppcor,
    RColorBrewer,
    stats,
    tensorflow,
    tibble,
    tidyr,
    utils,
    withr
Suggests:
    NMF,
    PLIER,
    sgnesR,
    testthat (>= 3.0.0),
    tidyestimate
Config/testthat/edition: 3
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.2
```

- [ ] **Step 5: Regenerate the documentation and run the test**

Run: `make document`
Expected: `NAMESPACE` is rewritten with no `import(PLIER)`, `import(sgnesR)`, `import(NMF)`, `import(dynutils)` or `useDynLib`. A warning about the `RoxygenNote` version is fine.

Run: `make test FILTER=load`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add DESCRIPTION NAMESPACE man R/Benchmarks.R R/simulation_functions.R R/GRACKLE.R R/Similarity.R tests/testthat.R tests/testthat/test-load.R
git commit -m "fix: package loads without GitHub-only or Bioconductor dependencies

PLIER, sgnesR and NMF move to Suggests behind requireNamespace checks;
drop the unused dynutils import and the orphaned useDynLib/RcppExports
(there is no src/ directory).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

(`git rm` already staged the deletion of `R/RcppExports.R`.)

---

### Task 3: GRACKLE loss, convergence and precision fixes

**Files:**
- Modify: `R/GRACKLE.R` (replace everything from the top of the file through the closing brace of `GRACKLE()`, which is the line before `#' colScale`; keep `colScale` and the commented `matMult` below it)
- Create: `tests/testthat/helper-grackle.R`, `tests/testthat/test-grackle.R`

**Interfaces:**
- Produces: `GRACKLE(Y, net_similarity, patient_similarity, diff_threshold = 1e-4, lambda_1 = .5, lambda_2 = .5, k = 5, verbose = FALSE, beta = 0, error_terms = FALSE, iterations = 100, svd_W, svd_H, dtype = c("float32", "float16", "float64"))`. It returns `list(W, H, svd_W, svd_H, loss_trace, n_iter, converged)`:
  - `loss_trace` has columns `iteration, reconstruction, sample_graph, gene_graph, total`.
  - `W` has dimnames `(rownames(Y), LV1..LVk)`; `H` has `(LV1..LVk, colnames(Y))`.
- Produces: `svdInit(Y, k)` returning `list(W, H)`.

- [ ] **Step 1: Write the test helpers**

`tests/testthat/helper-grackle.R`:

```r
skip_if_no_tf <- function() {
  testthat::skip_if_not(reticulate::py_module_available("tensorflow"), "TensorFlow not available")
}

toy_problem <- function(n = 20, g = 30, k = 3, seed = 1) {
  withr::with_seed(seed, {
    Y <- matrix(stats::runif(n * g), n, g,
                dimnames = list(paste0("s", seq_len(n)), paste0("g", seq_len(g))))
    S_p <- matrix(stats::rbinom(n * n, 1, 0.2), n, n)
    S_g <- matrix(stats::rbinom(g * g, 1, 0.1), g, g)
  })
  S_p <- pmax(S_p, t(S_p)); diag(S_p) <- 0
  S_g <- pmax(S_g, t(S_g)); diag(S_g) <- 0
  list(Y = Y, S_p = S_p, S_g = S_g, init = svdInit(Y, k), k = k)
}

mu_nmf <- function(Y, W, H, iterations, eps = 1e-10) {
  for (i in seq_len(iterations)) {
    H <- H * (t(W) %*% Y) / (t(W) %*% W %*% H + eps)
    W <- W * (Y %*% t(H)) / (W %*% H %*% t(H) + eps)
  }
  list(W = W, H = H)
}
```

(`mu_nmf` loops over iterations because the algorithm is iterative. That is not a data loop.)

- [ ] **Step 2: Write the failing tests**

`tests/testthat/test-grackle.R`:

```r
test_that("svdInit returns nonnegative factors scaled to [0, 1]", {
  p <- toy_problem()
  expect_equal(dim(p$init$W), c(20, 3))
  expect_equal(dim(p$init$H), c(3, 30))
  expect_true(all(p$init$W >= 0) && all(p$init$H >= 0))
  expect_equal(unname(apply(p$init$W, 2, max)), rep(1, 3))
  expect_equal(unname(apply(p$init$H, 1, max)), rep(1, 3))
})

test_that("the objective does not increase across iterations", {
  skip_if_no_tf()
  p <- toy_problem()
  traces <- lapply(c(0, 0.1), function(lam) {
    GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = lam, lambda_2 = lam, k = p$k,
            iterations = 50, svd_W = p$init$W, svd_H = p$init$H, dtype = "float64")$loss_trace$total
  })
  lapply(traces, function(tr) expect_true(all(diff(tr) <= 1e-9 * tr[-1])))
})

test_that("loss_trace reports the true objective parts", {
  skip_if_no_tf()
  p <- toy_problem()
  fit <- GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = 0.1, lambda_2 = 0.1, k = p$k,
                 iterations = 5, svd_W = p$init$W, svd_H = p$init$H, dtype = "float64")
  lt <- fit$loss_trace
  expect_named(lt, c("iteration", "reconstruction", "sample_graph", "gene_graph", "total"))
  expect_equal(lt$total, lt$reconstruction + lt$sample_graph + lt$gene_graph, tolerance = 1e-10)
  L_p <- diag(rowSums(p$S_p)) - p$S_p
  L_g <- diag(rowSums(p$S_g)) - p$S_g
  W <- unname(fit$W); H <- unname(fit$H)
  expect_equal(lt$reconstruction[5], sum((p$Y - W %*% H)^2), tolerance = 1e-8)
  expect_equal(lt$sample_graph[5], 0.1 * sum(diag(t(W) %*% L_p %*% W)), tolerance = 1e-8)
  expect_equal(lt$gene_graph[5], 0.1 * sum(diag(H %*% L_g %*% t(H))), tolerance = 1e-8)
})

test_that("iterations is honored and convergence is reported", {
  skip_if_no_tf()
  p <- toy_problem()
  fit <- GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = 0.1, lambda_2 = 0.1, k = p$k,
                 iterations = 7, svd_W = p$init$W, svd_H = p$init$H)
  expect_equal(fit$n_iter, 7L)
  expect_false(fit$converged)
  expect_equal(nrow(fit$loss_trace), 7L)
  fit2 <- GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 1, lambda_1 = 0.1, lambda_2 = 0.1, k = p$k,
                  iterations = 50, svd_W = p$init$W, svd_H = p$init$H)
  expect_equal(fit2$n_iter, 1L)
  expect_true(fit2$converged)
})

test_that("dtype is honored", {
  skip_if_no_tf()
  p <- toy_problem()
  run <- function(dt) GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = 0.1, lambda_2 = 0.1,
                              k = p$k, iterations = 20, svd_W = p$init$W, svd_H = p$init$H, dtype = dt)
  f64 <- run("float64"); f32 <- run("float32"); f16 <- run("float16")
  expect_equal(f32$W, f64$W, tolerance = 1e-3)
  expect_false(isTRUE(all.equal(f32$W, f64$W, tolerance = 1e-12)))
  expect_true(all(is.finite(f16$loss_trace$total)))
  expect_error(run("int8"))
})

test_that("with lambda = 0 GRACKLE equals plain multiplicative-update NMF", {
  skip_if_no_tf()
  p <- toy_problem()
  fit <- GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = 0, lambda_2 = 0, k = p$k,
                 iterations = 25, svd_W = p$init$W, svd_H = p$init$H, dtype = "float64")
  ref <- mu_nmf(p$Y, p$init$W, p$init$H, 25)
  expect_equal(unname(fit$W), unname(ref$W), tolerance = 1e-8)
  expect_equal(unname(fit$H), unname(ref$H), tolerance = 1e-8)
})

test_that("dimnames carry sample and gene names", {
  skip_if_no_tf()
  p <- toy_problem()
  fit <- GRACKLE(p$Y, p$S_g, p$S_p, k = p$k, iterations = 2, svd_W = p$init$W, svd_H = p$init$H)
  expect_equal(rownames(fit$W), rownames(p$Y))
  expect_equal(colnames(fit$H), colnames(p$Y))
  expect_equal(colnames(fit$W), paste0("LV", 1:3))
})
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `make test FILTER=grackle`
Expected: FAIL. `svdInit` is not found, and `loss_trace` is `NULL`.

- [ ] **Step 4: Implement**

In `R/GRACKLE.R`, replace everything from the top of the file through the closing brace of `GRACKLE()` with:

```r
## GRACKLE algorithm

#' GRACKLE
#'
#' Graph Regularized non-negative matrix factorization Across Contextual KnowLedgE.
#' Factorizes `Y` (samples x genes) as `W %*% H` with a sample-graph penalty on `W` and a
#' gene-graph penalty on `H`, using multiplicative updates. The objective is
#' `||Y - WH||_F^2 + lambda_1 tr(W' L_p W) + lambda_2 tr(H L_g H')` with `L = D - S` for each
#' similarity matrix `S`. Fitting stops when the relative change in that objective falls
#' below `diff_threshold`, or after `iterations` iterations.
#'
#' @importFrom Matrix Matrix
#' @import tensorflow
#' @param Y Numeric matrix, samples x genes, nonnegative.
#' @param net_similarity Gene similarity (genes x genes), symmetric and nonnegative.
#' @param patient_similarity Sample similarity (samples x samples), symmetric and nonnegative.
#' @param diff_threshold Relative change in the objective that counts as converged. Default 1e-4.
#' @param lambda_1 Sample-graph penalty weight. Default 0.5.
#' @param lambda_2 Gene-graph penalty weight. Default 0.5.
#' @param k Number of latent variables. Default 5.
#' @param verbose Print the objective at each iteration. Default FALSE.
#' @param beta Unused; kept for backward compatibility.
#' @param error_terms Unused; kept for backward compatibility.
#' @param iterations Maximum number of iterations. Default 100.
#' @param svd_W Initial W (samples x k), for example from [svdInit()].
#' @param svd_H Initial H (k x genes), for example from [svdInit()].
#' @param dtype TensorFlow precision: "float32" (default), "float16" or "float64".
#' @return A list with `W`, `H`, `svd_W`, `svd_H`, `loss_trace` (a data frame with columns
#'   iteration, reconstruction, sample_graph, gene_graph, total), `n_iter` and `converged`.
#' @export
GRACKLE <- function(Y, net_similarity, patient_similarity, diff_threshold = 1e-4,
                    lambda_1 = .5, lambda_2 = .5, k = 5, verbose = FALSE, beta = 0,
                    error_terms = FALSE, iterations = 100, svd_W, svd_H,
                    dtype = c("float32", "float16", "float64")) {
  dtype <- match.arg(dtype)
  Sys.setenv(TF_ENABLE_ONEDNN_OPTS = "0", TF_CPP_MIN_LOG_LEVEL = "3")
  gpus <- tf$config$list_physical_devices("GPU")
  if (length(gpus) > 0) tf$config$experimental$set_memory_growth(gpus[[1]], TRUE)

  tf_dtype <- switch(dtype, float32 = tf$float32, float16 = tf$float16, float64 = tf$float64)
  as_tf <- function(x) tf$constant(unname(as.matrix(x)), dtype = tf_dtype)
  ## Bookkeeping in at least float32 so float16 sums cannot overflow.
  wide <- function(x) if (dtype == "float16") tf$cast(x, tf$float32) else x

  Y_tf <- as_tf(Y)
  W <- as_tf(svd_W)
  H <- as_tf(svd_H)
  S_p <- as_tf(patient_similarity)
  S_g <- as_tf(net_similarity)
  d_p <- tf$reshape(tf$reduce_sum(S_p, axis = 1L), c(-1L, 1L))  # sample degrees, n x 1
  d_g <- tf$reshape(tf$reduce_sum(S_g, axis = 0L), c(1L, -1L))  # gene degrees, 1 x g
  l1 <- tf$constant(lambda_1, dtype = tf_dtype)
  l2 <- tf$constant(lambda_2, dtype = tf_dtype)
  ## Guards 0/0 in the multiplicative updates; far below any real denominator.
  eps <- tf$constant(if (dtype == "float16") 1e-7 else 1e-10, dtype = tf_dtype)

  objective <- function(W, H) {
    recon <- tf$reduce_sum(tf$square(wide(Y_tf - tf$matmul(W, H))))
    LpW <- wide(d_p * W - tf$matmul(S_p, W))   # L_p W without forming L_p
    HLg <- wide(H * d_g - tf$matmul(H, S_g))   # H L_g without forming L_g
    parts <- c(reconstruction = as.numeric(recon$numpy()),
               sample_graph = lambda_1 * as.numeric(tf$reduce_sum(wide(W) * LpW)$numpy()),
               gene_graph = lambda_2 * as.numeric(tf$reduce_sum(wide(H) * HLg)$numpy()))
    c(parts, total = sum(parts))
  }

  trace <- matrix(NA_real_, nrow = iterations, ncol = 4,
                  dimnames = list(NULL, c("reconstruction", "sample_graph", "gene_graph", "total")))
  prev <- objective(W, H)
  if (!is.finite(prev[["total"]])) stop("GRACKLE: the initial objective is not finite")
  converged <- FALSE
  n_iter <- 0L
  for (i in seq_len(iterations)) {
    H <- H * ((tf$matmul(W, Y_tf, transpose_a = TRUE) + l2 * tf$matmul(H, S_g)) /
              (tf$matmul(tf$matmul(W, W, transpose_a = TRUE), H) + l2 * (H * d_g) + eps))
    W <- W * ((tf$matmul(Y_tf, H, transpose_b = TRUE) + l1 * tf$matmul(S_p, W)) /
              (tf$matmul(W, tf$matmul(H, H, transpose_b = TRUE)) + l1 * (d_p * W) + eps))
    curr <- objective(W, H)
    if (!is.finite(curr[["total"]])) stop("GRACKLE: the objective became non-finite at iteration ", i)
    trace[i, ] <- curr
    n_iter <- i
    if (verbose) message(sprintf("iteration %d: total %.6g", i, curr[["total"]]))
    if (abs(prev[["total"]] - curr[["total"]]) / prev[["total"]] < diff_threshold) {
      converged <- TRUE
      break
    }
    prev <- curr
  }

  W <- W$numpy()
  H <- H$numpy()
  lv <- paste0("LV", seq_len(k))
  dimnames(W) <- list(rownames(Y), lv)
  dimnames(H) <- list(lv, colnames(Y))
  loss_trace <- data.frame(iteration = seq_len(n_iter), trace[seq_len(n_iter), , drop = FALSE])
  list(W = W, H = H, svd_W = svd_W, svd_H = svd_H, loss_trace = loss_trace,
       n_iter = n_iter, converged = converged)
}

#' svdInit
#'
#' SVD-based nonnegative starting values for [GRACKLE()], matching the HTP analysis scripts:
#' absolute values of the leading singular vectors (W scaled by the singular values), each
#' factor min-max scaled.
#'
#' @param Y Numeric matrix, samples x genes.
#' @param k Number of latent variables.
#' @return A list with `W` (samples x k) and `H` (k x genes).
#' @export
svdInit <- function(Y, k) {
  s <- svd(Y, nu = k, nv = k)
  W <- abs(s$u[, seq_len(k), drop = FALSE] %*% diag(s$d[seq_len(k)], k))
  H <- t(abs(s$v[, seq_len(k), drop = FALSE]))
  W <- min_max_scale(W, apply(W, 2, min), apply(W, 2, max))
  H <- t(min_max_scale(t(H), apply(H, 1, min), apply(H, 1, max)))
  list(W = W, H = H)
}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `make document && make test FILTER=grackle`
Expected: PASS (7 tests). Tests that need TensorFlow are not skipped, because the environment from Task 1 has it.

- [ ] **Step 6: Commit**

```bash
git add R/GRACKLE.R NAMESPACE man/GRACKLE.Rd man/svdInit.Rd tests/testthat/helper-grackle.R tests/testthat/test-grackle.R
git commit -m "fix: GRACKLE stopping rule follows the true objective

Loss uses the squared Frobenius norm and traces instead of sums over all
entries; iterations is honored instead of a hard-coded 292; convergence
and the per-iteration loss parts are returned; dtype argument (float32
default) replaces the global mixed_float16 policy; retry loops removed
so errors fail fast. Adds svdInit(). Updates are unchanged.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Folds

**Files:**
- Create: `R/folds.R`, `tests/testthat/test-folds.R`

**Interfaces:**
- Produces: `makeFolds(samples, n_folds = 5, seed = 20261015)`. It takes a data frame with `case_id` and `cancer_type` and returns it with an integer `fold` column added.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-folds.R`:

```r
make_cases <- function() {
  data.frame(profile_id = paste0("p", 1:40),
             case_id = rep(paste0("c", 1:20), each = 2),
             cancer_type = rep(rep(c("COAD", "PAAD", "BRCA", "SKCM"), times = c(8, 6, 4, 2)), each = 2),
             stringsAsFactors = FALSE)
}

test_that("every case lands in exactly one fold", {
  s <- makeFolds(make_cases())
  per_case <- tapply(s$fold, s$case_id, function(f) length(unique(f)))
  expect_true(all(per_case == 1))
  expect_true(all(s$fold %in% 1:5))
})

test_that("folds are balanced within each cancer type", {
  s <- makeFolds(make_cases())
  cases <- unique(s[, c("case_id", "cancer_type", "fold")])
  balance <- tapply(cases$fold, cases$cancer_type, function(f) diff(range(tabulate(f, nbins = 5))))
  expect_true(all(balance <= 1))
})

test_that("folds are reproducible with the seed and change with it", {
  expect_identical(makeFolds(make_cases())$fold, makeFolds(make_cases())$fold)
  expect_false(identical(makeFolds(make_cases())$fold, makeFolds(make_cases(), seed = 7)$fold))
})

test_that("a case with two cancer types is an error", {
  s <- make_cases(); s$cancer_type[2] <- "PAAD"
  expect_error(makeFolds(s), "more than one cancer_type")
})

test_that("makeFolds leaves the caller's RNG state alone", {
  set.seed(99); before <- runif(1)
  set.seed(99); invisible(makeFolds(make_cases())); after <- runif(1)
  expect_identical(before, after)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=folds`
Expected: FAIL with `could not find function "makeFolds"`.

- [ ] **Step 3: Implement**

`R/folds.R`:

```r
#' makeFolds
#'
#' Assigns every case to one of `n_folds` cross-validation folds, stratified by cancer type.
#' Cases are shuffled, sorted by cancer type, and dealt round-robin, so fold sizes within a
#' type differ by at most one. All profiles of a case share its fold.
#'
#' @param samples Data frame with columns `case_id` and `cancer_type`.
#' @param n_folds Number of folds. Default 5.
#' @param seed Random seed. Default 20261015. The caller's RNG state is restored.
#' @return `samples` with an integer `fold` column.
#' @export
makeFolds <- function(samples, n_folds = 5, seed = 20261015) {
  cases <- unique(samples[, c("case_id", "cancer_type")])
  if (anyDuplicated(cases$case_id)) {
    stop("makeFolds: a case has more than one cancer_type: ",
         paste(unique(cases$case_id[duplicated(cases$case_id)]), collapse = ", "))
  }
  cases <- withr::with_seed(seed, cases[sample.int(nrow(cases)), , drop = FALSE])
  cases <- cases[order(cases$cancer_type, method = "radix"), , drop = FALSE]
  cases$fold <- ((seq_len(nrow(cases)) - 1L) %% n_folds) + 1L
  samples$fold <- cases$fold[match(samples$case_id, cases$case_id)]
  samples
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `make document && make test FILTER=folds`
Expected: PASS (5 tests).

- [ ] **Step 5: Commit**

```bash
git add R/folds.R NAMESPACE man/makeFolds.Rd tests/testthat/test-folds.R
git commit -m "feat: makeFolds, case-grouped folds stratified by cancer type

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Sample graph

**Files:**
- Create: `R/sample_graph.R`, `tests/testthat/helper-samples.R`, `tests/testthat/test-sample-graph.R`

**Interfaces:**
- Consumes: a sample table with columns `profile_id, case_id, sample_id, parent_sample_id, role ("tumor"/"organoid"), cancer_type, tumor_descriptor, prior_treatment, age, sex, ancestry, purity_dna`. Organoid rows already carry their parent tumor's annotations; Task 17 does that inheritance.
- Produces:
  - `knnSparsify(S, k_nn, secondary = NULL)` returns a symmetric matrix with a zero diagonal.
  - `scaleMeanDegree(S)`.
  - `pairTable(samples)` returns `data.frame(tumor, organoid)` of profile ids.
  - `pairGraph(samples)`, `annotGraph(samples, vars, k_nn = 10)`, `purityGraph(samples, k_nn = 10)`.
  - `buildSampleGraph(samples, weights = c(pair = 1, annot = 1, purity = 1), k_nn = 10)` returns `list(S, summary)`, where `summary` has columns `part, n_edges, mean_degree, sparsity, weight`.

- [ ] **Step 1: Write the test fixture**

`tests/testthat/helper-samples.R`:

```r
## Six cases. c1: one tumor + organoid. c2: primary + metastasis, organoid from the metastasis.
## c3: two tumors, organoid with unknown parent. c4, c5: unpaired tumors. c6: tumor + organoid.
## Annotations include NA age, NA sex, a constant ancestry and one tumor without purity (t3b).
toy_samples <- function() {
  data.frame(
    profile_id = c("t1", "o1", "t2a", "t2b", "o2", "t3a", "t3b", "o3", "t4", "t5", "t6", "o6"),
    case_id = c("c1", "c1", "c2", "c2", "c2", "c3", "c3", "c3", "c4", "c5", "c6", "c6"),
    sample_id = c("c1-01", "c1-85", "c2-01", "c2-06", "c2-85", "c3-01", "c3-06", "c3-85",
                  "c4-01", "c5-01", "c6-01", "c6-85"),
    parent_sample_id = c(NA, "c1-01", NA, NA, "c2-06", NA, NA, NA, NA, NA, NA, "c6-01"),
    role = c("tumor", "organoid", "tumor", "tumor", "organoid", "tumor", "tumor", "organoid",
             "tumor", "tumor", "tumor", "organoid"),
    cancer_type = c("COAD", "COAD", "PAAD", "PAAD", "PAAD", "COAD", "COAD", "COAD", "BRCA", "COAD",
                    "PAAD", "PAAD"),
    tumor_descriptor = c("primary", "primary", "primary", "metastatic", "metastatic", "primary",
                         "metastatic", NA, "primary", "primary", "primary", "primary"),
    prior_treatment = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, NA, TRUE, FALSE, FALSE),
    age = c(60, 60, 55, 55, 55, NA, NA, NA, 70, 62, 48, 48),
    sex = c("F", "F", "M", "M", "M", "F", "F", "F", NA, "M", "F", "F"),
    ancestry = "EUR",
    purity_dna = c(0.5, NA, 0.7, 0.4, NA, 0.6, NA, NA, 0.8, 0.45, 0.65, NA),
    stringsAsFactors = FALSE)
}
```

- [ ] **Step 2: Write the failing tests**

`tests/testthat/test-sample-graph.R`:

```r
test_that("knnSparsify keeps each row's k largest positive similarities, then symmetrizes", {
  S <- matrix(c(0, .9, .8, .1,
                .9, 0, .2, 0,
                .8, .2, 0, .3,
                .1, 0, .3, 0), 4, byrow = TRUE, dimnames = list(letters[1:4], letters[1:4]))
  K <- knnSparsify(S, 1)
  expect_equal(K["a", "b"], .9); expect_equal(K["b", "a"], .9)
  expect_equal(K["a", "c"], .8); expect_equal(K["c", "d"], .3)
  expect_equal(K["a", "d"], 0)
  expect_true(isSymmetric(K)); expect_equal(unname(diag(K)), rep(0, 4))
})

test_that("pair edges link each organoid only to its own parent tumor", {
  P <- pairGraph(toy_samples())
  expect_equal(P["o1", "t1"], 1)
  expect_equal(P["o2", "t2b"], 1); expect_equal(P["o2", "t2a"], 0)
  expect_equal(P["o3", "t3a"], 1); expect_equal(P["o3", "t3b"], 1)
  expect_equal(P["o6", "t6"], 1)
  expect_equal(sum(P[upper.tri(P)]), 5)
  org <- toy_samples()$role == "organoid"
  expect_equal(sum(P[org, org]), 0)
})

test_that("pairTable lists the same pairs as pairGraph", {
  pt <- pairTable(toy_samples())
  expect_setequal(paste(pt$tumor, pt$organoid), c("t1 o1", "t2b o2", "t3a o3", "t3b o3", "t6 o6"))
})

test_that("purity edges never touch organoids or tumors without purity", {
  s <- toy_samples()
  Pu <- purityGraph(s, k_nn = 2)
  expect_equal(sum(Pu[s$role == "organoid", ]), 0)
  expect_equal(sum(Pu["t3b", ]), 0)
  expect_gt(sum(Pu), 0)
  expect_true(isSymmetric(Pu))
})

test_that("the annotation graph handles missing and constant annotations", {
  A <- annotGraph(toy_samples(), k_nn = 3)
  expect_false(anyNA(A))
  expect_true(isSymmetric(A))
  expect_equal(unname(diag(A)), rep(0, nrow(A)))
  expect_gt(sum(A), 0)
})

test_that("buildSampleGraph is symmetric, NA-free and scaled per part", {
  g <- buildSampleGraph(toy_samples(), k_nn = 3)
  expect_true(isSymmetric(g$S))
  expect_false(anyNA(g$S))
  expect_equal(unname(diag(g$S)), rep(0, 12))
  expect_equal(g$summary$part, c("pair", "annot", "purity"))
  expect_equal(g$summary$mean_degree, c(1, 1, 1))
})

test_that("the graph does not depend on row order", {
  s <- toy_samples()
  g1 <- buildSampleGraph(s, k_nn = 3)$S
  g2 <- buildSampleGraph(s[12:1, ], k_nn = 3)$S
  expect_equal(g2[rownames(g1), colnames(g1)], g1)
})

test_that("a graph built from a training subset contains only those samples", {
  s <- toy_samples()
  train <- s[!s$case_id %in% c("c2", "c6"), ]
  g <- buildSampleGraph(train, k_nn = 3)$S
  expect_equal(rownames(g), train$profile_id)
  expect_false(any(c("o2", "o6", "t2a") %in% rownames(g)))
})
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `make test FILTER=sample-graph`
Expected: FAIL with `could not find function "knnSparsify"`.

- [ ] **Step 4: Implement**

`R/sample_graph.R`:

```r
## Sample graph for GRACKLE: pairing, patient annotations and tumor purity (spec Section 6).

#' knnSparsify
#'
#' Keeps each row's `k_nn` largest positive similarities, then symmetrizes with `pmax`.
#' Ties are broken by `secondary` (smaller first, NA last) and then by row name, so the
#' result does not depend on row order.
#'
#' @param S Square similarity matrix with row and column names.
#' @param k_nn Neighbors kept per row.
#' @param secondary Optional square matrix of tie-break distances, aligned with `S`.
#' @return Symmetric matrix with a zero diagonal.
#' @export
knnSparsify <- function(S, k_nn, secondary = NULL) {
  n <- nrow(S)
  ids <- rownames(S)
  if (is.null(ids)) ids <- as.character(seq_len(n))
  id_rank <- rank(ids, ties.method = "first")
  pick <- lapply(seq_len(n), function(i) {
    s <- S[i, ]
    s[i] <- -Inf
    s[is.na(s)] <- -Inf
    tb <- if (is.null(secondary)) rep(0, n) else {
      x <- secondary[i, ]
      x[is.na(x)] <- Inf
      x
    }
    o <- order(-s, tb, id_rank)
    utils::head(o[s[o] > 0], k_nn)
  })
  out <- matrix(0, n, n, dimnames = dimnames(S))
  idx <- cbind(rep(seq_len(n), lengths(pick)), unlist(pick, use.names = FALSE))
  out[idx] <- S[idx]
  out <- pmax(out, t(out))
  diag(out) <- 0
  out
}

#' scaleMeanDegree
#'
#' Divides a similarity matrix by its mean degree over nodes with at least one edge.
#' @param S Square similarity matrix.
#' @return `S` scaled so that its mean nonzero degree is 1 (unchanged if it has no edges).
#' @export
scaleMeanDegree <- function(S) {
  deg <- rowSums(S)
  m <- mean(deg[deg > 0])
  if (is.na(m) || m == 0) return(S)
  S / m
}

#' pairTable
#'
#' Tumor-organoid pairs. An organoid pairs with the tumor whose `sample_id` equals its
#' `parent_sample_id`. If `parent_sample_id` is NA, it pairs with every tumor of its case.
#' @param samples Sample table (see [buildSampleGraph()]).
#' @return Data frame with columns `tumor` and `organoid` (profile ids).
#' @export
pairTable <- function(samples) {
  tum <- samples[samples$role %in% "tumor", c("profile_id", "case_id", "sample_id"), drop = FALSE]
  org <- samples[samples$role %in% "organoid", c("profile_id", "case_id", "parent_sample_id"), drop = FALSE]
  if (!nrow(tum) || !nrow(org)) return(data.frame(tumor = character(0), organoid = character(0)))
  e <- merge(org, tum, by = "case_id", suffixes = c("_org", "_tum"))
  e <- e[is.na(e$parent_sample_id) | e$parent_sample_id == e$sample_id, , drop = FALSE]
  data.frame(tumor = e$profile_id_tum, organoid = e$profile_id_org, stringsAsFactors = FALSE)
}

#' pairGraph
#'
#' Weight 1 between each organoid and its parent tumor(s) (see [pairTable()]).
#' @param samples Sample table.
#' @return Square matrix over `samples$profile_id`.
#' @export
pairGraph <- function(samples) {
  ids <- samples$profile_id
  S <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  pt <- pairTable(samples)
  if (!nrow(pt)) return(S)
  idx <- cbind(match(pt$organoid, ids), match(pt$tumor, ids))
  S[idx] <- 1
  S[idx[, 2:1, drop = FALSE]] <- 1
  S
}

#' annotGraph
#'
#' `1 - Gower distance` over patient annotations, sparsified to `k_nn` neighbors. Columns
#' that are constant or entirely missing are dropped; missing values are skipped by Gower.
#' Ties are broken by absolute age difference, then by profile id.
#' @param samples Sample table.
#' @param vars Annotation columns to use (those absent from `samples` are ignored).
#' @param k_nn Neighbors kept per sample.
#' @return Square matrix over `samples$profile_id`.
#' @export
annotGraph <- function(samples,
                       vars = c("cancer_type", "tumor_descriptor", "prior_treatment", "age", "sex", "ancestry"),
                       k_nn = 10) {
  ids <- samples$profile_id
  empty <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  df <- samples[, intersect(vars, names(samples)), drop = FALSE]
  df[] <- lapply(df, function(x) if (is.numeric(x)) x else factor(x))
  informative <- vapply(df, function(x) length(unique(stats::na.omit(x))) > 1, logical(1))
  df <- df[, informative, drop = FALSE]
  if (ncol(df) == 0) return(empty)
  S <- 1 - as.matrix(cluster::daisy(df, metric = "gower"))
  S[is.na(S)] <- 0
  dimnames(S) <- list(ids, ids)
  age_gap <- if ("age" %in% names(samples)) abs(outer(samples$age, samples$age, "-")) else NULL
  knnSparsify(S, k_nn, secondary = age_gap)
}

#' purityGraph
#'
#' Gaussian-kernel similarity of DNA-based purity between tumors only, sparsified to `k_nn`
#' neighbors. The bandwidth is the median distance from each tumor to its `k_nn`-th nearest
#' purity neighbor. Organoids and tumors with NA purity get no edges.
#' @param samples Sample table with `role` and `purity_dna`.
#' @param k_nn Neighbors kept per tumor.
#' @return Square matrix over `samples$profile_id`.
#' @export
purityGraph <- function(samples, k_nn = 10) {
  ids <- samples$profile_id
  S <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  t_idx <- which(samples$role %in% "tumor" & !is.na(samples$purity_dna))
  if (length(t_idx) < 2) return(S)
  p <- samples$purity_dna[t_idx]
  D <- abs(outer(p, p, "-"))
  kth <- apply(D, 1, function(r) sort(r)[min(k_nn + 1, length(r))])  # +1 skips the zero self-distance
  sigma <- stats::median(kth)
  if (sigma <= 0) sigma <- min(D[D > 0], 1)  # every purity tied
  K <- exp(-D^2 / (2 * sigma^2))
  dimnames(K) <- list(ids[t_idx], ids[t_idx])
  S[t_idx, t_idx] <- knnSparsify(K, k_nn)
  S
}

#' buildSampleGraph
#'
#' `S = w_pair * S_pair + w_annot * S_annot + w_pur * S_pur`, with each part first divided by
#' its mean degree (spec Section 6). Build it from training samples only.
#' @param samples Sample table with columns `profile_id, case_id, sample_id, parent_sample_id,
#'   role, cancer_type, tumor_descriptor, prior_treatment, age, sex, ancestry, purity_dna`.
#' @param weights Named weights for the parts `pair`, `annot`, `purity`.
#' @param k_nn Neighbors kept per sample in the annotation and purity parts.
#' @return A list with `S` (square matrix over `samples$profile_id`) and `summary`
#'   (columns part, n_edges, mean_degree, sparsity, weight).
#' @export
buildSampleGraph <- function(samples, weights = c(pair = 1, annot = 1, purity = 1), k_nn = 10) {
  parts <- list(pair = pairGraph(samples),
                annot = annotGraph(samples, k_nn = k_nn),
                purity = purityGraph(samples, k_nn = k_nn))
  scaled <- lapply(parts, scaleMeanDegree)
  S <- Reduce(`+`, Map(function(m, w) w * m, scaled, weights[names(scaled)]))
  diag(S) <- 0
  mean_deg <- function(m) {
    d <- rowSums(m)
    if (any(d > 0)) mean(d[d > 0]) else 0
  }
  summary <- data.frame(part = names(parts),
                        n_edges = vapply(parts, function(m) sum(m[upper.tri(m)] > 0), numeric(1)),
                        mean_degree = vapply(scaled, mean_deg, numeric(1)),
                        sparsity = vapply(parts, function(m) mean(m[upper.tri(m)] == 0), numeric(1)),
                        weight = unname(weights[names(parts)]),
                        row.names = NULL, stringsAsFactors = FALSE)
  list(S = S, summary = summary)
}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `make document && make test FILTER=sample-graph`
Expected: PASS (8 tests).

- [ ] **Step 6: Commit**

```bash
git add R/sample_graph.R NAMESPACE man tests/testthat/helper-samples.R tests/testthat/test-sample-graph.R
git commit -m "feat: sample graph from pairing, patient annotations and tumor purity

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: STRING gene graphs

**Files:**
- Create: `R/gene_graph.R` (STRING section), `tests/testthat/test-gene-graph-string.R`, `tests/testthat/fixtures/string_links_full_sample.txt`, `tests/testthat/fixtures/README.md`

**Interfaces:**
- Produces:
  - `readStringLinks(path)`: a data frame with STRING's full-links columns.
  - `readStringAliases(path)`: columns `string_protein_id, alias, source`.
  - `stringProteinToSymbol(aliases, hgnc)`: columns `string_protein_id, symbol`. `hgnc` needs columns `symbol, ensembl_gene_id`.
  - `stringRecombine(links, channels = c("experiments", "database"), prior = 0.041)`: numeric in [0, 1], and 0 where the selected channels have no evidence.
  - `stringCombinedAll(links)`: an integer in 0..1000 that reproduces `combined_score`.
  - `stringGraph(links, protein_map, genes, mode = c("mech", "all"))`: genes x genes, symmetric, zero diagonal.
  - `graphSummary(G)`: columns `n_genes_with_edges, n_edges, mean_degree`.

- [ ] **Step 1: Extract a real-data fixture from the local STRING file**

Run:

```bash
mkdir -p tests/testthat/fixtures
awk 'NR==1 || NR % 5000 == 0' GRACKLE_data/9606.protein.links.full.v12.0.txt > tests/testthat/fixtures/string_links_full_sample.txt
wc -l tests/testthat/fixtures/string_links_full_sample.txt
```

Expected: `2744` lines (a header plus every 5,000th edge).

Write `tests/testthat/fixtures/README.md`:

```markdown
# Test fixtures

- `string_links_full_sample.txt`: every 5,000th row (plus the header) of STRING v12
  `9606.protein.links.full.v12.0.txt`. STRING is distributed under CC BY 4.0
  (Szklarczyk et al., Nucleic Acids Res 2023). Used to check that `stringCombinedAll()`
  reproduces STRING's `combined_score`.
```

- [ ] **Step 2: Write the failing tests**

`tests/testthat/test-gene-graph-string.R`:

```r
toy_links <- function() {
  data.frame(protein1 = c("P1", "P1", "P2", "P3", "P1", "P5"),
             protein2 = c("P2", "P3", "P1", "P4", "P4", "P1"),
             neighborhood = 0, neighborhood_transferred = 0, fusion = 0, cooccurence = 0, homology = 0,
             coexpression = c(0, 500, 0, 0, 900, 0), coexpression_transferred = 0,
             experiments = c(400, 0, 400, 300, 0, 800), experiments_transferred = 0,
             database = c(0, 0, 0, 900, 0, 0), database_transferred = 0,
             textmining = c(200, 300, 200, 0, 0, 0), textmining_transferred = 0,
             combined_score = c(500, 600, 500, 930, 900, 800))
}
## P2 and P3 both map to B; P4 maps to C and C2; P5 maps to nothing.
toy_map <- data.frame(string_protein_id = c("P1", "P2", "P3", "P4", "P4"),
                      symbol = c("A", "B", "B", "C", "C2"), stringsAsFactors = FALSE)

test_that("stringCombinedAll reproduces STRING's combined_score", {
  x <- readStringLinks(test_path("fixtures", "string_links_full_sample.txt"))
  expect_gt(nrow(x), 2000)
  expect_lte(max(abs(stringCombinedAll(x) - x$combined_score)), 1)
})

test_that("the mechanistic score never exceeds combined_score and is 0 without evidence", {
  x <- readStringLinks(test_path("fixtures", "string_links_full_sample.txt"))
  m <- stringRecombine(x)
  expect_true(all(m * 1000 <= x$combined_score + 1))
  none <- x$experiments + x$experiments_transferred + x$database + x$database_transferred == 0
  expect_true(all(m[none] == 0))
})

test_that("stringRecombine follows STRING's prior-corrected combination", {
  m <- stringRecombine(toy_links())
  e <- (0.3 - 0.041) / 0.959; d <- (0.9 - 0.041) / 0.959
  expect_equal(m[4], (1 - (1 - e) * (1 - d)) * 0.959 + 0.041)
  expect_equal(m[1], 0.4)
  expect_equal(m[2], 0)
})

test_that("stringGraph maps proteins to genes, aggregates by max, and keeps unmapped genes empty", {
  genes <- c("A", "B", "C", "D")
  mech <- stringGraph(toy_links(), toy_map, genes, "mech")
  expect_equal(dim(mech), c(4, 4))
  expect_equal(mech["A", "B"], 0.4)
  expect_equal(mech["B", "C"], (1 - (1 - (0.3 - 0.041) / 0.959) * (1 - (0.9 - 0.041) / 0.959)) * 0.959 + 0.041)
  expect_equal(mech["A", "C"], 0)
  expect_equal(sum(mech["D", ]), 0)
  expect_true(isSymmetric(mech)); expect_false(anyNA(mech)); expect_equal(unname(diag(mech)), rep(0, 4))
  all_ch <- stringGraph(toy_links(), toy_map, genes, "all")
  expect_equal(all_ch["A", "B"], 0.6)
  expect_equal(all_ch["A", "C"], 0.9)
})

test_that("stringProteinToSymbol maps through Ensembl gene ids", {
  aliases <- data.frame(string_protein_id = c("9606.P1", "9606.P1", "9606.P2"),
                        alias = c("ENSG01", "TP53X", "ENSG02"),
                        source = c("Ensembl_gene", "Ensembl_HGNC_symbol", "Ensembl_gene"))
  hgnc <- data.frame(symbol = c("TP53", "KRAS"), ensembl_gene_id = c("ENSG01", "ENSG03"))
  m <- stringProteinToSymbol(aliases, hgnc)
  expect_equal(m$string_protein_id, "9606.P1")
  expect_equal(m$symbol, "TP53")
})

test_that("graphSummary counts genes, edges and mean degree", {
  G <- matrix(0, 3, 3); G[1, 2] <- G[2, 1] <- 0.5
  s <- graphSummary(G)
  expect_equal(s$n_genes_with_edges, 2); expect_equal(s$n_edges, 1); expect_equal(s$mean_degree, 0.5)
})
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `make test FILTER=gene-graph-string`
Expected: FAIL with `could not find function "readStringLinks"`.

- [ ] **Step 4: Implement**

`R/gene_graph.R` (create it with this STRING section; Task 7 appends to it):

```r
## Gene graphs for GRACKLE (spec Section 4.5): STRING arms and co-essentiality arms.

.datatable.aware <- TRUE
.string_prior <- 0.041

#' readStringLinks
#' @param path STRING `protein.links.full` file (plain or .gz).
#' @return Data frame with protein1, protein2, the channel columns and combined_score.
#' @export
readStringLinks <- function(path) {
  data.table::fread(path, sep = " ", header = TRUE, data.table = FALSE)
}

#' readStringAliases
#' @param path STRING `protein.aliases` file (plain or .gz).
#' @return Data frame with columns string_protein_id, alias, source.
#' @export
readStringAliases <- function(path) {
  x <- data.table::fread(path, sep = "\t", header = TRUE, quote = "", data.table = FALSE)
  names(x) <- c("string_protein_id", "alias", "source")
  x
}

#' stringProteinToSymbol
#'
#' Maps STRING proteins to HGNC symbols through their Ensembl gene ids, using the pinned HGNC table.
#' @param aliases Output of [readStringAliases()].
#' @param hgnc Data frame with columns `symbol` and `ensembl_gene_id`.
#' @return Unique data frame with columns string_protein_id and symbol.
#' @export
stringProteinToSymbol <- function(aliases, hgnc) {
  a <- aliases[aliases$source == "Ensembl_gene", c("string_protein_id", "alias")]
  a$symbol <- hgnc$symbol[match(a$alias, hgnc$ensembl_gene_id)]
  out <- unique(a[!is.na(a$symbol), c("string_protein_id", "symbol")])
  rownames(out) <- NULL
  out
}

priorAway <- function(score, prior = .string_prior) {
  (pmax(score, prior) - prior) / (1 - prior)
}

#' stringRecombine
#'
#' Recombines selected STRING evidence channels with STRING's procedure: remove the prior
#' from each channel, combine direct and transferred scores per channel, combine channels as
#' `1 - prod(1 - s)`, then add the prior back once.
#' @param links Data frame from [readStringLinks()] (channels on the 0-1000 scale).
#' @param channels Channel names; each `<channel>_transferred` column is included when present.
#' @param prior STRING's prior. Default 0.041.
#' @return Numeric vector in [0, 1]; 0 where no selected channel has evidence.
#' @export
stringRecombine <- function(links, channels = c("experiments", "database"), prior = .string_prior) {
  per_channel <- vapply(channels, function(ch) {
    direct <- priorAway(links[[ch]] / 1000, prior)
    tr <- paste0(ch, "_transferred")
    transferred <- if (tr %in% names(links)) priorAway(links[[tr]] / 1000, prior) else 0
    1 - (1 - direct) * (1 - transferred)
  }, numeric(nrow(links)))
  per_channel <- matrix(per_channel, nrow = nrow(links))
  combined <- (1 - matrixStats::rowProds(1 - per_channel)) * (1 - prior) + prior
  ifelse(rowSums(per_channel > 0) > 0, combined, 0)
}

#' stringCombinedAll
#'
#' Recomputes STRING's combined_score from all channels, including the homology correction
#' of the co-occurrence and text-mining channels. Used to validate [stringRecombine()].
#' @param links Data frame from [readStringLinks()].
#' @param prior STRING's prior. Default 0.041.
#' @return Integer scores on the 0-1000 scale.
#' @export
stringCombinedAll <- function(links, prior = .string_prior) {
  pa <- function(col) priorAway(links[[col]] / 1000, prior)
  both <- function(a, b) 1 - (1 - pa(a)) * (1 - pa(b))
  homology <- links$homology / 1000
  parts <- cbind(both("neighborhood", "neighborhood_transferred"),
                 pa("fusion"),
                 pa("cooccurence") * (1 - homology),
                 both("coexpression", "coexpression_transferred"),
                 both("experiments", "experiments_transferred"),
                 both("database", "database_transferred"),
                 both("textmining", "textmining_transferred") * (1 - homology))
  combined <- (1 - matrixStats::rowProds(1 - parts)) * (1 - prior) + prior
  floor(combined * 1000)
}

#' stringGraph
#'
#' Gene-by-gene STRING graph. `mode = "mech"` uses the experiments and database channels
#' recombined by [stringRecombine()]; `mode = "all"` uses combined_score / 1000. Proteins
#' that map to the same gene pair are aggregated by max. Genes without edges have zero rows.
#' @param links Data frame from [readStringLinks()].
#' @param protein_map Output of [stringProteinToSymbol()].
#' @param genes Character vector of HGNC symbols giving the graph's rows and columns.
#' @param mode "mech" or "all".
#' @return Symmetric genes x genes matrix with a zero diagonal.
#' @export
stringGraph <- function(links, protein_map, genes, mode = c("mech", "all")) {
  mode <- match.arg(mode)
  map <- protein_map[protein_map$symbol %in% genes, , drop = FALSE]
  keep <- links$protein1 %in% map$string_protein_id & links$protein2 %in% map$string_protein_id
  links <- links[keep, , drop = FALSE]
  G <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  if (!nrow(links)) return(G)
  w <- if (mode == "mech") stringRecombine(links) else links$combined_score / 1000
  e <- data.table::data.table(p1 = links$protein1, p2 = links$protein2, w = w)
  e <- e[e$w > 0]
  m <- data.table::data.table(string_protein_id = map$string_protein_id, symbol = map$symbol)
  e <- merge(e, m, by.x = "p1", by.y = "string_protein_id", allow.cartesian = TRUE)
  data.table::setnames(e, "symbol", "g1")
  e <- merge(e, m, by.x = "p2", by.y = "string_protein_id", allow.cartesian = TRUE)
  data.table::setnames(e, "symbol", "g2")
  e <- e[e$g1 != e$g2]
  e <- e[, list(w = max(w)), by = c("g1", "g2")]
  G[cbind(match(e$g1, genes), match(e$g2, genes))] <- e$w
  G <- pmax(G, t(G))
  diag(G) <- 0
  G
}

#' graphSummary
#' @param G Symmetric gene graph.
#' @return One-row data frame: n_genes_with_edges, n_edges, mean_degree (over genes with edges).
#' @export
graphSummary <- function(G) {
  deg <- rowSums(G)
  data.frame(n_genes_with_edges = sum(deg > 0),
             n_edges = sum(G[upper.tri(G)] > 0),
             mean_degree = if (any(deg > 0)) mean(deg[deg > 0]) else 0)
}

utils::globalVariables("w")
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `make document && make test FILTER=gene-graph-string`
Expected: PASS (6 tests).

If `stringCombinedAll` differs from `combined_score` by more than 1 on some rows, print the worst rows (`x[order(-abs(stringCombinedAll(x) - x$combined_score))[1:5], ]`). Compare them with STRING's published `combine_subscores.py` from the STRING download page and fix the formula. Do not loosen the test tolerance.

- [ ] **Step 6: Commit**

```bash
git add R/gene_graph.R NAMESPACE man tests/testthat/test-gene-graph-string.R tests/testthat/fixtures/string_links_full_sample.txt tests/testthat/fixtures/README.md
git commit -m "feat: STRING mechanistic and all-channel gene graphs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Co-essentiality gene graphs

**Files:**
- Modify: `R/gene_graph.R` (append the co-essentiality section)
- Create: `tests/testthat/test-gene-graph-coess.R`

**Interfaces:**
- Produces:
  - `geneArm(location)`: e.g. "17q", or NA.
  - `hcmiToken(x)`: the 17-character case id `HCM-XXXX-0000-C00`, or NA.
  - `informativeGenes(effect, nonessential, quantile = 0.95)`: a character vector.
  - `glsWhiten(effect, ridge = 1e-8)`: genes x models with attribute `"intercept"`.
  - `glsCorrelation(Ew)`: genes x genes.
  - `coessGraph(effect, genes, nonessential, arm, k_nn = 10)`: genes x genes.
  - `dropHeldOutModels(effect, model_case, heldout_cases)`.
- `effect` is always genes x models (rows are HGNC symbols, columns are model ids).

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-gene-graph-coess.R`:

```r
toy_effect <- function(seed = 1, n_models = 30) {
  withr::with_seed(seed, {
    driver <- stats::rnorm(n_models)
    module <- t(vapply(1:5, function(i) driver + stats::rnorm(n_models, sd = 0.3), numeric(n_models)))
    noise <- matrix(stats::rnorm(12 * n_models, sd = 0.05), 12)
    other <- matrix(stats::rnorm(183 * n_models), 183)
  })
  E <- rbind(module, noise, other)
  dimnames(E) <- list(c(paste0("G", 1:5), paste0("N", 1:12), paste0("V", 1:183)), paste0("M", seq_len(n_models)))
  E
}

test_that("geneArm parses HGNC cytogenetic locations", {
  expect_equal(geneArm(c("17q21.31", "1p36.33", "Xp22.33 and Yp11.2", "mitochondria", "13cen-q12", NA)),
               c("17q", "1p", "Xp", NA, NA, NA))
})

test_that("hcmiToken extracts the case-level HCMI id", {
  expect_equal(hcmiToken(c("HCM-BROD-0577-C15-85M", "x HCM-CSHL-0143-C20 y", "ACH-000001", NA)),
               c("HCM-BROD-0577-C15", "HCM-CSHL-0143-C20", NA, NA))
})

test_that("informativeGenes keeps variable genes and drops noise-level genes", {
  E <- toy_effect()
  inf <- informativeGenes(E, nonessential = paste0("N", 1:12))
  expect_true(all(paste0("G", 1:5) %in% inf))
  expect_false(any(paste0("N", 1:12) %in% inf))
  expect_error(informativeGenes(E, nonessential = "N1"), "fewer than 10")
})

test_that("glsWhiten makes the model covariance the identity", {
  withr::with_seed(2, {
    A <- matrix(stats::rnorm(20 * 20), 20)
    E <- matrix(stats::rnorm(500 * 20), 500) %*% A
  })
  Ew <- glsWhiten(E)
  expect_lt(max(abs(stats::cov(Ew) - diag(20))), 1e-6)
  expect_length(attr(Ew, "intercept"), 20)
})

test_that("GLS correlation removes each gene's mean effect", {
  withr::with_seed(3, {
    A <- matrix(stats::rnorm(25 * 25), 25)                 # correlated models
    E <- matrix(stats::rnorm(300 * 25), 300) %*% A
  })
  E[1:2, ] <- E[1:2, ] - 3   # two strongly essential genes whose profiles are otherwise independent
  rownames(E) <- paste0("g", 1:300)
  C <- glsCorrelation(glsWhiten(E))
  expect_lt(abs(C["g1", "g2"]), 0.5)
})

test_that("coessGraph links module genes, skips same-arm pairs and leaves other genes empty", {
  E <- toy_effect()
  arm <- c(G1 = "1p", G2 = "1p", G3 = "2q", G4 = "3q", G5 = "4p")
  genes <- c(paste0("G", 1:5), "N1", "V1", "Z")
  G <- coessGraph(E, genes, nonessential = paste0("N", 1:12), arm = arm, k_nn = 3)
  expect_equal(dim(G), c(8, 8))
  expect_true(isSymmetric(G)); expect_false(anyNA(G)); expect_true(all(G >= 0))
  expect_equal(unname(diag(G)), rep(0, 8))
  expect_equal(G["G1", "G2"], 0)
  expect_gt(G["G3", "G4"], 0)
  expect_equal(sum(G["N1", ]), 0)
  expect_equal(sum(G["Z", ]), 0)
})

test_that("dropHeldOutModels removes only the held-out cases' models", {
  E <- toy_effect()
  model_case <- stats::setNames(c("caseA", "caseB", rep(NA, 28)), colnames(E))
  out <- dropHeldOutModels(E, model_case, heldout_cases = "caseA")
  expect_false("M1" %in% colnames(out))
  expect_true(all(c("M2", "M3") %in% colnames(out)))
  expect_equal(ncol(out), 29)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=gene-graph-coess`
Expected: FAIL with `could not find function "geneArm"`.

- [ ] **Step 3: Implement (append to `R/gene_graph.R`)**

```r
## Co-essentiality arms ---------------------------------------------------------------

#' geneArm
#' @param location HGNC cytogenetic locations, for example "17q21.31".
#' @return Chromosome arm ("17q"), or NA when the location has no single arm.
#' @export
geneArm <- function(location) {
  m <- regmatches(location, regexec("^([0-9]{1,2}|X|Y)([pq])", location))
  vapply(m, function(x) if (length(x) == 3) paste0(x[2], x[3]) else NA_character_, character(1))
}

.hcmi_pattern <- "HCM-[A-Z]{4}-[0-9]{4}-[A-Z][0-9]{2}"

#' hcmiToken
#' @param x Character vector that may contain an HCMI identifier.
#' @return The case-level HCMI id (for example "HCM-BROD-0577-C15"), or NA.
#' @export
hcmiToken <- function(x) {
  out <- rep(NA_character_, length(x))
  hit <- !is.na(x) & grepl(.hcmi_pattern, x)
  out[hit] <- regmatches(x[hit], regexpr(.hcmi_pattern, x[hit]))
  out
}

#' informativeGenes
#'
#' Genes whose dependency varies across models: SD across models above the `quantile` of the
#' SDs of DepMap's non-essential control genes.
#' @param effect Genes x models matrix of gene effects (Chronos).
#' @param nonessential HGNC symbols of the non-essential controls.
#' @param quantile Control quantile used as the threshold. Default 0.95.
#' @return Character vector of informative genes.
#' @export
informativeGenes <- function(effect, nonessential, quantile = 0.95) {
  sds <- matrixStats::rowSds(effect, na.rm = TRUE)
  names(sds) <- rownames(effect)
  ctrl <- sds[intersect(nonessential, names(sds))]
  if (length(ctrl) < 10) stop("informativeGenes: fewer than 10 non-essential controls found in the effect matrix")
  thr <- stats::quantile(ctrl, quantile, na.rm = TRUE)
  names(sds)[!is.na(sds) & sds > thr]
}

#' glsWhiten
#'
#' Whitens the model dimension by `Sigma^(-1/2)`, where `Sigma` is the model-by-model
#' covariance estimated across all genes (Wainberg et al. 2021). Missing values are replaced
#' by the gene's mean. Eigenvalues are floored at `ridge * max`.
#' @param effect Genes x models matrix.
#' @param ridge Relative eigenvalue floor. Default 1e-8.
#' @return Whitened genes x models matrix with attribute `"intercept"`: the whitened
#'   all-ones vector (`colSums(Sigma^(-1/2))`), needed by [glsCorrelation()].
#' @export
glsWhiten <- function(effect, ridge = 1e-8) {
  E <- effect
  na <- is.na(E)
  if (any(na)) E[na] <- rowMeans(E, na.rm = TRUE)[row(E)[na]]
  ev <- eigen(stats::cov(E), symmetric = TRUE)
  vals <- pmax(ev$values, ridge * max(ev$values))
  W_half <- ev$vectors %*% (t(ev$vectors) / sqrt(vals))
  out <- E %*% W_half
  dimnames(out) <- dimnames(effect)
  attr(out, "intercept") <- colSums(W_half)
  out
}

#' glsCorrelation
#'
#' GLS correlation between genes: the cosine of whitened profiles after projecting out the
#' whitened intercept. This is the GLS analogue of a Pearson correlation.
#' @param Ew Output of [glsWhiten()] (optionally subset to some genes; keep the attribute).
#' @return Genes x genes similarity matrix.
#' @export
glsCorrelation <- function(Ew) {
  v <- attr(Ew, "intercept")
  R <- Ew - outer(as.vector(Ew %*% v) / sum(v^2), v)
  nr <- sqrt(rowSums(R^2))
  nr[nr == 0] <- 1
  C <- tcrossprod(R / nr)
  dimnames(C) <- list(rownames(Ew), rownames(Ew))
  C
}

#' coessGraph
#'
#' Co-essentiality gene graph (spec Section 4.5): informative genes only, GLS correlation,
#' same-arm pairs and negative values set to 0, `k_nn` nearest neighbors among `genes`.
#' @param effect Genes x models gene-effect matrix for the arm's models.
#' @param genes Final gene set giving the graph's rows and columns.
#' @param nonessential Non-essential control genes.
#' @param arm Named character vector: gene symbol to chromosome arm.
#' @param k_nn Neighbors kept per gene.
#' @return Symmetric genes x genes matrix with a zero diagonal.
#' @export
coessGraph <- function(effect, genes, nonessential, arm, k_nn = 10) {
  G <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  use <- intersect(genes, informativeGenes(effect, nonessential))
  if (length(use) < 2) return(G)
  Ew <- glsWhiten(effect)
  sub <- Ew[use, , drop = FALSE]
  attr(sub, "intercept") <- attr(Ew, "intercept")
  C <- glsCorrelation(sub)
  a <- unname(arm[use])
  same_arm <- outer(a, a, "==")
  same_arm[is.na(same_arm)] <- FALSE
  C[same_arm] <- 0
  C[is.na(C) | C < 0] <- 0
  diag(C) <- 0
  G[use, use] <- knnSparsify(C, k_nn)
  G
}

#' dropHeldOutModels
#' @param effect Genes x models matrix.
#' @param model_case Named character vector: model id (column of `effect`) to HCMI case id (NA if not HCMI).
#' @param heldout_cases Case ids held out in this fold.
#' @return `effect` without the columns of held-out cases' models.
#' @export
dropHeldOutModels <- function(effect, model_case, heldout_cases) {
  drop <- names(model_case)[!is.na(model_case) & model_case %in% heldout_cases]
  effect[, setdiff(colnames(effect), drop), drop = FALSE]
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `make document && make test FILTER=gene-graph-coess`
Expected: PASS (7 tests).

- [ ] **Step 5: Commit**

```bash
git add R/gene_graph.R NAMESPACE man tests/testthat/test-gene-graph-coess.R
git commit -m "feat: co-essentiality gene graphs with GLS correlation

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Rank by parallel analysis

**Files:**
- Create: `R/rank_selection.R`, `tests/testthat/test-rank.R`

**Interfaces:**
- Produces: `parallelAnalysisRank(Y, n_perm = 100, quantile = 0.95, nv = 200, seed = 20261015)`. It returns `list(k, observed, threshold, null)`, where `k` is an integer: the number of leading components whose singular value exceeds the null quantile at the same index.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-rank.R`:

```r
planted <- function(r = 3, n = 60, g = 200, seed = 4) {
  withr::with_seed(seed, {
    U <- matrix(stats::runif(n * r), n)
    V <- matrix(0, r, g)
    blocks <- split(seq_len(r * 50), rep(seq_len(r), each = 50))
    V[cbind(rep(seq_len(r), each = 50), unlist(blocks))] <- 5
    Y <- pmax(U %*% V + matrix(stats::rnorm(n * g, sd = 0.3), n), 0)
  })
  Y
}

test_that("parallel analysis recovers a planted rank", {
  pa <- parallelAnalysisRank(planted(3), n_perm = 30, nv = 20)
  expect_equal(pa$k, 3L)
  expect_length(pa$observed, 20)
  expect_length(pa$threshold, 20)
  expect_equal(dim(pa$null), c(30, 20))
})

test_that("parallel analysis is reproducible with the seed", {
  a <- parallelAnalysisRank(planted(2), n_perm = 10, nv = 10)
  b <- parallelAnalysisRank(planted(2), n_perm = 10, nv = 10)
  expect_identical(a$threshold, b$threshold)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=rank`
Expected: FAIL with `could not find function "parallelAnalysisRank"`.

- [ ] **Step 3: Implement**

`R/rank_selection.R`:

```r
topSingularValues <- function(A, nv) {
  if (nv >= min(dim(A)) / 2) return(svd(A, nu = 0, nv = 0)$d[seq_len(nv)])
  irlba::irlba(A, nv = nv)$d
}

#' parallelAnalysisRank
#'
#' Horn's parallel analysis (spec Section 7.1). Null matrices permute each gene (column)
#' independently. k is the number of leading components whose observed singular value exceeds
#' the null `quantile` at the same index. If k reaches `nv`, the analysis reruns with `2 * nv`.
#' @param Y Samples x genes matrix (uncentered, as factorized).
#' @param n_perm Number of null permutations. Default 100.
#' @param quantile Null quantile. Default 0.95.
#' @param nv Number of singular values computed. Default 200.
#' @param seed Random seed. Default 20261015.
#' @return A list with `k`, `observed`, `threshold` and `null` (n_perm x nv).
#' @export
parallelAnalysisRank <- function(Y, n_perm = 100, quantile = 0.95, nv = 200, seed = 20261015) {
  nv <- min(nv, min(dim(Y)) - 1)
  withr::with_seed(seed, {
    observed <- topSingularValues(Y, nv)
    null <- t(vapply(seq_len(n_perm), function(b) topSingularValues(apply(Y, 2, sample), nv), numeric(nv)))
  })
  threshold <- apply(null, 2, stats::quantile, probs = quantile)
  above <- observed > threshold
  k <- if (all(above)) nv else which(!above)[1] - 1L
  if (k == nv && nv < min(dim(Y)) - 1) return(parallelAnalysisRank(Y, n_perm, quantile, 2 * nv, seed))
  list(k = as.integer(k), observed = observed, threshold = unname(threshold), null = null)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `make document && make test FILTER=rank`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add R/rank_selection.R NAMESPACE man/parallelAnalysisRank.Rd tests/testthat/test-rank.R
git commit -m "feat: rank k from Horn's parallel analysis

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Lambda scaling and setting selection

**Files:**
- Create: `R/tuning.R`, `tests/testthat/test-tuning.R`

**Interfaces:**
- Produces:
  - `traceLaplacian(X, S, by = c("row", "col"))`. With `by = "row"` it returns `tr(X' L X)` for X with samples in rows; with `by = "col"` it returns `tr(X L X')` for X with genes in columns.
  - `lambdaFromAlpha(Y, W0, H0, S_p, S_g, alpha_1, alpha_2)` returns `c(lambda_1 = , lambda_2 = )`.
  - `selectSetting(tab)` returns the `setting_id` it selects. `tab` needs columns `setting_id, alpha_1, alpha_2, all_converged, recon_mean, recon_se, reproducibility`.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-tuning.R`:

```r
test_that("traceLaplacian equals the explicit trace", {
  withr::with_seed(5, {
    S <- matrix(stats::runif(36), 6); S <- (S + t(S)) / 2; diag(S) <- 0
    W <- matrix(stats::runif(18), 6); H <- matrix(stats::runif(18), 3)
  })
  L <- diag(rowSums(S)) - S
  expect_equal(traceLaplacian(W, S, "row"), sum(diag(t(W) %*% L %*% W)))
  expect_equal(traceLaplacian(H, S, "col"), sum(diag(H %*% L %*% t(H))))
})

test_that("lambdaFromAlpha makes each penalty alpha times the starting reconstruction", {
  withr::with_seed(6, {
    Y <- matrix(stats::runif(60), 6)
    S_p <- matrix(stats::runif(36), 6); S_p <- (S_p + t(S_p)) / 2; diag(S_p) <- 0
    S_g <- matrix(stats::runif(100), 10); S_g <- (S_g + t(S_g)) / 2; diag(S_g) <- 0
  })
  init <- svdInit(Y, 2)
  lam <- lambdaFromAlpha(Y, init$W, init$H, S_p, S_g, alpha_1 = 0.1, alpha_2 = 0.3)
  recon <- sum((Y - init$W %*% init$H)^2)
  expect_equal(lam[["lambda_1"]] * traceLaplacian(init$W, S_p, "row"), 0.1 * recon)
  expect_equal(lam[["lambda_2"]] * traceLaplacian(init$H, S_g, "col"), 0.3 * recon)
  zero <- lambdaFromAlpha(Y, init$W, init$H, S_p * 0, S_g, 0.1, 0.3)
  expect_equal(zero[["lambda_1"]], 0)
})

test_that("selectSetting applies the one-SE rule, then reproducibility, then the tie-break", {
  tab <- data.frame(
    setting_id = c("a", "b", "c", "d", "e"),
    alpha_1 = c(0, 0.1, 0.3, 0.1, 0.01), alpha_2 = c(0, 0.1, 0.3, 0.3, 0.1),
    all_converged = c(TRUE, TRUE, TRUE, FALSE, TRUE),
    recon_mean = c(0.100, 0.104, 0.120, 0.090, 0.104),
    recon_se = c(0.005, 0.005, 0.005, 0.005, 0.005),
    reproducibility = c(0.80, 0.90, 0.99, 0.99, 0.90))
  ## d is best but not converged; the band is recon <= 0.105 -> a, b, e; b and e tie on
  ## reproducibility; e has the smaller alpha_1.
  expect_equal(selectSetting(tab), "e")
  expect_error(selectSetting(transform(tab, all_converged = FALSE)), "no setting")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=tuning`
Expected: FAIL with `could not find function "traceLaplacian"`.

- [ ] **Step 3: Implement**

`R/tuning.R`:

```r
#' traceLaplacian
#'
#' `tr(X' L X)` (by = "row", X is samples x k) or `tr(X L X')` (by = "col", X is k x genes)
#' for `L = diag(rowSums(S)) - S`, without forming `L`.
#' @param X Factor matrix.
#' @param S Symmetric similarity matrix.
#' @param by "row" or "col".
#' @return Numeric scalar.
#' @export
traceLaplacian <- function(X, S, by = c("row", "col")) {
  by <- match.arg(by)
  d <- rowSums(S)
  if (by == "row") sum(d * rowSums(X^2)) - sum(X * (S %*% X)) else sum(d * colSums(X^2)) - sum(X * (X %*% S))
}

#' lambdaFromAlpha
#'
#' Dimensionless penalty weights (spec Section 7.2): each penalty equals `alpha` times the
#' reconstruction error at the starting point. A graph with no edges gives lambda = 0.
#' @param Y Training matrix.
#' @param W0,H0 Starting factors, for example from [svdInit()].
#' @param S_p Sample similarity. @param S_g Gene similarity.
#' @param alpha_1,alpha_2 Penalty shares.
#' @return Named numeric vector with `lambda_1` and `lambda_2`.
#' @export
lambdaFromAlpha <- function(Y, W0, H0, S_p, S_g, alpha_1, alpha_2) {
  recon <- sum((Y - W0 %*% H0)^2)
  g1 <- traceLaplacian(W0, S_p, "row")
  g2 <- traceLaplacian(H0, S_g, "col")
  c(lambda_1 = if (g1 > 0) alpha_1 * recon / g1 else 0,
    lambda_2 = if (g2 > 0) alpha_2 * recon / g2 else 0)
}

#' selectSetting
#'
#' Spec Section 7.3. Settings with a non-converged or missing fold are excluded. Among
#' settings with held-out reconstruction within one SE of the best, the most reproducible
#' wins; ties go to the smaller alpha_1, then the smaller alpha_2.
#' @param tab Data frame with setting_id, alpha_1, alpha_2, all_converged, recon_mean, recon_se, reproducibility.
#' @return The selected setting_id.
#' @export
selectSetting <- function(tab) {
  ok <- tab[tab$all_converged %in% TRUE, , drop = FALSE]
  if (!nrow(ok)) stop("selectSetting: no setting has all folds present and converged")
  best <- which.min(ok$recon_mean)
  band <- ok[ok$recon_mean <= ok$recon_mean[best] + ok$recon_se[best], , drop = FALSE]
  band <- band[order(-band$reproducibility, band$alpha_1, band$alpha_2), , drop = FALSE]
  band$setting_id[1]
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `make document && make test FILTER=tuning`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add R/tuning.R NAMESPACE man tests/testthat/test-tuning.R
git commit -m "feat: alpha-scaled lambdas and the one-SE selection rule

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 10: Program matching across folds

**Files:**
- Create: `R/program_matching.R`, `tests/testthat/test-matching.R`

**Interfaces:**
- Produces:
  - `cosineSim(A, B)`: compares the rows of A with the rows of B.
  - `assignMax(sim)`: an integer vector giving the column assigned to each row.
  - `matchPrograms(H_ref, H_fold, tau)`: `data.frame(ref, fold_prog, cosine, matched)`.
  - `matchNull(H_ref, H_fold, n_perm = 1000, seed = 1)`: numeric, with length `nrow(H_fold) * n_perm`.
  - `programStability(matches, min_folds = 4)`: `data.frame(ref, n_matched, mean_cosine, stable)`. `matches` has columns `fold, ref, cosine, matched`.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-matching.R`:

```r
test_that("assignMax solves the case where greedy matching fails", {
  sim <- matrix(c(0.90, 0.80,
                  0.85, 0.10), 2, byrow = TRUE)   # rows: ref A, B; columns: fold programs 1, 2
  expect_equal(assignMax(sim), c(2L, 1L))
})

test_that("matchPrograms recovers a permutation and flags weak matches", {
  withr::with_seed(7, H_ref <- matrix(stats::runif(4 * 50), 4))
  perm <- c(3, 1, 4, 2)
  H_fold <- H_ref[perm, ]
  m <- matchPrograms(H_ref, H_fold, tau = 0.5)
  expect_equal(m$fold_prog, match(1:4, perm))
  expect_true(all(m$matched))
  expect_equal(m$cosine, rep(1, 4))
  m2 <- matchPrograms(H_ref, H_fold, tau = 1.01)
  expect_false(any(m2$matched))
})

test_that("the gene-permutation null sits below true matches", {
  withr::with_seed(8, H_ref <- matrix(stats::rexp(4 * 200), 4))
  null <- matchNull(H_ref, H_ref + 0.01, n_perm = 50, seed = 1)
  expect_length(null, 4 * 50)
  expect_true(all(null >= 0 & null <= 1))
  expect_lt(stats::quantile(null, 0.99), 0.99)
})

test_that("programStability counts folds matched", {
  matches <- data.frame(fold = rep(1:5, each = 2), ref = rep(1:2, 5),
                        cosine = c(.9, .5, .9, .5, .9, .9, .9, .5, .8, .5),
                        matched = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE))
  st <- programStability(matches)
  expect_equal(st$n_matched, c(5L, 1L))
  expect_equal(st$stable, c(TRUE, FALSE))
  expect_equal(st$mean_cosine, c(mean(c(.9, .9, .9, .9, .8)), mean(c(.5, .5, .9, .5, .5))))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=matching`
Expected: FAIL with `could not find function "assignMax"`.

- [ ] **Step 3: Implement**

`R/program_matching.R`:

```r
#' cosineSim
#' @param A,B Matrices with the same number of columns.
#' @return nrow(A) x nrow(B) matrix of cosine similarities between rows.
#' @export
cosineSim <- function(A, B) {
  na <- sqrt(rowSums(A^2)); nb <- sqrt(rowSums(B^2))
  na[na == 0] <- 1; nb[nb == 0] <- 1
  (A %*% t(B)) / outer(na, nb)
}

#' assignMax
#'
#' One-to-one assignment maximizing total similarity (Hungarian algorithm via `clue`).
#' @param sim Square nonnegative similarity matrix.
#' @return Integer vector: the column assigned to each row.
#' @export
assignMax <- function(sim) {
  sim[is.na(sim) | sim < 0] <- 0
  as.integer(clue::solve_LSAP(sim, maximum = TRUE))
}

#' matchPrograms
#'
#' Matches a fold's programs to the reference programs by gene-loading cosine (spec Section 8.2).
#' @param H_ref Reference H (k x genes). @param H_fold Fold H (k x genes).
#' @param tau Cosine below which a forced match counts as unmatched.
#' @return Data frame with ref, fold_prog, cosine, matched.
#' @export
matchPrograms <- function(H_ref, H_fold, tau) {
  sim <- cosineSim(H_ref, H_fold)
  a <- assignMax(sim)
  cos <- sim[cbind(seq_len(nrow(sim)), a)]
  data.frame(ref = seq_len(nrow(sim)), fold_prog = a, cosine = cos, matched = cos >= tau)
}

#' matchNull
#'
#' Best-match cosine of each fold program against all reference programs after permuting the
#' fold's gene loadings. Pool across folds and take the 99th percentile for `tau`.
#' @param H_ref,H_fold k x genes matrices.
#' @param n_perm Number of permutations. @param seed Random seed.
#' @return Numeric vector of length nrow(H_fold) * n_perm.
#' @export
matchNull <- function(H_ref, H_fold, n_perm = 1000, seed = 1) {
  withr::with_seed(seed, as.vector(vapply(seq_len(n_perm), function(b) {
    Hp <- H_fold[, sample.int(ncol(H_fold)), drop = FALSE]
    apply(cosineSim(Hp, H_ref), 1, max)
  }, numeric(nrow(H_fold)))))
}

#' programStability
#' @param matches Data frame with fold, ref, cosine, matched (one row per fold x reference program).
#' @param min_folds Folds a program must match in to count as stable. Default 4.
#' @return Data frame with ref, n_matched, mean_cosine, stable.
#' @export
programStability <- function(matches, min_folds = 4) {
  n_matched <- tapply(matches$matched, matches$ref, sum)
  mean_cos <- tapply(matches$cosine, matches$ref, mean)
  data.frame(ref = as.integer(names(n_matched)), n_matched = as.integer(n_matched),
             mean_cosine = as.numeric(mean_cos), stable = as.integer(n_matched) >= min_folds)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `make document && make test FILTER=matching`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add R/program_matching.R NAMESPACE man tests/testthat/test-matching.R
git commit -m "feat: Hungarian program matching with a permutation threshold

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 11: Persistence scores

**Files:**
- Create: `R/persistence.R`, `tests/testthat/test-persistence.R`

**Interfaces:**
- Consumes: a `pairs` data frame with columns `tumor_act, organoid_act, cancer_type, fold`. There is one row per held-out pair for a single program.
- Produces:
  - `normRankWithin(x, group)`, `concordanceStat(zt, zo)`, `permuteWithin(strata)`, `resampleWithin(strata)`.
  - `persistenceConcordance(pairs, n_perm = 10000, n_boot = 2000, seed = 1)`: a one-row data frame with columns `conc_rho, conc_p, conc_ci_low, conc_ci_high, detect_limit, n_pairs`.
  - `persistenceRetention(pairs, n_perm = 10000, n_boot = 2000, seed = 1)`: a one-row data frame with columns `ret_median_log2, ret_p, ret_ci_low, ret_ci_high, eps`.
  - `classifyPrograms(conc_fdr, ret_fdr, ret_median_log2, alpha = 0.05)`: a character vector with values in `drops_out, reduced_concordant, persists, present_not_patient_specific`.
  - `concordancePower(strata, rhos = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), n_sims = 200, n_perm = 999, alpha = 0.05, seed = 1)`: a data frame with `rho, power, n_pairs`. `strata` has `cancer_type, fold`, one row per held-out pair. This is the power curve of spec Section 8.3, simulated at the run's real pair structure.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-persistence.R`:

```r
make_pairs <- function(n = 60, shared = TRUE, seed = 9) {
  withr::with_seed(seed, {
    a <- stats::rgamma(n, 2, 1)
    t_act <- a + stats::rnorm(n, sd = 0.1)
    o_act <- if (shared) a + stats::rnorm(n, sd = 0.1) else stats::rgamma(n, 2, 1)
  })
  data.frame(tumor_act = pmax(t_act, 0), organoid_act = pmax(o_act, 0),
             cancer_type = rep(c("COAD", "PAAD", "BRCA"), length.out = n),
             fold = rep(1:5, each = 1, length.out = n))
}

test_that("shared activity is concordant and independent activity is not", {
  hit <- persistenceConcordance(make_pairs(shared = TRUE), n_perm = 999, n_boot = 100)
  miss <- persistenceConcordance(make_pairs(shared = FALSE), n_perm = 999, n_boot = 100)
  expect_lt(hit$conc_p, 0.01); expect_gt(hit$conc_rho, 0.8)
  expect_gt(miss$conc_p, 0.01)
  expect_true(hit$conc_ci_low < hit$conc_rho && hit$conc_rho < hit$conc_ci_high)
  expect_lt(hit$detect_limit, 0.5)
})

test_that("concordance p-values are calibrated under the null", {
  p <- vapply(1:200, function(s) persistenceConcordance(make_pairs(shared = FALSE, seed = s),
                                                        n_perm = 199, n_boot = 2, seed = s)$conc_p, numeric(1))
  expect_true(mean(p < 0.05) > 0.01 && mean(p < 0.05) < 0.10)
})

test_that("all-zero organoid activity gives rho 0 and p 1, never NA", {          # Review Focus 1
  pr <- make_pairs(); pr$organoid_act <- 0
  r <- persistenceConcordance(pr, n_perm = 99, n_boot = 20)
  expect_equal(r$conc_rho, 0); expect_equal(r$conc_p, 1)
  expect_false(anyNA(unlist(r)))
})

test_that("a stratum of size one stays fixed and shuffles never cross strata", { # Review Focus 2
  strata <- c("a", "b", "b", "c", "b")
  idx <- withr::with_seed(1, permuteWithin(strata))
  expect_equal(idx[c(1, 4)], c(1L, 4L))
  expect_setequal(idx[c(2, 3, 5)], c(2L, 3L, 5L))
  big <- sample(letters[1:6], 300, replace = TRUE)
  expect_identical(big[withr::with_seed(2, permuteWithin(big))], big)
  pr <- make_pairs(); pr$cancer_type[1] <- "UNIQUE"
  expect_false(anyNA(unlist(persistenceConcordance(pr, n_perm = 99, n_boot = 20))))
})

test_that("retention detects dropout and is calibrated when there is none", {
  pr <- make_pairs(); pr$organoid_act <- pr$organoid_act * 0.1
  r <- persistenceRetention(pr, n_perm = 999, n_boot = 100)
  expect_lt(r$ret_median_log2, -2); expect_lt(r$ret_p, 0.01)
  p <- vapply(1:200, function(s) persistenceRetention(make_pairs(shared = FALSE, seed = s),
                                                      n_perm = 199, n_boot = 2, seed = s)$ret_p, numeric(1))
  expect_true(mean(p < 0.05) > 0.01 && mean(p < 0.05) < 0.10)
})

test_that("concordance power is near alpha at rho 0 and high for strong concordance", {
  strata <- data.frame(cancer_type = rep(c("COAD", "PAAD", "BRCA"), length.out = 60), fold = rep(1:5, length.out = 60))
  pw <- concordancePower(strata, rhos = c(0, 0.6), n_sims = 40, n_perm = 99)
  expect_equal(pw$rho, c(0, 0.6))
  expect_lt(pw$power[1], 0.2)
  expect_gt(pw$power[2], 0.8)
  expect_equal(unique(pw$n_pairs), 60L)
})

test_that("classifyPrograms maps the two calls to the four groups", {
  expect_equal(classifyPrograms(conc_fdr = c(0.5, 0.01, 0.01, 0.5),
                                ret_fdr = c(0.01, 0.01, 0.5, 0.5),
                                ret_median_log2 = c(-2, -2, 0, 0)),
               c("drops_out", "reduced_concordant", "persists", "present_not_patient_specific"))
  expect_equal(classifyPrograms(0.5, 0.01, -0.5), "present_not_patient_specific")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=persistence`
Expected: FAIL with `could not find function "persistenceConcordance"`.

- [ ] **Step 3: Implement**

`R/persistence.R`:

```r
## Persistence scoring (spec Section 8.3).

#' normRankWithin
#' @param x Numeric vector. @param group Grouping vector (cancer type).
#' @return `(rank - (n + 1) / 2) / n` within each group, so every group is weighted by its size.
#' @export
normRankWithin <- function(x, group) {
  stats::ave(x, group, FUN = function(v) (rank(v) - (length(v) + 1) / 2) / length(v))
}

#' concordanceStat
#' @param zt,zo Normalized within-type ranks of tumor and organoid activity.
#' @return Pearson correlation, or 0 when either side has no variance.
#' @export
concordanceStat <- function(zt, zo) {
  if (length(zt) < 3 || !isTRUE(stats::sd(zt) > 0) || !isTRUE(stats::sd(zo) > 0)) return(0)
  stats::cor(zt, zo)
}

#' permuteWithin
#' @param strata Grouping vector.
#' @return Permuted indices that only move positions within a stratum. A size-one stratum stays fixed.
#' @export
permuteWithin <- function(strata) {
  groups <- split(seq_along(strata), strata)
  out <- seq_along(strata)
  out[unlist(groups, use.names = FALSE)] <- unlist(lapply(groups, function(i) i[sample.int(length(i))]), use.names = FALSE)
  out
}

#' resampleWithin
#' @param strata Grouping vector.
#' @return Bootstrap indices drawn with replacement within each stratum.
#' @export
resampleWithin <- function(strata) {
  groups <- split(seq_along(strata), strata)
  unlist(lapply(groups, function(i) i[sample.int(length(i), replace = TRUE)]), use.names = FALSE)
}

#' persistenceConcordance
#'
#' Stratified Spearman concordance between held-out tumor and organoid activity. The null
#' shuffles organoids among tumors within cancer type x fold.
#' @param pairs Data frame with tumor_act, organoid_act, cancer_type, fold.
#' @param n_perm Permutations. @param n_boot Bootstrap resamples. @param seed Random seed.
#' @return One-row data frame: conc_rho, conc_p (one-sided), conc_ci_low, conc_ci_high,
#'   detect_limit (95th percentile of the null), n_pairs.
#' @export
persistenceConcordance <- function(pairs, n_perm = 10000, n_boot = 2000, seed = 1) {
  zt <- normRankWithin(pairs$tumor_act, pairs$cancer_type)
  zo <- normRankWithin(pairs$organoid_act, pairs$cancer_type)
  obs <- concordanceStat(zt, zo)
  strata <- interaction(pairs$cancer_type, pairs$fold, drop = TRUE)
  draws <- withr::with_seed(seed, list(
    null = vapply(seq_len(n_perm), function(b) concordanceStat(zt, zo[permuteWithin(strata)]), numeric(1)),
    boot = vapply(seq_len(n_boot), function(b) {
      i <- resampleWithin(strata)
      concordanceStat(normRankWithin(pairs$tumor_act[i], pairs$cancer_type[i]),
                      normRankWithin(pairs$organoid_act[i], pairs$cancer_type[i]))
    }, numeric(1))))
  data.frame(conc_rho = obs,
             conc_p = (1 + sum(draws$null >= obs)) / (1 + n_perm),
             conc_ci_low = unname(stats::quantile(draws$boot, 0.025)),
             conc_ci_high = unname(stats::quantile(draws$boot, 0.975)),
             detect_limit = unname(stats::quantile(draws$null, 0.95)),
             n_pairs = nrow(pairs))
}

#' persistenceRetention
#'
#' Median over pairs of `log2((organoid + eps) / (tumor + eps))`, where eps is the 1st
#' percentile of the program's nonzero held-out activities. The one-sided test for reduction
#' flips signs within pairs.
#' @param pairs Data frame with tumor_act, organoid_act, cancer_type, fold.
#' @param n_perm Sign-flip draws. @param n_boot Bootstrap resamples. @param seed Random seed.
#' @return One-row data frame: ret_median_log2, ret_p, ret_ci_low, ret_ci_high, eps.
#' @export
persistenceRetention <- function(pairs, n_perm = 10000, n_boot = 2000, seed = 1) {
  act <- c(pairs$tumor_act, pairs$organoid_act)
  nz <- act[act > 0]
  eps <- if (length(nz)) unname(stats::quantile(nz, 0.01)) else 1e-8
  lr <- log2((pairs$organoid_act + eps) / (pairs$tumor_act + eps))
  obs <- stats::median(lr)
  strata <- interaction(pairs$cancer_type, pairs$fold, drop = TRUE)
  draws <- withr::with_seed(seed, list(
    null = vapply(seq_len(n_perm), function(b) stats::median(lr * sample(c(-1, 1), length(lr), replace = TRUE)), numeric(1)),
    boot = vapply(seq_len(n_boot), function(b) stats::median(lr[resampleWithin(strata)]), numeric(1))))
  data.frame(ret_median_log2 = obs,
             ret_p = (1 + sum(draws$null <= obs)) / (1 + n_perm),
             ret_ci_low = unname(stats::quantile(draws$boot, 0.025)),
             ret_ci_high = unname(stats::quantile(draws$boot, 0.975)),
             eps = eps)
}

#' classifyPrograms
#'
#' Spec Section 8.3 groups. Reduced: retention FDR < alpha and median log2 ratio <= -1.
#' Concordant: concordance FDR < alpha.
#' @param conc_fdr,ret_fdr,ret_median_log2 Vectors over programs. @param alpha FDR level.
#' @return Character vector of group labels.
#' @export
classifyPrograms <- function(conc_fdr, ret_fdr, ret_median_log2, alpha = 0.05) {
  reduced <- ret_fdr < alpha & ret_median_log2 <= -1
  concordant <- conc_fdr < alpha
  lookup <- c("TRUE FALSE" = "drops_out", "TRUE TRUE" = "reduced_concordant",
              "FALSE TRUE" = "persists", "FALSE FALSE" = "present_not_patient_specific")
  unname(lookup[paste(reduced %in% TRUE, concordant %in% TRUE)])
}

#' concordancePower
#'
#' Power of the concordance test at a given pair structure (spec Section 8.3). Within each
#' simulated pair, tumor and organoid activity are standard normal with correlation `rho`. The
#' test is [persistenceConcordance()] at one-sided level `alpha` (nominal, before FDR).
#' @param strata Data frame with cancer_type and fold, one row per held-out pair.
#' @param rhos True correlations. @param n_sims Simulations per rho. @param n_perm Permutations per test.
#' @param alpha One-sided level. @param seed Random seed.
#' @return Data frame with rho, power, n_pairs.
#' @export
concordancePower <- function(strata, rhos = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), n_sims = 200, n_perm = 999,
                             alpha = 0.05, seed = 1) {
  n <- nrow(strata)
  grid <- expand.grid(sim = seq_len(n_sims), rho = rhos)
  hit <- withr::with_seed(seed, vapply(seq_len(nrow(grid)), function(i) {
    z1 <- stats::rnorm(n)
    z2 <- grid$rho[i] * z1 + sqrt(1 - grid$rho[i]^2) * stats::rnorm(n)
    persistenceConcordance(data.frame(tumor_act = z1, organoid_act = z2, cancer_type = strata$cancer_type,
                                      fold = strata$fold), n_perm = n_perm, n_boot = 2, seed = i)$conc_p < alpha
  }, logical(1)))
  data.frame(rho = rhos, power = as.numeric(tapply(hit, grid$rho, mean)[as.character(rhos)]), n_pairs = n)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `make document && make test FILTER=persistence`
Expected: PASS (7 tests). The calibration and power tests take about 30 seconds.

- [ ] **Step 5: Commit**

```bash
git add R/persistence.R NAMESPACE man tests/testthat/test-persistence.R
git commit -m "feat: stratified concordance and retention with permutation nulls

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 12: Program annotation

**Files:**
- Create: `R/annotation.R`, `tests/testthat/test-annotation.R`

**Interfaces:**
- Produces:
  - `topGenes(h, frac = 0.05)`.
  - `enrichHypergeom(H, gene_sets, frac = 0.05)`: a data frame with `program, gene_set, overlap, set_size, n_top, p, fdr`.
  - `programAssociation(act, covariate, label)`: a data frame with `program, covariate, rho, p, fdr`.
  - `estimateScores(expr)`: a data frame with `profile_id, stromal, immune, estimate`. `expr` is genes x samples with HGNC rownames.
  - `mutationMatrix(mut, row_ids, row_sample, genes, profiled_samples)`: a logical matrix, NA for samples that were not sequenced.
  - `genotypeAssociation(act, mut, min_mut = 10)`: a data frame with `program, gene, n_mut, delta, p, fdr`.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-annotation.R`:

```r
test_that("enrichHypergeom finds the planted gene set", {
  genes <- sprintf("g%03d", 1:200)
  H <- matrix(0.01, 2, 200, dimnames = list(c("LV1", "LV2"), genes))
  H["LV1", 1:10] <- 1
  H["LV2", 101:110] <- 1
  sets <- list(SET_A = genes[1:20], SET_B = genes[101:120], SET_C = genes[150:170])
  e <- enrichHypergeom(H, sets, frac = 0.05)
  expect_named(e, c("program", "gene_set", "overlap", "set_size", "n_top", "p", "fdr"))
  top <- e[order(e$p), ]
  expect_equal(top$gene_set[top$program == "LV1"][1], "SET_A")
  expect_equal(top$gene_set[top$program == "LV2"][1], "SET_B")
  expect_equal(e$overlap[e$program == "LV1" & e$gene_set == "SET_A"], 10)
})

test_that("programAssociation is NA-safe and finds a negative purity association", {
  withr::with_seed(10, {
    purity <- stats::runif(50, 0.2, 0.9)
    act <- cbind(LV1 = 1 - purity + stats::rnorm(50, sd = 0.05), LV2 = stats::rnorm(50), LV3 = 0)
  })
  act[1, "LV2"] <- NA
  a <- programAssociation(act, purity, "purity_dna")
  expect_lt(a$rho[1], -0.8); expect_lt(a$fdr[1], 0.01)
  expect_equal(a$rho[3], 0); expect_equal(a$p[3], 1)
  expect_false(anyNA(a$p))
})

test_that("mutationMatrix keeps non-silent calls and marks unsequenced samples NA", {
  mut <- data.frame(sample_id = c("s1", "s1", "s2", "s3"),
                    Hugo_Symbol = c("KRAS", "TP53", "KRAS", "KRAS"),
                    Variant_Classification = c("Missense_Mutation", "Silent", "Nonsense_Mutation", "Missense_Mutation"))
  M <- mutationMatrix(mut, row_ids = c("p1", "p2", "p3", "p4"), row_sample = c("s1", "s2", "s3", "s4"),
                      genes = c("KRAS", "TP53"), profiled_samples = c("s1", "s2", "s3"))
  expect_equal(M[, "KRAS"], c(p1 = TRUE, p2 = TRUE, p3 = TRUE, p4 = NA))
  expect_equal(M[, "TP53"], c(p1 = FALSE, p2 = FALSE, p3 = FALSE, p4 = NA))
})

test_that("genotypeAssociation tests only genes with enough mutated and wild-type tumors", {
  withr::with_seed(11, {
    m <- c(rep(TRUE, 15), rep(FALSE, 35))
    act <- cbind(LV1 = ifelse(m, 2, 0) + stats::rnorm(50, sd = 0.2), LV2 = stats::rnorm(50))
  })
  mut <- cbind(KRAS = m, RARE = c(rep(TRUE, 3), rep(FALSE, 47)))
  g <- genotypeAssociation(act, mut, min_mut = 10)
  expect_equal(unique(g$gene), "KRAS")
  expect_lt(g$fdr[g$program == "LV1"], 0.01)
  expect_gt(g$delta[g$program == "LV1"], 1)
})

test_that("estimateScores returns stromal and immune scores per sample", {
  skip_if_not_installed("tidyestimate")
  genes <- unique(tidyestimate::common_genes$hgnc_symbol)
  withr::with_seed(12, expr <- matrix(stats::rnorm(length(genes) * 4, 5), length(genes),
                                     dimnames = list(genes, paste0("s", 1:4))))
  sc <- estimateScores(expr)
  expect_named(sc, c("profile_id", "stromal", "immune", "estimate"))
  expect_equal(sc$profile_id, paste0("s", 1:4))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=annotation`
Expected: FAIL with `could not find function "enrichHypergeom"`.

- [ ] **Step 3: Implement**

`R/annotation.R`:

```r
## Post-hoc annotation of programs (spec Section 9).

.nonsilent <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site")

#' topGenes
#' @param h Named loading vector. @param frac Fraction kept. Default 0.05.
#' @return Names of the top `floor(frac * length(h))` genes (at least one).
#' @export
topGenes <- function(h, frac = 0.05) {
  names(sort(h, decreasing = TRUE))[seq_len(max(1, floor(frac * length(h))))]
}

#' enrichHypergeom
#'
#' Hypergeometric enrichment of each program's top genes, with the program's gene set
#' (colnames(H)) as background. Vectorized over programs and gene sets.
#' @param H Programs x genes loadings. @param gene_sets Named list of gene symbol vectors.
#' @param frac Fraction of genes used as each program's top set.
#' @return Data frame with program, gene_set, overlap, set_size, n_top, p, fdr (BH over all tests).
#' @export
enrichHypergeom <- function(H, gene_sets, frac = 0.05) {
  bg <- colnames(H)
  sets <- lapply(gene_sets, intersect, bg)
  sets <- sets[lengths(sets) > 0]
  M <- vapply(sets, function(s) bg %in% s, logical(length(bg)))                 # genes x sets
  M <- matrix(M, nrow = length(bg), dimnames = list(bg, names(sets)))
  Tm <- t(apply(H, 1, function(h) bg %in% topGenes(stats::setNames(h, bg), frac)))   # programs x genes
  Tm <- matrix(Tm, nrow = nrow(H), dimnames = list(rownames(H), bg))
  overlap <- Tm %*% M
  size <- matrix(colSums(M), nrow(overlap), ncol(overlap), byrow = TRUE)
  n_top <- rowSums(Tm)
  p <- stats::phyper(overlap - 1, size, length(bg) - size, n_top, lower.tail = FALSE)
  out <- data.frame(program = rep(rownames(H), times = ncol(M)),
                    gene_set = rep(colnames(M), each = nrow(H)),
                    overlap = as.vector(overlap), set_size = as.vector(size),
                    n_top = rep(n_top, times = ncol(M)), p = as.vector(p), stringsAsFactors = FALSE)
  out$fdr <- stats::p.adjust(out$p, "BH")
  out
}

#' programAssociation
#'
#' Spearman correlation of each program's activity with a covariate. Columns with too few
#' complete values or no variance get rho = 0 and p = 1.
#' @param act Samples x programs activity. @param covariate Numeric vector aligned with rows.
#' @param label Name recorded in the `covariate` column.
#' @return Data frame with program, covariate, rho, p, fdr.
#' @export
programAssociation <- function(act, covariate, label) {
  res <- vapply(seq_len(ncol(act)), function(j) {
    a <- act[, j]
    ok <- !is.na(a) & !is.na(covariate)
    if (sum(ok) < 5 || !isTRUE(stats::sd(a[ok]) > 0) || !isTRUE(stats::sd(covariate[ok]) > 0)) return(c(0, 1))
    ct <- suppressWarnings(stats::cor.test(a[ok], covariate[ok], method = "spearman", exact = FALSE))
    c(unname(ct$estimate), ct$p.value)
  }, numeric(2))
  data.frame(program = colnames(act), covariate = label, rho = res[1, ], p = res[2, ],
             fdr = stats::p.adjust(res[2, ], "BH"), stringsAsFactors = FALSE)
}

#' estimateScores
#'
#' ESTIMATE stromal and immune scores (tidyestimate), used for annotation only.
#' @param expr Genes x samples expression (log2 TPM) with HGNC symbols as rownames.
#' @return Data frame with profile_id, stromal, immune, estimate.
#' @export
estimateScores <- function(expr) {
  if (!requireNamespace("tidyestimate", quietly = TRUE)) stop("estimateScores needs the tidyestimate package")
  df <- as.data.frame(expr)
  common <- tidyestimate::filter_common_genes(df, id = "hgnc_symbol", tidy = FALSE,
                                              tell_missing = FALSE, find_alias = FALSE)
  sc <- as.data.frame(tidyestimate::estimate_score(common, is_affymetrix = FALSE))
  data.frame(profile_id = sc$sample, stromal = sc$stromal, immune = sc$immune, estimate = sc$estimate,
             stringsAsFactors = FALSE)
}

#' mutationMatrix
#' @param mut Data frame with sample_id, Hugo_Symbol, Variant_Classification.
#' @param row_ids Row names of the result (profile ids). @param row_sample Sample id of each row.
#' @param genes Columns (for example OncoKB genes). @param profiled_samples Sequenced sample ids.
#' @return Logical matrix: TRUE for a non-silent mutation; NA for rows that were not sequenced.
#' @export
mutationMatrix <- function(mut, row_ids, row_sample, genes, profiled_samples) {
  m <- mut[mut$Variant_Classification %in% .nonsilent & mut$Hugo_Symbol %in% genes &
             mut$sample_id %in% row_sample, , drop = FALSE]
  M <- matrix(FALSE, length(row_ids), length(genes), dimnames = list(row_ids, genes))
  M[cbind(match(m$sample_id, row_sample), match(m$Hugo_Symbol, genes))] <- TRUE
  M[!(row_sample %in% profiled_samples), ] <- NA
  M
}

#' genotypeAssociation
#'
#' Wilcoxon test of program activity in mutated vs wild-type tumors, for genes with at least
#' `min_mut` mutated and at least 2 wild-type tumors.
#' @param act Tumors x programs. @param mut Tumors x genes logical (rows aligned with `act`).
#' @param min_mut Minimum number of mutated tumors. Default 10.
#' @return Data frame with program, gene, n_mut, delta (median difference), p, fdr.
#' @export
genotypeAssociation <- function(act, mut, min_mut = 10) {
  n_mut <- colSums(mut, na.rm = TRUE)
  n_wt <- colSums(!mut, na.rm = TRUE)
  mut <- mut[, n_mut >= min_mut & n_wt >= 2, drop = FALSE]
  if (!ncol(mut)) return(data.frame(program = character(0), gene = character(0), n_mut = integer(0),
                                    delta = numeric(0), p = numeric(0), fdr = numeric(0)))
  grid <- expand.grid(program = colnames(act), gene = colnames(mut), stringsAsFactors = FALSE)
  res <- mapply(function(p, g) {
    m <- mut[, g]; a <- act[, p]; ok <- !is.na(m) & !is.na(a)
    w <- suppressWarnings(stats::wilcox.test(a[ok & m], a[ok & !m], exact = FALSE))
    c(stats::median(a[ok & m]) - stats::median(a[ok & !m]), w$p.value)
  }, grid$program, grid$gene)
  data.frame(grid, n_mut = as.integer(colSums(mut, na.rm = TRUE)[grid$gene]), delta = res[1, ], p = res[2, ],
             fdr = stats::p.adjust(res[2, ], "BH"), row.names = NULL, stringsAsFactors = FALSE)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `make document && make test FILTER=annotation`
Expected: PASS (5 tests).

If the `estimateScores` test fails because `tidyestimate`'s function or dataset names differ, read `help(package = "tidyestimate")` in the environment. Adapt the two calls in `estimateScores` (not the test's expected output columns) to the installed version's API.

- [ ] **Step 5: Commit**

```bash
git add R/annotation.R NAMESPACE man tests/testthat/test-annotation.R
git commit -m "feat: enrichment, purity/ESTIMATE and genotype annotation

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 13: One cross-fitting task

**Files:**
- Create: `R/cross_fitting.R`, `tests/testthat/test-cross-fitting.R`

**Interfaces:**
- Consumes: `buildSampleGraph`, `svdInit`, `lambdaFromAlpha`, `GRACKLE`, and the existing `project_W(test_data, H_train)` in `R/Evaluation.R`.
- Produces:
  - `normalizeRows(H)`.
  - `fitTask(Y, samples, gene_graph, k, alpha_1, alpha_2, fold, train_role = c("all", "tumor"), test_role = c("all", "organoid"), iterations = 2000, diff_threshold = 1e-6, dtype = "float32", k_nn = 10)`. It returns `list(H, W_test, lambda, n_iter, converged, final_loss, recon_test, train_ids, test_ids)`:
    - `fold = 0` is the reference fit on all samples, with `W_test = NULL` and `recon_test = NA`.
    - `H` rows have unit L2 norm and names `LV1..LVk`.
    - `W_test` has held-out profile ids as rows.
    - `recon_test` is the mean squared error per entry of the held-out reconstruction.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-cross-fitting.R`:

```r
small_fixture <- function() {
  s <- toy_samples()
  s <- s[rep(seq_len(nrow(s)), 3), ]
  s$profile_id <- paste0(s$profile_id, "_", rep(1:3, each = 12))
  s$case_id <- paste0(s$case_id, "_", rep(1:3, each = 12))
  s$sample_id <- paste0(s$sample_id, "_", rep(1:3, each = 12))
  s$parent_sample_id <- ifelse(is.na(s$parent_sample_id), NA, paste0(s$parent_sample_id, "_", rep(1:3, each = 12)))
  s <- makeFolds(s, seed = 1)
  withr::with_seed(13, Y <- matrix(stats::runif(nrow(s) * 40), nrow(s),
                                   dimnames = list(s$profile_id, paste0("g", 1:40))))
  G <- matrix(0, 40, 40, dimnames = list(colnames(Y), colnames(Y)))
  G[cbind(1:39, 2:40)] <- 1; G <- pmax(G, t(G))
  list(Y = Y, samples = s, G = G)
}

test_that("fitTask holds out the fold and projects held-out samples", {
  skip_if_no_tf()
  fx <- small_fixture()
  r <- fitTask(fx$Y, fx$samples, fx$G, k = 2, alpha_1 = 0.1, alpha_2 = 0.1, fold = 1, iterations = 30)
  held <- fx$samples$profile_id[fx$samples$fold == 1]
  expect_setequal(r$test_ids, held)
  expect_length(intersect(r$train_ids, held), 0)
  expect_equal(rownames(r$W_test), held)
  expect_equal(unname(sqrt(rowSums(r$H^2))), c(1, 1))
  expect_equal(rownames(r$H), c("LV1", "LV2"))
  expect_true(is.finite(r$recon_test) && r$recon_test > 0)
  expect_named(r$lambda, c("lambda_1", "lambda_2"))
})

test_that("the reference fit uses every sample and the transfer roles are respected", {
  skip_if_no_tf()
  fx <- small_fixture()
  ref <- fitTask(fx$Y, fx$samples, fx$G, k = 2, alpha_1 = 0.1, alpha_2 = 0.1, fold = 0, iterations = 10)
  expect_setequal(ref$train_ids, fx$samples$profile_id)
  expect_null(ref$W_test)
  tr <- fitTask(fx$Y, fx$samples, fx$G, k = 2, alpha_1 = 0.1, alpha_2 = 0.1, fold = 2,
                train_role = "tumor", test_role = "organoid", iterations = 10)
  role <- stats::setNames(fx$samples$role, fx$samples$profile_id)
  expect_true(all(role[tr$train_ids] == "tumor"))
  expect_true(all(role[tr$test_ids] == "organoid"))
})

test_that("misaligned inputs are rejected", {
  fx <- small_fixture()
  expect_error(fitTask(fx$Y[, 40:1], fx$samples, fx$G, 2, 0.1, 0.1, 1), "aligned")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=cross-fitting`
Expected: FAIL with `could not find function "fitTask"`.

- [ ] **Step 3: Implement**

`R/cross_fitting.R`:

```r
#' normalizeRows
#' @param H Matrix. @return `H` with each row scaled to unit L2 norm (zero rows unchanged).
#' @export
normalizeRows <- function(H) {
  n <- sqrt(rowSums(H^2))
  n[n == 0] <- 1
  H / n
}

#' fitTask
#'
#' One cross-fitting task (spec Section 8.1). It builds the sample graph from the training
#' samples, scales lambda from alpha, fits GRACKLE, normalizes the rows of H, and projects the
#' held-out samples onto H with NNLS.
#' @param Y Samples x genes matrix; rownames are profile ids.
#' @param samples Sample table whose `profile_id` equals `rownames(Y)`, with a `fold` column.
#' @param gene_graph Genes x genes matrix aligned with `colnames(Y)`.
#' @param k Rank. @param alpha_1,alpha_2 Penalty shares (spec Section 7.2).
#' @param fold Held-out fold; 0 is the reference fit on all samples.
#' @param train_role "all", or "tumor" (transfer test). @param test_role "all", or "organoid".
#' @param iterations,diff_threshold,dtype Passed to [GRACKLE()].
#' @param k_nn Neighbors in the sample graph.
#' @return A list with H, W_test, lambda, n_iter, converged, final_loss, recon_test, train_ids, test_ids.
#' @export
fitTask <- function(Y, samples, gene_graph, k, alpha_1, alpha_2, fold,
                    train_role = c("all", "tumor"), test_role = c("all", "organoid"),
                    iterations = 2000, diff_threshold = 1e-6, dtype = "float32", k_nn = 10) {
  train_role <- match.arg(train_role)
  test_role <- match.arg(test_role)
  if (!identical(rownames(Y), samples$profile_id) || !identical(colnames(Y), rownames(gene_graph))) {
    stop("fitTask: Y, samples and gene_graph must be aligned (rownames(Y) == samples$profile_id, colnames(Y) == rownames(gene_graph))")
  }
  train <- if (fold == 0) rep(TRUE, nrow(samples)) else samples$fold != fold
  test <- if (fold == 0) rep(FALSE, nrow(samples)) else samples$fold == fold
  if (train_role == "tumor") train <- train & samples$role == "tumor"
  if (test_role == "organoid") test <- test & samples$role == "organoid"
  Ytr <- Y[train, , drop = FALSE]
  S <- buildSampleGraph(samples[train, , drop = FALSE], k_nn = k_nn)$S
  init <- svdInit(Ytr, k)
  lam <- lambdaFromAlpha(Ytr, init$W, init$H, S, gene_graph, alpha_1, alpha_2)
  fit <- GRACKLE(Ytr, net_similarity = gene_graph, patient_similarity = S,
                 diff_threshold = diff_threshold, lambda_1 = lam[["lambda_1"]], lambda_2 = lam[["lambda_2"]],
                 k = k, iterations = iterations, svd_W = init$W, svd_H = init$H, dtype = dtype)
  H <- normalizeRows(fit$H)
  W_test <- NULL
  recon_test <- NA_real_
  if (any(test)) {
    Yte <- Y[test, , drop = FALSE]
    W_test <- project_W(Yte, H)
    recon_test <- sum((Yte - W_test %*% H)^2) / length(Yte)
  }
  list(H = H, W_test = W_test, lambda = lam, n_iter = fit$n_iter, converged = fit$converged,
       final_loss = utils::tail(fit$loss_trace, 1), recon_test = recon_test,
       train_ids = samples$profile_id[train], test_ids = samples$profile_id[test])
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `make document && make test FILTER=cross-fitting`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add R/cross_fitting.R NAMESPACE man tests/testthat/test-cross-fitting.R
git commit -m "feat: fitTask, one leakage-free cross-fitting task

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 14: Combine one gene-graph arm

**Files:**
- Create: `R/combine.R`, `tests/testthat/helper-combine.R`, `tests/testthat/test-combine.R`

**Interfaces:**
- Consumes: task results as written by `05_fit.R` (Task 21). Each one is `list(meta = list(task_id, arm, mode, setting_id, alpha_1, alpha_2, fold, k_mult, variant, k), fit = <fitTask output or NULL>, error = <NULL or message>)`.
- Produces:
  - `settingTable(results, expected, n_folds = 5)`. `expected` has `setting_id, alpha_1, alpha_2`. The result has `setting_id, alpha_1, alpha_2, complete, all_converged, recon_mean, recon_se, missing_or_failed`.
  - `matchSetting(fits, n_perm = 1000, seed = 1)`: `list(tau, matches, reproducibility, stability)`.
  - `assembleOutOfFold(folds, matches)`: held-out profiles x `LV<ref>` columns, NA where a program is unmatched.
  - `pairActivities(oof, pairs, samples, program)`.
  - `scoreSetting(oof, pairs, samples, stability, n_perm, n_boot, seed)`: one row per stable program, with the Task 11 columns plus `program, conc_fdr, ret_fdr, group`.
  - `perTypeConcordance(oof, pairs, samples, programs, min_pairs = 10)`: `program, cancer_type, n_pairs, rho`.
  - `annotatePrograms(H_ref, oof, samples, gene_sets, mut = NULL, estimate = NULL)`: `list(programs, enrichment, genotype)`.
  - `combineArm(results, expected, samples, pairs, gene_sets = list(), mut = NULL, estimate = NULL, n_perm = 10000, n_boot = 2000, n_perm_match = 1000, seed = 20261015, cores = 1, n_folds = 5)`: `list(selection, selected, tau_selected, programs, maps_long, matches, per_type, enrichment, genotype, ref_H, oof, convergence)`.

- [ ] **Step 1: Write the test helper (synthetic fit results; no TensorFlow needed)**

`tests/testthat/helper-combine.R`:

```r
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
```

- [ ] **Step 2: Write the failing tests**

`tests/testthat/test-combine.R`:

```r
build_case <- function() {
  w <- fake_world()
  res <- c(fake_results(w, "s1", 0.01, noise = 0.08, recon = c(.100, .104, .096, .102, .098), seed = 1),
           fake_results(w, "s2", 0.10, noise = 0.02, recon = c(.101, .103, .099, .102, .100), seed = 2),
           fake_results(w, "s3", 0.30, noise = 0.02, recon = rep(.090, 5), converged = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE), seed = 3),
           fake_results(w, "s4", 1.00, noise = 0.02, recon = rep(.080, 5), missing_folds = 3L, seed = 4))
  expected <- data.frame(setting_id = c("s1", "s2", "s3", "s4", "s5"), alpha_1 = c(0.01, 0.1, 0.3, 1, 3), alpha_2 = 0)
  list(w = w, res = res, expected = expected)
}

test_that("settingTable flags missing folds, failed folds and non-convergence", {      # Review Focus 5
  b <- build_case()
  failed <- b$res[[1]]; failed$fit <- NULL; failed$error <- "boom"; failed$meta$setting_id <- "s5"
  tab <- settingTable(c(b$res, list(failed)), b$expected)
  expect_equal(tab$complete, c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(tab$all_converged, c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(tab$missing_or_failed[tab$setting_id == "s4"], "3")
  expect_true(is.na(tab$recon_mean[tab$setting_id == "s4"]))
})

test_that("assembleOutOfFold undoes each fold's program order", {
  b <- build_case()
  fits <- Filter(function(r) r$meta$setting_id == "s2", b$res)
  m <- matchSetting(fits, n_perm = 100)
  oof <- assembleOutOfFold(Filter(function(r) r$meta$fold > 0, fits), m$matches)
  expect_equal(oof[rownames(oof), "LV1"], b$w$act[rownames(oof), 1])
  expect_equal(oof[rownames(oof), "LV2"], b$w$act[rownames(oof), 2])
})

test_that("combineArm selects, scores and excludes incomplete settings without erroring", {
  b <- build_case()
  out <- combineArm(b$res, b$expected, b$w$samples, b$w$pairs, gene_sets = list(),
                    n_perm = 199, n_boot = 50, n_perm_match = 100, seed = 1)
  expect_equal(out$selected, "s2")
  expect_false(out$selection$all_converged[out$selection$setting_id == "s3"])
  expect_false(out$selection$complete[out$selection$setting_id == "s4"])
  grp <- stats::setNames(out$programs$group, out$programs$program)
  expect_equal(grp[["LV1"]], "persists")
  expect_equal(grp[["LV2"]], "drops_out")
  expect_true(all(c("conc_rho", "ret_median_log2", "conc_fdr", "stable", "purity_rho", "microenv") %in% names(out$programs)))
  expect_setequal(unique(out$maps_long$setting_id), c("s1", "s2", "s3"))
  expect_true(all(out$per_type$n_pairs >= 10))
  expect_equal(nrow(out$convergence), length(b$res))
})
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `make test FILTER=combine`
Expected: FAIL with `could not find function "settingTable"`.

- [ ] **Step 4: Implement**

`R/combine.R`:

```r
## Combine the cross-fitting results of one gene-graph arm (spec Sections 7.3, 8.2, 8.3, 9).

#' settingTable
#' @param results Task results (see Task 14 interfaces). @param expected Data frame with
#'   setting_id, alpha_1, alpha_2 for every setting that should exist.
#' @param n_folds Number of cross-validation folds (fold 0 is the reference fit).
#' @return One row per expected setting: complete, all_converged, recon_mean, recon_se, missing_or_failed.
#' @export
settingTable <- function(results, expected, n_folds = 5) {
  meta <- do.call(rbind, lapply(results, function(r) data.frame(
    setting_id = r$meta$setting_id, fold = r$meta$fold,
    ok = is.null(r$error) && !is.null(r$fit),
    converged = !is.null(r$fit) && isTRUE(r$fit$converged),
    recon_test = if (is.null(r$fit)) NA_real_ else r$fit$recon_test, stringsAsFactors = FALSE)))
  if (is.null(meta)) meta <- data.frame(setting_id = character(0), fold = integer(0), ok = logical(0),
                                        converged = logical(0), recon_test = numeric(0))
  do.call(rbind, lapply(seq_len(nrow(expected)), function(i) {
    m <- meta[meta$setting_id == expected$setting_id[i], , drop = FALSE]
    present <- m$fold[m$ok]
    complete <- all(0:n_folds %in% present)
    cv <- m$recon_test[m$ok & m$fold > 0]
    data.frame(setting_id = expected$setting_id[i], alpha_1 = expected$alpha_1[i], alpha_2 = expected$alpha_2[i],
               complete = complete, all_converged = complete && all(m$converged[m$ok]),
               recon_mean = if (complete) mean(cv) else NA_real_,
               recon_se = if (complete) stats::sd(cv) / sqrt(length(cv)) else NA_real_,
               missing_or_failed = paste(setdiff(0:n_folds, present), collapse = ","),
               stringsAsFactors = FALSE)
  }))
}

#' matchSetting
#' @param fits The successful results of one setting (reference fold 0 plus folds 1..5).
#' @param n_perm Gene-permutation draws per fold. @param seed Random seed.
#' @return A list with tau, matches (fold, ref, fold_prog, cosine, matched), reproducibility, stability.
#' @export
matchSetting <- function(fits, n_perm = 1000, seed = 1) {
  ref <- Filter(function(r) r$meta$fold == 0, fits)[[1]]$fit$H
  folds <- Filter(function(r) r$meta$fold > 0, fits)
  null <- unlist(lapply(folds, function(r) matchNull(ref, r$fit$H, n_perm, seed + r$meta$fold)))
  tau <- unname(stats::quantile(null, 0.99))
  matches <- do.call(rbind, lapply(folds, function(r) cbind(fold = r$meta$fold, matchPrograms(ref, r$fit$H, tau))))
  list(tau = tau, matches = matches, reproducibility = mean(matches$cosine), stability = programStability(matches))
}

#' assembleOutOfFold
#' @param folds Results for folds 1..5 of one setting. @param matches From [matchSetting()].
#' @return Held-out profiles x reference programs (`LV<ref>`), NA where a program is unmatched in that fold.
#' @export
assembleOutOfFold <- function(folds, matches) {
  do.call(rbind, lapply(folds, function(fr) {
    m <- matches[matches$fold == fr$meta$fold, , drop = FALSE]
    m <- m[order(m$ref), , drop = FALSE]
    A <- fr$fit$W_test[, m$fold_prog, drop = FALSE]
    A[, !m$matched] <- NA
    colnames(A) <- paste0("LV", m$ref)
    A
  }))
}

#' pairActivities
#' @param oof Out-of-fold activity. @param pairs From [pairTable()]. @param samples Sample table.
#' @param program Column of `oof`.
#' @return Data frame with tumor_act, organoid_act, cancer_type, fold (complete rows only).
#' @export
pairActivities <- function(oof, pairs, samples, program) {
  p <- pairs[pairs$tumor %in% rownames(oof) & pairs$organoid %in% rownames(oof), , drop = FALSE]
  out <- data.frame(tumor_act = oof[p$tumor, program], organoid_act = oof[p$organoid, program],
                    cancer_type = samples$cancer_type[match(p$tumor, samples$profile_id)],
                    fold = samples$fold[match(p$tumor, samples$profile_id)], stringsAsFactors = FALSE)
  out[stats::complete.cases(out), , drop = FALSE]
}

#' scoreSetting
#' @param oof,pairs,samples See [pairActivities()]. @param stability From [matchSetting()].
#' @param n_perm,n_boot,seed Passed to the persistence tests.
#' @return One row per stable program with the persistence columns, conc_fdr, ret_fdr and group.
#' @export
scoreSetting <- function(oof, pairs, samples, stability, n_perm = 10000, n_boot = 2000, seed = 1) {
  stable <- paste0("LV", stability$ref[stability$stable])
  if (!length(stable)) return(data.frame())
  res <- do.call(rbind, lapply(stable, function(lv) {
    pa <- pairActivities(oof, pairs, samples, lv)
    cbind(program = lv, persistenceConcordance(pa, n_perm, n_boot, seed),
          persistenceRetention(pa, n_perm, n_boot, seed), stringsAsFactors = FALSE)
  }))
  res$conc_fdr <- stats::p.adjust(res$conc_p, "BH")
  res$ret_fdr <- stats::p.adjust(res$ret_p, "BH")
  res$group <- classifyPrograms(res$conc_fdr, res$ret_fdr, res$ret_median_log2)
  res
}

#' perTypeConcordance
#' @param oof,pairs,samples See [pairActivities()]. @param programs Programs to report.
#' @param min_pairs Minimum number of held-out pairs for a cancer type to be shown. Default 10.
#' @return Data frame with program, cancer_type, n_pairs, rho.
#' @export
perTypeConcordance <- function(oof, pairs, samples, programs, min_pairs = 10) {
  if (!length(programs)) return(data.frame())
  long <- do.call(rbind, lapply(programs, function(lv) cbind(program = lv, pairActivities(oof, pairs, samples, lv))))
  groups <- split(long, list(long$program, long$cancer_type), drop = TRUE)
  res <- do.call(rbind, lapply(groups, function(s) data.frame(
    program = s$program[1], cancer_type = s$cancer_type[1], n_pairs = nrow(s),
    rho = concordanceStat(normRankWithin(s$tumor_act, s$cancer_type), normRankWithin(s$organoid_act, s$cancer_type)),
    stringsAsFactors = FALSE)))
  res <- res[res$n_pairs >= min_pairs, , drop = FALSE]
  rownames(res) <- NULL
  res
}

#' annotatePrograms
#' @param H_ref Reference H. @param oof Out-of-fold activity. @param samples Sample table.
#' @param gene_sets Named list of gene sets (empty list to skip enrichment).
#' @param mut Tumor x gene mutation matrix, or NULL. @param estimate Output of [estimateScores()], or NULL.
#' @return A list with programs (program, top_enrichment, purity_rho, purity_fdr, microenv,
#'   and stromal_rho, immune_rho, top_genotype when available), enrichment and genotype.
#' @export
annotatePrograms <- function(H_ref, oof, samples, gene_sets, mut = NULL, estimate = NULL, alpha = 0.05) {
  programs <- rownames(H_ref)
  out <- data.frame(program = programs, top_enrichment = NA_character_, stringsAsFactors = FALSE)
  enrich <- if (length(gene_sets)) enrichHypergeom(H_ref, gene_sets) else NULL
  if (!is.null(enrich)) {
    sig <- enrich[enrich$fdr < alpha, , drop = FALSE]
    sig <- sig[order(sig$program, sig$p), , drop = FALSE]
    top <- sig[!duplicated(sig$program), , drop = FALSE]
    out$top_enrichment <- top$gene_set[match(programs, top$program)]
  }
  tum <- samples$profile_id[samples$role == "tumor" & samples$profile_id %in% rownames(oof)]
  A <- oof[tum, programs, drop = FALSE]
  pur <- programAssociation(A, samples$purity_dna[match(tum, samples$profile_id)], "purity_dna")
  out$purity_rho <- pur$rho
  out$purity_fdr <- pur$fdr
  out$microenv <- pur$fdr < alpha & pur$rho < 0
  if (!is.null(estimate)) {
    out$stromal_rho <- programAssociation(A, estimate$stromal[match(tum, estimate$profile_id)], "stromal")$rho
    out$immune_rho <- programAssociation(A, estimate$immune[match(tum, estimate$profile_id)], "immune")$rho
  }
  geno <- NULL
  if (!is.null(mut)) {
    rows <- intersect(tum, rownames(mut))
    geno <- genotypeAssociation(A[rows, , drop = FALSE], mut[rows, , drop = FALSE])
    gs <- geno[geno$fdr < alpha, , drop = FALSE]
    gs <- gs[order(gs$program, gs$p), , drop = FALSE]
    gt <- gs[!duplicated(gs$program), , drop = FALSE]
    i <- match(programs, gt$program)
    out$top_genotype <- ifelse(is.na(i), NA_character_, paste0(gt$gene[i], ifelse(gt$delta[i] > 0, "(+)", "(-)")))
  }
  list(programs = out, enrichment = enrich, genotype = geno)
}

#' combineArm
#'
#' Selection, matching, persistence scoring and annotation for one gene-graph arm. Settings
#' with a missing, failed or non-converged fold are listed in `selection` and never scored.
#' @param results Task results for the arm. @param expected Expected settings.
#' @param samples Sample table (included samples, with fold). @param pairs From [pairTable()].
#' @param gene_sets,mut,estimate Annotation inputs (see [annotatePrograms()]).
#' @param n_perm,n_boot Persistence draws. @param n_perm_match Matching-null draws.
#' @param seed Random seed. @param cores Cores used to score settings in parallel.
#' @param n_folds Number of folds.
#' @return A list with selection, selected, tau_selected, programs, maps_long, matches,
#'   per_type, enrichment, genotype, ref_H, oof, convergence.
#' @export
combineArm <- function(results, expected, samples, pairs, gene_sets = list(), mut = NULL, estimate = NULL,
                       n_perm = 10000, n_boot = 2000, n_perm_match = 1000, seed = 20261015, cores = 1,
                       n_folds = 5) {
  ok <- Filter(function(r) is.null(r$error) && !is.null(r$fit), results)
  tab <- settingTable(results, expected, n_folds)
  by_setting <- split(ok, vapply(ok, function(r) r$meta$setting_id, character(1)))
  ids <- tab$setting_id[tab$complete]
  matched <- stats::setNames(lapply(ids, function(id) matchSetting(by_setting[[id]], n_perm_match, seed)), ids)
  tab$reproducibility <- NA_real_
  tab$tau <- NA_real_
  tab$reproducibility[match(ids, tab$setting_id)] <- vapply(matched, function(m) m$reproducibility, numeric(1))
  tab$tau[match(ids, tab$setting_id)] <- vapply(matched, function(m) m$tau, numeric(1))
  selected <- selectSetting(tab)
  tab$selected <- tab$setting_id == selected
  scored <- parallel::mclapply(ids, function(id) {
    folds <- Filter(function(r) r$meta$fold > 0, by_setting[[id]])
    oof <- assembleOutOfFold(folds, matched[[id]]$matches)
    list(oof = oof, scores = scoreSetting(oof, pairs, samples, matched[[id]]$stability, n_perm, n_boot, seed))
  }, mc.cores = cores)
  bad <- vapply(scored, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("combineArm: scoring failed for ", paste(ids[bad], collapse = ", "), ": ", as.character(scored[[which(bad)[1]]]))
  names(scored) <- ids
  ref_H <- Filter(function(r) r$meta$fold == 0, by_setting[[selected]])[[1]]$fit$H
  oof <- scored[[selected]]$oof
  ann <- annotatePrograms(ref_H, oof, samples, gene_sets, mut, estimate)
  st <- matched[[selected]]$stability
  st$program <- paste0("LV", st$ref)
  programs <- merge(st, scored[[selected]]$scores, by = "program", all.x = TRUE)
  programs <- merge(programs, ann$programs, by = "program", all.x = TRUE)
  programs <- programs[order(programs$ref), , drop = FALSE]
  maps_long <- do.call(rbind, lapply(ids, function(id) {
    s <- scored[[id]]$scores
    if (nrow(s)) cbind(setting_id = id, s, stringsAsFactors = FALSE) else NULL
  }))
  field <- function(r, f) if (is.null(r$fit)) NA else r$fit[[f]]
  convergence <- do.call(rbind, lapply(results, function(r) data.frame(
    task_id = if (is.null(r$meta$task_id)) NA_character_ else r$meta$task_id,
    setting_id = r$meta$setting_id, fold = r$meta$fold,
    n_iter = field(r, "n_iter"), converged = field(r, "converged"),
    final_total = if (is.null(r$fit)) NA_real_ else r$fit$final_loss$total,
    error = if (is.null(r$error)) NA_character_ else r$error, stringsAsFactors = FALSE)))
  list(selection = tab, selected = selected, tau_selected = matched[[selected]]$tau, programs = programs,
       maps_long = maps_long, matches = matched[[selected]]$matches,
       per_type = perTypeConcordance(oof, pairs, samples, st$program[st$stable]),
       enrichment = ann$enrichment, genotype = ann$genotype, ref_H = ref_H, oof = oof, convergence = convergence)
}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `make document && make test FILTER=combine`
Expected: PASS (3 tests).

- [ ] **Step 6: Commit**

```bash
git add R/combine.R NAMESPACE man tests/testthat/helper-combine.R tests/testthat/test-combine.R
git commit -m "feat: combineArm, selection, matching, scoring and annotation for one arm

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 15: Compare gene priors

**Files:**
- Create: `R/compare.R`, `tests/testthat/test-compare.R`

**Interfaces:**
- Produces:
  - `transferTable(results)`. Each result's meta has `arm, variant ("selected"/"baseline"), fold`, and its fit has `recon_test`. The table has columns `arm, selected_mean, selected_se, baseline_mean, baseline_se, diff_mean, diff_se, n_folds`.
  - `cohenKappa(a, b)`.
  - `mapAgreement(H_a, H_b, programs_a, programs_b, tau)` returns `list(pairs, summary)`. `programs_*` have `program, stable, group, conc_rho`. The summary has columns `n_matched_stable, agreement, kappa, concordance_cor, unmatched_a, unmatched_b`.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-compare.R`:

```r
test_that("transferTable pairs folds within each arm", {
  mk <- function(arm, variant, fold, recon) list(meta = list(arm = arm, variant = variant, fold = fold),
                                                 fit = list(recon_test = recon), error = NULL)
  res <- c(lapply(1:5, function(f) mk("a", "selected", f, 0.10 + f / 1000)),
           lapply(1:5, function(f) mk("a", "baseline", f, 0.12 + f / 1000)),
           lapply(1:5, function(f) mk("b", "selected", f, 0.20)),
           lapply(1:5, function(f) mk("b", "baseline", f, 0.20)))
  tt <- transferTable(res)
  expect_equal(tt$arm, c("a", "b"))
  expect_equal(tt$diff_mean, c(-0.02, 0))
  expect_equal(tt$n_folds, c(5L, 5L))
})

test_that("cohenKappa matches a hand calculation", {
  expect_equal(cohenKappa(c("x", "x", "y", "y"), c("x", "y", "y", "y")), 0.5)
  expect_equal(cohenKappa(c("x", "y"), c("x", "y")), 1)
})

test_that("mapAgreement matches programs across arms and scores agreement", {
  withr::with_seed(15, H <- normalizeRows(matrix(stats::rexp(3 * 80), 3, dimnames = list(paste0("LV", 1:3), NULL))))
  pa <- data.frame(program = paste0("LV", 1:3), stable = TRUE, group = c("persists", "drops_out", "persists"), conc_rho = c(.8, .1, .6))
  Hb <- H[c(2, 3, 1), ]; rownames(Hb) <- paste0("LV", 1:3)
  pb <- data.frame(program = paste0("LV", 1:3), stable = TRUE, group = c("drops_out", "persists", "persists"), conc_rho = c(.1, .6, .8))
  ag <- mapAgreement(H, Hb, pa, pb, tau = 0.5)
  expect_equal(ag$summary$n_matched_stable, 3L)
  expect_equal(ag$summary$agreement, 1)
  expect_equal(ag$summary$kappa, 1)
  expect_equal(ag$pairs$program_b, c("LV3", "LV1", "LV2"))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=compare`
Expected: FAIL with `could not find function "transferTable"`.

- [ ] **Step 3: Implement**

`R/compare.R`:

```r
#' transferTable
#'
#' Tumor-to-organoid transfer reconstruction per arm (spec Section 8.4). The gene prior's
#' effect is the paired per-fold difference between the selected setting and alpha_2 = 0.
#' @param results Transfer-task results (meta: arm, variant, fold; fit: recon_test).
#' @return Data frame with arm, selected_mean/se, baseline_mean/se, diff_mean/se, n_folds.
#' @export
transferTable <- function(results) {
  d <- do.call(rbind, lapply(results, function(r) data.frame(
    arm = r$meta$arm, variant = r$meta$variant, fold = r$meta$fold,
    recon = if (is.null(r$fit)) NA_real_ else r$fit$recon_test, stringsAsFactors = FALSE)))
  paired <- merge(d[d$variant == "selected", ], d[d$variant == "baseline", ],
                  by = c("arm", "fold"), suffixes = c("_selected", "_baseline"))
  paired$diff <- paired$recon_selected - paired$recon_baseline
  se <- function(x) stats::sd(x) / sqrt(length(x))
  arms <- sort(unique(paired$arm))
  agg <- function(x, f) as.numeric(tapply(x, paired$arm, f)[arms])
  data.frame(arm = arms,
             selected_mean = agg(paired$recon_selected, mean), selected_se = agg(paired$recon_selected, se),
             baseline_mean = agg(paired$recon_baseline, mean), baseline_se = agg(paired$recon_baseline, se),
             diff_mean = agg(paired$diff, mean), diff_se = agg(paired$diff, se),
             n_folds = as.integer(table(paired$arm)[arms]), stringsAsFactors = FALSE)
}

#' cohenKappa
#' @param a,b Equal-length label vectors. @return Cohen's kappa (NA for empty input).
#' @export
cohenKappa <- function(a, b) {
  if (!length(a)) return(NA_real_)
  lv <- union(a, b)
  tab <- table(factor(a, lv), factor(b, lv))
  n <- sum(tab)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  if (pe == 1) return(1)
  (po - pe) / (1 - pe)
}

#' mapAgreement
#'
#' Matches two arms' reference programs (Hungarian method on H cosine, threshold `tau`) and
#' compares group calls and concordance among programs that are matched and stable in both.
#' @param H_a,H_b Reference H of each arm (same genes, same k).
#' @param programs_a,programs_b Program tables with program, stable, group, conc_rho.
#' @param tau Matching threshold.
#' @return A list with pairs (one row per program of arm a) and summary.
#' @export
mapAgreement <- function(H_a, H_b, programs_a, programs_b, tau) {
  m <- matchPrograms(H_a, H_b, tau)
  m$program_a <- rownames(H_a)[m$ref]
  m$program_b <- rownames(H_b)[m$fold_prog]
  ga <- programs_a[match(m$program_a, programs_a$program), , drop = FALSE]
  gb <- programs_b[match(m$program_b, programs_b$program), , drop = FALSE]
  pairs <- data.frame(m, stable_a = ga$stable %in% TRUE, stable_b = gb$stable %in% TRUE,
                      group_a = ga$group, group_b = gb$group, conc_a = ga$conc_rho, conc_b = gb$conc_rho,
                      stringsAsFactors = FALSE)
  sub <- pairs[pairs$matched & pairs$stable_a & pairs$stable_b & !is.na(pairs$group_a) & !is.na(pairs$group_b), , drop = FALSE]
  list(pairs = pairs,
       summary = data.frame(n_matched_stable = nrow(sub),
                            agreement = if (nrow(sub)) mean(sub$group_a == sub$group_b) else NA_real_,
                            kappa = cohenKappa(sub$group_a, sub$group_b),
                            concordance_cor = if (nrow(sub) > 2) stats::cor(sub$conc_a, sub$conc_b) else NA_real_,
                            unmatched_a = paste(m$program_a[!m$matched], collapse = ";"),
                            unmatched_b = paste(setdiff(rownames(H_b), m$program_b[m$matched]), collapse = ";"),
                            stringsAsFactors = FALSE))
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `make document && make test FILTER=compare`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add R/compare.R NAMESPACE man tests/testthat/test-compare.R
git commit -m "feat: transfer test table and 2D/3D map agreement

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 16: Synthetic fixture and the end-to-end test

**Files:**
- Create: `R/fixture.R`, `tests/testthat/test-e2e-fixture.R`

**Interfaces:**
- Produces:
  - `simulatePersistenceFixture(n_pairs_per_type = 20, n_types = 3, n_unpaired = 15, n_genes = 500, noise_sd = 0.2, seed = 20261015)` returns `list(Y, samples, gene_graph, blocks, k, expected)`.
  - `fixtureCalls(combined, fixture)` returns a named character vector of groups for P1..P4. Task 21 reuses both.

- [ ] **Step 1: Write the failing test**

`tests/testthat/test-e2e-fixture.R`:

```r
test_that("the fixture has the planted structure", {
  fx <- simulatePersistenceFixture()
  expect_equal(dim(fx$Y), c(75 + 60, 500))
  expect_equal(sum(fx$samples$role == "organoid"), 60)
  expect_true(all(fx$Y >= 0 & fx$Y <= 1))
  expect_equal(nrow(pairTable(fx$samples)), 60)
})

test_that("fixture programs land in their expected persistence groups", {
  skip_if_not(Sys.getenv("GRACKLE_RUN_E2E") == "true", "set GRACKLE_RUN_E2E=true (make e2e) to run")
  skip_if_no_tf()
  fx <- simulatePersistenceFixture()
  grid <- expand.grid(alpha_2 = c(0, 0.1), alpha_1 = c(0, 0.1))[, c("alpha_1", "alpha_2")]
  grid$setting_id <- sprintf("a%02d", seq_len(nrow(grid)))
  tasks <- merge(grid, data.frame(fold = 0:5))
  results <- lapply(seq_len(nrow(tasks)), function(i) {
    t <- tasks[i, ]
    list(meta = list(task_id = paste(t$setting_id, t$fold), arm = "fixture", mode = "main",
                     setting_id = t$setting_id, alpha_1 = t$alpha_1, alpha_2 = t$alpha_2, fold = t$fold),
         fit = fitTask(fx$Y, fx$samples, fx$gene_graph, fx$k, t$alpha_1, t$alpha_2, t$fold),
         error = NULL)
  })
  out <- combineArm(results, grid, fx$samples, pairTable(fx$samples),
                    gene_sets = list(BLOCK_P1 = colnames(fx$Y)[fx$blocks$P1]),
                    n_perm = 999, n_boot = 200, n_perm_match = 200)
  expect_equal(fixtureCalls(out, fx), fx$expected)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=e2e`
Expected: FAIL with `could not find function "simulatePersistenceFixture"`.

- [ ] **Step 3: Implement**

`R/fixture.R`:

```r
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
#' @param combined Output of [combineArm()]. @param fixture Output of [simulatePersistenceFixture()].
#' @return Named character vector of groups for P1..P4.
#' @export
fixtureCalls <- function(combined, fixture) {
  H <- combined$ref_H
  prog <- vapply(fixture$blocks, function(b) rownames(H)[which.max(rowMeans(H[, b, drop = FALSE]))], character(1))
  stats::setNames(combined$programs$group[match(prog, combined$programs$program)], names(fixture$blocks))
}
```

- [ ] **Step 4: Run the quick test, then the end-to-end test**

Run: `make document && make test FILTER=e2e`
Expected: the structure test passes and the end-to-end test is skipped.

Run: `make e2e`
Expected: PASS (24 GRACKLE fits; about 5-10 minutes).

If a planted program lands in the wrong group:
1. Print `out$programs` and `fixtureCalls(out, fx)`. Check whether that program is unstable (group `NA`) or mis-scored.
2. If it is unstable at k = 5, check `out$matches` for the program's per-fold cosines.
3. Do not change the expected groups or the thresholds. Find the cause in the code under test and report it.

- [ ] **Step 5: Commit**

```bash
git add R/fixture.R NAMESPACE man tests/testthat/test-e2e-fixture.R
git commit -m "test: synthetic fixture with planted programs and end-to-end check

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 17: Parse HCMI, cBioPortal and DepMap inputs

**Files:**
- Create: `R/hcmi_data.R`, `R/depmap_data.R`, `tests/testthat/test-hcmi-data.R`, `tests/testthat/test-depmap-data.R`

**Interfaces:**
- Produces (`R/hcmi_data.R`):
  - `readStarCounts(path)`: a data frame with `gene_id, gene_name, gene_type, tpm_unstranded`.
  - `buildExpressionMatrix(files, hgnc)`: profiles x HGNC symbols, log2(TPM + 1). `files` is a named character vector (names are profile ids); `hgnc` has `symbol, ensembl_gene_id`.
  - `selectGenes(logexpr, n = 5000)`.
  - `buildSampleTable(gdc, clinical)` returns `list(samples, filter_counts)`:
    - `gdc` columns: `profile_id, case_id, case_submitter_id, sample_submitter_id, sample_type, tumor_descriptor, preservation_method, passage_count`.
    - `clinical` columns: `sample_id, growth_pattern, culture_media, treatment_status, cancer_type, age, sex, race_category, purity_absolute, purity_ascat, purity_hatchet, purity_purple, consensus_purity` (all character).
    - `samples` gains `role, sample_id, model_type, medium, passage, derivation_center, model_id, parent_sample_id, pair_flag, cancer_type, prior_treatment, age, sex, ancestry, purity_dna, consensus_purity, included, exclusion_reason`.
- Produces (`R/depmap_data.R`): `depmapSymbols(cols, hgnc)`, `readDepmapEffect(path, hgnc, genes_in_rows = FALSE)` (genes x models), `readNonessential(path, hgnc)`, `hcmiTokenFromRow(df)`, `combine3dSources(E_dm, E_sg = NULL, min_cor = 0.8)` (returns `list(effect, decision)`). `hgnc` also needs `entrez_id`.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-hcmi-data.R`:

```r
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
```

`tests/testthat/test-depmap-data.R`:

```r
hgnc_dm <- data.frame(symbol = c("A1BG", "KRAS"), entrez_id = c(1, 3845), ensembl_gene_id = c("E1", "E2"))

test_that("DepMap gene columns map to HGNC symbols", {
  expect_equal(depmapSymbols(c("A1BG (1)", "KRAS (3845)", "FAKE (999999)", "NOSYM"), hgnc_dm),
               c("A1BG", "KRAS", "FAKE", "NOSYM"))
})

test_that("readDepmapEffect returns genes x models", {
  f <- tempfile(fileext = ".csv")
  writeLines(c(",A1BG (1),KRAS (3845)", "ACH-1,0.1,-1.2", "ACH-2,0.2,-0.8"), f)
  E <- readDepmapEffect(f, hgnc_dm)
  expect_equal(dim(E), c(2, 2))
  expect_equal(rownames(E), c("A1BG", "KRAS"))
  expect_equal(colnames(E), c("ACH-1", "ACH-2"))
  expect_equal(E["KRAS", "ACH-2"], -0.8)
  g <- tempfile(fileext = ".csv")
  writeLines(c("gene,ORG-1,ORG-2", "A1BG,0.1,0.2", "KRAS,-1.2,-0.8"), g)
  E2 <- readDepmapEffect(g, hgnc_dm, genes_in_rows = TRUE)
  expect_equal(rownames(E2), c("A1BG", "KRAS")); expect_equal(colnames(E2), c("ORG-1", "ORG-2"))
  expect_equal(E2["KRAS", "ORG-2"], -0.8)
})

test_that("readNonessential parses the control gene list", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("Gene", "A1BG (1)", "KRAS (3845)"), f)
  expect_equal(readNonessential(f, hgnc_dm), c("A1BG", "KRAS"))
})

test_that("hcmiTokenFromRow finds HCMI ids in any text column", {
  df <- data.frame(ModelID = c("ACH-1", "ACH-2", "ACH-3"),
                   CellLineName = c("HCM-BROD-0416-C71", "HT29", NA),
                   Alias = c(NA, NA, "HCM-SANG-0270-C20-A"), stringsAsFactors = FALSE)
  expect_equal(hcmiTokenFromRow(df), c("HCM-BROD-0416-C71", NA, "HCM-SANG-0270-C20"))
})

test_that("combine3dSources combines agreeing sources and otherwise keeps the larger one", {
  withr::with_seed(16, {
    base <- stats::rnorm(100)
    E1 <- matrix(base + stats::rnorm(100 * 4, sd = 0.2), 100, dimnames = list(paste0("g", 1:100), paste0("dm", 1:4)))
    E2 <- matrix(base + stats::rnorm(100 * 3, sd = 0.2), 100, dimnames = list(paste0("g", 1:100), paste0("sg", 1:3)))
    E3 <- matrix(stats::rnorm(100 * 6), 100, dimnames = list(paste0("g", 1:100), paste0("sg", 1:6)))
  })
  a <- combine3dSources(E1, E2)
  expect_equal(ncol(a$effect), 7); expect_match(a$decision, "combined")
  b <- combine3dSources(E1, E3)
  expect_equal(colnames(b$effect), colnames(E3)); expect_match(b$decision, "Sanger alone")
  expect_equal(combine3dSources(E1)$decision, "DepMap 3D only")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `make test FILTER=hcmi-data && make test FILTER=depmap-data`
Expected: FAIL with `could not find function "readStarCounts"` and `"depmapSymbols"`.

- [ ] **Step 3: Implement `R/hcmi_data.R`**

```r
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
#' @param logexpr Profiles x genes matrix. @param n Number of genes. Default 5000.
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
#' partner out. Organoids inherit their parent tumor's annotations; tumors get DNA purity as
#' the median of the four DNA callers.
#' @param gdc GDC file metadata (see Task 17 interfaces). @param clinical cBioPortal fields per sample.
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
  ann <- data.frame(cancer_type = cl$cancer_type,
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
```

- [ ] **Step 4: Implement `R/depmap_data.R`**

```r
## DepMap and Sanger CRISPR inputs (spec Section 4.5).

#' depmapSymbols
#' @param cols DepMap gene column names such as "KRAS (3845)". @param hgnc Data frame with symbol, entrez_id.
#' @return HGNC symbols from the Entrez id, falling back to the text before " (".
#' @export
depmapSymbols <- function(cols, hgnc) {
  entrez <- ifelse(grepl("\\([0-9]+\\)$", cols), sub("^.*\\(([0-9]+)\\)$", "\\1", cols), NA_character_)
  sym <- hgnc$symbol[match(entrez, as.character(hgnc$entrez_id), incomparables = NA)]
  ifelse(is.na(sym), trimws(sub("\\s*\\(.*$", "", cols)), sym)
}

#' readDepmapEffect
#' @param path Gene-effect CSV whose first column holds ids: model ids (models in rows, the DepMap
#'   layout) or gene names (`genes_in_rows = TRUE`). @param hgnc HGNC table.
#' @param genes_in_rows TRUE when genes are rows and models are columns.
#' @return Genes x models matrix (duplicate symbols: the first is kept).
#' @export
readDepmapEffect <- function(path, hgnc, genes_in_rows = FALSE) {
  x <- data.table::fread(path, data.table = FALSE)
  m <- as.matrix(x[, -1, drop = FALSE])
  rownames(m) <- x[[1]]
  if (genes_in_rows) m <- t(m)
  sym <- depmapSymbols(colnames(m), hgnc)
  keep <- !is.na(sym) & !duplicated(sym)
  out <- t(m[, keep, drop = FALSE])
  rownames(out) <- sym[keep]
  out
}

#' readNonessential
#' @param path DepMap non-essential control list (first column: "SYMBOL (ENTREZ)"). @param hgnc HGNC table.
#' @return HGNC symbols.
#' @export
readNonessential <- function(path, hgnc) {
  x <- data.table::fread(path, data.table = FALSE)
  depmapSymbols(as.character(x[[1]]), hgnc)
}

#' hcmiTokenFromRow
#' @param df Data frame of model metadata. @return The first HCMI case id found in each row's text columns, or NA.
#' @export
hcmiTokenFromRow <- function(df) {
  chr <- vapply(df, is.character, logical(1))
  tok <- do.call(cbind, lapply(df[chr], hcmiToken))
  apply(tok, 1, function(r) {
    r <- r[!is.na(r)]
    if (length(r)) r[1] else NA_character_
  })
}

#' combine3dSources
#'
#' Combines DepMap NextGen and Sanger 3D gene effects when they agree: the correlation of
#' per-gene means over shared genes is at least `min_cor`. In that case each model is z-scored
#' within its source. Otherwise the source with more models is used alone.
#' @param E_dm,E_sg Genes x models matrices (E_sg may be NULL). @param min_cor Agreement threshold.
#' @return A list with effect and decision (a sentence recorded in gene_graph_summary.csv).
#' @export
combine3dSources <- function(E_dm, E_sg = NULL, min_cor = 0.8) {
  if (is.null(E_sg) || !ncol(E_sg)) return(list(effect = E_dm, decision = "DepMap 3D only"))
  genes <- intersect(rownames(E_dm), rownames(E_sg))
  r <- stats::cor(rowMeans(E_dm[genes, , drop = FALSE], na.rm = TRUE),
                  rowMeans(E_sg[genes, , drop = FALSE], na.rm = TRUE), use = "complete.obs")
  if (is.finite(r) && r >= min_cor) {
    z <- function(E) {
      out <- scale(E)
      attr(out, "scaled:center") <- NULL
      attr(out, "scaled:scale") <- NULL
      out
    }
    return(list(effect = cbind(z(E_dm[genes, , drop = FALSE]), z(E_sg[genes, , drop = FALSE])),
                decision = sprintf("combined (gene-mean r = %.2f)", r)))
  }
  if (ncol(E_sg) > ncol(E_dm)) {
    list(effect = E_sg, decision = sprintf("Sanger alone (gene-mean r = %.2f < %.2f)", r, min_cor))
  } else {
    list(effect = E_dm, decision = sprintf("DepMap alone (gene-mean r = %.2f < %.2f)", r, min_cor))
  }
}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `make document && make test FILTER=hcmi-data && make test FILTER=depmap-data`
Expected: PASS (6 + 5 tests).

- [ ] **Step 6: Run the whole suite**

Run: `make test`
Expected: every test passes. The end-to-end test is skipped unless `make e2e` is used.

- [ ] **Step 7: Commit**

```bash
git add R/hcmi_data.R R/depmap_data.R NAMESPACE man tests/testthat/test-hcmi-data.R tests/testthat/test-depmap-data.R
git commit -m "feat: parse GDC STAR counts, the HCMI sample table and DepMap inputs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 18: Alpine access through `ralpine`, then push and create the worktree

**Files:**
- Create: `scripts/alpine/ralpine` (copied from fm-pm-evaluator, then edited)
- Create: `tests/testthat/test-ralpine-boundary.R`

**Interfaces:**
- Produces: `scripts/alpine/ralpine`, with the same verbs as in fm-pm-evaluator (`status, sq, sacct, ls, cat, tail, log, rev, du, run, pull, update, switch, submit, cancel, jobinfo`). Its root is the GRACKLE worktree on Alpine. Every later Alpine step goes through it.

- [ ] **Step 1: Copy and adapt `ralpine`**

```bash
mkdir -p scripts/alpine
cp ~/Repositories/fm-pdo-evaluator/scripts/alpine/ralpine scripts/alpine/ralpine
chmod +x scripts/alpine/ralpine
```

Edit `scripts/alpine/ralpine` in three places:
1. Replace the line `REMOTE_ROOT="${ALPINE_ROOT:-/projects/lgillenwater@xsede.org/repositories/fm-pdo-evaluator}"` with:

```bash
REMOTE_ROOT="${ALPINE_ROOT:-/projects/lgillenwater@xsede.org/repositories/GRACKLE-organoid-transfer}"
```

2. Directly under the first header line `# Alpine access for the harness workflow: read-only, plus job submission.`, add:

```bash
# GRACKLE copy: taken from fm-pm-evaluator (lagillenwater/fm-pm-evaluator, commit dff524f). The
# root is the organoid-transfer worktree; the shared checkout at repositories/GRACKLE (detached
# HEAD, uncommitted HTP work) is never touched.
```

3. In the `submit` verb, change `git -C "$LOCAL_REPO" status --porcelain -- src scripts` to `git -C "$LOCAL_REPO" status --porcelain -- R analysis scripts`. In the warning line just below it, change `src/ or scripts/` to `R/, analysis/ or scripts/`.

- [ ] **Step 2: Write the boundary tests (ported from fm-pm-evaluator's `test_ralpine_boundary.py`)**

`tests/testthat/test-ralpine-boundary.R`:

```r
ralpine_path <- function() testthat::test_path("..", "..", "scripts", "alpine", "ralpine")
ralpine_text <- function() paste(readLines(ralpine_path(), warn = FALSE), collapse = "\n")
first_group <- function(pattern, text) {
  m <- regmatches(text, regexec(pattern, text, perl = TRUE))[[1]]
  if (length(m) < 2) NA_character_ else m[2]
}
case_block <- function(verb) first_group(sprintf("(?s)\\n  %s\\)\\n(.*?)\\n    ;;\\n", verb), ralpine_text())

test_that("READ_ONLY excludes every escape-capable command", {
  skip_if_not(file.exists(ralpine_path()))
  block <- first_group("(?s)READ_ONLY=\\(\\n(.*?)\\n\\)", ralpine_text())
  expect_false(is.na(block))
  tokens <- strsplit(trimws(block), "\\s+")[[1]]
  escape <- c("find", "scontrol", "file", "nvidia-smi", "git", "bash", "sh", "xargs", "python", "perl", "R", "Rscript")
  expect_equal(intersect(escape, tokens), character(0))
})

test_that("log validates its pattern before calling remote_fixed", {
  skip_if_not(file.exists(ralpine_path()))
  b <- case_block("log")
  expect_false(is.na(b))
  v <- regexpr("^[A-Za-z0-9._-]*$", b, fixed = TRUE)
  r <- regexpr('remote_fixed "', b, fixed = TRUE)
  expect_true(v > 0 && r > 0 && v < r)
})

test_that("jobinfo validates the job id before calling remote_fixed and only shows the job", {
  skip_if_not(file.exists(ralpine_path()))
  b <- case_block("jobinfo")
  expect_false(is.na(b))
  v <- regexpr("^[0-9]+$", b, fixed = TRUE)
  r <- regexpr('remote_fixed "', b, fixed = TRUE)
  expect_true(v > 0 && r > 0 && v < r)
  expect_true(grepl("scontrol show job", b, fixed = TRUE))
})

test_that("remote() enforces the allowlist and rejects metacharacters", {
  skip_if_not(file.exists(ralpine_path()))
  body <- first_group("(?s)\\nremote\\(\\) \\{(.*?)\\n\\}\\n", ralpine_text())
  expect_false(is.na(body))
  expect_true(grepl("is_read_only", body, fixed = TRUE))
  expect_true(grepl("reject_metacharacters", body, fixed = TRUE))
})

test_that("the default root is the organoid-transfer worktree, not the shared checkout", {
  skip_if_not(file.exists(ralpine_path()))
  expect_true(grepl("repositories/GRACKLE-organoid-transfer}", ralpine_text(), fixed = TRUE))
  expect_false(grepl("repositories/GRACKLE}", ralpine_text(), fixed = TRUE))
})
```

- [ ] **Step 3: Run the tests**

Run: `make test FILTER=ralpine`
Expected: PASS (5 tests).

- [ ] **Step 4: Commit and push the branch**

```bash
git add scripts/alpine/ralpine tests/testthat/test-ralpine-boundary.R
git commit -m "build: ralpine for GRACKLE, rooted at the organoid-transfer worktree

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push -u origin organoid-transfer
```

Expected: the branch is published to `origin` (`lagillenwater/GRACKLE`). Pushing to this branch falls under the standing permission.

- [ ] **Step 5: HUMAN STEP: open the Alpine session and create the worktree**

Ask the user to run these two commands at the prompt (the `!` prefix runs them in this session):

```
! ssh -fN alpine
! ssh alpine 'git -C /projects/lgillenwater@xsede.org/repositories/GRACKLE fetch origin && git -C /projects/lgillenwater@xsede.org/repositories/GRACKLE worktree add ../GRACKLE-organoid-transfer organoid-transfer'
```

The first opens the Duo-authenticated master connection; skip it if `scripts/alpine/ralpine status` already says the connection is up. The second creates a worktree next to the shared checkout. It does not modify that checkout's working tree.

- [ ] **Step 6: Verify the worktree**

Run: `scripts/alpine/ralpine rev`
Expected: the first line equals the output of `git rev-parse HEAD` run locally, and the status shows `## organoid-transfer...origin/organoid-transfer` with nothing modified.

---

### Task 19: The environment on Alpine

**Files:**
- Create: `scripts/alpine/env.sh`, `scripts/alpine/setup_env.sbatch`, `scripts/alpine/test.sbatch`
- Create (pulled from Alpine): `environment.lock.linux-64.txt`

**Interfaces:**
- Produces: `/projects/$USER/software/envs/grackle` on Alpine, plus `scripts/alpine/env.sh`, which every job sources. It exports `GRACKLE_ENV`, `RETICULATE_PYTHON`, `GRACKLE_DATA` and `GRACKLE_SCRATCH`, and puts the environment's `bin` first on `PATH`.

- [ ] **Step 1: Write `scripts/alpine/env.sh`**

```bash
# Sourced by every GRACKLE Alpine job. Puts the micromamba environment on PATH: `module load
# anaconda` fails inside batch jobs (see fm-pm-evaluator's ralpine notes), and putting an
# environment's bin first on PATH is the pattern proven to work. Also sets the data and scratch
# roots read by analysis/hcmi_persistence/config.R.
export GRACKLE_ENV="/projects/$USER/software/envs/grackle"
export PATH="$GRACKLE_ENV/bin:$PATH"
export RETICULATE_PYTHON="$GRACKLE_ENV/bin/python"
export R_LIBS_USER="$GRACKLE_ENV/lib/R/library"
export GRACKLE_DATA="/projects/$USER/grackle_data/hcmi"
export GRACKLE_SCRATCH="/scratch/alpine/$USER/grackle/hcmi_persistence"
export XDG_CACHE_HOME="/scratch/alpine/$USER/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export TF_NUM_INTRAOP_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export TF_NUM_INTEROP_THREADS=1
export TF_CPP_MIN_LOG_LEVEL=3
mkdir -p logs "$XDG_CACHE_HOME"
```

- [ ] **Step 2: Write `scripts/alpine/setup_env.sbatch`**

```bash
#!/bin/bash
#SBATCH --job-name=grackle-setup-env
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Builds /projects/$USER/software/envs/grackle from environment.yml with a standalone micromamba,
# then writes environment.lock.linux-64.txt into this checkout. Pull it back and commit it.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; mkdir -p logs
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
export XDG_CACHE_HOME="/scratch/alpine/$USER/cache"; mkdir -p "$XDG_CACHE_HOME"
scripts/dev/setup_env.sh "/projects/$USER/software/envs/grackle" "/projects/$USER/software/bin"
```

- [ ] **Step 3: Write `scripts/alpine/test.sbatch`**

```bash
#!/bin/bash
#SBATCH --job-name=grackle-test
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=7G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Runs the testthat suite on an Alpine compute node in the pinned environment. Set E2E=true to
# include the end-to-end fixture test.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
if [[ "${E2E:-false}" == "true" ]]; then make e2e GRACKLE_ENV="$GRACKLE_ENV"; fi
make test GRACKLE_ENV="$GRACKLE_ENV"
```

Run: `bash -n scripts/alpine/setup_env.sbatch && bash -n scripts/alpine/test.sbatch && bash -n scripts/alpine/env.sh`
Expected: no output.

- [ ] **Step 4: Commit, push, deploy and build the environment**

```bash
git add scripts/alpine/env.sh scripts/alpine/setup_env.sbatch scripts/alpine/test.sbatch
git commit -m "build: Alpine environment and test jobs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push
scripts/alpine/ralpine update
scripts/alpine/ralpine submit scripts/alpine/setup_env.sbatch
```

Expected: `update` prints the new short HEAD; `submit` prints `Submitted batch job <id>`.

- [ ] **Step 5: Wait for the build and check it**

Check with `scripts/alpine/ralpine sacct <id>` every few minutes; the build takes 20 to 60 minutes. When the state is `COMPLETED`, run `scripts/alpine/ralpine log grackle-setup-env`.
Expected: the log ends with `tensorflow 2.21.<x> ok` and `environment ready: /projects/lgillenwater@xsede.org/software/envs/grackle (linux-64)`.

If the log shows the micromamba download or solve failing, stop and report the error to the user. The fallback (a pinned Apptainer image, spec Section 14) is a separate decision.

- [ ] **Step 6: Pull the Linux lockfile, commit it, and run the tests on Alpine**

```bash
scripts/alpine/ralpine pull /projects/lgillenwater@xsede.org/repositories/GRACKLE-organoid-transfer/environment.lock.linux-64.txt environment.lock.linux-64.txt
git add environment.lock.linux-64.txt
git commit -m "build: Linux lockfile from the Alpine environment build

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push
scripts/alpine/ralpine update
scripts/alpine/ralpine submit scripts/alpine/test.sbatch --export=ALL,E2E=true
```

`update` moves the untracked lockfile on Alpine aside before fast-forwarding, because the incoming commit adds the same path.

When the test job finishes, run: `scripts/alpine/ralpine log grackle-test`
Expected: `[ FAIL 0 | WARN <n> | SKIP 0 | PASS <n> ]` for the end-to-end run and for the full suite.

---

### Task 20: Pipeline scripts, part 1: configuration and data steps 01 to 04

**Files:**
- Create: `analysis/hcmi_persistence/config.R`
- Create: `analysis/hcmi_persistence/01_download.R`, `02_sample_table.R`, `03_expression.R`, `03b_gene_graphs.R`, `04_rank.R`
- Create: `scripts/alpine/01_download.sbatch`, `02_sample_table.sbatch`, `03_expression.sbatch`, `03b_gene_graphs.sbatch`, `04_rank.sbatch`

**Interfaces:**
- Consumes: the package functions from Tasks 3 to 17.
- Produces, on Alpine, under `$GRACKLE_DATA`:
  - `raw/`: the downloads plus `raw/manifest.csv`.
  - `processed/samples_included.rds`: the included sample table, with `fold`.
  - `processed/sample_table.csv`, `folds.csv`, `filter_counts.csv`, `field_coverage.csv`.
  - `processed/mutations.rds`.
  - `processed/expression_log2tpm.rds`, `Y.rds`, `genes.txt`, `estimate.csv`.
  - `processed/gene_graphs/{string_mech,string_all,coess_2d,coess_3d_fold0..5}.rds`.
  - `processed/gene_graph_summary.csv`.
  - `processed/k.rds` and `k_parallel_analysis.csv`.
- Reads two committed tables that Task 22 creates: `analysis/hcmi_persistence/depmap_files.csv` (columns `role, filename, release`; roles `gene_effect`, `models`, `nonessential`, and optionally `gene_effect_3d`) and `analysis/hcmi_persistence/depmap_growth_patterns.csv` (columns `value, class`, where class is `2D`, `3D` or `exclude`). `analysis/hcmi_persistence/sanger_files.csv` (same columns plus an optional `genes_in_rows` TRUE/FALSE column; role `gene_effect_3d`) is optional.

- [ ] **Step 1: Write `analysis/hcmi_persistence/config.R`**

```r
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
```

- [ ] **Step 2: Write `analysis/hcmi_persistence/01_download.R`**

```r
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
muts <- request(sprintf("%s/molecular-profiles/pancan_hcmi_2025_mutations/mutations", cb)) |>
  req_url_query(sampleListId = "pancan_hcmi_2025_sequenced", projection = "DETAILED") |>
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

## 5. Manual downloads: DepMap (required) and Sanger (optional) -------------------------------
dm_tab <- file.path(cfg$analysis, "depmap_files.csv")
if (!file.exists(dm_tab)) stop("DepMap files are a manual download (the portal blocks scripts): see plan Task 22")
dm <- as.data.frame(fread(dm_tab))
missing_roles <- setdiff(c("gene_effect", "models", "nonessential"), dm$role)
if (length(missing_roles)) stop("depmap_files.csv lacks roles: ", paste(missing_roles, collapse = ", "))
dm_paths <- file.path(cfg$manual, "depmap", dm$filename)
if (!all(file.exists(dm_paths))) stop("missing DepMap files on Alpine: ", paste(dm$filename[!file.exists(dm_paths)], collapse = ", "))
invisible(lapply(dm_paths, note, source = paste("manual DepMap download", unique(dm$release))))
models <- as.data.frame(fread(dm_paths[dm$role == "models"]))
writeLines(names(models), file.path(cfg$manual, "depmap", "model_columns.txt"))
if ("GrowthPattern" %in% names(models)) {
  fwrite(as.data.frame(table(value = models$GrowthPattern, useNA = "ifany")),
         file.path(cfg$manual, "depmap", "growth_pattern_values.csv"))
}
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
```

- [ ] **Step 3: Write `analysis/hcmi_persistence/02_sample_table.R`**

```r
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
coverage$purity_vs_consensus_r <- c(suppressWarnings(stats::cor(inc$purity_dna, inc$consensus_purity, use = "complete.obs")), NA, NA)
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
```

- [ ] **Step 4: Write `analysis/hcmi_persistence/03_expression.R`**

```r
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
```

- [ ] **Step 5: Write `analysis/hcmi_persistence/03b_gene_graphs.R`**

```r
## Step 03b: the four gene-graph arms (spec Section 4.5). coess_3d is rebuilt per fold.
source("analysis/hcmi_persistence/config.R")
genes <- readLines(file.path(cfg$processed, "genes.txt"))
samples <- readRDS(file.path(cfg$processed, "samples_included.rds"))
hgnc <- as.data.frame(fread(file.path(cfg$raw, "annotation", "hgnc_complete_set_2024-10-01.txt"), quote = ""))
gg_dir <- file.path(cfg$processed, "gene_graphs"); dir.create(gg_dir, recursive = TRUE, showWarnings = FALSE)
summ <- list()
record <- function(arm, fold, G, n_models = NA, n_informative = NA, decision = NA) {
  summ[[length(summ) + 1L]] <<- cbind(data.frame(arm = arm, fold = fold, n_models = n_models,
                                                 n_informative = n_informative, decision = decision), graphSummary(G))
}

## STRING arms
links <- readStringLinks(file.path(cfg$raw, "string", "9606.protein.links.full.v12.0.txt.gz"))
pmap <- stringProteinToSymbol(readStringAliases(file.path(cfg$raw, "string", "9606.protein.aliases.v12.0.txt.gz")), hgnc)
invisible(lapply(c("mech", "all"), function(mode) {
  G <- stringGraph(links, pmap, genes, mode)
  saveRDS(G, file.path(gg_dir, sprintf("string_%s.rds", mode)))
  record(paste0("string_", mode), NA, G)
}))
rm(links); invisible(gc())

## Co-essentiality arms
dm <- as.data.frame(fread(file.path(cfg$analysis, "depmap_files.csv")))
gp_map <- file.path(cfg$analysis, "depmap_growth_patterns.csv")
if (!file.exists(gp_map)) stop("depmap_growth_patterns.csv is missing: see plan Task 22")
gp <- as.data.frame(fread(gp_map))
dm_path <- function(role) file.path(cfg$manual, "depmap", dm$filename[dm$role == role])
effect <- readDepmapEffect(dm_path("gene_effect"), hgnc)
if ("gene_effect_3d" %in% dm$role) {
  e3 <- readDepmapEffect(dm_path("gene_effect_3d"), hgnc)
  shared <- intersect(rownames(effect), rownames(e3))
  effect <- cbind(effect[shared, , drop = FALSE], e3[shared, setdiff(colnames(e3), colnames(effect)), drop = FALSE])
}
models <- as.data.frame(fread(dm_path("models")))
models$class <- gp$class[match(models$GrowthPattern, gp$value)]
models$hcmi <- hcmiTokenFromRow(models)
nonessential <- readNonessential(dm_path("nonessential"), hgnc)
arm <- stats::setNames(geneArm(hgnc$location), hgnc$symbol)

ids_2d <- models$ModelID[models$class %in% "2D" & is.na(models$hcmi)]
E2 <- effect[, intersect(colnames(effect), ids_2d), drop = FALSE]
G2 <- coessGraph(E2, genes, nonessential, arm)
saveRDS(G2, file.path(gg_dir, "coess_2d.rds"))
record("coess_2d", NA, G2, ncol(E2), length(informativeGenes(E2, nonessential)))

ids_3d <- models$ModelID[models$class %in% "3D"]
E3_dm <- effect[, intersect(colnames(effect), ids_3d), drop = FALSE]
E3_sg <- NULL
sg_tab <- file.path(cfg$analysis, "sanger_files.csv")
if (file.exists(sg_tab)) {
  sg <- as.data.frame(fread(sg_tab))
  row <- sg$role == "gene_effect_3d"
  E3_sg <- readDepmapEffect(file.path(cfg$manual, "sanger", sg$filename[row]), hgnc,
                            genes_in_rows = "genes_in_rows" %in% names(sg) && isTRUE(as.logical(sg$genes_in_rows[row])))
  colnames(E3_sg) <- paste0("SANGER:", colnames(E3_sg))
}
combined_3d <- combine3dSources(E3_dm, E3_sg)
E3 <- combined_3d$effect
model_token <- c(stats::setNames(models$hcmi, models$ModelID)[colnames(E3_dm)],
                 if (!is.null(E3_sg)) stats::setNames(hcmiToken(sub("^SANGER:", "", colnames(E3_sg))), colnames(E3_sg)))
model_case <- stats::setNames(samples$case_id[match(model_token, samples$case_submitter_id, incomparables = NA)], names(model_token))
model_case <- model_case[colnames(E3)]
invisible(lapply(0:5, function(f) {
  held <- if (f == 0) character(0) else unique(samples$case_id[samples$fold == f])
  Ef <- dropHeldOutModels(E3, model_case, held)
  G3 <- coessGraph(Ef, genes, nonessential, arm)
  saveRDS(G3, file.path(gg_dir, sprintf("coess_3d_fold%d.rds", f)))
  record("coess_3d", f, G3, ncol(Ef), length(informativeGenes(Ef, nonessential)), combined_3d$decision)
}))
fwrite(rbindlist(summ, fill = TRUE), file.path(cfg$processed, "gene_graph_summary.csv"))
writeProvenance(file.path(cfg$processed, "gene_graph_summary.csv"),
                inputs = c(file.path(cfg$processed, "genes.txt"), dm_path("gene_effect"), dm_path("models")))
print(rbindlist(summ, fill = TRUE))
```

- [ ] **Step 6: Write `analysis/hcmi_persistence/04_rank.R`**

```r
## Step 04: rank k by Horn's parallel analysis (spec Section 7.1).
source("analysis/hcmi_persistence/config.R")
Y <- readRDS(file.path(cfg$processed, "Y.rds"))
pa <- parallelAnalysisRank(Y, n_perm = 100, quantile = 0.95, nv = 200, seed = cfg$seed)
saveRDS(pa$k, file.path(cfg$processed, "k.rds"))
fwrite(data.frame(component = seq_along(pa$observed), observed = pa$observed, threshold = pa$threshold),
       file.path(cfg$processed, "k_parallel_analysis.csv"))
writeProvenance(file.path(cfg$processed, "k.rds"), inputs = file.path(cfg$processed, "Y.rds"), extra = list(k = pa$k))
message("k = ", pa$k)
```

- [ ] **Step 7: Write the job scripts**

Each script has the same header. Only the job name, resources and command differ.

`scripts/alpine/01_download.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=grackle-01-download
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
Rscript analysis/hcmi_persistence/01_download.R
```

`scripts/alpine/02_sample_table.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=grackle-02-samples
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3840M
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
Rscript analysis/hcmi_persistence/02_sample_table.R
```

`scripts/alpine/03_expression.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=grackle-03-expression
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
Rscript analysis/hcmi_persistence/03_expression.R
```

`scripts/alpine/03b_gene_graphs.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=grackle-03b-genegraphs
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
Rscript analysis/hcmi_persistence/03b_gene_graphs.R
```

`scripts/alpine/04_rank.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=grackle-04-rank
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
Rscript analysis/hcmi_persistence/04_rank.R
```

- [ ] **Step 8: Check syntax**

Run:

```bash
for f in config 01_download 02_sample_table 03_expression 03b_gene_graphs 04_rank; do
  ~/.local/envs/grackle/bin/Rscript -e "invisible(parse('analysis/hcmi_persistence/$f.R'))" || echo "PARSE ERROR in $f"
done
for f in 01_download 02_sample_table 03_expression 03b_gene_graphs 04_rank; do bash -n scripts/alpine/$f.sbatch; done
```

Expected: no output.

- [ ] **Step 9: Commit**

```bash
git add analysis/hcmi_persistence/config.R analysis/hcmi_persistence/01_download.R analysis/hcmi_persistence/02_sample_table.R analysis/hcmi_persistence/03_expression.R analysis/hcmi_persistence/03b_gene_graphs.R analysis/hcmi_persistence/04_rank.R scripts/alpine/01_download.sbatch scripts/alpine/02_sample_table.sbatch scripts/alpine/03_expression.sbatch scripts/alpine/03b_gene_graphs.sbatch scripts/alpine/04_rank.sbatch
git commit -m "feat: HCMI pipeline steps 01-04 and their Alpine jobs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 21: Pipeline scripts, part 2: fits, combine, figures and the job chain

**Files:**
- Create: `analysis/hcmi_persistence/make_tasks.R`, `analysis/hcmi_persistence/tasks/main_<arm>.csv` (4 files, generated)
- Create: `analysis/hcmi_persistence/05_fit.R`, `06_combine.R`, `07_figures.R`
- Create: `analysis/hcmi_persistence/fixture/run_fixture.R`
- Create: `scripts/alpine/fit_array.sbatch`, `scripts/alpine/combine.sbatch`, `scripts/alpine/figures.sbatch`, `scripts/alpine/submit_chain.sh`

**Interfaces:**
- Consumes: `processed/*` from Task 20, and `fitTask`, `combineArm`, `transferTable`, `mapAgreement`, `simulatePersistenceFixture`, `fixtureCalls`.
- Produces:
  - Task tables with columns `task_id, arm, mode, setting_id, alpha_1, alpha_2, fold, k_mult, variant`.
  - Fit results in `<run>/fits/<arm>/<mode>/<task_id>.rds`.
  - Results in `<run>/results/<arm>/{combined.rds, programs.csv, selection.csv, convergence.csv, maps_all_settings.csv, per_type.csv, enrichment.csv, genotype.csv, sensitivity_k0.5/, sensitivity_k2/}`, `<run>/results/compare/{transfer.csv, map_agreement.csv, map_agreement_pairs.csv, power_curve.csv}` and `<run>/results/figures/`.

- [ ] **Step 1: Write and run `make_tasks.R`**

`analysis/hcmi_persistence/make_tasks.R`:

```r
## Writes the main task tables: 36 alpha settings x (reference + 5 folds) for each arm (216 rows each).
arms <- c("coess_2d", "coess_3d", "string_mech", "string_all")
alphas <- c(0, 0.01, 0.03, 0.1, 0.3, 1)
grid <- expand.grid(alpha_2 = alphas, alpha_1 = alphas)[, c("alpha_1", "alpha_2")]
grid$setting_id <- sprintf("a%02d", seq_len(nrow(grid)))
dir.create("analysis/hcmi_persistence/tasks", showWarnings = FALSE)
invisible(lapply(arms, function(arm) {
  t <- merge(grid, data.frame(fold = 0:5))
  t <- t[order(t$setting_id, t$fold), ]
  out <- data.frame(task_id = sprintf("%s_%s_f%d", arm, t$setting_id, t$fold), arm = arm, mode = "main",
                    setting_id = t$setting_id, alpha_1 = t$alpha_1, alpha_2 = t$alpha_2, fold = t$fold,
                    k_mult = 1, variant = NA)
  utils::write.csv(out, sprintf("analysis/hcmi_persistence/tasks/main_%s.csv", arm), row.names = FALSE, na = "")
}))
```

Run: `~/.local/envs/grackle/bin/Rscript analysis/hcmi_persistence/make_tasks.R && wc -l analysis/hcmi_persistence/tasks/*.csv`
Expected: each of the four files has 217 lines (a header plus 216 tasks).

- [ ] **Step 2: Write `05_fit.R`**

```r
## Step 05: one cross-fitting task (spec Section 8.1); one Slurm array index per task.
##   Rscript analysis/hcmi_persistence/05_fit.R --tasks <tasks.csv> --index <1-based row>
## A task whose result exists without an error is skipped, so failed indices can be resubmitted alone.
source("analysis/hcmi_persistence/config.R")
args <- parseArgs(list(tasks = "", index = as.character(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")) + 1L)))
tasks <- as.data.frame(fread(args$tasks))
t <- tasks[as.integer(args$index), , drop = FALSE]
if (!nrow(t) || is.na(t$task_id)) stop("no task at index ", args$index, " of ", args$tasks)
out <- file.path(cfg$run, "fits", t$arm, t$mode, paste0(t$task_id, ".rds"))
if (file.exists(out) && is.null(readRDS(out)$error)) {
  message("already done: ", out)
  quit(save = "no", status = 0)
}
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
Y <- readRDS(file.path(cfg$processed, "Y.rds"))
samples <- readRDS(file.path(cfg$processed, "samples_included.rds"))
samples <- samples[match(rownames(Y), samples$profile_id), , drop = FALSE]
k <- as.integer(round(readRDS(file.path(cfg$processed, "k.rds")) * t$k_mult))
graph_file <- file.path(cfg$processed, "gene_graphs",
                        if (t$arm == "coess_3d") sprintf("coess_3d_fold%d.rds", t$fold) else paste0(t$arm, ".rds"))
G <- readRDS(graph_file)[colnames(Y), colnames(Y)]
transfer <- t$mode == "transfer"
meta <- c(as.list(t), k = k)
res <- tryCatch(
  list(fit = fitTask(Y, samples, G, k, t$alpha_1, t$alpha_2, t$fold,
                     train_role = if (transfer) "tumor" else "all",
                     test_role = if (transfer) "organoid" else "all"),
       error = NULL),
  error = function(e) list(fit = NULL, error = conditionMessage(e)))
saveRDS(c(list(meta = meta), res), out)
writeProvenance(out, inputs = c(file.path(cfg$processed, c("Y.rds", "samples_included.rds", "k.rds")), graph_file),
                extra = list(task = meta))
if (!is.null(res$error)) {
  message("fit failed: ", res$error)
  quit(save = "no", status = 1)
}
message(sprintf("%s: n_iter %d, converged %s", t$task_id, res$fit$n_iter, res$fit$converged))
```

- [ ] **Step 3: Write `06_combine.R`**

```r
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
```

- [ ] **Step 4: Write `07_figures.R`**

REQUIRED: load the `dataviz` skill before writing this file, and apply its palette and mark guidance to the colors and themes below. Keep the panel structure; it is spec Section 10.

```r
## Step 07: main and supplementary figures (spec Section 10).
source("analysis/hcmi_persistence/config.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
res_dir <- file.path(cfg$run, "results")
fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
group_shapes <- c(persists = 16, drops_out = 17, reduced_concordant = 15, present_not_patient_specific = 1)
compartment_cols <- c(`tumor-intrinsic` = "#2F6690", microenvironment = "#C2593F")
theme_set(theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank()))
read_csv_df <- function(...) as.data.frame(fread(file.path(...)))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
mapPanel <- function(dir, title) {
  p <- read_csv_df(dir, "programs.csv")
  p <- p[p$stable %in% TRUE, , drop = FALSE]
  p$compartment <- ifelse(p$microenv %in% TRUE, "microenvironment", "tumor-intrinsic")
  p$label <- ifelse(is.na(p$top_enrichment) | p$top_enrichment == "", p$program, sub("^HALLMARK_", "", p$top_enrichment))
  ggplot(p, aes(ret_median_log2, conc_rho)) +
    geom_vline(xintercept = -1, linetype = "dashed", linewidth = 0.3, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
    geom_errorbar(aes(ymin = conc_ci_low, ymax = conc_ci_high), width = 0, linewidth = 0.25, colour = "grey60") +
    geom_point(aes(colour = compartment, shape = group), size = 2) +
    geom_text(aes(label = label), size = 2.2, vjust = -0.9, check_overlap = TRUE) +
    scale_colour_manual(values = compartment_cols) + scale_shape_manual(values = group_shapes, drop = FALSE) +
    labs(title = title, x = "Retention: median log2(organoid / tumor)", y = "Concordance (stratified Spearman)")
}
agreementPanel <- function() {
  a <- read_csv_df(res_dir, "compare", "map_agreement_pairs.csv")
  s <- read_csv_df(res_dir, "compare", "map_agreement.csv")
  a <- a[a$matched & a$stable_a & a$stable_b, , drop = FALSE]
  a$call <- ifelse(a$group_a == a$group_b, "same group", "different group")
  ggplot(a, aes(conc_a, conc_b, colour = call)) + geom_abline(linewidth = 0.3, colour = "grey75") + geom_point(size = 2) +
    scale_colour_manual(values = c(`same group` = "#2F6690", `different group` = "#C2593F")) +
    labs(title = sprintf("2D vs 3D co-essentiality: kappa %.2f (n = %d)", s$kappa, s$n_matched_stable),
         x = "Concordance, coess_2d", y = "Concordance, coess_3d", colour = NULL)
}
heatmapPanel <- function(arm) {
  d <- read_csv_df(res_dir, arm, "per_type.csv")
  ggplot(d, aes(cancer_type, program, fill = rho)) + geom_tile() +
    scale_fill_gradient2(low = "#C2593F", mid = "white", high = "#2F6690", midpoint = 0, limits = c(-1, 1)) +
    labs(title = paste(arm, "concordance by cancer type"), x = NULL, y = NULL)
}

save_fig((mapPanel(file.path(res_dir, "coess_2d"), "coess_2d") | mapPanel(file.path(res_dir, "coess_3d"), "coess_3d") |
            agreementPanel()) + plot_layout(guides = "collect"), "main_persistence_maps", 13, 4.5)
save_fig(heatmapPanel("coess_2d") | heatmapPanel("coess_3d"), "main_per_type_heatmaps", 10, 6)
save_fig(mapPanel(file.path(res_dir, "string_mech"), "string_mech") | mapPanel(file.path(res_dir, "string_all"), "string_all"),
         "supp_string_maps", 9, 4.5)
tr <- read_csv_df(res_dir, "compare", "transfer.csv")
save_fig(ggplot(tr, aes(arm, diff_mean)) + geom_hline(yintercept = 0, colour = "grey75") +
           geom_pointrange(aes(ymin = diff_mean - 1.96 * diff_se, ymax = diff_mean + 1.96 * diff_se)) +
           labs(x = NULL, y = "Tumor-to-organoid error: selected minus no gene prior"), "supp_transfer", 5, 3.5)
invisible(lapply(cfg$arms, function(arm) {
  m <- read_csv_df(res_dir, arm, "maps_all_settings.csv")
  save_fig(ggplot(m, aes(ret_median_log2, conc_rho, shape = group)) + geom_point(size = 1.1) + facet_wrap(~setting_id) +
             scale_shape_manual(values = group_shapes) + labs(title = paste(arm, ": every alpha setting")),
           paste0("supp_alpha_grid_", arm), 12, 10)
  sens <- file.path(res_dir, arm, c("sensitivity_k0.5", "sensitivity_k2"))
  sens <- sens[file.exists(file.path(sens, "programs.csv"))]
  if (length(sens)) save_fig(Reduce(`|`, lapply(sens, function(d) mapPanel(d, paste(arm, basename(d))))),
                             paste0("supp_rank_sensitivity_", arm), 4.5 * length(sens), 4.5)
}))
pw <- read_csv_df(res_dir, "compare", "power_curve.csv")
save_fig(ggplot(pw, aes(rho, power)) + geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey60") +
           geom_line() + geom_point() + labs(x = "True concordance (rho)", y = "Power, one-sided alpha 0.05",
                                             title = sprintf("Concordance test power at %d held-out pairs", pw$n_pairs[1])),
         "supp_power_curve", 5, 3.5)
invisible(lapply(cfg$arms, function(arm) {
  p <- read_csv_df(res_dir, arm, "programs.csv")
  p$program <- factor(p$program, levels = p$program[order(p$ref)])
  save_fig(ggplot(p, aes(program, n_matched, fill = stable)) + geom_col() + geom_hline(yintercept = 4, linetype = "dashed") +
             scale_fill_manual(values = c(`TRUE` = "#2F6690", `FALSE` = "#BBBBBB")) +
             labs(title = paste(arm, ": folds matched per reference program"), x = NULL, y = "Folds matched (of 5)") +
             theme(axis.text.x = element_text(angle = 90, size = 5)),
           paste0("supp_stability_", arm), 10, 3.5)
}))
pa <- read_csv_df(cfg$processed, "k_parallel_analysis.csv")
save_fig(ggplot(pa, aes(component)) + geom_line(aes(y = observed, colour = "observed")) +
           geom_line(aes(y = threshold, colour = "null, 95th percentile")) + scale_y_log10() +
           labs(y = "Singular value", colour = NULL,
                title = sprintf("Parallel analysis: k = %d", readRDS(file.path(cfg$processed, "k.rds")))),
         "supp_parallel_analysis", 6, 4)
message("figures written to ", fig_dir)
```

- [ ] **Step 5: Write the fixture runner `analysis/hcmi_persistence/fixture/run_fixture.R`**

```r
## Runs the real 05_fit.R and 06_combine.R scripts end to end on the synthetic fixture (make fixture).
## Exits with status 1 unless every planted program lands in its expected group.
root <- normalizePath(file.path("analysis", "hcmi_persistence", "fixture", "out"), mustWork = FALSE)
unlink(root, recursive = TRUE)
Sys.setenv(GRACKLE_DATA = file.path(root, "data"), GRACKLE_SCRATCH = file.path(root, "scratch"), GRACKLE_RUN_ID = "fixture")
source("analysis/hcmi_persistence/config.R")
dir.create(file.path(cfg$processed, "gene_graphs"), recursive = TRUE, showWarnings = FALSE)
fx <- simulatePersistenceFixture()
saveRDS(fx$Y, file.path(cfg$processed, "Y.rds"))
saveRDS(fx$samples, file.path(cfg$processed, "samples_included.rds"))
saveRDS(fx$k, file.path(cfg$processed, "k.rds"))
saveRDS(fx$gene_graph, file.path(cfg$processed, "gene_graphs", "coess_2d.rds"))
saveRDS(list(BLOCK_P1 = colnames(fx$Y)[fx$blocks$P1]), file.path(cfg$processed, "hallmark.rds"))
grid <- expand.grid(alpha_2 = c(0, 0.1), alpha_1 = c(0, 0.1))[, c("alpha_1", "alpha_2")]
grid$setting_id <- sprintf("a%02d", seq_len(nrow(grid)))
t <- merge(grid, data.frame(fold = 0:5))
t <- t[order(t$setting_id, t$fold), ]
tasks <- data.frame(task_id = sprintf("coess_2d_%s_f%d", t$setting_id, t$fold), arm = "coess_2d", mode = "main",
                    setting_id = t$setting_id, alpha_1 = t$alpha_1, alpha_2 = t$alpha_2, fold = t$fold,
                    k_mult = 1, variant = NA)
tasks_file <- file.path(root, "tasks_main_coess_2d.csv")
utils::write.csv(tasks, tasks_file, row.names = FALSE, na = "")
rscript <- file.path(R.home("bin"), "Rscript")
status <- vapply(seq_len(nrow(tasks)), function(i)
  system2(rscript, c("analysis/hcmi_persistence/05_fit.R", "--tasks", tasks_file, "--index", i)), integer(1))
if (any(status != 0)) stop("05_fit.R failed for tasks: ", paste(which(status != 0), collapse = ", "))
if (system2(rscript, c("analysis/hcmi_persistence/06_combine.R", "--arm", "coess_2d", "--stage", "main",
                       "--tasks", tasks_file, "--n_perm", "999", "--n_boot", "200", "--n_perm_match", "200")) != 0) {
  stop("06_combine.R failed")
}
out <- readRDS(file.path(cfg$run, "results", "coess_2d", "combined.rds"))
calls <- fixtureCalls(out, fx)
print(rbind(expected = fx$expected, observed = calls[names(fx$expected)]))
if (!identical(unname(calls[names(fx$expected)]), unname(fx$expected))) quit(save = "no", status = 1)
cat("fixture PASS\n")
```

- [ ] **Step 6: Run the fixture through the real scripts**

Run: `make fixture`
Expected: a two-row table of expected and observed groups that match, then `fixture PASS`. This runs 24 `05_fit.R` processes and one `06_combine.R`, about 10-20 minutes.

This is the check required before any Alpine submission (spec Section 13). If it fails, debug locally with superpowers:systematic-debugging. Do not submit to Alpine until it passes.

- [ ] **Step 7: Write the array, combine and figure job scripts**

`scripts/alpine/fit_array.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=grackle-fit
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=01:30:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.out
#
# One cross-fitting task per array index. ARM and KIND (main | sensitivity | transfer) choose the
# task table; GRACKLE_RUN_ID names the run. The array size is set at submission: --array=0-215%50
# for main, 0-11 for sensitivity, 0-9 for transfer. A completed task exits at once, so a failed
# index is resubmitted alone with --array=<index>.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
: "${ARM:?ARM not set}" "${KIND:?KIND not set}" "${GRACKLE_RUN_ID:?GRACKLE_RUN_ID not set}"
echo "commit=$(git rev-parse --short HEAD) job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} arm=$ARM kind=$KIND run=$GRACKLE_RUN_ID"
if [[ "$KIND" == "main" ]]; then
  TASKS="analysis/hcmi_persistence/tasks/main_${ARM}.csv"
else
  TASKS="$GRACKLE_SCRATCH/$GRACKLE_RUN_ID/tasks/${KIND}_${ARM}.csv"
fi
Rscript analysis/hcmi_persistence/05_fit.R --tasks "$TASKS" --index "$((SLURM_ARRAY_TASK_ID + 1))"
```

`scripts/alpine/combine.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=grackle-combine
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Runs 06_combine.R with ARM and STAGE (main | sensitivity | compare). Submitted with
# --dependency=afterany:<array>, so it runs even if some fits failed; those settings are excluded
# and listed in selection.csv.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
: "${ARM:?ARM not set}" "${STAGE:?STAGE not set}" "${GRACKLE_RUN_ID:?GRACKLE_RUN_ID not set}"
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID arm=$ARM stage=$STAGE run=$GRACKLE_RUN_ID"
Rscript analysis/hcmi_persistence/06_combine.R --arm "$ARM" --stage "$STAGE" --cores "$SLURM_CPUS_PER_TASK"
```

`scripts/alpine/figures.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=grackle-figures
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=7G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
: "${GRACKLE_RUN_ID:?GRACKLE_RUN_ID not set}"
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID run=$GRACKLE_RUN_ID"
Rscript analysis/hcmi_persistence/07_figures.R
```

- [ ] **Step 8: Write `scripts/alpine/submit_chain.sh` and make it executable**

```bash
#!/usr/bin/env bash
# Submits the HCMI persistence pipeline as a Slurm dependency chain through ralpine.
# Run it from the repository root on the laptop. Each dependency is read back with
# `ralpine jobinfo`, because a chain that silently drops a dependency runs out of order.
#
#   scripts/alpine/submit_chain.sh --from prep       # 02 -> 03 -> 03b -> 04 (data already downloaded)
#   scripts/alpine/submit_chain.sh --from fit        # fit arrays, combine, sensitivity, transfer, compare, figures
#   RUN_ID=<id> scripts/alpine/submit_chain.sh --from compare   # compare and figures for an existing run
set -euo pipefail
RALPINE="$(dirname "${BASH_SOURCE[0]}")/ralpine"
[[ "${1:-}" == "--from" ]] || { echo "usage: $0 --from prep|fit|compare" >&2; exit 1; }
FROM="${2:?--from needs prep, fit or compare}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%d)_$(git rev-parse --short HEAD)}"
EXP="--export=ALL,GRACKLE_RUN_ID=$RUN_ID"
ARMS=(coess_2d coess_3d string_mech string_all)

job_id() { grep -oE '[0-9]+$' <<<"$1" | tail -1; }
submit() { local out; out="$("$RALPINE" submit "$@")"; echo "$out" >&2; job_id "$out"; }
check_dependency() {
  local id="$1" want="$2" dep
  dep="$("$RALPINE" jobinfo "$id" | grep -oE 'Dependency=[^ ]+' || true)"
  if [[ "$dep" != *"$want"* ]]; then
    echo "job $id was submitted without its dependency ($dep, wanted $want); cancel it with: $RALPINE cancel $id" >&2
    exit 1
  fi
}
echo "run_id: $RUN_ID"

if [[ "$FROM" == "prep" ]]; then
  J2=$(submit scripts/alpine/02_sample_table.sbatch)
  J3=$(submit scripts/alpine/03_expression.sbatch "--dependency=afterok:$J2"); check_dependency "$J3" "afterok:$J2"
  J3b=$(submit scripts/alpine/03b_gene_graphs.sbatch "--dependency=afterok:$J3"); check_dependency "$J3b" "afterok:$J3"
  J4=$(submit scripts/alpine/04_rank.sbatch "--dependency=afterok:$J3"); check_dependency "$J4" "afterok:$J3"
  echo "prep chain: 02=$J2 03=$J3 03b=$J3b 04=$J4"
  exit 0
fi

if [[ "$FROM" == "fit" ]]; then
  WAIT=()
  for ARM in "${ARMS[@]}"; do   # four arms, each its own chain
    F=$(submit scripts/alpine/fit_array.sbatch --array=0-215%50 "$EXP,ARM=$ARM,KIND=main")
    C=$(submit scripts/alpine/combine.sbatch "--dependency=afterany:$F" "$EXP,ARM=$ARM,STAGE=main"); check_dependency "$C" "afterany:$F"
    S=$(submit scripts/alpine/fit_array.sbatch --array=0-11 "--dependency=afterok:$C" "$EXP,ARM=$ARM,KIND=sensitivity"); check_dependency "$S" "afterok:$C"
    T=$(submit scripts/alpine/fit_array.sbatch --array=0-9 "--dependency=afterok:$C" "$EXP,ARM=$ARM,KIND=transfer"); check_dependency "$T" "afterok:$C"
    CS=$(submit scripts/alpine/combine.sbatch "--dependency=afterany:$S" "$EXP,ARM=$ARM,STAGE=sensitivity"); check_dependency "$CS" "afterany:$S"
    echo "$ARM: fits=$F combine=$C sensitivity=$S transfer=$T combine_sensitivity=$CS"
    WAIT+=("$C" "$T" "$CS")
  done
  DEP="afterany:$(IFS=:; echo "${WAIT[*]}")"
fi

if [[ "$FROM" == "fit" || "$FROM" == "compare" ]]; then
  CMP=$(submit scripts/alpine/combine.sbatch ${DEP:+"--dependency=$DEP"} "$EXP,ARM=all,STAGE=compare")
  [[ -n "${DEP:-}" ]] && check_dependency "$CMP" "afterany:"
  FIG=$(submit scripts/alpine/figures.sbatch "--dependency=afterok:$CMP" "$EXP"); check_dependency "$FIG" "afterok:$CMP"
  echo "compare=$CMP figures=$FIG"
fi
echo "watch with: $RALPINE sq"
```

Run: `chmod +x scripts/alpine/submit_chain.sh && bash -n scripts/alpine/submit_chain.sh && for f in fit_array combine figures; do bash -n scripts/alpine/$f.sbatch; done`
Expected: no output.

- [ ] **Step 9: Commit**

```bash
git add analysis/hcmi_persistence/make_tasks.R analysis/hcmi_persistence/tasks analysis/hcmi_persistence/05_fit.R analysis/hcmi_persistence/06_combine.R analysis/hcmi_persistence/07_figures.R analysis/hcmi_persistence/fixture/run_fixture.R scripts/alpine/fit_array.sbatch scripts/alpine/combine.sbatch scripts/alpine/figures.sbatch scripts/alpine/submit_chain.sh
git commit -m "feat: fit, combine and figure steps, the job chain and the fixture runner

make fixture runs the real 05/06 scripts on the synthetic fixture and
passes before any Alpine submission.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 22: Get the data on Alpine

**Files:**
- Create: `analysis/hcmi_persistence/depmap_files.csv`, `analysis/hcmi_persistence/depmap_growth_patterns.csv`
- Create (only if Sanger scores are obtained): `analysis/hcmi_persistence/sanger_files.csv`
- Create: `analysis/hcmi_persistence/data_sources.md`

**Interfaces:**
- Produces: a complete `raw/` on Alpine, verified by `01_download.R`, and the committed tables that `01_download.R` and `03b_gene_graphs.R` read (see Task 20 interfaces). This task also settles spec Section 15 items 5 and 6.

- [ ] **Step 1: HUMAN STEP: download the DepMap files in a browser**

Ask the user to download these from https://depmap.org/portal/data_page/ into `~/Downloads/depmap/`. Use the newest "DepMap Public" release whose release notes include the NextGen organoid and spheroid models (the 2026 expansion, Neiswender et al.):
- `CRISPRGeneEffect.csv`: Chronos gene effect. If the NextGen 3D screens are released as a separate gene-effect file, download that too.
- `Model.csv`: model metadata.
- `AchillesNonessentialControls.csv`: the non-essential control genes.

Optionally, Sanger HCMI organoid CRISPR screens (Herranz-Ors et al. 2026): only if the paper's data availability statement points to gene-level scores per model, saved to `~/Downloads/sanger/`.

The user should tell you the exact file names and the release name.

- [ ] **Step 2: HUMAN STEP: copy them to Alpine**

Ask the user to run:

```
! ssh alpine 'mkdir -p /projects/lgillenwater@xsede.org/grackle_data/hcmi/raw/manual/depmap /projects/lgillenwater@xsede.org/grackle_data/hcmi/raw/manual/sanger'
! rsync -av ~/Downloads/depmap/ alpine-dtn:/projects/lgillenwater@xsede.org/grackle_data/hcmi/raw/manual/depmap/
```

Add `! rsync -av ~/Downloads/sanger/ alpine-dtn:/projects/lgillenwater@xsede.org/grackle_data/hcmi/raw/manual/sanger/` if Sanger files were downloaded.

Then run: `scripts/alpine/ralpine ls /projects/lgillenwater@xsede.org/grackle_data/hcmi/raw/manual/depmap`
Expected: the downloaded files are listed with nonzero sizes.

- [ ] **Step 3: Write `depmap_files.csv` (and `sanger_files.csv` if used)**

`analysis/hcmi_persistence/depmap_files.csv`. Use the real file names and release from Step 1. If the NextGen screens came as a separate file, add a `gene_effect_3d` row:

```
role,filename,release
gene_effect,CRISPRGeneEffect.csv,<release name from Step 1>
models,Model.csv,<release name from Step 1>
nonessential,AchillesNonessentialControls.csv,<release name from Step 1>
```

If Sanger scores were obtained, check their layout first: `scripts/alpine/ralpine run head -c 400 /projects/lgillenwater@xsede.org/grackle_data/hcmi/raw/manual/sanger/<file>`. Then write `analysis/hcmi_persistence/sanger_files.csv`:

```
role,filename,release,genes_in_rows
gene_effect_3d,<file>,Herranz-Ors et al. 2026,<TRUE if the first column holds genes, FALSE if it holds models>
```

- [ ] **Step 4: Commit, deploy and run step 01**

```bash
git add analysis/hcmi_persistence/depmap_files.csv
git add analysis/hcmi_persistence/sanger_files.csv 2>/dev/null || true
git commit -m "data: DepMap (and Sanger) file roles for the manual downloads

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push
scripts/alpine/ralpine update
scripts/alpine/ralpine submit scripts/alpine/01_download.sbatch
```

When the job is `COMPLETED` (`scripts/alpine/ralpine sacct <id>`; expect 1-4 hours), run: `scripts/alpine/ralpine log grackle-01-download`
Expected: lines `GDC: <n> STAR files, md5 verified, release <release>` (n was 1153 on 2026-10-15) and `01_download complete: <m> files in raw/manifest.csv`. If the job failed, read the error at the end of the log. Network errors from a single source can be rerun; the script skips files it already has.

- [ ] **Step 5: Map DepMap growth patterns to 2D and 3D**

Run:

```bash
scripts/alpine/ralpine cat /projects/lgillenwater@xsede.org/grackle_data/hcmi/raw/manual/depmap/model_columns.txt
scripts/alpine/ralpine cat /projects/lgillenwater@xsede.org/grackle_data/hcmi/raw/manual/depmap/growth_pattern_values.csv
```

Check that `ModelID` and `GrowthPattern` appear in `model_columns.txt`; `03b_gene_graphs.R` reads them. If either is missing, stop and report to the user with the column list.

Write `analysis/hcmi_persistence/depmap_growth_patterns.csv` with one row per value in `growth_pattern_values.csv`:
- adherent, suspension and mixed-2D values map to `2D`;
- organoid, spheroid and neurosphere values map to `3D` (neurospheres are spheroid cultures; spec Section 4.5 defines 3D as organoid and spheroid screens);
- empty or unknown values map to `exclude`.

For example:

```
value,class
Adherent,2D
Suspension,2D
Mixed,2D
Organoid,3D
Spheroid,3D
Neurosphere,3D
Unknown,exclude
```

- [ ] **Step 6: Record the data sources**

Write `analysis/hcmi_persistence/data_sources.md` from `raw/manifest.csv` (`scripts/alpine/ralpine cat /projects/lgillenwater@xsede.org/grackle_data/hcmi/raw/manifest.csv`). Record:
- the GDC data release and STAR file count;
- the cBioPortal access date;
- STRING v12;
- the HGNC file (2024-10-01);
- the OncoKB download date and sha256;
- the msigdbr version;
- the DepMap release, files and sha256s;
- the growth-pattern mapping and why;
- whether Sanger screens are used.

It also records the outcome of spec Section 15 items 5 and 6.

```bash
git add analysis/hcmi_persistence/depmap_growth_patterns.csv analysis/hcmi_persistence/data_sources.md
git commit -m "data: growth-pattern map and data-source record

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push
scripts/alpine/ralpine update
```

---

### Task 23: Prep run (steps 02 to 04) and checkpoint

**Files:**
- Create (pulled from Alpine): `analysis/hcmi_persistence/sample_table.csv`, `analysis/hcmi_persistence/folds.csv`, `analysis/hcmi_persistence/filter_counts.csv`

**Interfaces:**
- Produces: everything in `processed/` (Task 20 interfaces) on Alpine, and the three committed tables spec Section 4.2 requires.

- [ ] **Step 1: Submit the prep chain**

Run: `scripts/alpine/ralpine rev`. The commit must equal local `git rev-parse HEAD`.
Run: `scripts/alpine/submit_chain.sh --from prep`
Expected: `prep chain: 02=<id> 03=<id> 03b=<id> 04=<id>`, and no dependency error.

- [ ] **Step 2: Check step 02 and apply the purity stop rule**

When `02` has finished, run: `scripts/alpine/ralpine log grackle-02-samples`
Expected:
- the coverage table, with `purity_dna, included tumors` having `n_available > 0` and `purity_vs_consensus_r` printed;
- the filter-count table;
- `pairs: <n>; flagged organoids: <m>`.

If the job failed with `no DNA-based tumor purity`, **stop and ask the user** (spec Section 4.6). Do not continue the chain; cancel the dependent jobs with `scripts/alpine/ralpine cancel <id>`, one id at a time.

- [ ] **Step 3: Check steps 03, 03b and 04**

When all have finished:

```bash
scripts/alpine/ralpine log grackle-03-expression
scripts/alpine/ralpine log grackle-03b-genegraphs
scripts/alpine/ralpine log grackle-04-rank
```

Expected:
- `Y: <n> profiles x 5000 genes`.
- A gene-graph summary with rows for `string_mech`, `string_all`, `coess_2d` and `coess_3d` folds 0-5, each with `n_genes_with_edges > 0`. The `coess_3d` rows carry the 3D-source decision.
- `k = <k>`.

- [ ] **Step 4: Pull the committed tables and the checkpoint numbers**

```bash
P=/projects/lgillenwater@xsede.org/grackle_data/hcmi/processed
for f in sample_table.csv folds.csv filter_counts.csv; do scripts/alpine/ralpine pull "$P/$f" "analysis/hcmi_persistence/$f"; done
mkdir -p analysis/hcmi_persistence/results/checkpoint
for f in field_coverage.csv gene_graph_summary.csv k_parallel_analysis.csv; do scripts/alpine/ralpine pull "$P/$f" "analysis/hcmi_persistence/results/checkpoint/$f"; done
git add analysis/hcmi_persistence/sample_table.csv analysis/hcmi_persistence/folds.csv analysis/hcmi_persistence/filter_counts.csv
git commit -m "data: HCMI sample table, folds and filter counts

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push
scripts/alpine/ralpine update
```

- [ ] **Step 5: CHECKPOINT: report to the user and wait**

Send the user one message containing:
- **Filter counts:** the final number of 3D-organoid pairs and tumors. Compare it with the 421 / 297 / 248 numbers in spec Section 2 for the grant text.
- **Field coverage:** purity and growth-pattern coverage, and the purity agreement with `CONSENSUS_PURITY`.
- **Gene graphs:** for each arm, the genes with edges and the informative genes, plus the 3D-source decision.
- **k:** the value from parallel analysis. The user asked to see it.

Ask for a go-ahead before Task 24; it uses about 65-190 core-hours.

---

### Task 24: Full run: fits, combine, sensitivity, transfer, compare, figures

**Files:** none created locally; the outputs go to Alpine scratch.

**Interfaces:**
- Produces: `/scratch/alpine/lgillenwater@xsede.org/grackle/hcmi_persistence/<RUN_ID>/results/` (Task 21 interfaces).

- [ ] **Step 1: Submit the fit chain**

Run: `scripts/alpine/ralpine rev`. The commit must equal local `git rev-parse HEAD`.
Run: `scripts/alpine/submit_chain.sh --from fit`
Expected: `run_id: <YYYYMMDD>_<hash>`, then one line per arm with five job ids, then `compare=<id> figures=<id>`. Record `RUN_ID` and every job id in the task notes.

- [ ] **Step 2: Watch the fit arrays**

Check `scripts/alpine/ralpine sq` periodically; each arm's 216 tasks run at most 50 at a time. When an arm's array is done, run `scripts/alpine/ralpine sacct <array id>` and count `FAILED` and `TIMEOUT` tasks.
- **No failures:** nothing to do. The combine job starts automatically.
- **TIMEOUT:**
  1. Raise `#SBATCH --time` in `scripts/alpine/fit_array.sbatch` (for example to `03:00:00`), then commit, push and `scripts/alpine/ralpine update`.
  2. Resubmit only the failed indices: `scripts/alpine/ralpine submit scripts/alpine/fit_array.sbatch --array=<i,j,...> --export=ALL,GRACKLE_RUN_ID=<RUN_ID>,ARM=<arm>,KIND=main`.
  3. Rerun that arm's combine, then its sensitivity and transfer steps, with the Step 3 commands.
- **FAILED:**
  1. Read the task log: `scripts/alpine/ralpine tail /projects/lgillenwater@xsede.org/repositories/GRACKLE-organoid-transfer/logs/grackle-fit-<array id>_<index>.out 60`.
  2. If it is a code error, use superpowers:systematic-debugging locally, add a failing test, fix, commit, push, update, then resubmit the failed indices as above.
  3. If it is a node or file-system error, resubmit the indices unchanged.

- [ ] **Step 3: Rerun an arm's tail after resubmitting fits (only when Step 2 resubmitted something)**

```bash
C=$(scripts/alpine/ralpine submit scripts/alpine/combine.sbatch --dependency=afterany:<resubmitted array id> --export=ALL,GRACKLE_RUN_ID=<RUN_ID>,ARM=<arm>,STAGE=main | grep -oE '[0-9]+$')
S=$(scripts/alpine/ralpine submit scripts/alpine/fit_array.sbatch --array=0-11 --dependency=afterok:$C --export=ALL,GRACKLE_RUN_ID=<RUN_ID>,ARM=<arm>,KIND=sensitivity | grep -oE '[0-9]+$')
T=$(scripts/alpine/ralpine submit scripts/alpine/fit_array.sbatch --array=0-9 --dependency=afterok:$C --export=ALL,GRACKLE_RUN_ID=<RUN_ID>,ARM=<arm>,KIND=transfer | grep -oE '[0-9]+$')
scripts/alpine/ralpine submit scripts/alpine/combine.sbatch --dependency=afterany:$S --export=ALL,GRACKLE_RUN_ID=<RUN_ID>,ARM=<arm>,STAGE=sensitivity
```

Then rerun compare and figures for the run: `RUN_ID=<RUN_ID> scripts/alpine/submit_chain.sh --from compare`. Cancel the original compare and figures jobs from Step 1 if they are still pending (`scripts/alpine/ralpine cancel <id>`).

Sensitivity and transfer tasks that already completed are skipped, because `05_fit.R` skips finished tasks. The fold-0 sensitivity reference and the transfer fits depend on the selected setting. If the rerun of combine selects a different setting, delete nothing: the new task tables have new `task_id`s, so the old results are left unused.

- [ ] **Step 4: Confirm the run finished**

Run: `scripts/alpine/ralpine log grackle-figures`
Expected: `figures written to /scratch/alpine/lgillenwater@xsede.org/grackle/hcmi_persistence/<RUN_ID>/results/figures`.

---

### Task 25: Pull, verify and summarize for review

**Files:**
- Create (gitignored until the user approves): `analysis/hcmi_persistence/results/<RUN_ID>/` (results, figures, logs and `SUMMARY.md`)

**Interfaces:**
- Produces: a summary the user reviews before any result is committed (spec Section 10; the results policy in Global Constraints).

- [ ] **Step 1: Pull the results and logs**

```bash
R=<RUN_ID>
scripts/alpine/ralpine pull /scratch/alpine/lgillenwater@xsede.org/grackle/hcmi_persistence/$R/results/ analysis/hcmi_persistence/results/$R/
scripts/alpine/ralpine pull /projects/lgillenwater@xsede.org/repositories/GRACKLE-organoid-transfer/logs/ analysis/hcmi_persistence/results/$R/logs/
```

- [ ] **Step 2: Run the verification checks**

REQUIRED: use superpowers:verification-before-completion. Every claim in the summary must come from this output.

Run (replace `<RUN_ID>`):

```bash
~/.local/envs/grackle/bin/Rscript -e '
r <- "analysis/hcmi_persistence/results/<RUN_ID>"
arms <- c("coess_2d", "coess_3d", "string_mech", "string_all")
for (a in arms) {
  sel <- read.csv(file.path(r, a, "selection.csv")); conv <- read.csv(file.path(r, a, "convergence.csv"))
  p <- read.csv(file.path(r, a, "programs.csv"))
  cat("\n==", a, "\n")
  cat("settings complete:", sum(sel$complete), "/", nrow(sel), "| all converged:", sum(sel$all_converged), "\n")
  cat("fits converged:", sum(conv$converged %in% TRUE), "/", nrow(conv), "| errors:", sum(!is.na(conv$error) & conv$error != ""), "\n")
  cat("selected:", sel$setting_id[sel$selected %in% TRUE], "alpha_1 =", sel$alpha_1[sel$selected %in% TRUE], "alpha_2 =", sel$alpha_2[sel$selected %in% TRUE], "\n")
  cat("stable programs:", sum(p$stable %in% TRUE), "/", nrow(p), "\n")
  print(table(group = p$group[p$stable %in% TRUE], microenv = p$microenv[p$stable %in% TRUE], useNA = "ifany"))
  for (k in c("sensitivity_k0.5", "sensitivity_k2")) {
    f <- file.path(r, a, k, "programs.csv")
    cat(k, if (file.exists(f)) paste("stable:", sum(read.csv(f)$stable %in% TRUE)) else "missing (see ERROR.txt)", "\n")
  }
}
cat("\n== transfer\n"); print(read.csv(file.path(r, "compare", "transfer.csv")))
cat("\n== 2D vs 3D agreement\n"); print(read.csv(file.path(r, "compare", "map_agreement.csv")))
cat("\n== concordance power at the real pair structure\n"); print(read.csv(file.path(r, "compare", "power_curve.csv")))
'
```

Check, and report each item honestly whatever the outcome:
1. **Convergence:** the fraction of fits converged per arm, and any errors.
2. **Positive control** (spec Section 9): how many microenvironment-associated programs are `drops_out` or `reduced_concordant`, against how many are not.
3. **Groups:** the counts of stable programs per group in the two co-primary arms.
4. **Sensitivity:** whether stable-program counts and groups shift at k/2 and 2k.
5. **Transfer:** each arm's paired difference with its SE. A negative `diff_mean` means the gene prior improved tumor-to-organoid transfer.
6. **Agreement:** kappa and the concordance correlation between 2D and 3D.
7. **Power:** the smallest true concordance detected with at least 80% power at the real number of pairs. This answers the Research Plan's "NEEDS A POWER ANALYSIS" note.

- [ ] **Step 3: Write the summary and show the figures**

Write `analysis/hcmi_persistence/results/<RUN_ID>/SUMMARY.md` with:
- the seven checks above, numbers copied from the Step 2 output;
- the paths to `figures/main_persistence_maps.png`, `figures/main_per_type_heatmaps.png` and the supplementary figures;
- anything that departed from the spec during the run (resubmissions, timeouts, the 3D-source decision).

Read `figures/main_persistence_maps.png` and `figures/main_per_type_heatmaps.png` with the Read tool, so the summary describes what the figures actually show.

- [ ] **Step 4: Hand off for review; nothing gets committed yet**

Tell the user where `SUMMARY.md` and the figures are, and give the headline numbers. Ask which artifacts to commit. After they approve, commit only those, with `git add -f` (for example `analysis/hcmi_persistence/results/<RUN_ID>/figures`, `.../coess_2d/programs.csv`, `.../coess_3d/programs.csv`, `.../compare/*.csv` and `SUMMARY.md`), using a `results:` commit message.

Then use superpowers:finishing-a-development-branch. The pull request's base branch (`alpine` or `main`) is the user's decision.

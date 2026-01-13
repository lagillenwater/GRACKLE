# GRACKLE Repository Improvement Recommendations

This document outlines recommended improvements for the GRACKLE package based on a comprehensive code review conducted January 2025.

## High Priority

### 1. Package Metadata (DESCRIPTION file)

The DESCRIPTION file contains placeholder values that must be updated:

```r
# Current placeholders needing replacement:
Title: What the Package Does (One Line, Title Case)
Authors@R: person("First", "Last", , "first.last@example.com", ...)
Description: What the package does (one paragraph).
License: `use_mit_license()`, `use_gpl3_license()` or friends to pick a license
```

**Actions:**

- Add proper title (e.g., "Graph Regularization Across Contextual KnowLedgE for NMF")
- Add author information with ORCID
- Write a one-paragraph description
- Specify license (BSD-2-Clause Plus Patent per Greene Lab standards)
- Update version from 0.0.0.9000 for release

### 2. Missing Package Dependencies

The NAMESPACE imports packages not declared in DESCRIPTION Imports. This will cause installation failures.

**Missing from DESCRIPTION:**

- tensorflow
- matrixStats
- PLIER
- tibble

**Additionally used in analysis scripts but not declared:**

- enrichR
- msigdbr
- aricode
- SNFtool
- NEMO
- SummarizedExperiment

### 3. Fix Broken Code Paths in GRACKLE.R

The `error_terms = TRUE` parameter path references undefined variables and will cause runtime errors:

```r
# Line 131 references variables that are never defined:
out <- list(
    residual = (Y - W %*% H),
    H = H,
    W = W,
    error = reconstruction_error_vec,
    H_diff = H_diff_vec,
    W_diff = W_diff_vec,
    pat_sim_error_vec = pat_sim_error_vec,
    grn_error_vec = grn_error_vec
)
```

**Actions:**

- Either implement the error tracking code (currently commented out at lines 137-153)
- Or remove the `error_terms` parameter entirely

### 4. Resolve C++ Compilation Issue

The NAMESPACE declares `useDynLib(GRACKLE)` and `R/RcppExports.R` references C++ functions, but no `src/` directory exists.

**Actions:**

- Add the missing `src/` directory with C++ source files
- Or remove `useDynLib(GRACKLE)` from NAMESPACE and `LinkingTo: Rcpp` from DESCRIPTION if C++ is no longer used

---

## Medium Priority

### 5. Branch Reconciliation

The `alpine` branch contains significant additional work not in `main`.

**Alpine branch has:**

- k=25 HTP analysis results (April 2025)
- Additional test scripts and data files
- Enrichment analysis data (enrichments_25_2025-04-09.RData, sig_conditions_25_2025-04-09.RData)
- Approximately 70,000 lines of additions vs main

**Actions:**

- Decide on branch strategy (merge alpine to main, or vice versa)
- Clean up backup files (*.R~ files) before merging
- Remove temporary/lock files (.#HTP_GNMF.R)

### 6. Add Formal Unit Testing

Currently no testthat framework exists. The `tests/` directory contains analysis scripts, not unit tests.

**Actions:**

- Add `tests/testthat/` directory structure
- Create basic unit tests:
  - `test-GRACKLE.R` - verify output dimensions, non-negativity constraints
  - `test-similarity.R` - test similarity matrix calculations
  - `test-split_data.R` - verify train/test proportions
  - `test-evaluation.R` - test evaluation metrics
- Add `tests/testthat.R` runner script

### 7. Re-enable Convergence Criterion

The documented convergence threshold (relative change < 1e-4) is commented out. Currently runs fixed 100 iterations regardless of convergence.

**Actions:**

- Re-enable convergence checking
- Or document that fixed iterations is intentional behavior
- Consider exposing `diff_threshold` behavior to users

### 8. Add Publication Figure Scripts

The dual-panel heatmap figure code (Condition Enrichments + MSigDB Hallmark) is missing from the repository.

**Actions:**

- Locate code on HPC and add to repository
- Create `figures/` or `scripts/figures/` directory for publication figure generation
- Ensure figures are reproducible from committed code and data

---

## Lower Priority

### 9. Code Cleanup

**In R/GRACKLE.R:**

- Remove commented-out code blocks (lines 91-119, 137-153, 193-201)
- Clean up inline comments that are no longer relevant

**In tests/ directory:**

- Organize scripts into logical subdirectories
- Remove or archive old simulation versions
- Clean up backup files (*.R~)

### 10. Hardcoded Random Seed

The GRACKLE function has hardcoded `set.seed(42)` for W and H initialization:

```r
set.seed(42)  # Lines 61, 63
```

**Actions:**

- Consider adding a `seed` parameter to allow user control
- Or document that seed is fixed for reproducibility

### 11. Export Missing Functions

- `euclideanDist` has a man page but is not exported in NAMESPACE
- `colScale` is defined but not exported (verify if intentional)

### 12. Documentation Improvements

**Actions:**

- Add a package-level vignette demonstrating typical workflow
- Add examples to man pages that currently lack them
- Consider adding a `docs/` website using pkgdown

### 13. Data Management

**Actions:**

- Document the data sources and processing pipelines
- Consider whether `data/breast_igraph_prob_1_cor_0_05.RData` should include provenance information
- Add data dictionary or metadata for included datasets

### 14. Analysis Script Improvements

Several scripts have hardcoded paths that will not work for other users:

```r
# Example from HTP_karyotype.R:
use_virtualenv("/mnt/grackle_env")
setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
```

**Actions:**

- Use relative paths or configurable paths
- Add setup instructions for reproducing analysis environment
- Consider using `here::here()` for path management

---

## Summary Checklist

| Priority | Item | Status |
|----------|------|--------|
| High | Update DESCRIPTION metadata | Pending |
| High | Add missing dependencies | Pending |
| High | Fix error_terms code path | Pending |
| High | Resolve src/ directory issue | Pending |
| Medium | Reconcile alpine/main branches | Pending |
| Medium | Add testthat unit tests | Pending |
| Medium | Re-enable convergence criterion | Pending |
| Medium | Add publication figure scripts | Pending |
| Lower | Code cleanup | Pending |
| Lower | Parameterize random seed | Pending |
| Lower | Export missing functions | Pending |
| Lower | Add vignette | Pending |
| Lower | Document data provenance | Pending |
| Lower | Fix hardcoded paths in scripts | Pending |

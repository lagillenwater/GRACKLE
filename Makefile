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

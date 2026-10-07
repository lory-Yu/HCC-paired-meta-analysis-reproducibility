SHELL := /usr/bin/env bash
PROJECT_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PYTHON ?= python3
RSCRIPT ?= Rscript

.PHONY: verify-analysis verify-sensitivity verify-cbc verify-identifier verify-targeted-revision checksums verify

verify-analysis:
	bash "$(PROJECT_ROOT)/workflow/scripts/verify_analysis.sh" "$(PROJECT_ROOT)"

verify-sensitivity:
	bash "$(PROJECT_ROOT)/workflow/scripts/verify_sensitivity.sh" "$(PROJECT_ROOT)"

verify-cbc:
	$(PYTHON) "$(PROJECT_ROOT)/workflow/scripts/verify_cbc_external_validation_v1.py"

verify-identifier:
	$(PYTHON) "$(PROJECT_ROOT)/workflow/scripts/verify_cbc_icgc_identifier_harmonization_v1.py"

verify-targeted-revision:
	@test -s "$(PROJECT_ROOT)/results/validation/cbc_targeted_revision_audits_v1/allain2016_overlap_summary.tsv"
	@test -s "$(PROJECT_ROOT)/results/validation/cbc_targeted_revision_audits_v1/discovery_external_sample_identifier_overlap_audit.tsv"
	@echo "targeted revision audit inputs present"

checksums:
	@sha256sum -c "$(PROJECT_ROOT)/SHA256SUMS.sha256"

verify: checksums verify-analysis verify-sensitivity verify-cbc verify-identifier verify-targeted-revision

# --- Added in v0.1.1: figures, Supplementary Table S1 and GEO metadata audit ---
# The submitted figures were built with matplotlib 3.11.2; cbc-figures refuses other versions.
# cbc-supplementary-s1 rewrites the S1 workbook (xlsx bytes carry a timestamp), and
# cbc-geo-metadata-audit needs re-downloaded GEO SOFT exports (set GEO_SOFT_DIR).
FIG_MPL_VERSION ?= 3.11.2
.PHONY: cbc-figures cbc-supplementary-s1 cbc-geo-metadata-audit

cbc-figures:
	$(PYTHON) -c 'import matplotlib, sys; v = matplotlib.__version__; sys.exit(0 if v == "$(FIG_MPL_VERSION)" else f"matplotlib {v} != $(FIG_MPL_VERSION); refusing to rebuild figures")'
	$(PYTHON) "$(PROJECT_ROOT)/workflow/scripts/build_cbc_manuscript_figures.py"

cbc-supplementary-s1:
	$(PYTHON) "$(PROJECT_ROOT)/workflow/scripts/build_cbc_supplementary_table_s1.py"

cbc-geo-metadata-audit:
	$(PYTHON) "$(PROJECT_ROOT)/workflow/scripts/audit_geo_characteristics_fields.py"

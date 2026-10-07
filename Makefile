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

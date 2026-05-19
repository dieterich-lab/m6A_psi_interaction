NOW = `date --utc +%y%m%d.%H%M%SZ`

include makefile.vars

KD = samples_KD.txt
OE = samples_OE.txt
M3I = samples_M3I.txt
IVT = samples_IVT.txt

# default 
DRY ?= 0

export DRY

.PHONY: run-xcalls-selected-% run-cooccurence-% logs

# pre-processing ---

logs-xcalls:
	mkdir -p logs/xcalls

# modkit extract "selected" calls
# make run-xcalls-selected-KD DRY=1
run-xcalls-selected-%: logs-xcalls
	@INPUT=$($*); \
	DRY_FLAG=""; \
	if [ "$$DRY" -eq 1 ]; then \
	  DRY_FLAG="--dry-run"; \
	fi; \
	N=$$(wc -l < $$INPUT); \
	sbatch \
	  --job-name=xcalls-$* \
	  --array=1-$$N \
	  run_modkit_extract_selected_calls \
	    $$INPUT \
	    $(MODKIT_BINARY) \
	    $(XCALLS_SELECTED_OPTS) \
	    $$DRY_FLAG

# analyses ---

# uses output of extract selected calls (motifs)
# co-occurence (no plots)
# make run-cooccurence-KD DRY=1
run-cooccurence-%: logs-xcalls
	@INPUT=$*; \
	DRY_FLAG=""; \
	if [ "$$DRY" -eq 1 ]; then \
	  DRY_FLAG="--dry-run"; \
	fi; \
	bash \
	  create_run_scripts_cooccurence_motif \
	    $$INPUT \
	    $(COOC_OPTS) \
	    $$DRY_FLAG
 

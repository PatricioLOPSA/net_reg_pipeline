#!/usr/bin/env bash
set -euo pipefail

PY_SCRIPT="/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/run_sisana_kfold.py"
KFOLD_DIR="/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort"
TEMPLATE="/storage/kuijjerarea/plos/sisana160/example_inputs/params.yml"
OUT_ROOT="/storage/kuijjerarea/plos/sisana_runs_multicohort"

# Optional settings
COMPUTE="cpu"                  # cpu | gpu
INCLUDE_TEST="--include-test"  # leave empty to disable
RUN_FLAG=""               
#RUN_FLAG="--run"               # leave empty to only write params.yml files
CONDA_ENV="sisana_env"         # required if RUN_FLAG is --run
TRAIN_FILENAME="Train_Set_Log2norm_Counts.tsv" # Filename of the training expression matrix
TEST_FILENAME="Test_Set_Log2norm_Counts.tsv"   # Filename of the test expression matrix

python "$PY_SCRIPT" \
  --kfold-dir "$KFOLD_DIR" \
  --template "$TEMPLATE" \
  --out-root "$OUT_ROOT" \
  --compute "$COMPUTE" \
  --train-filename "$TRAIN_FILENAME" \
  --test-filename "$TEST_FILENAME" \
  $INCLUDE_TEST \
  $RUN_FLAG \
  --conda-env "$CONDA_ENV"

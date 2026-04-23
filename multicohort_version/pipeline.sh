#!/bin/bash
#========================================================================================
# Master Pipeline Execution Script for BeatAML dataset
#
# Purpose:
#   Executes the entire data preprocessing pipeline: mutation formatting,
#   clinical/count data cleaning, K-fold generation, and preprocessing.
#
# Usage:
#   bash pipeline.sh /path/to/raw_data /path/to/output_folds_dir <true|false>
#
# Arguments:
#   1: RAW_DATA_DIR - Path to the directory containing raw input data files [required]
#   2: KFOLD_OUTDIR - Directory name or path for storing the generated K-Folds [required]
#   3: FILTER_PROT  - Logical (true/false) whether to filter for protein-coding genes [required]
#========================================================================================

set -e

if [ "$#" -ne 3 ]; then
    echo "Error: Missing required arguments."
    echo "Usage: bash pipeline.sh /path/to/raw_data /path/to/output_folds_dir <true|false>"
    exit 1
fi

RAW_DATA_DIR=$1
KFOLD_OUTDIR=$2
FILTER_PROT=$3
PROTEIN_CODING_GENES_FILE="protein_coding_genes.tsv"
SEED=1111
FOLDS=5
#Defaults for seed = 1111 and folds = 5



echo "==================================================="
echo "Starting BeatAML Data Processing Pipeline"
echo "Raw Data Directory: $RAW_DATA_DIR"
echo "K-Fold Output Directory: $KFOLD_OUTDIR"
echo "Filter Protein Coding Genes: $FILTER_PROT"
echo "==================================================="

echo ""
echo "[1/4] Formatting Mutation Data..."
Rscript mutation_data_clean.R --raw_data_dir "$RAW_DATA_DIR"

echo ""
if [ "$FILTER_PROT" = "true" ]; then
    echo "[2/4] Cleaning Clinical, Drug, and Count Data (Filtering for Protein Coding Genes)..."
    Rscript Data_cleaning.R --raw_data_dir "$RAW_DATA_DIR" --biomart_file "$PROTEIN_CODING_GENES_FILE"
else
    echo "[2/4] Cleaning Clinical, Drug, and Count Data (No Filtering)..."
    Rscript Data_cleaning.R --raw_data_dir "$RAW_DATA_DIR"
fi

echo ""
echo "[3/4] Generating $FOLDS-Fold Cross Validation Splits..."
Rscript Generate_Folds.R --input_file "cleaned_data/training/training_RawCountData.tsv" --folds "$FOLDS" --outdir "$KFOLD_OUTDIR" --seed "$SEED"

echo ""
echo "[4/4] Running Preprocessing on All Sets of Folds..."
for FOLD_DIR in "$KFOLD_OUTDIR"/Set_*; do
    if [ -d "$FOLD_DIR" ]; then
        echo "  -> Preprocessing $FOLD_DIR"
        Rscript preprocess_merge.R --fold_dir "$FOLD_DIR"
    fi
done

echo ""
echo "==================================================="
echo "Preprocesing completed successfully!"
echo "==================================================="

#========================================================================================
# Generate K-Folds for BeatAML dataset
#
# Purpose:
#   Generates K-fold cross-validation splits from the BeatAML training data,
#   stratified by clinical center ID.
#
# Usage:
#   Rscript Generate_Folds.R --input_file /path/to/training_RawCountData.tsv
#
# Options:
#   -i, --input_file  Character string. Path to the training count data [default: cleaned_data/training/training_RawCountData.tsv]
#   -s, --seed        Integer. Random seed for reproducibility [default: 1111]
#   -k, --folds       Integer. Number of folds for K-fold CV [default: 5]
#   -o, --outdir      Character string. Output directory for storing splits [default: kFold_splits]
#   -h, --help        Show help message and exit
#========================================================================================
library(data.table)
library(rsample)
library(tidyverse)
library(optparse)

option_list <- list(
  make_option(
    c("-i", "--input_file"),
    type = "character",
    default = "cleaned_data/training/training_RawCountData.tsv",
    help = "Path to the training count data [default: %default]"
  ),
  make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 1111,
    help = "Random seed for reproducibility [default: %default]"
  ),
  make_option(
    c("-k", "--folds"),
    type = "integer",
    default = 5,
    help = "Number of folds for K-fold CV [default: %default]"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "kFold_splits",
    help = "Output directory for storing splits [default: %default]"
  )
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list
)

opts <- parse_args(parser)

matrix_file <- opts$input_file
seed <- opts$seed
v_folds <- opts$folds
outdir <- opts$outdir

if (!file.exists(matrix_file)) {
  stop(paste("The file", matrix_file, "does not exist."), call. = FALSE)
}

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

#==============================
#==========Read data
#==============================

data <- fread(matrix_file, sep='\t') %>%
as.data.frame() %>%
rename(sample_id = 1) %>%
column_to_rownames('sample_id') %>%
mutate(centerID = as.factor(centerID))

#==============================
#==========split data
#==============================

#TODO: add option to specify stratification variable. For now, use centerID as default stratification variable.
#TODO: test that parameter k should only be a psitive integer less or equal to n.

#Generate  v train/test splits. Use CenterID (or any categorical var) as strat. factor

set.seed(seed)

#count_split <- initial_split(data = data, prop = split_size, strata = centerID)

vfold_split  <- vfold_cv(data, v = v_folds, repeats = 1, strata = centerID)

for (i in 1:v_folds) {

fold_dir <- file.path(outdir, paste0('Set_', i))
dir.create(fold_dir, showWarnings = FALSE)

train <- analysis(vfold_split$splits[[i]])
test <- assessment(vfold_split$splits[[i]])

train_fold_id  <- file.path(fold_dir, 'Train_Set_Counts.tsv')
test_fold_id  <- file.path(fold_dir, 'Test_Set_Counts.tsv')

#print intersections of sample IDs in train and test sets to confirm no overlap
train_samp_ids <- rownames(train)
test_samp_ids <- rownames(test)
intersection <- intersect(train_samp_ids, test_samp_ids)

ifelse(length(intersection) > 0,
    warning(paste0("Overlap detected in fold ", i, ": ", paste(intersection, collapse = ", "))), 
    print(paste0("No overlap detected in fold ", i)))

fwrite(train, file = train_fold_id, sep = '\t', row.names = T)
fwrite(test, file = test_fold_id, sep = '\t', row.names = T)

#save ids of samples in train and test as metadata for downstream use in model training and evaluation
train_ids_file <- file.path(fold_dir, 'Train_Sample_IDs.tsv')
test_ids_file <- file.path(fold_dir, 'Test_Sample_IDs.tsv')

fwrite(data.frame(sample_id = train_samp_ids), file = train_ids_file, sep = '\t', row.names = F)
fwrite(data.frame(sample_id = test_samp_ids), file = test_ids_file, sep = '\t', row.names = F)

}





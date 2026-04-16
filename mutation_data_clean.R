#========================================================================================
# Data formatting of mutation data for BeatAML dataset
#
# Purpose:
#   Formats and cleans mutation data from the BEAT AML WES dataset.
#   Outputs a binary matrix of samples (rows) and genes (columns) where
#   1 indicates a mutation and 0 indicates no mutation.
#
# Usage:
#   Rscript mutation_data_clean.R --raw_data_dir /path/to/raw_data_dir
#
# Options:
#   -r, --raw_data_dir  Character string. Path to directory containing raw input data files
#   -o, --output_file   Character string. Path to output file for the binary mutation matrix [default: mutation_binary_matrix_wv1to4.tsv]
#   -h, --help          Show help message and exit
#========================================================================================
library(dplyr)
library(data.table)
library(tidyverse)
library(optparse)

#################
#Read arguments
#################

option_list <- list(
  make_option(
    c("-r", "--raw_data_dir"),
    type = "character",
    default = NULL,
    help = "Path to directory containing raw input data files [required]"
  ),
  make_option(
    c("-o", "--output_file"),
    type = "character",
    default = "mutation_binary_matrix_wv1to4.tsv",
    help = "Path to output file for the binary mutation matrix [default: %default]"
  )
)

parser <- OptionParser(
  usage = "%prog --raw_data_dir /path/to/raw_data_dir [--output_file /path/to/output.tsv]",
  option_list = option_list
)

opts <- parse_args(parser)

if (is.null(opts$raw_data_dir)) {
  print_help(parser)
  stop("Missing required argument: --raw_data_dir", call. = FALSE)
}

if (!dir.exists(opts$raw_data_dir)) {
  stop(paste("The directory", opts$raw_data_dir, "does not exist."), call. = FALSE)
}

raw_data_dir <- opts$raw_data_dir
output_file <- opts$output_file

mutation_data_file <- file.path(raw_data_dir, 'beataml_wes_wv1to4_mutations_dbgap.txt')
mutation_data <- fread(mutation_data_file, sep = '\t') %>%
as.data.frame()

####################################
# --- Format mutation data
###################################

# Replace empty strings and '""' with NA
mutation_data <- mutation_data %>%
  mutate(across(everything(), ~ ifelse(. == '""' | . == "", NA_character_, .)))

mut_binary <- mutation_data %>%
    distinct(dbgap_sample_id, symbol) %>%
    mutate(mutated = 1L) %>%
    pivot_wider(
      names_from  = symbol,
      values_from = mutated,
      values_fill = 0L,
      names_prefix = ""
    ) %>%
    column_to_rownames("dbgap_sample_id")
  
  colnames(mut_binary) <- paste0(colnames(mut_binary), "_mutated")

#Eliminate D flag from mutation sample names
mut_samples <- rownames(mut_binary) %>% str_replace('D$', '')
rownames(mut_binary) <- mut_samples

# Save the binary mutation matrix to a file
fwrite(mut_binary, file = output_file, sep = "\t", row.names = TRUE)

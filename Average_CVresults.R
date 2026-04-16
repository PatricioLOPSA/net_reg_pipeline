#========================================================================================
# Average_CVresults.R
#
# Purpose:
#   This script takes the output from the regression benchmark pipeline, averages the results across all CV folds, and saves the averaged results to a new CSV file.
#
# Usage:
#   Rscript Average_CVresults.R --input_dir <file> --output_file <file>
#
# Options:
#   --input_dir     Character string. Path to directory containing CV results [required]
#   --output_file_avg   Character string. Path to output file [default: "Average_CV_results.csv"]
#   --output_file_se   Character string. Path to output file for standard errors [default: "StErr_CV_results.csv"]
#   -h, --help      Show help message and exit
#========================================================================================

library(data.table)
library(optparse)


option_list <- list(
  make_option(
    c("--input_dir"),
    type = "character",
    default = NULL,
    help = "Path to directory containing CV results [required]"
  ),
  make_option(
    c("--output_file"),
    type = "character",
    default = "Average_CV_results.csv",
    help = "Path to output file [default: \"Average_CV_results.csv\"]"
  )
)

parser <- OptionParser(
  usage = "%prog --input_dir <file> --output_file <file>",
  option_list = option_list
)
opts <- parse_args(parser)
if (is.null(opts$input_dir)) {
  print_help(parser)
  stop("Missing required argument: --input_dir", call. = FALSE)
}
input_dir <- opts$input_dir
output_file <- file.path(input_dir, opts$output_file)
output_file_se <- file.path(input_dir, "StErr_CV_results.csv")

# List all CSV files in the input directory
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) stop("No CSV files found in the input directory.")

# Read and combine all CSV files into a single data frame
combined_data <- rbindlist(lapply(csv_files, fread), fill = TRUE)
combined_data[is.na(combined_data)] <- 0

# Identify grouping columns (non-numeric, e.g. V1, drug) and numeric columns
group_cols <- names(combined_data)[!sapply(combined_data, is.numeric)]
numeric_cols <- names(combined_data)[sapply(combined_data, is.numeric)]

# Compute the mean for all numeric columns across the 5 sets
averaged_data <- combined_data[, lapply(.SD, mean, na.rm = TRUE), by = group_cols, .SDcols = numeric_cols]
StErr_data <- combined_data[, lapply(.SD, function(x) sd(x, na.rm = TRUE)/sqrt(length(x))), by = group_cols, .SDcols = numeric_cols]
# Save the averaged results to the target output file
fwrite(StErr_data, output_file_se, sep = ",")
fwrite(averaged_data, output_file, sep = ",")


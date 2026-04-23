#========================================================================================
# Preprocessing and Merging of TCGA and BeatAML count data
#
# Purpose:
#   Merges TCGA and BeatAML train/test sets, filters low count genes, 
#   applies ComBat-seq for batch correction, and normalizes (Log2).
#
# Usage:
#   Rscript preprocess_merge.R --fold_dir /path/to/kFold_splits/Set_i
#
# Options:
#   -d, --fold_dir  Character string. Path to the directory containing the fold train/test count data [required]
#   -h, --help      Show help message and exit
#========================================================================================

#Limit number of cores for PCAs
library(RhpcBLASctl)
blas_set_num_threads(16)
omp_set_num_threads(16)
omp_get_max_threads()

library(data.table)
library(gridExtra)
library(ggfortify)
library(janitor)
library(tidyverse)
library(SummarizedExperiment)
library(optparse)


library(sva)

# --- Parse Command Line Arguments ---
option_list <- list(
  make_option(
    c("-d", "--fold_dir"),
    type = "character",
    default = NULL,
    help = "Path to the directory containing the fold train/test count data [required]"
  )
)

parser <- OptionParser(
  usage = "%prog --fold_dir /path/to/kFold_splits/Set_i",
  option_list = option_list
)

opts <- parse_args(parser)

if (is.null(opts$fold_dir)) {
  print_help(parser)
  stop("Missing required argument: --fold_dir", call. = FALSE)
}

fold_dir <- opts$fold_dir

if (!dir.exists(fold_dir)) {
  stop(paste("The directory", fold_dir, "does not exist."), call. = FALSE)
}

# --- Setup Output Directories ---
out_dir <- file.path(fold_dir, 'preprocessing_output')
fig_dir <- file.path(out_dir, 'figures')

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Define input file paths ---
path_helper_funcs <- "Helper_functions.R"
path_tcga_counts <- "/storage/kuijjerarea/plos/multicohort/uncorrected_counts.csv"
path_metadata <- "/storage/kuijjerarea/plos/multicohort/metadata_multicohort.csv"
path_beat_train_counts <- file.path(fold_dir, "Train_Set_Counts.tsv")
path_beat_test_counts <- file.path(fold_dir, "Test_Set_Counts.tsv")
path_beat_clinical <- "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/raw_data/beataml_wv1to4_clinical.xlsx"
# -------------------------------

source(path_helper_funcs)

# 1. Load and process TCGA
tcga <- fread(path_tcga_counts, sep = ",") %>% 
  as.data.frame() %>% column_to_rownames('V1')
tcga <- tcga[,grepl("TCGA", colnames(tcga))]
tcga$rownames <- sapply(strsplit(rownames(tcga), split = "\\."), function(x) x[1])
mat <- as.matrix(tcga[, setdiff(names(tcga), "rownames"), drop = FALSE])
tcga <- rowsum(mat, tcga$rownames) %>% as.data.frame()

# 2. Load Metadata and TCGA Clinical
metadata <- fread(path_metadata, sep = ',') %>% as.data.frame()
tcga_clinical <- metadata %>% filter(study == 'TCGA_AML') %>% select(sample_id, tissue, sex)

# 3. Processing Function for BEAT-AML Splits
process_split <- function(beat_counts_path, beat_clin_path, split_name) {
  # Load BEAT-AML (Handling structure from Generate_Folds.R: samples as rows, genes + centerID as columns)
  beat_raw <- fread(beat_counts_path, sep = "\t") %>% as.data.frame() %>% column_to_rownames('V1')
  if ("centerID" %in% colnames(beat_raw)) {
    beat_raw$centerID <- NULL
  }
  beat <- beat_raw %>% t() %>% as.data.frame()
  beat <- mutate_all(beat, function(x) as.numeric(as.character(x)))
  
  # Load and clean BEAT Clinical from raw data
  beat_clinical <- readxl::read_xlsx(beat_clin_path) %>% 
    filter(dbgap_rnaseq_sample %in% colnames(beat)) %>% select(dbgap_rnaseq_sample, specimenType, consensus_sex)
  colnames(beat_clinical) <- c('sample_id','tissue','sex')
  beat_clinical$tissue <- gsub("Bone Marrow Aspirate", "bone marrow", beat_clinical$tissue)
  beat_clinical$tissue <- gsub('Peripheral Blood', 'blood', beat_clinical$tissue)
  beat_clinical$tissue <- gsub('Leukapheresis', 'leukapheresis', beat_clinical$tissue)
  
  # Merge genes
  common_genes <- intersect(rownames(tcga), rownames(beat))
  tcga_beat <- cbind.data.frame(tcga[common_genes, ], beat[common_genes, ])
  
  # Filter zeros
  zero_filtered <- filter.zeros.means(tcga_beat, prctg_threshold = 0.20, mean_reads = 10)
  
  # Merge clinical data
  all_clinical <- rbind(tcga_clinical, beat_clinical)
  zero_filtered <- zero_filtered[, all_clinical$sample_id]
  
  # Covariates and Batch Correction
  batch <- rep(c('TCGA_AML', 'BEAT_AML'), times = c(nrow(tcga_clinical), nrow(beat_clinical)))
  covar_mat <- cbind(tissue = as.factor(all_clinical$tissue), sex = as.factor(all_clinical$sex))
  
  message(paste("Running ComBat-seq for", split_name, "..."))
  corrected_counts <- ComBat_seq(as.matrix(zero_filtered), batch=batch, group=NULL, covar_mod=covar_mat)
  
  # Apply log2 normalization
  message(paste("Running log2 normalization for", split_name, "..."))
  col_data <- all_clinical %>% column_to_rownames('sample_id')
  col_data$batch <- as.factor(batch)
  
  log2_matrix <- log2(corrected_counts + 1)
  
  # Create SummarizedExperiment
  se <- SummarizedExperiment(assays = list(counts = corrected_counts, log2 = log2_matrix), colData = col_data)
  
  # PCA Visualization on log2 data
  pca_corr <- prcomp(t(log2_matrix), scale. = TRUE)
  var_exp <- (pca_corr$sdev^2) / sum(pca_corr$sdev^2) * 100
  pca_df <- data.frame(
    Sample = rownames(pca_corr$x),
    PC1 = pca_corr$x[, 1],
    PC2 = pca_corr$x[, 2],
    Cohort = ifelse(grepl("BA", rownames(pca_corr$x)), "BEAT-AML", "TCGA-AML")
  )
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cohort)) +
    geom_point(alpha = 0.7, size = 2) +
    labs(title = paste("PCA of Log2 Normalized Corrected Data -", split_name),
         x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
         y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
    theme_minimal() +
    scale_color_manual(values = c("TCGA-AML" = "blue", "BEAT-AML" = "red"))
  
  return(list(se = se, plot = p, corrected_counts = corrected_counts, log2_matrix = log2_matrix))
}

# 4. Run Pipeline for Train and Test Splits
train_results <- process_split(path_beat_train_counts, path_beat_clinical, "Train")
test_results <- process_split(path_beat_test_counts, path_beat_clinical, "Test")

# 5. Output Outputs
plots_list <- list(train_results$plot, test_results$plot)
multi_page_plot <- marrangeGrob(plots_list, nrow = 1, ncol = 1, top = NULL)

ggsave(file.path(fig_dir, "batch_correction_pcas.pdf"), 
       plot = multi_page_plot, width = 7, height = 6)


saveRDS(train_results$se, file.path(out_dir, "train_corrected_se.rds"))
saveRDS(test_results$se, file.path(out_dir, "test_corrected_se.rds"))

# Save Log2 matrices for network generation
log2_export_train <- train_results$log2_matrix %>% as.data.frame() %>% rownames_to_column('Target')
fwrite(log2_export_train, file.path(out_dir, "Train_Set_Log2norm_Counts.tsv"), sep = '\t', row.names = FALSE)

log2_export_test <- test_results$log2_matrix %>% as.data.frame() %>% rownames_to_column('Target')
fwrite(log2_export_test, file.path(out_dir, "Test_Set_Log2norm_Counts.tsv"), sep = '\t', row.names = FALSE)

message("Preprocessing and integration complete.")

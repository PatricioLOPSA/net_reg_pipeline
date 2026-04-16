#========================================================================================
# Preprocessing of count data for network generation
#
# Purpose:
#   Filters low count genes and normalizes count data (VST, Log2) using DESeq2
#   for train and test sets independently within a specified fold directory.
#
# Usage:
#   Rscript preprocessing.R --fold_dir /path/to/kFold_splits/Set_i
#
# Options:
#   -d, --fold_dir  Character string. Path to the directory containing the fold train/test count data [required]
#   -h, --help      Show help message and exit
#========================================================================================

#preprocessing and data exploration of count data
library(RhpcBLASctl)
blas_set_num_threads(16)
omp_set_num_threads(16)

library(data.table)
library(DESeq2)
library(gridExtra)
library(ggplot2)
library(viridis)
library(tidyverse)
library(optparse)
source("Helper_functions.R")

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

out_dir <- file.path(fold_dir, 'preprocessing_output')
fig_dir <- file.path(out_dir, 'figures')

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

train_count_file  <- file.path(fold_dir, 'Train_Set_Counts.tsv')
test_count_file  <- file.path(fold_dir, 'Test_Set_Counts.tsv')

output_matrix_file_tr <- file.path(out_dir, 'Train_Set_VSTnorm_Counts.tsv')
output_fig_file_tr <- file.path(fig_dir, 'Train_Set_VSTnorm_Counts.pdf')

output_matrix_file_te <- file.path(out_dir, 'Test_Set_VSTnorm_Counts.tsv')
output_fig_file_te <- file.path(fig_dir, 'Test_Set_VSTnorm_Counts.pdf')


# count_data_file <- 'test_countdata.tsv'
# output_matrix_file <- 'test_vst_normdata.tsv'
# output_fig_file <- 'test_vst_normdata.pdf'

transform_and_explore <- function(count_file, output_matrix_file, output_fig_file) {

#----------------------------
#------Read data
#----------------------------


count_data <- fread(count_file, sep = '\t') %>%
as.data.frame() %>%
rename(sample_id = 1) %>%
column_to_rownames('sample_id') %>%
mutate(centerID = as.factor(centerID))


#-------------------------------
#------Explore complete tr. data
#-------------------------------

#1. filter genes with very low counts
count_matrix <- count_data %>% select(!c(centerID)) %>%
t() %>%
filter.zeros.means(.,prctg_threshold = 0.20, mean_reads = 10)


#2. Norm. using deseq norm factors and use vst and/or log2 vst

#Create coldata object with center ids. TODO generalize for any categ. var.
coldata <- data.frame(sample=row.names(count_data), center = as.factor(count_data$centerID))

#normalize for sequencing depth and rna composition
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ center) %>% estimateSizeFactors()
norm_matrix <- counts(dds, normalized = T)

vst_transform <- vst(count_matrix)
vst_nonblind <- vst(dds,blind = F)
log2_transform <- log2(count_matrix + 1)
log2norm_transform <- log2(norm_matrix + 1)


#explore mean-variance relationship
# Compute mean and variance for each transformation
raw_df <- compute_mean_sd(count_matrix)
log2_df <- compute_mean_sd(log2_transform)
log2_norm_df <- compute_mean_sd(log2norm_transform)
vst_df <- compute_mean_sd(vst_transform)
vst_nonblind_df <- compute_mean_sd(assay(vst_nonblind))

p1 <- plot_mean_sd(raw_df, "Raw Counts",sample_size = 10000)
p2 <- plot_mean_sd(log2_df, "Log2 Counts",sample_size = 10000)
p3 <- plot_mean_sd(log2_norm_df, "Log2 Normalized Counts",sample_size = 10000)
p4 <- plot_mean_sd(vst_df, "VST Counts",sample_size = 10000)
p5 <- plot_mean_sd(vst_nonblind_df, 'VST non-blind Counts', sample_size = 10000)

# 3. plot PCs
vst_dds <- vst(dds)
p6 <- plotPCA(vst_dds, intgroup = 'center', ntop=NULL) + labs(title='VST normalized')

log2_dds <- normTransform(dds)
p7 <- plotPCA(log2_dds, intgroup = 'center', ntop=NULL) + labs(title='log2 normalized')

p8 <- plotPCA(vst_nonblind, intgroup = 'center', ntop= NULL) + labs(title = 'VST non-blind normalized')

# Arrange in 2x2 grid
 grid.arrange(p1, p2, p3, p4, p5, p6,p7, p8,ncol = 2)

#--------------------------------------
#------Export transformed data and figs
#--------------------------------------
norm_export <- vst_transform %>%
 as.data.frame() %>%
  rownames_to_column('Target')

fwrite(norm_export, output_matrix_file, sep = '\t', row.names = F)

pdf(output_fig_file)
grid.arrange(p1, p2, p3, p4, p5, p6,p7, p8,ncol = 2)
dev.off()

}

transform_and_explore(train_count_file, output_matrix_file_tr, output_fig_file_tr)

transform_and_explore(test_count_file, output_matrix_file_te, output_fig_file_te)

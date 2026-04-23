#========================================================================================
# Data Cleaning for BeatAML dataset
#
# Purpose:
#   Cleans and splits the BeatAML raw data into training (Waves 1+2 and Both) and 
#   testing (Waves 3+4) subsets.
#
# Usage:
#   Rscript Data_cleaning.R --raw_data_dir /path/to/raw_data_dir
#
# Options:
#   -r, --raw_data_dir  Character string. Path to directory containing raw input data files [required]
#   -b, --biomart_file  Character string. Path to protein coding genes mapping file (.tsv). If not provided, data will not be filtered for protein coding genes.
#   -h, --help          Show help message and exit
#========================================================================================
library(data.table)
library(tidyverse)
library(readxl)
library(optparse)
source('Helper_functions.R')

option_list <- list(
  make_option(
    c("-r", "--raw_data_dir"),
    type = "character",
    default = NULL,
    help = "Path to the directory containing the raw input data files [required]"
  ),
  make_option(
    c("-b", "--biomart_file"),
    type = "character",
    default = NULL,
    help = "Path to the file containing the list of protein coding genes. If not provided, no filtering is applied."
  )
)

parser <- OptionParser(
  usage = "%prog --raw_data_dir /path/to/raw_data",
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
biomart_prot_cod_file <- opts$biomart_file

protein_coding_genes <- NULL
if (!is.null(biomart_prot_cod_file) && file.exists(biomart_prot_cod_file)) {
  message("Biomart file provided. Filtering for protein coding genes...")
  protein_coding_genes <- fread(biomart_prot_cod_file) %>% as.data.frame()
} else {
  message("No valid biomart file provided. Continuing without filtering for protein coding genes.")
}

#read expression data, clinical data, sample map data and drug data
count_matrix_file <- file.path(raw_data_dir, 'beataml_waves1to4_counts_dbgap.txt')
drugdata_file <- file.path(raw_data_dir, 'beataml_probit_curve_fits_v4_dbgap.txt')
clinic_file <- file.path(raw_data_dir, 'beataml_wv1to4_clinical.xlsx')
clinic_mapping_file <- file.path(raw_data_dir, 'beataml_waves1to4_sample_mapping.xlsx')

count_matrix <- fread(count_matrix_file, sep= '\t') %>% as.data.frame()
drug_curves <- fread(drugdata_file) %>% as.data.frame()
clinical_data <- readxl::read_xlsx(clinic_file) %>% as.data.frame()
sample_map_data  <- readxl::read_xlsx(clinic_mapping_file)  %>% as.data.frame()

#####################
#------Count Data
#####################

count_mat_sampleid  <- colnames(count_matrix[,-c(1:4)])
clin_rnaseq_sampleid  <- clinical_data$dbgap_rnaseq_sample %>% .[!is.na(.)]
ids_not_present_in_clin  <- count_mat_sampleid[!(count_mat_sampleid %in% clin_rnaseq_sampleid)]
control_ids <- sample_map_data %>%
    filter(., .$rna_control != 'No') %>%
    pull(., dbgap_rnaseq_sample)

#Check that ids not present in clinical file are all healthy controls
intersect(ids_not_present_in_clin, control_ids)

#Filter waves 1+2 (training) and 3+4 (validation)
wv_12_clinical <- clinical_data %>% filter(cohort == 'Waves1+2') %>% filter(rnaSeq == 'y')
wv_34_clinical <- clinical_data %>% filter(cohort == 'Waves3+4') %>% filter(rnaSeq == 'y')
wv_both_clinical <- clinical_data %>% filter(cohort == 'Both') %>% filter(rnaSeq == 'y')

Complete_training_clinical <- rbind.data.frame(wv_12_clinical, wv_both_clinical)

wv_12_count <- count_matrix %>%
 select(c('stable_id'),colnames(count_matrix[,colnames(count_matrix) %in% wv_12_clinical$dbgap_rnaseq_sample])) %>%
  column_to_rownames('stable_id')

wv_34_count <- count_matrix %>%
 select(c('stable_id'),colnames(count_matrix[,colnames(count_matrix) %in% wv_34_clinical$dbgap_rnaseq_sample])) %>%
  column_to_rownames('stable_id')

wv_both_count <- count_matrix %>%
 select(c('stable_id'),colnames(count_matrix[,colnames(count_matrix) %in% wv_both_clinical$dbgap_rnaseq_sample])) %>%
  column_to_rownames('stable_id')

healthy_controls <- count_matrix %>% 
 select(c('stable_id'),colnames(count_matrix[,colnames(count_matrix) %in% control_ids])) %>%
  column_to_rownames('stable_id')

if (!is.null(protein_coding_genes)) {
  wv_12_count <- wv_12_count %>% filter(rownames(.) %in% protein_coding_genes$ensembl_gene_id)
  wv_34_count <- wv_34_count %>% filter(rownames(.) %in% protein_coding_genes$ensembl_gene_id)
  wv_both_count <- wv_both_count %>% filter(rownames(.) %in% protein_coding_genes$ensembl_gene_id)
  healthy_controls <- healthy_controls %>% filter(rownames(.) %in% protein_coding_genes$ensembl_gene_id)
}

#####################
#------Drug Data
#####################
drug_to_sample <- drug_curves %>% filter(dbgap_rnaseq_sample != '')


#Create drug - sample matrix usuing both auc and IC50
# create matrix of AUC values
AUC_sample_matrix <- drug_to_sample %>%
  select(inhibitor, dbgap_rnaseq_sample, auc) %>%
  pivot_wider(names_from = dbgap_rnaseq_sample, values_from = auc) %>%
  column_to_rownames('inhibitor')

IC50_sample_matrix <- drug_to_sample %>%
  select(inhibitor, dbgap_rnaseq_sample, ic50) %>%
  pivot_wider(names_from = dbgap_rnaseq_sample, values_from = ic50) %>%
  column_to_rownames('inhibitor')

wv_12_AUC <- AUC_sample_matrix[, colnames(AUC_sample_matrix) %in% colnames(wv_12_count)]
wv_34_AUC <- AUC_sample_matrix[, colnames(AUC_sample_matrix) %in% colnames(wv_34_count)]
wv_both_AUC <- AUC_sample_matrix[, colnames(AUC_sample_matrix) %in% colnames(wv_both_count)]



###################################################
#------Get center id for splitting
###################################################

#Set colnames from count matrix we want to keep
tmp <- data.frame(centerID = c(NA))
rownames(tmp) <- c('stable_id')

#pull center ID from clinical data

center_ids_wv12 <- wv_12_clinical %>%
 select(dbgap_rnaseq_sample,centerID) %>%
  column_to_rownames('dbgap_rnaseq_sample') %>% rbind(.,tmp) 

center_ids_wv34 <- wv_34_clinical %>%
 select(dbgap_rnaseq_sample,centerID) %>%
  column_to_rownames('dbgap_rnaseq_sample') %>% rbind(.,tmp)

center_ids_both <- wv_both_clinical %>%
 select(dbgap_rnaseq_sample,centerID) %>%
  column_to_rownames('dbgap_rnaseq_sample') %>% rbind(.,tmp)


wv12_CountCenter <- merge_rows(center_ids_wv12, t(wv_12_count))
wv34_CountCenter <- merge_rows(center_ids_wv34, t(wv_34_count))
wvboth_CountCenter <- merge_rows(center_ids_both, t(wv_both_count))

Complete_training_CountCenter <- rbind.data.frame(wv12_CountCenter,wvboth_CountCenter)

#Uncomment if merging entire response matrix with feature matrix + Center IDs is desired
# training_center_count_drug <- merge(Complete_training_CountCenter, t(AUC_sample_matrix), by = "row.names", all.x = TRUE, sort = FALSE)
# validation_center_count_drug <- merge(wv34_CountCenter, t(AUC_sample_matrix), by = "row.names", all.x = TRUE, sort = FALSE)


#####################
#------Export
#####################

train_dir <- file.path('cleaned_data', 'training')
test_dir <- file.path('cleaned_data', 'testing')

dir.create('cleaned_data', showWarnings = FALSE)
dir.create(train_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)

#Save separate raw count data for each Wave ID
fwrite(wv_12_count, file.path(train_dir, 'waves_12_rnaseq_count_data.tsv'), sep='\t', row.names = T)
fwrite(wv_both_count, file.path(train_dir, 'waves_both_count_data.tsv'), sep='\t', row.names = T)
fwrite(wv_34_count, file.path(test_dir, 'waves_34_rnaseq_count_data.tsv'), sep='\t', row.names = T)

#Save training (wv1+2,both) and validation (wv3+4) data count matrix ready for splitting by Center ID
fwrite(Complete_training_CountCenter, file.path(train_dir, 'training_RawCountData.tsv'), sep = '\t', row.names = T)
fwrite(wv34_CountCenter, file.path(test_dir, 'validation_RawCountData.tsv'), sep='\t', row.names = T)

#Save separate clinical data for each Wave ID
fwrite(wv_12_clinical, file.path(train_dir, 'waves_12_rnaseq_clinical_data.tsv'), sep='\t', row.names = T)
fwrite(wv_both_clinical, file.path(train_dir, 'waves_both_clinical_data.tsv'), sep='\t', row.names = T)
fwrite(Complete_training_clinical, file.path(train_dir, 'training_clinical_data.tsv'), sep='\t', row.names = T)
fwrite(wv_34_clinical, file.path(test_dir, 'waves_34_rnaseq_clinical_data.tsv'), sep='\t', row.names = T)

#Save Healthy controls count data
fwrite(healthy_controls, file.path('cleaned_data', 'healthy_controls_counts.tsv'), sep='\t', row.names = T)

#Save Drug response matrices (Shared/Complete matrices in common root)
fwrite(AUC_sample_matrix, file.path('cleaned_data', 'AUC_matrix_complete.tsv'), sep='\t', row.names = T)
fwrite(IC50_sample_matrix, file.path('cleaned_data', 'IC50_matrix_complete.tsv'), sep='\t', row.names = T)




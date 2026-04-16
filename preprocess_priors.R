library(data.table)
library(tidyverse)

#----------------------------
#------PARAMS
#----------------------------
#' Preprocess priors 
#'
#' Takes an input set path containing train and test normalized data along with
#' motif and PPI prior files, then prepares priors for downstream analysis.
#'
#' @param set_i_path Path to the Set_i directory with train/test normalized data.
#' @param motif_prior_file Path to the motif prior file.
#' @param ppi_prior_file Path to the PPI prior file.
#'
#' @examples
#' # Rscript preprocess_priors.R /path/to/Set_i motif_prior_file ppi_prior_file

args <- commandArgs(trailingOnly = TRUE)
set_i_path <- args[1]
motif_prior_file <- args[2]
ppi_prior_file <- args[3]

#read motif and ppi priors
#motif_prior_file <- 'motif_prior_ids_2024.tsv'
#ppi_prior_file <- 'ppi_prior_2024.tsv'

motif_prior <- fread(motif_prior_file, sep = "\t") %>% as.data.frame() %>% set_names(c("tf","gene", "edge"))
ppi_prior <- fread(ppi_prior_file, sep = "\t") %>% as.data.frame() %>% set_names(c("tf1","tf2","edge"))

#set_i_path <- 'kFold_splits/Set_1/preprocessing_output' #for testing. Comment out for general use.
setwd(set_i_path)

tr_exp_matrix_file <- 'Train_Set_VSTnorm_Counts.tsv'
te_exp_matrix_file <- 'Test_Set_VSTnorm_Counts.tsv'



#--Read expression data
tr_exp_matrix <- fread(tr_exp_matrix_file, sep = "\t") %>% as.data.frame()
te_exp_matrix <- fread(te_exp_matrix_file, sep = "\t") %>% as.data.frame()




preprocess_priors <- function(exp_matrix, motif_prior, ppi_prior, split_type) {

split <- split_type

filt_exp_matrix_file_header <- paste0('filt_w_header_',split,'_exp_matrix.tsv')
filt_exp_matrix_file <- paste0('filt_no_header_',split,'_exp_matrix.tsv')
filt_motif_prior_file <- paste0('filt_',split,'_motif_prior.tsv')
filt_ppi_prior_file <- paste0('filt_',split,'_ppi_prior.tsv')


#Ensure correspondence between motif prior and expression
exp_genes <- exp_matrix$Target
motif_genes <- motif_prior$gene %>% unique()
common_genes <- intersect(exp_genes,motif_genes)
exp_filt <- exp_matrix %>% filter(Target %in% common_genes)
motif_gene_filt <- motif_prior %>% filter(gene %in% common_genes)

#Ensure all TFs from PPI prior are in the motif prior
all_tfs <- union(ppi_prior$tf1, ppi_prior$tf2)
motif_prior_tfs <- motif_gene_filt$tf %>% unique()
extra_tfs <- setdiff(all_tfs,motif_prior_tfs)
ppi_prior_filt <- ppi_prior %>%
  filter(!(tf1 %in% extra_tfs) & !(tf2 %in% extra_tfs))

#Save everything properly
fwrite(exp_filt,filt_exp_matrix_file_header, sep='\t', col.names = T)
fwrite(exp_filt,filt_exp_matrix_file, sep='\t', col.names = F)
fwrite(motif_gene_filt, filt_motif_prior_file, sep = '\t', col.names = F)
fwrite(ppi_prior_filt, filt_ppi_prior_file,sep = '\t', col.names = F )

}

preprocess_priors(tr_exp_matrix, motif_prior, ppi_prior, 'Train')
preprocess_priors(te_exp_matrix, motif_prior, ppi_prior, 'Test')

#Reference in python https://github.com/ladislav-hovan/lioness_sparsity_scripts/blob/main/filter_expression_and_priors.py

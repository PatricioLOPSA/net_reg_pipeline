#========================================================================================
# Regression Benchmark Parallel Script
#
# Purpose:
#   Compare multiple methods of regularized regression (Lasso, elastic net, adaptive lasso, 
#   IPF-Lasso, priority-lasso, and 2-step refitting procedures) using indegree, outdegree, 
#   and expression data to predict drug response.
#
# Usage:
#   Rscript Regression_benchmark_par.R --train_indeg <file> --test_indeg <file> ...
#
# Options:
#   --train_indeg   Character string. Path to train indegree network data [required]
#   --test_indeg    Character string. Path to test indegree network data [required]
#   --train_outdeg  Character string. Path to train outdegree network data [required]
#   --test_outdeg   Character string. Path to test outdegree network data [required]
#   --train_expr    Character string. Path to train expression data [required]
#   --test_expr     Character string. Path to test expression data [required]
#   --drug_data     Character string. Path to drug response matrix [required]
#   -o, --out_file  Character string. Output results file [default: regression_benchmark_results.csv]
#   -c, --cores     Integer. Number of parallel workers/cores to use [default: 32]
#   -h, --help      Show help message and exit
#========================================================================================

#-----------------------------
#   Load Libraries and functions
#----------------------------- 
library(ipflasso)
library(prioritylasso)
library(c060)
library(tidyverse)
library(tidymodels)
library(glmnet)
library(glmnetUtils)
library(data.table)
library(janitor)
library(parallel)
library(foreach)
library(doParallel)
library(tictoc)
library(optparse)

option_list <- list(
  make_option("--train_indeg", type = "character", default = NULL, help = "Path to train indegree network data [required]"),
  make_option("--test_indeg", type = "character", default = NULL, help = "Path to test indegree network data [required]"),
  make_option("--train_outdeg", type = "character", default = NULL, help = "Path to train outdegree network data [required]"),
  make_option("--test_outdeg", type = "character", default = NULL, help = "Path to test outdegree network data [required]"),
  make_option("--train_expr", type = "character", default = NULL, help = "Path to train expression matrix [required]"),
  make_option("--test_expr", type = "character", default = NULL, help = "Path to test expression matrix [required]"),
  make_option("--drug_data", type = "character", default = NULL, help = "Path to drug response matrix [required]"),
  make_option(c("-o", "--out_file"), type = "character", default = "regression_benchmark_results.csv", help = "Output results file [default: %default]"),
  make_option(c("-c", "--cores"), type = "integer", default = 32, help = "Number of parallel workers/cores to use [default: %default]")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
opts <- parse_args(parser)

required_args <- c("train_indeg", "test_indeg", "train_outdeg", "test_outdeg", "train_expr", "test_expr", "drug_data")
for (arg in required_args) {
  if (is.null(opts[[arg]])) {
    print_help(parser)
    stop(paste("Missing required argument:", arg), call. = FALSE)
  }
}

n_workers <- opts$cores

source("Helper_functions.R")

#Train data with integrated TCGA+BeatAML (or other mapped train data)
train_indeg <- read_expr_data(opts$train_indeg)

#Test data with integrated TCGA+BeatAML (or other mapped test data)
test_indeg <- read_expr_data(opts$test_indeg)

flag <- ".in"
colnames(train_indeg) <- paste0(colnames(train_indeg),flag)
colnames(test_indeg) <- paste0(colnames(test_indeg),flag)

#Train outdegree data
train_outdeg <- read_TF_data(opts$train_outdeg)

#test outdegree data
test_outdeg <- read_TF_data(opts$test_outdeg)

#Expression data
train_expr <- read_expr_data(opts$train_expr)
test_expr <- read_expr_data(opts$test_expr)

all_sample_ids  <- union(rownames(train_expr), rownames(test_expr))
#load drug data a filter drug responses only on train and validation samples
drug.data <- read_drug_response(opts$drug_data) %>% clean_names()  %>% filter(rownames(.) %in% all_sample_ids)

#set threshold of number of samples per drug to create models
nas <- drug.data %>% is.na() %>% `!` %>% colSums() 
to_keep <- names(nas[nas > 100])
drug.data <- drug.data %>% select(all_of(to_keep))

n_train_val<- ncol(drug.data)

output_file <- opts$out_file

#-----------------------------
# Prepare output table (column schema):
#----------------------------- 
output_cols <- c("drug","n_samples_train","n_samples_test",
                 "p_nzero_indeg_lasso","p_nzero_outdeg_lasso","p_nzero_expr_lasso",
                 "p_nzero_indeg_enet","p_nzero_outdeg_enet","p_nzero_expr_enet",
                 "p_nzero_ipflasso","p_nzero_indeg_ipflasso","p_nzero_outdeg_ipflasso","p_nzero_expr_ipflasso",
                 "p_stb_indeg","p_stb_outdeg","p_stb_expr",
                 "prs_indeg_lasso","prs_outdeg_lasso","prs_expr_lasso",
                 "prs_indeg_enet","prs_outdeg_enet","prs_expr_enet",
                 "prs_indeg_adalasso","prs_outdeg_adalasso","prs_expr_adalasso",
                 "prs_ipflasso",'prs_prioritylasso',
                 "prs_sep_inout","prs_sep_outexpr","prs_sep_inexpr","prs_sep_allmods",
                 "prs_stb_inout","prs_stb_outexpr","prs_stb_inexpr","prs_stb_allmods",
                 "prs_stb_indeg","prs_stb_outdeg","prs_stb_expr",
                 "mse_indeg_lasso","mse_outdeg_lasso","mse_expr_lasso",
                 "mse_indeg_enet","mse_outdeg_enet","mse_expr_enet",
                 "mse_indeg_adalasso","mse_outdeg_adalasso","mse_expr_adalasso",
                 "mse_ipflasso","mse_prioritylasso",
                 "mse_sep_inout","mse_sep_outexpr","mse_sep_inexpr","mse_sep_allmods",
                 "mse_stb_inout","mse_stb_outexpr","mse_stb_inexpr","mse_stb_allmods",
                 "mse_stb_indeg","mse_stb_outdeg","mse_stb_expr",
                 "rmse_indeg_lasso","rmse_outdeg_lasso","rmse_expr_lasso",
                 "rmse_indeg_enet","rmse_outdeg_enet","rmse_expr_enet",
                 "rmse_indeg_adalasso","rmse_outdeg_adalasso","rmse_expr_adalasso",
                 "rmse_ipflasso","rmse_prioritylasso",
                 "rmse_sep_inout","rmse_sep_outexpr","rmse_sep_inexpr","rmse_sep_allmods",
                 "rmse_stb_inout","rmse_stb_outexpr","rmse_stb_inexpr","rmse_stb_allmods",
                 "rmse_stb_indeg","rmse_stb_outdeg","rmse_stb_expr")

#-----------------------------
# Parallel setup
#-----------------------------

#n_workers <- max(1, parallel::detectCores() - 1)

cl <- parallel::makeCluster(n_workers)
doParallel::registerDoParallel(cl)
parallel::clusterEvalQ(cl, {
  library(ipflasso); library(prioritylasso); library(c060)
  library(glmnet); library(glmnetUtils)
  library(data.table); library(janitor)
  library(tidyverse); library(tidymodels)
  source("Helper_functions.R")
  NULL
})
parallel::clusterExport(cl, varlist = c(
  "train_indeg","test_indeg",
  "train_outdeg","test_outdeg",
  "train_expr","test_expr",
  "drug.data"
), envir = environment())
parallel::clusterSetRNGStream(cl, 1234)

# choose which drugs to process (was 1:5 in your script)
idx <- 1:n_train_val
res_list <- foreach::foreach(i = idx, .combine = rbind, .inorder = TRUE) %dopar% {
  #-----------------------------
  # Per-drug work (this is your for-loop body)
  #-----------------------------
  # set drug of interest
  drug <- colnames(drug.data)[i]

  # Make seeds reproducible per drug
  set.seed(1234)

  # Prepare data for regression modeling:
  indegree.data <- prep_data(train_data = train_indeg, test_data = test_indeg, drug_data = drug.data, target_drug = drug)
  outdegree.data <- prep_data(train_data = train_outdeg, test_data = test_outdeg, drug_data = drug.data, target_drug = drug)
  expr.data     <- prep_data(train_data = train_expr,  test_data = test_expr,  drug_data = drug.data, target_drug = drug)

  # - Merged features for Multi-modal frameworks (IPF-lasso, etc.)
  in_out_deg_train <- merge_rows(indegree.data$x_train, outdegree.data$x_train)
  in_out_deg_test  <- merge_rows(indegree.data$x_test,  outdegree.data$x_test)
  out_expr_train   <- merge_rows(outdegree.data$x_train, expr.data$x_train)
  out_expr_test    <- merge_rows(outdegree.data$x_test,  expr.data$x_test)
  in_expr_train    <- merge_rows(indegree.data$x_train, expr.data$x_train)
  in_expr_test     <- merge_rows(indegree.data$x_test,  expr.data$x_test)
  all_merged_train <- merge_rows(in_out_deg_train, expr.data$x_train)
  all_merged_test  <- merge_rows(in_out_deg_test,  expr.data$x_test)

  # - Create blocks (indices for different modalities in training data)
  n.mod1       <- ncol(indegree.data$x_train)
  n.combmods   <- ncol(in_out_deg_train) # indegree + outdegree
  n.inoutexpr  <- ncol(all_merged_train) # indegree + outdegree + expression
  blocks_all   <- list(block1 = 1:n.mod1, block2 = (n.mod1+1):n.combmods, block3 = (n.combmods+1):n.inoutexpr)

  # PFs for ipf-lasso CV
  pflist <- generate_PFlist(3, 4)

  # impose hierarchy or blocks for priority lasso
  #blocks_expr_tf_in <- list(block1 = (n.combmods+1):n.inoutexpr , block2 = (n.mod1+1):n.combmods, block3 = 1:n.mod1)
  #blocks_tf_expr_in <- list(block1 =(n.mod1+1):n.combmods, block2 = (n.combmods+1):n.inoutexpr, block3 = 1:n.mod1)
  blocks_expr_in_tf <- list(block1 = (n.combmods+1):n.inoutexpr, block2 = 1:n.mod1, block3 =  (n.mod1+1):n.combmods)

  # Get fold ids for internal 10 fold cv and benchmark between models (per-drug due to possible NA patterns)
  set.seed(1234 + i)
  foldid <- sample(rep(seq(10), length = nrow(indegree.data$x_train)))

  # Define a null model and mock data for prediction 
  n_samples_train <- length(expr.data$y_train)
  n_samples_test  <- length(expr.data$y_test)
  mock_x_train <- matrix(rnorm(2 * n_samples_train), n_samples_train, 2)
  mock_x_test  <- matrix(rnorm(2 * n_samples_test),  n_samples_test, 2)
  nullmodel <- glmnet(x = mock_x_train, y = expr.data$y_train, family = "gaussian", alpha = 1, lambda = Inf)

  #-----------------------------
  # Model fitting: Lasso
  #----------------------------- 
  lasso.indeg <- cv.glmnet(x = indegree.data$x_train, y = indegree.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian")
  lasso.outdeg <- cv.glmnet(x = outdegree.data$x_train, y = outdegree.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian")
  lasso.expr <- cv.glmnet(x = expr.data$x_train, y = expr.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian")

  coefs.lasso.indeg <- get_coefs(lasso.indeg, "lambda.min")
  coefs.lasso.outdeg <- get_coefs(lasso.outdeg, "lambda.min")
  coefs.lasso.expr <- get_coefs(lasso.expr, "lambda.min")
  p_lasso_indeg <- length(coefs.lasso.indeg[-1,])
  p_lasso_outdeg <- length(coefs.lasso.outdeg[-1,])
  p_lasso_expr <- length(coefs.lasso.expr[-1,])

  #-----------------------------
  # Model fitting: E-net
  #----------------------------- 
  enet.indeg <- cv.glmnet(x = indegree.data$x_train, y = indegree.data$y_train, foldid = foldid, alpha = 0.5, standardize = TRUE, type.measure = "mse", family = "gaussian")
  enet.outdeg <- cv.glmnet(x = outdegree.data$x_train, y = outdegree.data$y_train, foldid = foldid, alpha = 0.5, standardize = TRUE, type.measure = "mse", family = "gaussian")
  enet.expr <- cv.glmnet(x = expr.data$x_train, y = expr.data$y_train, foldid = foldid, alpha = 0.5, standardize = TRUE, type.measure = "mse", family = "gaussian")

  coefs.enet.indeg <- get_coefs(enet.indeg, "lambda.min")
  coefs.enet.outdeg <- get_coefs(enet.outdeg, "lambda.min")
  coefs.enet.expr <- get_coefs(enet.expr, "lambda.min")
  p_enet_indeg <- length(coefs.enet.indeg[-1,])
  p_enet_outdeg <- length(coefs.enet.outdeg[-1,])
  p_enet_expr <- length(coefs.enet.expr[-1,])

  #-----------------------------
  # Model fitting: Ridge
  #----------------------------- 
  ridge.indeg <- cv.glmnet(x = indegree.data$x_train, y = indegree.data$y_train, foldid = foldid, alpha = 0, standardize = TRUE, type.measure = "mse", family = "gaussian")
  ridge.outdeg <- cv.glmnet(x = outdegree.data$x_train, y = outdegree.data$y_train, foldid = foldid, alpha = 0, standardize = TRUE, type.measure = "mse", family = "gaussian")
  ridge.expr <- cv.glmnet(x = expr.data$x_train, y = expr.data$y_train, foldid = foldid, alpha = 0, standardize = TRUE, type.measure = "mse", family = "gaussian")

  #-----------------------------
  # Model fitting: Adaptive Lasso (ridge -> lasso)
  #----------------------------- 
  coefs.ridge.indeg <- get_coefs(ridge.indeg, lambda = "lambda.min") %>% .[-1,,drop = FALSE]
  coefs.ridge.outdeg <- get_coefs(ridge.outdeg, lambda = "lambda.min") %>% .[-1,,drop = FALSE]
  coefs.ridge.expr <- get_coefs(ridge.expr, lambda = "lambda.min") %>% .[-1,,drop = FALSE]

  adalasso.indeg <- cv.glmnet(x = indegree.data$x_train, y = indegree.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian", penalty.factor = 1/abs(coefs.ridge.indeg))
  adalasso.outdeg <- cv.glmnet(x = outdegree.data$x_train, y = outdegree.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian", penalty.factor = 1/abs(coefs.ridge.outdeg))
  adalasso.expr <- cv.glmnet(x = expr.data$x_train, y = expr.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian", penalty.factor = 1/abs(coefs.ridge.expr))

  #-----------------------------
  # Model fitting: IPF-Lasso
  #----------------------------- 
  ipflasso.allmods <- cvr2.ipflasso(
    X = as.matrix(all_merged_train), Y = indegree.data$y_train, family = "gaussian", type.measure = "mse", standardize = TRUE,
    blocks = blocks_all, pflist = pflist, nfolds = 5, ncv = 1, plot = FALSE
  )

  ipf_ind_bestlambda <- ipflasso.allmods$ind.bestlambda
  coefs_ipf <- ipflasso.allmods$coeff[-1,ipf_ind_bestlambda]
  coefs_ipf <- coefs_ipf[coefs_ipf != 0]
  p_nonzero_ipflasso <- length(coefs_ipf)
  indeg_in_ipflasso <- intersect(colnames(indegree.data$x_train), names(coefs_ipf)) %>% length()
  outdeg_in_ipflasso <- intersect(colnames(outdegree.data$x_train), names(coefs_ipf)) %>% length()
  expr_in_ipflasso <- intersect(colnames(expr.data$x_train), names(coefs_ipf)) %>% length()

#-----------------------------
# Model fitting: Priority-Lasso
#----------------------------- 
 pl_eio <- prioritylasso(X = as.matrix(all_merged_train), Y = indegree.data$y_train, family = "gaussian", type.measure = "mse", 
                      blocks = blocks_expr_in_tf, standardize = TRUE, cvoffset = T, foldid = foldid)

  #-----------------------------
  # Model fitting: Separate models (then ridge combine)
  #----------------------------- 
  indeg_active <- indegree.data$x_train[, colnames(indegree.data$x_train) %in% rownames(coefs.lasso.indeg), drop = FALSE]
  outdeg_active <- outdegree.data$x_train[, colnames(outdegree.data$x_train) %in% rownames(coefs.lasso.outdeg), drop = FALSE]
  expr_active <- expr.data$x_train[, colnames(expr.data$x_train) %in% rownames(coefs.lasso.expr), drop = FALSE]

  in_out_active <- merge_rows(indeg_active, outdeg_active)
  out_expr_active <- merge_rows(outdeg_active, expr_active)
  in_expr_active <- merge_rows(indeg_active, expr_active)
  in_out_expr_active <- merge_rows(in_out_active, expr_active)

  in_out_active.data <- prep_data(train_data = in_out_active, test_data = in_out_deg_test, drug_data = drug.data, target_drug = drug)
  out_expr_active.data <- prep_data(train_data = out_expr_active, test_data = out_expr_test, drug_data = drug.data, target_drug = drug)
  in_expr_active.data <- prep_data(train_data = in_expr_active, test_data = in_expr_test, drug_data = drug.data, target_drug = drug)
  all_mods_active.data <- prep_data(train_data = in_out_expr_active, test_data = all_merged_test, drug_data = drug.data, target_drug = drug)

  sep_model.inout <- if (ncol(in_out_active.data$x_train) >= 2) cv.glmnet(x = in_out_active.data$x_train, y = in_out_active.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(sep_model.inout, nullmodel)) { in_out_active.data$x_train <- mock_x_train; in_out_active.data$x_test <- mock_x_test }

  sep_model.outexpr <- if (ncol(out_expr_active.data$x_train) >= 2) cv.glmnet(x = out_expr_active.data$x_train, y = out_expr_active.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(sep_model.outexpr, nullmodel)) { out_expr_active.data$x_train <- mock_x_train; out_expr_active.data$x_test <- mock_x_test }

  sep_model.inexpr <- if (ncol(in_expr_active.data$x_train) >= 2) cv.glmnet(x = in_expr_active.data$x_train, y = in_expr_active.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(sep_model.inexpr, nullmodel)) { in_expr_active.data$x_train <- mock_x_train; in_expr_active.data$x_test <- mock_x_test }

  sep_model.all <- if (ncol(all_mods_active.data$x_train) >= 2) cv.glmnet(x = all_mods_active.data$x_train, y = all_mods_active.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(sep_model.all, nullmodel)) { all_mods_active.data$x_train <- mock_x_train; all_mods_active.data$x_test <- mock_x_test }

  #-----------------------------
  # Model fitting: Separate stable models (mc.cores = 1 here)
  #----------------------------- 
  spath_indegree <- stabpath(y = indegree.data$y_train, x = indegree.data$x_train, mc.cores = 1, family = "gaussian", weakness = .8, alpha = 0.3)
  spath_outdegree <- stabpath(y = outdegree.data$y_train, x = outdegree.data$x_train, mc.cores = 1, family = "gaussian", weakness = .8, alpha = 0.3)
  spath_expr <- stabpath(y = expr.data$y_train, x = expr.data$x_train, mc.cores = 1, family = "gaussian", weakness = .8, alpha = 0.3)

  stable_predictors_indegree <- stabsel(spath_indegree, error = 0.05, type = "pcer", pi_thr = 0.6)$stable
  stable_predictors_outdegree <- stabsel(spath_outdegree, error = 0.05, type = "pcer", pi_thr = 0.6)$stable
  stable_predictors_expr <- stabsel(spath_expr, error = 0.05, type = "pcer", pi_thr = 0.6)$stable

  p_stb_indeg <- length(stable_predictors_indegree)
  p_stb_outdeg <- length(stable_predictors_outdegree)
  p_stb_expr <- length(stable_predictors_expr)

  indeg_stable.train <- indegree.data$x_train[, colnames(indegree.data$x_train) %in% names(stable_predictors_indegree), drop = FALSE]
  outdeg_stable.train <- outdegree.data$x_train[, colnames(outdegree.data$x_train) %in% names(stable_predictors_outdegree), drop = FALSE]
  expr_stable.train <- expr.data$x_train[, colnames(expr.data$x_train) %in% names(stable_predictors_expr), drop = FALSE]

  in_out_stable <- merge_rows(indeg_stable.train, outdeg_stable.train)
  out_expr_stable <- merge_rows(outdeg_stable.train, expr_stable.train)
  in_expr_stable <- merge_rows(indeg_stable.train, expr_stable.train)
  in_out_expr_stable <- merge_rows(in_out_stable, expr_stable.train)

  in_out_stable.data <- prep_data(train_data = in_out_stable, test_data = all_merged_test, drug_data = drug.data, target_drug = drug)
  out_expr_stable.data <- prep_data(train_data = out_expr_stable, test_data = out_expr_test, drug_data = drug.data, target_drug = drug)
  in_expr_stable.data <- prep_data(train_data = in_expr_stable, test_data = in_expr_test, drug_data = drug.data, target_drug = drug)
  all_mods_stable.data <- prep_data(train_data = in_out_expr_stable, test_data = all_merged_test, drug_data = drug.data, target_drug = drug)

  in_stable.data <- prep_data(train_data = indeg_stable.train, test_data = indegree.data$x_test, drug_data = drug.data, target_drug = drug)
  out_stable.data <- prep_data(train_data = outdeg_stable.train, test_data = outdegree.data$x_test, drug_data = drug.data, target_drug = drug)
  expr_stable.data <- prep_data(train_data = expr_stable.train, test_data = expr.data$x_test, drug_data = drug.data, target_drug = drug)

  sep_model_inout.stb <- if (ncol(in_out_stable.data$x_train) >= 2) cv.glmnet(x = in_out_stable.data$x_train, y = in_out_stable.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(sep_model_inout.stb, nullmodel)) { in_out_stable.data$x_train <- mock_x_train; in_out_stable.data$x_test <- mock_x_test }

  sep_model_outexpr.stb <- if (ncol(out_expr_stable.data$x_train) >= 2) cv.glmnet(x = out_expr_stable.data$x_train, y = out_expr_stable.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(sep_model_outexpr.stb, nullmodel)) { out_expr_stable.data$x_train <- mock_x_train; out_expr_stable.data$x_test <- mock_x_test }

  sep_model_inexpr.stb <- if (ncol(in_expr_stable.data$x_train) >= 2) cv.glmnet(x = in_expr_stable.data$x_train, y = in_expr_stable.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(sep_model_inexpr.stb, nullmodel)) { in_expr_stable.data$x_train <- mock_x_train; in_expr_stable.data$x_test <- mock_x_test }

  sep_model_all.stb <- if (ncol(all_mods_stable.data$x_train) >= 2) cv.glmnet(x = all_mods_stable.data$x_train, y = all_mods_stable.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(sep_model_all.stb, nullmodel)) { all_mods_stable.data$x_train <- mock_x_train; all_mods_stable.data$x_test <- mock_x_test }

  indeg_model.stb <- if (ncol(in_stable.data$x_train) >= 2) cv.glmnet(x = in_stable.data$x_train, y = in_stable.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(indeg_model.stb, nullmodel)) { in_stable.data$x_train <- mock_x_train; in_stable.data$x_test <- mock_x_test }

  out_model.stb <- if (ncol(out_stable.data$x_train) >= 2) cv.glmnet(x = out_stable.data$x_train, y = out_stable.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(out_model.stb, nullmodel)) { out_stable.data$x_train <- mock_x_train; out_stable.data$x_test <- mock_x_test }

  expr_model.stb <- if (ncol(expr_stable.data$x_train) >= 2) cv.glmnet(x = expr_stable.data$x_train, y = expr_stable.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
  if (identical(expr_model.stb, nullmodel)) { expr_stable.data$x_train <- mock_x_train; expr_stable.data$x_test <- mock_x_test }

  #-----------------------------
  # Predictions
  #-----------------------------
  pred.lasso.indeg <- predict(lasso.indeg, indegree.data$x_test, s = "lambda.min")
  pred.lasso.outdeg <- predict(lasso.outdeg, outdegree.data$x_test, s = "lambda.min")
  pred.lasso.expr <- predict(lasso.expr, expr.data$x_test, s = "lambda.min")

  pred.enet.indeg <- predict(enet.indeg, indegree.data$x_test, s = "lambda.min")
  pred.enet.outdeg <- predict(enet.outdeg, outdegree.data$x_test, s = "lambda.min")
  pred.enet.expr <- predict(enet.expr, expr.data$x_test, s = "lambda.min")

  pred.adalasso.indeg <- predict(adalasso.indeg, indegree.data$x_test, s = "lambda.min")
  pred.adalasso.outdeg <- predict(adalasso.outdeg, outdegree.data$x_test, s = "lambda.min")
  pred.adalasso.expr <- predict(adalasso.expr, expr.data$x_test, s = "lambda.min")

  pred.ipflasso <- ipflasso.predict(object = ipflasso.allmods, Xtest = as.matrix(all_merged_test))$linpredtest
  pred.prioritylasso <- predict(pl_eio, newdata = all_merged_test[,names(pl_eio$coefficients)], type = 'response', include.allintercepts = FALSE)

  pred.sep_inout <- predict(sep_model.inout, in_out_active.data$x_test, s = "lambda.min")
  pred.sep_outexpr <- predict(sep_model.outexpr, out_expr_active.data$x_test, s = "lambda.min")
  pred.sep_inexpr <- predict(sep_model.inexpr, in_expr_active.data$x_test, s = "lambda.min")
  pred.sep_all <- predict(sep_model.all, all_mods_active.data$x_test, s = "lambda.min") 

  pred.sep_model_inout.stb <- predict(sep_model_inout.stb, in_out_stable.data$x_test, s = "lambda.min")
  pred.sep_model_outexpr.stb <- predict(sep_model_outexpr.stb, out_expr_stable.data$x_test, s = "lambda.min")
  pred.sep_model_inexpr.stb <- predict(sep_model_inexpr.stb, in_expr_stable.data$x_test, s = "lambda.min")
  pred.sep_model_all.stb <- predict(sep_model_all.stb, all_mods_stable.data$x_test, s = "lambda.min")
  pred.indeg_model.stb <- predict(indeg_model.stb, in_stable.data$x_test, s = "lambda.min")
  pred.outdeg_model.stb <- predict(out_model.stb, out_stable.data$x_test, s = "lambda.min" )
  pred.expr_model.stb  <- predict(expr_model.stb, expr_stable.data$x_test, s = "lambda.min" )

  # Suppress warnings from cor() when standard deviation is zero (returns NA instead of error)
  safe_cor <- function(x, y) suppressWarnings(cor(x, y))
  calc_mse <- function(x, y) mean((x - y)^2, na.rm = TRUE)
  calc_rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))

  prs_indeg_lasso <- safe_cor(pred.lasso.indeg, indegree.data$y_test)
  prs_outdeg_lasso <- safe_cor(pred.lasso.outdeg, outdegree.data$y_test)
  prs_expr_lasso <- safe_cor(pred.lasso.expr, expr.data$y_test)

  prs_indeg_enet <- safe_cor(pred.enet.indeg, indegree.data$y_test)
  prs_outdeg_enet <- safe_cor(pred.enet.outdeg, outdegree.data$y_test)
  prs_expr_enet <- safe_cor(pred.enet.expr, expr.data$y_test)

  prs_indeg_adalasso <- safe_cor(pred.adalasso.indeg, indegree.data$y_test)
  prs_outdeg_adalasso <- safe_cor(pred.adalasso.outdeg, outdegree.data$y_test)
  prs_expr_adalasso <- safe_cor(pred.adalasso.expr, expr.data$y_test)

  prs_ipflasso <- safe_cor(pred.ipflasso, indegree.data$y_test)
  prs_prioritylasso <- safe_cor(pred.prioritylasso, indegree.data$y_test)

  prs_sep_inout <- safe_cor(pred.sep_inout, in_out_active.data$y_test)
  prs_sep_outexpr <- safe_cor(pred.sep_outexpr, out_expr_active.data$y_test)
  prs_sep_inexpr <- safe_cor(pred.sep_inexpr, in_expr_active.data$y_test)
  prs_sep_allmods <- safe_cor(pred.sep_all, all_mods_active.data$y_test)

  prs_stb_inout <- safe_cor(pred.sep_model_inout.stb, in_out_stable.data$y_test)
  prs_stb_outexpr <- safe_cor(pred.sep_model_outexpr.stb, out_expr_stable.data$y_test)
  prs_stb_inexpr <- safe_cor(pred.sep_model_inexpr.stb, in_expr_stable.data$y_test)
  prs_stb_allmods <- safe_cor(pred.sep_model_all.stb, all_mods_stable.data$y_test)
  prs_stb_indeg <- safe_cor(pred.indeg_model.stb ,  in_stable.data$y_test)
  prs_stb_outdeg <- safe_cor(pred.outdeg_model.stb ,  out_stable.data$y_test)
  prs_stb_expr <- safe_cor(pred.expr_model.stb ,  expr_stable.data$y_test)

  mse_indeg_lasso <- calc_mse(pred.lasso.indeg, indegree.data$y_test)
  mse_outdeg_lasso <- calc_mse(pred.lasso.outdeg, outdegree.data$y_test)
  mse_expr_lasso <- calc_mse(pred.lasso.expr, expr.data$y_test)
  mse_indeg_enet <- calc_mse(pred.enet.indeg, indegree.data$y_test)
  mse_outdeg_enet <- calc_mse(pred.enet.outdeg, outdegree.data$y_test)
  mse_expr_enet <- calc_mse(pred.enet.expr, expr.data$y_test)
  mse_indeg_adalasso <- calc_mse(pred.adalasso.indeg, indegree.data$y_test)
  mse_outdeg_adalasso <- calc_mse(pred.adalasso.outdeg, outdegree.data$y_test)
  mse_expr_adalasso <- calc_mse(pred.adalasso.expr, expr.data$y_test)
  mse_ipflasso <- calc_mse(pred.ipflasso, indegree.data$y_test)
  mse_prioritylasso <- calc_mse(pred.prioritylasso, indegree.data$y_test)
  mse_sep_inout <- calc_mse(pred.sep_inout, in_out_active.data$y_test)
  mse_sep_outexpr <- calc_mse(pred.sep_outexpr, out_expr_active.data$y_test)
  mse_sep_inexpr <- calc_mse(pred.sep_inexpr, in_expr_active.data$y_test)
  mse_sep_allmods <- calc_mse(pred.sep_all, all_mods_active.data$y_test)
  mse_stb_inout <- calc_mse(pred.sep_model_inout.stb, in_out_stable.data$y_test)
  mse_stb_outexpr <- calc_mse(pred.sep_model_outexpr.stb, out_expr_stable.data$y_test)
  mse_stb_inexpr <- calc_mse(pred.sep_model_inexpr.stb, in_expr_stable.data$y_test)
  mse_stb_allmods <- calc_mse(pred.sep_model_all.stb, all_mods_stable.data$y_test)
  mse_stb_indeg <- calc_mse(pred.indeg_model.stb ,  in_stable.data$y_test)
  mse_stb_outdeg <- calc_mse(pred.outdeg_model.stb ,  out_stable.data$y_test)
  mse_stb_expr <- calc_mse(pred.expr_model.stb ,  expr_stable.data$y_test)

  rmse_indeg_lasso <- calc_rmse(pred.lasso.indeg, indegree.data$y_test)
  rmse_outdeg_lasso <- calc_rmse(pred.lasso.outdeg, outdegree.data$y_test)
  rmse_expr_lasso <- calc_rmse(pred.lasso.expr, expr.data$y_test)
  rmse_indeg_enet <- calc_rmse(pred.enet.indeg, indegree.data$y_test)
  rmse_outdeg_enet <- calc_rmse(pred.enet.outdeg, outdegree.data$y_test)
  rmse_expr_enet <- calc_rmse(pred.enet.expr, expr.data$y_test)
  rmse_indeg_adalasso <- calc_rmse(pred.adalasso.indeg, indegree.data$y_test)
  rmse_outdeg_adalasso <- calc_rmse(pred.adalasso.outdeg, outdegree.data$y_test)
  rmse_expr_adalasso <- calc_rmse(pred.adalasso.expr, expr.data$y_test)
  rmse_ipflasso <- calc_rmse(pred.ipflasso, indegree.data$y_test)
  rmse_prioritylasso <- calc_rmse(pred.prioritylasso, indegree.data$y_test)
  rmse_sep_inout <- calc_rmse(pred.sep_inout, in_out_active.data$y_test)
  rmse_sep_outexpr <- calc_rmse(pred.sep_outexpr, out_expr_active.data$y_test)
  rmse_sep_inexpr <- calc_rmse(pred.sep_inexpr, in_expr_active.data$y_test)
  rmse_sep_allmods <- calc_rmse(pred.sep_all, all_mods_active.data$y_test)
  rmse_stb_inout <- calc_rmse(pred.sep_model_inout.stb, in_out_stable.data$y_test)
  rmse_stb_outexpr <- calc_rmse(pred.sep_model_outexpr.stb, out_expr_stable.data$y_test)
  rmse_stb_inexpr <- calc_rmse(pred.sep_model_inexpr.stb, in_expr_stable.data$y_test)
  rmse_stb_allmods <- calc_rmse(pred.sep_model_all.stb, all_mods_stable.data$y_test)
  rmse_stb_indeg <- calc_rmse(pred.indeg_model.stb ,  in_stable.data$y_test)
  rmse_stb_outdeg <- calc_rmse(pred.outdeg_model.stb ,  out_stable.data$y_test)
  rmse_stb_expr <- calc_rmse(pred.expr_model.stb ,  expr_stable.data$y_test)

  # Return one-row data.frame
  data.frame(
    drug = drug,
    n_samples_train = n_samples_train,
    n_samples_test = n_samples_test,
    p_nzero_indeg_lasso = p_lasso_indeg,
    p_nzero_outdeg_lasso = p_lasso_outdeg,
    p_nzero_expr_lasso = p_lasso_expr,
    p_nzero_indeg_enet = p_enet_indeg,
    p_nzero_outdeg_enet = p_enet_outdeg,   # fixed
    p_nzero_expr_enet = p_enet_expr,
    p_nzero_ipflasso = p_nonzero_ipflasso,
    p_nzero_indeg_ipflasso = indeg_in_ipflasso,
    p_nzero_outdeg_ipflasso = outdeg_in_ipflasso,
    p_nzero_expr_ipflasso = expr_in_ipflasso,
    p_stb_indeg = p_stb_indeg,
    p_stb_outdeg = p_stb_outdeg,
    p_stb_expr = p_stb_expr,
    prs_indeg_lasso = prs_indeg_lasso,
    prs_outdeg_lasso = prs_outdeg_lasso,
    prs_expr_lasso = prs_expr_lasso,
    prs_indeg_enet = prs_indeg_enet,
    prs_outdeg_enet = prs_outdeg_enet,
    prs_expr_enet = prs_expr_enet,
    prs_indeg_adalasso = prs_indeg_adalasso,
    prs_outdeg_adalasso = prs_outdeg_adalasso,
    prs_expr_adalasso = prs_expr_adalasso,
    prs_ipflasso = prs_ipflasso,
    prs_prioritylasso = prs_prioritylasso,
    prs_sep_inout = prs_sep_inout,
    prs_sep_outexpr = prs_sep_outexpr,
    prs_sep_inexpr = prs_sep_inexpr,
    prs_sep_allmods = prs_sep_allmods,
    prs_stb_inout = prs_stb_inout,
    prs_stb_outexpr = prs_stb_outexpr,
    prs_stb_inexpr = prs_stb_inexpr,
    prs_stb_allmods = prs_stb_allmods,
    prs_stb_indeg = prs_stb_indeg,
    prs_stb_outdeg = prs_stb_outdeg,
    prs_stb_expr = prs_stb_expr,
    mse_indeg_lasso = mse_indeg_lasso,
    mse_outdeg_lasso = mse_outdeg_lasso,
    mse_expr_lasso = mse_expr_lasso,
    mse_indeg_enet = mse_indeg_enet,
    mse_outdeg_enet = mse_outdeg_enet,
    mse_expr_enet = mse_expr_enet,
    mse_indeg_adalasso = mse_indeg_adalasso,
    mse_outdeg_adalasso = mse_outdeg_adalasso,
    mse_expr_adalasso = mse_expr_adalasso,
    mse_ipflasso = mse_ipflasso,
    mse_prioritylasso = mse_prioritylasso,
    mse_sep_inout = mse_sep_inout,
    mse_sep_outexpr = mse_sep_outexpr,
    mse_sep_inexpr = mse_sep_inexpr,
    mse_sep_allmods = mse_sep_allmods,
    mse_stb_inout = mse_stb_inout,
    mse_stb_outexpr = mse_stb_outexpr,
    mse_stb_inexpr = mse_stb_inexpr,
    mse_stb_allmods = mse_stb_allmods,
    mse_stb_indeg = mse_stb_indeg,
    mse_stb_outdeg = mse_stb_outdeg,
    mse_stb_expr = mse_stb_expr,
    rmse_indeg_lasso = rmse_indeg_lasso,
    rmse_outdeg_lasso = rmse_outdeg_lasso,
    rmse_expr_lasso = rmse_expr_lasso,
    rmse_indeg_enet = rmse_indeg_enet,
    rmse_outdeg_enet = rmse_outdeg_enet,
    rmse_expr_enet = rmse_expr_enet,
    rmse_indeg_adalasso = rmse_indeg_adalasso,
    rmse_outdeg_adalasso = rmse_outdeg_adalasso,
    rmse_expr_adalasso = rmse_expr_adalasso,
    rmse_ipflasso = rmse_ipflasso,
    rmse_prioritylasso = rmse_prioritylasso,
    rmse_sep_inout = rmse_sep_inout,
    rmse_sep_outexpr = rmse_sep_outexpr,
    rmse_sep_inexpr = rmse_sep_inexpr,
    rmse_sep_allmods = rmse_sep_allmods,
    rmse_stb_inout = rmse_stb_inout,
    rmse_stb_outexpr = rmse_stb_outexpr,
    rmse_stb_inexpr = rmse_stb_inexpr,
    rmse_stb_allmods = rmse_stb_allmods,
    rmse_stb_indeg = rmse_stb_indeg,
    rmse_stb_outdeg = rmse_stb_outdeg,
    rmse_stb_expr = rmse_stb_expr,
    check.names = FALSE
  )
}

# Convert to data.frame with expected columns
output_df <- as.data.frame(res_list, stringsAsFactors = FALSE)
# Optional: enforce column order
output_df <- output_df[, output_cols]

# Clean up cluster
parallel::stopCluster(cl)

# Write results
fwrite(output_df, file = output_file, sep = ",", row.names = TRUE, col.names = TRUE)







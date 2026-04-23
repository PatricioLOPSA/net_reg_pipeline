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

source("Helper_functions.R")

#-----------------------------
# Hardcoded paths from run_regression_benchmark.sh (Set 1)
#-----------------------------
train_indeg_path  <- "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort_nets/Set_1/train/lioness_indegree.csv"
test_indeg_path   <- "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort_nets/Set_1/test/lioness_indegree.csv"
train_outdeg_path <- "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort_nets/Set_1/train/lioness_outdegree.csv"
test_outdeg_path  <- "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort_nets/Set_1/test/lioness_outdegree.csv"
train_expr_path   <- "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort/Set_1/preprocessing_output/Train_Set_Log2norm_Counts.tsv"
test_expr_path    <- "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort/Set_1/preprocessing_output/Test_Set_Log2norm_Counts.tsv"
drug_data_path    <- "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/cleaned_data/AUC_matrix_complete.tsv"

#-----------------------------
# Load Data
#-----------------------------
train_indeg <- read_expr_data(train_indeg_path)
test_indeg  <- read_expr_data(test_indeg_path)

flag <- ".in"
colnames(train_indeg) <- paste0(colnames(train_indeg), flag)
colnames(test_indeg)  <- paste0(colnames(test_indeg), flag)

train_outdeg <- read_TF_data(train_outdeg_path)
test_outdeg  <- read_TF_data(test_outdeg_path)

train_expr <- read_expr_data(train_expr_path)
test_expr  <- read_expr_data(test_expr_path)

all_sample_ids  <- union(rownames(train_expr), rownames(test_expr))

# load drug data and filter drug responses only on train and validation samples
drug.data <- read_drug_response(drug_data_path) %>% clean_names() %>% filter(rownames(.) %in% all_sample_ids)

# set threshold of number of samples per drug to create models
nas <- drug.data %>% is.na() %>% `!` %>% colSums() 
to_keep <- names(nas[nas > 100])
drug.data <- drug.data %>% select(all_of(to_keep))

n_train_val <- ncol(drug.data)

#-----------------------------
# Setup for a single drug execution (simulating the first loop iteration)
#-----------------------------
i <- 26
drug <- colnames(drug.data)[i]

print(paste("============================================================"))
print(paste("Debugging drug:", drug, "(", i, "/", n_train_val, ")"))
print(paste("============================================================"))

set.seed(1234)

# Prepare data for regression modeling:
indegree.data  <- prep_data(train_data = train_indeg,  test_data = test_indeg,  drug_data = drug.data, target_drug = drug)
outdegree.data <- prep_data(train_data = train_outdeg, test_data = test_outdeg, drug_data = drug.data, target_drug = drug)
expr.data      <- prep_data(train_data = train_expr,   test_data = test_expr,   drug_data = drug.data, target_drug = drug)

# Get fold ids for internal 10 fold cv
set.seed(1234 + i)
foldid <- sample(rep(seq(10), length = nrow(indegree.data$x_train)))

#-----------------------------
# Model fitting: Lasso
#----------------------------- 
print(paste(" -> [", drug, "] Fitting Lasso..."))
lasso.indeg  <- cv.glmnet(x = indegree.data$x_train,  y = indegree.data$y_train,  foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian")
lasso.outdeg <- cv.glmnet(x = outdegree.data$x_train, y = outdegree.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian")
lasso.expr   <- cv.glmnet(x = expr.data$x_train,      y = expr.data$y_train,      foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian")

#-----------------------------
# Model fitting: E-net
#----------------------------- 
print(paste(" -> [", drug, "] Fitting Elastic-Net..."))
enet.indeg  <- cv.glmnet(x = indegree.data$x_train,  y = indegree.data$y_train,  foldid = foldid, alpha = 0.5, standardize = TRUE, type.measure = "mse", family = "gaussian")
enet.outdeg <- cv.glmnet(x = outdegree.data$x_train, y = outdegree.data$y_train, foldid = foldid, alpha = 0.5, standardize = TRUE, type.measure = "mse", family = "gaussian")
enet.expr   <- cv.glmnet(x = expr.data$x_train,      y = expr.data$y_train,      foldid = foldid, alpha = 0.5, standardize = TRUE, type.measure = "mse", family = "gaussian")

#-----------------------------
# Model fitting: Ridge
#----------------------------- 
print(paste(" -> [", drug, "] Fitting Ridge (for Adaptive Lasso)..."))
ridge.indeg  <- cv.glmnet(x = indegree.data$x_train,  y = indegree.data$y_train,  foldid = foldid, alpha = 0, standardize = TRUE, type.measure = "mse", family = "gaussian")
ridge.outdeg <- cv.glmnet(x = outdegree.data$x_train, y = outdegree.data$y_train, foldid = foldid, alpha = 0, standardize = TRUE, type.measure = "mse", family = "gaussian")
ridge.expr   <- cv.glmnet(x = expr.data$x_train,      y = expr.data$y_train,      foldid = foldid, alpha = 0, standardize = TRUE, type.measure = "mse", family = "gaussian")

#-----------------------------
# Model fitting: Adaptive Lasso (ridge -> lasso)
#----------------------------- 
print(paste(" -> [", drug, "] Fitting Adaptive Lasso..."))
coefs.ridge.indeg  <- get_coefs(ridge.indeg,  lambda = "lambda.min") %>% .[-1,,drop = FALSE]
coefs.ridge.outdeg <- get_coefs(ridge.outdeg, lambda = "lambda.min") %>% .[-1,,drop = FALSE]
coefs.ridge.expr   <- get_coefs(ridge.expr,   lambda = "lambda.min") %>% .[-1,,drop = FALSE]

# >>> THIS IS LIKELY WHERE YOU WANT TO DEBUG <<<
adalasso.indeg  <- cv.glmnet(x = indegree.data$x_train,  y = indegree.data$y_train,  foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian", penalty.factor = 1/abs(coefs.ridge.indeg))
adalasso.outdeg <- cv.glmnet(x = outdegree.data$x_train, y = outdegree.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian", penalty.factor = 1/abs(coefs.ridge.outdeg))
adalasso.expr   <- cv.glmnet(x = expr.data$x_train,      y = expr.data$y_train,      foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian", penalty.factor = 1/abs(coefs.ridge.expr))

#-----------------------------
# Setup for Multi-modal models
#-----------------------------
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
blocks_expr_in_tf <- list(block1 = (n.combmods+1):n.inoutexpr, block2 = 1:n.mod1, block3 =  (n.mod1+1):n.combmods)

# Define a null model and mock data for prediction 
n_samples_train <- length(expr.data$y_train)
n_samples_test  <- length(expr.data$y_test)
mock_x_train <- matrix(rnorm(2 * n_samples_train), n_samples_train, 2)
mock_x_test  <- matrix(rnorm(2 * n_samples_test),  n_samples_test, 2)
nullmodel <- glmnet(x = mock_x_train, y = expr.data$y_train, family = "gaussian", alpha = 1, lambda = Inf)

#-----------------------------
# Model fitting: IPF-Lasso
#----------------------------- 
print(paste(" -> [", drug, "] Fitting IPF-Lasso..."))
ipflasso.allmods <- cvr2.ipflasso(
  X = as.matrix(all_merged_train), Y = indegree.data$y_train, family = "gaussian", type.measure = "mse", standardize = TRUE,
  blocks = blocks_all, pflist = pflist, nfolds = 5, ncv = 1, plot = FALSE
)

#-----------------------------
# Model fitting: Priority-Lasso
#----------------------------- 
print(paste(" -> [", drug, "] Fitting Priority-Lasso..."))
pl_eio <- prioritylasso(X = as.matrix(all_merged_train), Y = indegree.data$y_train, family = "gaussian", type.measure = "mse", 
                    blocks = blocks_expr_in_tf, standardize = TRUE, cvoffset = T, foldid = foldid)

#-----------------------------
# Model fitting: Separate models (then ridge combine)
#----------------------------- 
coefs.lasso.indeg <- get_coefs(lasso.indeg, "lambda.min")
coefs.lasso.outdeg <- get_coefs(lasso.outdeg, "lambda.min")
coefs.lasso.expr <- get_coefs(lasso.expr, "lambda.min")

print(paste(" -> [", drug, "] Fitting Separate Models (Ridge refine)..."))
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
print(paste(" -> [", drug, "] Fitting Subsampling/Stability Models..."))
spath_indegree <- stabpath(y = indegree.data$y_train, x = indegree.data$x_train, mc.cores = 1, family = "gaussian", weakness = .8, alpha = 0.3)
spath_outdegree <- stabpath(y = outdegree.data$y_train, x = outdegree.data$x_train, mc.cores = 1, family = "gaussian", weakness = .8, alpha = 0.3)
spath_expr <- stabpath(y = expr.data$y_train, x = expr.data$x_train, mc.cores = 1, family = "gaussian", weakness = .8, alpha = 0.3)

stable_predictors_indegree <- stabsel(spath_indegree, error = 0.05, type = "pcer", pi_thr = 0.6)$stable
stable_predictors_outdegree <- stabsel(spath_outdegree, error = 0.05, type = "pcer", pi_thr = 0.6)$stable
stable_predictors_expr <- stabsel(spath_expr, error = 0.05, type = "pcer", pi_thr = 0.6)$stable

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
print(paste(" -> [", drug, "] Running Predictions & Calculating Metrics..."))
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

print(paste(" -> [", drug, "] Finished Debugging Calculations!"))

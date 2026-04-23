#script to remove any TCGA sample from data for training with BeatAML only
library(data.table)
library(tidyverse)
library(ipflasso)
library(prioritylasso)
library(c060)
library(glmnet)
library(glmnetUtils)
library(janitor)

source("/storage/kuijjerarea/plos/net_reg_pipeline/Helper_functions.R")

dir <- 'kFold_splits_ProtCod_nets'
set.tr <- paste0(dir, '/Set_1/train')
set.te <- paste0(dir, '/Set_1/test')

train_indeg <- read_expr_data(paste0(set.tr, '/lioness_indegree.csv'))
test_indeg <- read_expr_data(paste0(set.te, '/lioness_indegree.csv'))
train_outdeg <- read_TF_data(paste0(set.tr, '/lioness_outdegree.csv'))
test_outdeg <- read_TF_data(paste0(set.te, '/lioness_outdegree.csv'))

#remove TCGA samples
tcga_samples <- rownames(train_indeg)[grepl('TCGA', rownames(train_indeg))]
train_indeg <- train_indeg[!rownames(train_indeg) %in% tcga_samples, ]
test_indeg <- test_indeg[!rownames(test_indeg) %in% tcga_samples, ]
train_outdeg <- train_outdeg[!rownames(train_outdeg) %in% tcga_samples, ]
test_outdeg <- test_outdeg[!rownames(test_outdeg) %in% tcga_samples, ]

#test, read the expression data 
train_expr <- read_expr_data('/storage/kuijjerarea/plos/net_reg_pipeline/kFold_splits_ProtCod/Set_1/preprocessing_output/Train_Set_VSTnorm_Counts.tsv')
test_expr <- read_expr_data('/storage/kuijjerarea/plos/net_reg_pipeline/kFold_splits_ProtCod/Set_1/preprocessing_output/Test_Set_VSTnorm_Counts.tsv')

#remove TCGA samples from expression data
train_expr <- train_expr[!rownames(train_expr) %in% tcga_samples, ]
test_expr <- test_expr[!rownames(test_expr) %in% tcga_samples, ]


# --- DEBUG PIPELINE RUN ---

flag <- ".in"
colnames(train_indeg) <- paste0(colnames(train_indeg), flag)
colnames(test_indeg) <- paste0(colnames(test_indeg), flag)

all_sample_ids  <- union(rownames(train_expr), rownames(test_expr))

# load drug data
drug.data <- read_drug_response("/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/cleaned_data/AUC_matrix_complete.tsv") %>% 
  clean_names() %>% 
  filter(rownames(.) %in% all_sample_ids)

nas <- drug.data %>% is.na() %>% `!` %>% colSums() 
to_keep <- names(nas[nas > 100])
drug.data <- drug.data %>% select(all_of(to_keep))

# Debug run for the first drug
drug <- colnames(drug.data)[1]
set.seed(1234)

# Prepare data for regression modeling:
indegree.data <- prep_data(train_data = train_indeg, test_data = test_indeg, drug_data = drug.data, target_drug = drug)
outdegree.data <- prep_data(train_data = train_outdeg, test_data = test_outdeg, drug_data = drug.data, target_drug = drug)
expr.data     <- prep_data(train_data = train_expr,  test_data = test_expr,  drug_data = drug.data, target_drug = drug)

indegree.data$x_train[1:5, 1:5]
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

# Get fold ids for internal 10 fold cv
set.seed(1234 + 1)
foldid <- sample(rep(seq(10), length = nrow(indegree.data$x_train)))

# Define a null model and mock data for prediction (as in multicohort script)
n_samples_train <- length(expr.data$y_train)
n_samples_test  <- length(expr.data$y_test)
mock_x_train <- matrix(rnorm(2 * n_samples_train), n_samples_train, 2)
mock_x_test  <- matrix(rnorm(2 * n_samples_test),  n_samples_test, 2)
nullmodel <- glmnet(x = mock_x_train, y = expr.data$y_train, family = "gaussian", alpha = 1, lambda = Inf)

#-----------------------------
# Model fitting: Lasso
#----------------------------- 
print(paste("Fitting Lasso models for:", drug))
lasso.indeg <- cv.glmnet(x = indegree.data$x_train, y = indegree.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian")
lasso.outdeg <- cv.glmnet(x = outdegree.data$x_train, y = outdegree.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian")
lasso.expr <- cv.glmnet(x = expr.data$x_train, y = expr.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse", family = "gaussian")

#-----------------------------
# Model fitting: Priority-Lasso
#----------------------------- 
print(paste("Fitting PriorityLasso for:", drug))
pl_eio <- prioritylasso(X = as.matrix(all_merged_train), Y = indegree.data$y_train, family = "gaussian", type.measure = "mse", 
                      blocks = blocks_expr_in_tf, standardize = TRUE, cvoffset = T, foldid = foldid)

print("Debug setup loaded successfully. Features and models are ready in your environment.")









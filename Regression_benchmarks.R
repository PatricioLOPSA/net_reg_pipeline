# Script to compare multiple methods of regularized regression
# Lasso, elastic net, adaptive lasso, IPF-Lasso, Tree-Lasso, as well as 2-step refitting procedures

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

#tic()
source("Helper_functions.R")

#-----------------------------
# Load training and testing data for different modalities:
# TargetGene-Indegrees, TF-Outdegrees, Gene expression  
#----------------------------- 
train_indeg <- read_expr_data('/storage/kuijjerarea/plos/sisana/train/output/network/lioness_indegree.csv')
test_indeg <- read_expr_data('/storage/kuijjerarea/plos/sisana/testsplit/output/network/lioness_indegree.csv')
flag <- ".in"
colnames(train_indeg) <- paste0(colnames(train_indeg),flag)
colnames(test_indeg) <- paste0(colnames(test_indeg),flag)

train_outdeg <- read_TF_data('/storage/kuijjerarea/plos/sisana/train/output/network/lioness_outdegree.csv')
test_outdeg <- read_TF_data('/storage/kuijjerarea/plos/sisana/testsplit/output/network/lioness_outdegree.csv')

train_expr <- read_expr_data('/storage/kuijjerarea/plos/BeatAML/beataml2.0_data-main/filt_w_header_training_vst_normdata.tsv')
test_expr <- read_expr_data('/storage/kuijjerarea/plos/BeatAML/beataml2.0_data-main/filt_w_header_test_vst_normdata.tsv')

drug.data <- read_drug_response('AUC_matrix_complete.tsv') %>% clean_names()

#set threshold of number of samples per drug to create models
nas <- drug.data %>% is.na() %>% `!` %>% colSums() 
to_keep <- names(nas[nas > 80])
drug.data <- drug.data %>% select(all_of(to_keep))

n_train_val<- ncol(drug.data)

#-----------------------------
# Prepare output table:
#----------------------------- 
output_df <- cbind.data.frame(drug = NA,
                              n_samples_train = NA,
                              n_samples_test = NA,
                              p_nzero_indeg_lasso = NA,
                              p_nzero_outdeg_lasso = NA,
                              p_nzero_expr_lasso = NA,
                              p_nzero_indeg_enet = NA,
                              p_nzero_outdeg_enet = NA,
                              p_nzero_expr_enet = NA,
                              p_nzero_ipflasso = NA,
                              p_nzero_indeg_ipflasso = NA,
                              p_nzero_outdeg_ipflasso = NA,
                              p_nzero_expr_ipflasso = NA,
                              p_stb_indeg = NA,
                              p_stb_outdeg = NA,
                              p_stb_expr = NA,
                              prs_indeg_lasso = NA,
                              prs_outdeg_lasso = NA,
                              prs_expr_lasso = NA,
                              prs_indeg_enet = NA,
                              prs_outdeg_enet = NA,
                              prs_expr_enet = NA,
                              prs_indeg_adalasso = NA,
                              prs_outdeg_adalasso = NA,
                              prs_expr_adalasso = NA,
                              prs_ipflasso = NA,
                              prs_sep_inout = NA,
                              prs_sep_outexpr = NA,
                              prs_sep_inexpr = NA,
                              prs_sep_allmods = NA,
                              prs_stb_inout = NA,
                              prs_stb_outexpr = NA,
                              prs_stb_inexpr = NA,
                              prs_stb_allmods = NA,
                              prs_stb_indeg = NA,
                              prs_stb_outdeg = NA,
                              prs_stb_expr = NA)


for (i in 1:n_train_val) {
   



#-----------------------------
# Prepare data for regression modeling:
#----------------------------- 

#set drug of interest
drug <- colnames(drug.data)[i]

indegree.data <- prep_data(train_data = train_indeg,
 test_data = test_indeg,
 drug_data = drug.data,
 target_drug = drug)

outdegree.data <- prep_data(train_data = train_outdeg,
 test_data = test_outdeg,
 drug_data = drug.data,
 target_drug = drug)

expr.data <- prep_data(train_data = train_expr,
 test_data = test_expr,
 drug_data = drug.data,
 target_drug = drug)



# - Merged features for Multi-modal frameworks (IPF-lasso, priority lasso, TODO: SGL)
in_out_deg_train <- merge_rows(indegree.data$x_train, outdegree.data$x_train)
in_out_deg_test <- merge_rows(indegree.data$x_test, outdegree.data$x_test)

out_expr_train <- merge_rows(outdegree.data$x_train, expr.data$x_train)
out_expr_test <- merge_rows(outdegree.data$x_test, expr.data$x_test)

in_expr_train <- merge_rows(indegree.data$x_train, expr.data$x_train)
in_expr_test <- merge_rows(indegree.data$x_test, expr.data$x_test)

all_merged_train <- merge_rows(in_out_deg_train, expr.data$x_train)
all_merged_test <- merge_rows(in_out_deg_test, expr.data$x_test)

# - Create blocks (indices for different modalities in training data)
n.mod1 <- ncol(indegree.data$x_train)
n.combmods <- ncol(in_out_deg_train) #indegree and outdegree
n.inoutexpr <- ncol(all_merged_train) #indegree, outdegree and expression

blocks_all <- list(block1 = 1:n.mod1, block2 = (n.mod1+1):n.combmods, block3 = (n.combmods+1):n.inoutexpr) #indegree - outdegree - expression.

#Create list of PFs for ipf-lasso CV
pflist <- generate_PFlist(3, 4)

# impose hierarchy or blocks for priority lasso
blocks_expr_tf_in <- list(block1 = (n.combmods+1):n.inoutexpr , block2 = (n.mod1+1):n.combmods, block3 = 1:n.mod1) #expr - outdegree - indegree
blocks_tf_expr_in <- list(block1 =(n.mod1+1):n.combmods, block2 = (n.combmods+1):n.inoutexpr, block3 = 1:n.mod1) #outdegree - expr - indegree
blocks_expr_in_tf <- list(block1 = (n.combmods+1):n.inoutexpr, block2 = 1:n.mod1, block3 =  (n.mod1+1):n.combmods) #expr - indegree - outdegree

#Check number of samples match
expr.data$x_train %>% dim()
indegree.data$x_train %>% dim()
outdegree.data$x_train %>% dim()

#Check order of samples match
samples_expr <- expr.data$x_train %>% rownames()
samples_in <- indegree.data$x_train %>% rownames()
samples_out <- outdegree.data$x_train %>% rownames() 

identical(samples_expr, samples_in)
identical(samples_expr, samples_out)

# Get fold ids for internal 10 fold cv and benchmark between models
set.seed(1234)
foldid <- sample(rep(seq(10),length=dim(indegree.data$x_train)[1])) 

#Define a null model as intercept only = mean response and mock data for prediction 
n_samples_train <- length(expr.data$y_train)     #LOG
n_samples_test <- length(expr.data$y_test)       #LOG
mock_x_train <- matrix(rnorm(10),n_samples_train, 2)
mock_x_test <- matrix(rnorm(10),n_samples_test, 2)

nullmodel <- glmnet(x = mock_x_train, y = expr.data$y_train, family = "gaussian", alpha = 1, lambda = Inf)



#-----------------------------
# Model fitting: Lasso
#----------------------------- 
set.seed(1234)
 lasso.indeg <- cv.glmnet(x = indegree.data$x_train,
                          y = indegree.data$y_train,
                          foldid = foldid,
                          alpha = 1,
                          standardize = T,
                          type.measure = "mse",
                          family = "gaussian")

set.seed(1234)
 lasso.outdeg <- cv.glmnet(x = outdegree.data$x_train,
                           y = outdegree.data$y_train,
                           foldid = foldid,
                           alpha = 1,
                           standardize = T,
                           type.measure = "mse",
                           family = "gaussian")

set.seed(1234)
 lasso.expr <- cv.glmnet(x = expr.data$x_train,
                         y = expr.data$y_train,
                         foldid = foldid,
                         alpha = 1,
                         standardize = T,
                         type.measure = "mse",
                         family = "gaussian")

#Get active coefficients from individual lasso model
coefs.lasso.indeg <- get_coefs(lasso.indeg, "lambda.min")
coefs.lasso.outdeg <- get_coefs(lasso.outdeg, "lambda.min")
coefs.lasso.expr <- get_coefs(lasso.expr, "lambda.min")

p_lasso_indeg <- length(coefs.lasso.indeg[-1,])
p_lasso_outdeg <- length(coefs.lasso.outdeg[-1,])
p_lasso_expr <- length(coefs.lasso.expr[-1,])

#-----------------------------
# Model fitting: E-net
#----------------------------- 

set.seed(1234)
enet.indeg <- cv.glmnet(x = indegree.data$x_train,
                          y = indegree.data$y_train,
                          foldid = foldid,
                          alpha = 0.5,
                          standardize = T,
                          type.measure = "mse",
                          family = "gaussian")

set.seed(1234)
enet.outdeg <- cv.glmnet(x = outdegree.data$x_train,
                           y = outdegree.data$y_train,
                           foldid = foldid,
                           alpha = 0.5,
                           standardize = T,
                           type.measure = "mse",
                           family = "gaussian")

set.seed(1234)
enet.expr <- cv.glmnet(x = expr.data$x_train,
                         y = expr.data$y_train,
                         foldid = foldid,
                         alpha = 0.5,
                         standardize = T,
                         type.measure = "mse",
                         family = "gaussian")

#Get active coefficients from individual enet model
coefs.enet.indeg <- get_coefs(enet.indeg, "lambda.min")
coefs.enet.outdeg <- get_coefs(enet.outdeg, "lambda.min")
coefs.enet.expr <- get_coefs(enet.expr, "lambda.min")

p_enet_indeg <- length(coefs.enet.indeg[-1,])
p_enet_outdeg <- length(coefs.enet.outdeg[-1,])
p_enet_expr <- length(coefs.enet.expr[-1,])

#-----------------------------
# Model fitting: Ridge
#----------------------------- 

set.seed(1234)
 ridge.indeg <- cv.glmnet(x = indegree.data$x_train,
                         y = indegree.data$y_train,
                         foldid = foldid,
                         alpha = 0,
                         standardize = T,
                         type.measure = "mse",
                         family = "gaussian")

set.seed(1234)
 ridge.outdeg <- cv.glmnet(x = outdegree.data$x_train,
                         y = outdegree.data$y_train,
                         foldid = foldid,
                         alpha = 0,
                         standardize = T,
                         type.measure = "mse",
                         family = "gaussian")

set.seed(1234)
 ridge.expr <- cv.glmnet(x = expr.data$x_train,
                         y = expr.data$y_train,
                         foldid = foldid,
                         alpha = 0,
                         standardize = T,
                         type.measure = "mse",
                         family = "gaussian")

#-----------------------------
# Model fitting: Adaptive Lasso (vanilla version: ridge -> lasso)
#----------------------------- 

# get coefficiets from ridge regression
coefs.ridge.indeg <- get_coefs(ridge.indeg, lambda = "lambda.min") %>% .[-1,,drop = F]
coefs.ridge.outdeg <- get_coefs(ridge.outdeg, lambda = "lambda.min") %>% .[-1,,drop = F]
coefs.ridge.expr <- get_coefs(ridge.expr, lambda = "lambda.min") %>% .[-1,,drop = F]

adalasso.indeg <- cv.glmnet(x = indegree.data$x_train,
                          y = indegree.data$y_train,
                          foldid = foldid,
                          alpha = 1,
                          standardize = T,
                          type.measure = "mse",
                          family = "gaussian",
                          penalty.factor = 1/abs(coefs.ridge.indeg))

adalasso.outdeg <- cv.glmnet(x = outdegree.data$x_train,
                          y = outdegree.data$y_train,
                          foldid = foldid,
                          alpha = 1,
                          standardize = T,
                          type.measure = "mse",
                          family = "gaussian",
                          penalty.factor = 1/abs(coefs.ridge.outdeg))

adalasso.expr <- cv.glmnet(x = expr.data$x_train,
                          y = expr.data$y_train,
                          foldid = foldid,
                          alpha = 1,
                          standardize = T,
                          type.measure = "mse",
                          family = "gaussian",
                          penalty.factor = 1/abs(coefs.ridge.expr))

#-----------------------------
# Model fitting: IPF-Lasso
#----------------------------- 

#note, ipflasso doesnt allow control of folds. set nfolds to 5
set.seed(1234)
ipflasso.allmods <- cvr2.ipflasso(X=as.matrix(all_merged_train),Y=indegree.data$y_train,family="gaussian",type.measure="mse",standardize=T,
                                 blocks=blocks_all,
                                 pflist= pflist,
                                 nfolds = 5,
                                 ncv=1,
                                 plot = F)

ipf_ind_bestlambda <- ipflasso.allmods$ind.bestlambda
ipf_ind_bestpf <- ipflasso.allmods$ind.bestpf

coefs_ipf <- ipflasso.allmods$coeff[-1,ipf_ind_bestlambda] #no intercept
coefs_ipf <- coefs_ipf[coefs_ipf != 0]

#Check presence of active modalities in final cv ipflasso
p_nonzero_ipflasso <- length(coefs_ipf)
indeg_in_ipflasso <- indegree.data$x_train %>% colnames() %>% intersect(.,names(coefs_ipf)) %>% length()
outdeg_in_ipflasso <- outdegree.data$x_train %>% colnames() %>% intersect(.,names(coefs_ipf)) %>% length()
expr_in_ipflasso <- expr.data$x_train %>% colnames() %>% intersect(.,names(coefs_ipf)) %>% length() 

#-----------------------------
# Model fitting: Priority-Lasso
#----------------------------- 
# set.seed(1234)
# pl1 <- prioritylasso(X = as.matrix(all_merged_train), Y = indegree.data$y_train, family = "gaussian", type.measure = "mse", 
#                      blocks = blocks_expr_tf_in, standardize = T, cvoffset = T, foldid = foldid)

# set.seed(1234)
# pl2 <- prioritylasso(X = as.matrix(all_merged_train), Y = indegree.data$y_train, family = "gaussian", type.measure = "mse", 
#                      blocks = blocks_expr_in_tf, standardize = TRUE, cvoffset = T, foldid = foldid)

# set.seed(1234)
# pl3 <- prioritylasso(X = as.matrix(all_merged_train), Y = indegree.data$y_train, family = "gaussian", type.measure = "mse", 
#                      blocks = blocks_tf_expr_in, standardize = TRUE, cvoffset = T, foldid = foldid)

#-----------------------------
# Model fitting: Seperate models
#----------------------------- 


#filter training data sets uto only active coefs
indeg_active <- indegree.data$x_train[,colnames(indegree.data$x_train) %in% rownames(coefs.lasso.indeg)]
outdeg_active <- outdegree.data$x_train[,colnames(outdegree.data$x_train) %in% rownames(coefs.lasso.outdeg)]
expr_active <- expr.data$x_train[,colnames(expr.data$x_train) %in% rownames(coefs.lasso.expr)]

#Combine active coefs from separate modalities and prep data
in_out_active <- merge_rows(indeg_active, outdeg_active)
out_expr_active <- merge_rows(outdeg_active, expr_active)
in_expr_active <- merge_rows(indeg_active, expr_active)
in_out_expr_active <- merge_rows(in_out_active, expr_active)

in_out_active.data <- prep_data(train_data = in_out_active,
 test_data = in_out_deg_test,
 drug_data = drug.data,
 target_drug = drug)

out_expr_active.data <- prep_data(train_data = out_expr_active,
 test_data = out_expr_test,
 drug_data = drug.data,
 target_drug = drug)

in_expr_active.data <- prep_data(train_data = in_expr_active,
 test_data = in_expr_test,
 drug_data = drug.data,
 target_drug = drug)

all_mods_active.data <- prep_data(train_data = in_out_expr_active,
 test_data = all_merged_test,
 drug_data = drug.data,
 target_drug = drug)

#combine using ridge. Control if no predictors were selected in previous selection step .
if (ncol(in_out_active.data$x_train) >= 2) {
  set.seed(1234)
  sep_model.inout <- cv.glmnet(x = in_out_active.data$x_train, y = in_out_active.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid)
} else {
  #if no predictors, fill with mock data
   sep_model.inout <- nullmodel
   in_out_active.data$x_train <- mock_x_train
   in_out_active.data$x_test <- mock_x_test
}

if (ncol(out_expr_active.data$x_train) >= 2) {
  set.seed(1234)
  sep_model.outexpr <- cv.glmnet(x = out_expr_active.data$x_train, y = out_expr_active.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid)
} else {
  sep_model.outexpr <- nullmodel
   out_expr_active.data$x_train <- mock_x_train
   out_expr_active.data$x_test <- mock_x_test
}

if (ncol(in_expr_active.data$x_train) >= 2) {
  set.seed(1234)
  sep_model.inexpr <- cv.glmnet(x = in_expr_active.data$x_train, y = in_expr_active.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid) 
} else {
  sep_model.inexpr <- nullmodel
   in_expr_active.data$x_train <- mock_x_train
   in_expr_active.data$x_test <- mock_x_test
}

if (ncol(all_mods_active.data$x_train) >= 2) {
  set.seed(1234)
  sep_model.all <- cv.glmnet(x = all_mods_active.data$x_train, y = all_mods_active.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid) 
} else {
  sep_model.all <- nullmodel
   all_mods_active.data$x_train <- mock_x_train
   all_mods_active.data$x_test <- mock_x_test
}

#-----------------------------
# Model fitting: Seperate stable models
#----------------------------- 
#Get stability paths of an elastic net model 
set.seed(1234)
spath_indegree <-  stabpath(y = indegree.data$y_train, x = indegree.data$x_train, mc.cores = 8,
  family = "gaussian", weakness = .8, alpha=0.3)

set.seed(1234)
spath_outdegree <- stabpath(y = outdegree.data$y_train, x = outdegree.data$x_train, mc.cores = 8,
  family = "gaussian", weakness = .8, alpha=0.3)

set.seed(1234)
spath_expr <- stabpath(y = expr.data$y_train, x = expr.data$x_train, mc.cores = 8,
  family = "gaussian", weakness = .8, alpha=0.3)

#Get stable predictors at fixed probability threshold
stable_predictors_indegree <- stabsel(spath_indegree, error = 0.05, type = "pcer", pi_thr = 0.6)$stable
stable_predictors_outdegree <- stabsel(spath_outdegree, error = 0.05, type = "pcer", pi_thr = 0.6)$stable
stable_predictors_expr <- stabsel(spath_expr, error = 0.05, type = "pcer", pi_thr = 0.6)$stable

p_stb_indeg <- length(stable_predictors_indegree)
p_stb_outdeg <- length(stable_predictors_outdegree)
p_stb_expr <- length(stable_predictors_expr)


#filter training data sets to only stable coefs
indeg_stable.train <- indegree.data$x_train[,colnames(indegree.data$x_train) %in% names(stable_predictors_indegree)]
outdeg_stable.train <- outdegree.data$x_train[,colnames(outdegree.data$x_train) %in% names(stable_predictors_outdegree)]
expr_stable.train <- expr.data$x_train[,colnames(expr.data$x_train) %in% names(stable_predictors_expr)]

#Combine and prep data
in_out_stable <- merge_rows(indeg_stable.train, outdeg_stable.train)
out_expr_stable <- merge_rows(outdeg_stable.train, expr_stable.train)
in_expr_stable <- merge_rows(indeg_stable.train, expr_stable.train)
in_out_expr_stable <- merge_rows(in_out_stable, expr_stable.train)

in_out_stable.data <- prep_data(train_data = in_out_stable,
 test_data = all_merged_test,
 drug_data = drug.data,
 target_drug = drug)

out_expr_stable.data <- prep_data(train_data = out_expr_stable,
 test_data = out_expr_test,
 drug_data = drug.data,
 target_drug = drug)

in_expr_stable.data <- prep_data(train_data = in_expr_stable,
                                 test_data = in_expr_test,
                                 drug_data = drug.data,
                                 target_drug = drug)

all_mods_stable.data <- prep_data(
                                  train_data = in_out_expr_stable,
                                  test_data = all_merged_test,
                                  drug_data = drug.data,
                                  target_drug = drug)

#unimodal stable predictors only model
 
in_stable.data <- prep_data(train_data = indeg_stable.train,
 test_data = indegree.data$x_test,
 drug_data = drug.data,
 target_drug = drug)

 out_stable.data <- prep_data(train_data = outdeg_stable.train,
 test_data = outdegree.data$x_test,
 drug_data = drug.data,
 target_drug = drug)

 expr_stable.data <- prep_data(train_data = expr_stable.train,
 test_data = expr.data$x_test,
 drug_data = drug.data,
 target_drug = drug)


if (ncol(in_out_stable.data$x_train) >= 2) {   
set.seed(1234)
sep_model_inout.stb <- cv.glmnet(x = in_out_stable.data$x_train, y = in_out_stable.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid)
} else {
   sep_model.inout.stb <- nullmodel
   in_out_stable.data$x_train <- mock_x_train
   in_out_stable.data$x_test <- mock_x_test
}

if (ncol(out_expr_stable.data$x_train) >= 2) {
set.seed(1234)
sep_model_outexpr.stb <- cv.glmnet(x = out_expr_stable.data$x_train, y = out_expr_stable.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid)   
} else {
   sep_model_outexpr.stb <- nullmodel
   out_expr_stable.data$x_train <- mock_x_train
   out_expr_stable.data$x_test <- mock_x_test
}

if (ncol(in_expr_stable.data$x_train) >= 2) {
set.seed(1234)
sep_model_inexpr.stb <- cv.glmnet(x = in_expr_stable.data$x_train, y = in_expr_stable.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid)   
} else {
   sep_model_inexpr.stb <- nullmodel
   in_expr_stable.data$x_train <- mock_x_train
   in_expr_stable.data$x_test <- mock_x_test
}

if (ncol(all_mods_stable.data$x_train) >= 2) {
set.seed(1234)
sep_model_all.stb <- cv.glmnet(x = all_mods_stable.data$x_train, y = all_mods_stable.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid)  
} else {
   sep_model_all.stb <- nullmodel
   all_mods_stable.data$x_train <- mock_x_train
   all_mods_stable.data$x_test <- mock_x_test
}

if (ncol(in_stable.data$x_train) >= 2) {
set.seed(1234)
indeg_model.stb <- cv.glmnet(x = in_stable.data$x_train, y = in_stable.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid)   
} else {
   indeg_model.stb <- nullmodel
   in_stable.data$x_train <- mock_x_train
   in_stable.data$x_test <- mock_x_test
}

if (ncol(out_stable.data$x_train) >= 2) {
set.seed(1234)
out_model.stb <- cv.glmnet(x = out_stable.data$x_train, y = out_stable.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid)  
} else {
   out_model.stb <- nullmodel
   out_stable.data$x_train <- mock_x_train
   out_stable.data$x_test <- mock_x_test
}

if (ncol(expr_stable.data$x_train) >= 2) {
set.seed(1234)
expr_model.stb <- cv.glmnet(x = expr_stable.data$x_train, y = expr_stable.data$y_train, alpha=0, family = "gaussian", type.measure = "mse", foldid = foldid)
} else {
   expr_model.stb <- nullmodel
   expr_stable.data$x_train <- mock_x_train
   expr_stable.data$x_test <- mock_x_test
}


#-----------------------------
# Predictions
#-----------------------------

# lasso
pred.lasso.indeg <- predict(lasso.indeg, indegree.data$x_test, s = "lambda.min")
pred.lasso.outdeg <- predict(lasso.outdeg, outdegree.data$x_test, s = "lambda.min")
pred.lasso.expr <- predict(lasso.expr, expr.data$x_test, s = "lambda.min")

#elastic net
pred.enet.indeg <- predict(enet.indeg, indegree.data$x_test, s = "lambda.min")
pred.enet.outdeg <- predict(enet.outdeg, outdegree.data$x_test, s = "lambda.min")
pred.enet.expr <- predict(enet.expr, expr.data$x_test, s = "lambda.min")

#adalasso
pred.adalasso.indeg <- predict(adalasso.indeg, indegree.data$x_test, s="lambda.min")
pred.adalasso.outdeg <- predict(adalasso.outdeg, outdegree.data$x_test, s = "lambda.min")
pred.adalasso.expr <- predict(adalasso.expr, expr.data$x_test, s = "lambda.min")

#Ipf-lasso
pred.ipflasso <- ipflasso.predict(object=ipflasso.allmods,Xtest=as.matrix(all_merged_test))$linpredtest

#Priority-lasso
# pred.pl1 <- predict(pl1, newdata = all_merged_test[,names(pl1$coefficients)], type = 'response', include.allintercepts = FALSE)
# pred.pl2 <- predict(pl2, newdata = all_merged_test[,names(pl2$coefficients)], type = 'response', include.allintercepts = FALSE)
# pred.pl3 <- predict(pl3, newdata = all_merged_test[,names(pl3$coefficients)], type = 'response', include.allintercepts = FALSE)


#Separate models
pred.sep_inout <- predict(sep_model.inout, in_out_active.data$x_test, s = "lambda.min")
pred.sep_outexpr <- predict(sep_model.outexpr, out_expr_active.data$x_test, s = "lambda.min")
pred.sep_inexpr <- predict(sep_model.inexpr, in_expr_active.data$x_test, s = "lambda.min")
pred.sep_all <- predict(sep_model.all, all_mods_active.data$x_test, s = "lambda.min") 

#Separate stable models
pred.sep_model_inout.stb <-predict(sep_model_inout.stb, in_out_stable.data$x_test, s = "lambda.min")
pred.sep_model_outexpr.stb <- predict(sep_model_outexpr.stb, out_expr_stable.data$x_test, s = "lambda.min")
pred.sep_model_inexpr.stb <- predict(sep_model_inexpr.stb, in_expr_stable.data$x_test, s = "lambda.min")
pred.sep_model_all.stb <- predict(sep_model_all.stb, all_mods_stable.data$x_test, s = "lambda.min")
pred.indeg_model.stb <- predict(indeg_model.stb, in_stable.data$x_test, s = "lambda.min")
pred.outdeg_model.stb <- predict(out_model.stb, out_stable.data$x_test, s = "lambda.min" )
pred.expr_model.stb  <- predict(expr_model.stb, expr_stable.data$x_test, s = "lambda.min" )



prs_indeg_lasso <- cor(pred.lasso.indeg, indegree.data$y_test)
prs_outdeg_lasso <- cor(pred.lasso.outdeg, outdegree.data$y_test)
prs_expr_lasso <- cor(pred.lasso.expr, expr.data$y_test)

prs_indeg_enet <- cor(pred.enet.indeg, indegree.data$y_test)
prs_outdeg_enet <- cor(pred.enet.outdeg, outdegree.data$y_test)
prs_expr_enet <- cor(pred.enet.expr, expr.data$y_test)

prs_indeg_adalasso <- cor(pred.adalasso.indeg, indegree.data$y_test)
prs_outdeg_adalasso <- cor(pred.adalasso.outdeg, outdegree.data$y_test)
prs_expr_adalasso <- cor(pred.adalasso.expr, expr.data$y_test)

prs_ipflasso <- cor(pred.ipflasso, indegree.data$y_test)

prs_sep_inout <- cor(pred.sep_inout, in_out_active.data$y_test)
prs_sep_outexpr <- cor(pred.sep_outexpr, out_expr_active.data$y_test)
prs_sep_inexpr <- cor(pred.sep_inexpr, in_expr_active.data$y_test)
prs_sep_allmods <- cor(pred.sep_all, all_mods_active.data$y_test)

prs_stb_inout <- cor(pred.sep_model_inout.stb, in_out_stable.data$y_test)
prs_stb_outexpr <- cor(pred.sep_model_outexpr.stb, out_expr_stable.data$y_test)
prs_stb_inexpr <- cor(pred.sep_model_inexpr.stb, in_expr_stable.data$y_test)
prs_stb_allmods <- cor(pred.sep_model_all.stb, all_mods_stable.data$y_test)
prs_stb_indeg <- cor(pred.indeg_model.stb ,  in_stable.data$y_test)
prs_stb_outdeg <- cor(pred.outdeg_model.stb ,  out_stable.data$y_test)
prs_stb_expr <- cor(pred.expr_model.stb ,  expr_stable.data$y_test)

#-----------------------------
# ouput
#-----------------------------
output_df[i,1] <- drug
output_df[i,2] <- n_samples_train
output_df[i,3] <- n_samples_test
output_df[i,4] <- p_lasso_indeg
output_df[i,5] <- p_lasso_outdeg
output_df[i,6] <- p_lasso_expr
output_df[i,7] <- p_enet_indeg
output_df[i,8] <- p_enet_indeg
output_df[i,9] <- p_enet_expr
output_df[i,10] <- p_nonzero_ipflasso
output_df[i,11] <- indeg_in_ipflasso
output_df[i,12] <- outdeg_in_ipflasso
output_df[i,13] <- expr_in_ipflasso
output_df[i,14] <- p_stb_indeg
output_df[i,15] <- p_stb_outdeg
output_df[i,16] <- p_stb_expr
output_df[i,17] <- prs_indeg_lasso
output_df[i,18] <- prs_outdeg_lasso
output_df[i,19] <- prs_expr_lasso
output_df[i,20] <- prs_indeg_enet
output_df[i,21] <- prs_outdeg_enet
output_df[i,22] <- prs_expr_enet
output_df[i,23] <- prs_indeg_adalasso
output_df[i,24] <- prs_outdeg_adalasso
output_df[i,25] <- prs_expr_adalasso 
output_df[i,26] <- prs_ipflasso
output_df[i,27] <- prs_sep_inout
output_df[i,28] <- prs_sep_outexpr
output_df[i,29] <- prs_sep_inexpr
output_df[i,30] <- prs_sep_allmods
output_df[i,31] <- prs_stb_inout
output_df[i,32] <- prs_stb_outexpr
output_df[i,33] <- prs_stb_inexpr
output_df[i,34] <- prs_stb_allmods
output_df[i,35] <- prs_stb_indeg
output_df[i,36] <- prs_stb_outdeg
output_df[i,37] <- prs_stb_outdeg

gc()
}

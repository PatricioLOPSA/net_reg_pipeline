#this script will test if adding random gaussian noise gives similar results overall
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
library(ggplot2)
library(ggrepel)

setwd('net_reg_pipeline_multicohort')
source('Helper_functions.R')

# Setup output directory for sets
out_results_dir <- "random_noise_results"
if (!dir.exists(out_results_dir)) dir.create(out_results_dir, recursive = TRUE)

drug_data_path <- '/storage/kuijjerarea/plos/net_reg_pipeline/cleaned_data/AUC_matrix_complete.tsv'
drug_data_full <- read_drug_response(drug_data_path) %>% clean_names()

#-----------------------------
# Output Setup
#-----------------------------
output_cols <- c(
  "drug", "n_samples_train", "n_samples_test",
  "prs_expr_lasso", "prs_noise_lasso", "prs_noisyexpr_lasso", "prs_concat_lasso", "prs_concat2_lasso", "prs_ipflasso", "prs_ipflasso2", "prs_sep_mods"
)

#-----------------------------
# Parallel setup
#-----------------------------
n_workers <- 24 # Adjust as needed
cl <- parallel::makeCluster(n_workers)
doParallel::registerDoParallel(cl)
parallel::clusterEvalQ(cl, {
  library(ipflasso); library(c060)
  library(glmnet); library(glmnetUtils)
  library(data.table); library(janitor)
  library(tidyverse)
  source("Helper_functions.R")
  NULL
})
parallel::clusterSetRNGStream(cl, 1234)

#-----------------------------
# Set Loop
#-----------------------------
for (set_idx in 1:5) {
  print(paste("Processing Set:", set_idx))
  
  #paths
  expr_set_tr_path <- sprintf('/storage/kuijjerarea/plos/net_reg_pipeline/kFold_splits_ProtCod/Set_%d/preprocessing_output/Train_Set_VSTnorm_Counts.tsv', set_idx)
  expr_set_te_path <- sprintf('/storage/kuijjerarea/plos/net_reg_pipeline/kFold_splits_ProtCod/Set_%d/preprocessing_output/Test_Set_VSTnorm_Counts.tsv', set_idx)
  
  #Expression data
  train_expr <- read_expr_data(expr_set_tr_path)
  test_expr <- read_expr_data(expr_set_te_path)
  
  all_sample_ids  <- union(rownames(train_expr), rownames(test_expr))
  #load drug data a filter drug responses only on train and validation samples
  drug.data <- drug_data_full %>% filter(rownames(.) %in% all_sample_ids)
  
  #set threshold of number of samples per drug to create models
  nas <- drug.data %>% is.na() %>% `!` %>% colSums() 
  to_keep <- names(nas[nas > 100])
  drug.data <- drug.data %>% select(all_of(to_keep))
  
  n_train_val<- ncol(drug.data)
  
  #Create expression data with random noise
  set.seed(123 + set_idx) # change seed per set for varying noise
  noise_sd <- 0.1 * sd(as.matrix(train_expr)) # Set noise level as 10% of the standard deviation of the data
  noise_matrix_tr <- matrix(rnorm(nrow(train_expr) * ncol(train_expr), mean = 0, sd = noise_sd), nrow = nrow(train_expr), ncol = ncol(train_expr))
  noise_matrix_te <- matrix(rnorm(nrow(test_expr) * ncol(test_expr), mean = 0, sd = noise_sd), nrow = nrow(test_expr), ncol = ncol(test_expr))
  
  # Assign dimension names so prep_data works correctly
  rownames(noise_matrix_tr) <- rownames(train_expr)
  colnames(noise_matrix_tr) <- paste0("noise_", 1:ncol(train_expr))
  rownames(noise_matrix_te) <- rownames(test_expr)
  colnames(noise_matrix_te) <- paste0("noise_", 1:ncol(test_expr))
  
  train_expr_noisy <- train_expr + noise_matrix_tr
  test_expr_noisy <- test_expr + noise_matrix_te
  colnames(train_expr_noisy) <- paste0("noisy_", colnames(train_expr))
  colnames(test_expr_noisy) <- paste0("noisy_", colnames(test_expr))

  parallel::clusterExport(cl, varlist = c(
    "train_expr", "test_expr", "noise_matrix_tr", "noise_matrix_te", 
    "train_expr_noisy", "test_expr_noisy", "drug.data"
  ), envir = environment())

  #-----------------------------
  # Modeling Loop
  #-----------------------------
  idx <- 1:n_train_val
  res_list <- foreach::foreach(i = idx, .combine = rbind, .inorder = TRUE) %dopar% {
    drug <- colnames(drug.data)[i]
    set.seed(1234)
    
    # Prepare data 
    expr.data       <- prep_data(train_data = train_expr, test_data = test_expr, drug_data = drug.data, target_drug = drug)
    noise.data      <- prep_data(train_data = noise_matrix_tr, test_data = noise_matrix_te, drug_data = drug.data, target_drug = drug)
    noisy_expr.data <- prep_data(train_data = train_expr_noisy, test_data = test_expr_noisy, drug_data = drug.data, target_drug = drug)

    n_samples_train <- length(expr.data$y_train)
    n_samples_test  <- length(expr.data$y_test)
    
    # Merged explicit modalities
    concat_train <- merge_rows(expr.data$x_train, noise.data$x_train)
    concat_test  <- merge_rows(expr.data$x_test, noise.data$x_test)
    concat.data  <- prep_data(train_data = concat_train, test_data = concat_test, drug_data = drug.data, target_drug = drug)
    
    # Merged extra modalities (expression + noisy expression)
    concat2_train <- merge_rows(expr.data$x_train, noisy_expr.data$x_train)
    concat2_test  <- merge_rows(expr.data$x_test, noisy_expr.data$x_test)
    concat2.data  <- prep_data(train_data = concat2_train, test_data = concat2_test, drug_data = drug.data, target_drug = drug)
    
    # Blocks for IPF-lasso
    n.mod1 <- ncol(expr.data$x_train)
    n.combmods <- ncol(concat_train)
    blocks_all <- list(block1 = 1:n.mod1, block2 = (n.mod1+1):n.combmods)
    pflist <- generate_PFlist(2, 4)

    # Blocks for IPF-lasso 2 (expr + noisy_expr)
    n.mod1_2 <- ncol(expr.data$x_train)
    n.combmods_2 <- ncol(concat2_train)
    blocks_all2 <- list(block1 = 1:n.mod1_2, block2 = (n.mod1_2+1):n.combmods_2)

    set.seed(1234 + i)
    foldid <- sample(rep(seq(10), length = n_samples_train))

    # Mock null models
    mock_x_train <- matrix(rnorm(2 * n_samples_train), n_samples_train, 2)
    mock_x_test  <- matrix(rnorm(2 * n_samples_test),  n_samples_test, 2)
    nullmodel <- glmnet(x = mock_x_train, y = expr.data$y_train, family = "gaussian", alpha = 1, lambda = Inf)

    # Standard Lasso
    lasso.expr       <- cv.glmnet(x = expr.data$x_train, y = expr.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse")
    lasso.noise      <- cv.glmnet(x = noise.data$x_train, y = noise.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse")
    lasso.noisyexpr  <- cv.glmnet(x = noisy_expr.data$x_train, y = noisy_expr.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse")
    lasso.concat     <- cv.glmnet(x = concat.data$x_train, y = concat.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse")
    lasso.concat2    <- cv.glmnet(x = concat2.data$x_train, y = concat2.data$y_train, foldid = foldid, alpha = 1, standardize = TRUE, type.measure = "mse")

    coefs.expr       <- get_coefs(lasso.expr, "lambda.min")
    coefs.noise      <- get_coefs(lasso.noise, "lambda.min")
    coefs.noisyexpr  <- get_coefs(lasso.noisyexpr, "lambda.min")
    coefs.concat     <- get_coefs(lasso.concat, "lambda.min")
    coefs.concat2    <- get_coefs(lasso.concat2, "lambda.min")

    # IPF-Lasso
    ipflasso.mod <- cvr2.ipflasso(
      X = as.matrix(concat_train), Y = expr.data$y_train, family = "gaussian", type.measure = "mse", standardize = TRUE,
      blocks = blocks_all, pflist = pflist, nfolds = 5, ncv = 1, plot = FALSE
    )
    ipf_ind_bestlambda <- ipflasso.mod$ind.bestlambda
    coefs_ipf <- ipflasso.mod$coeff[-1,ipf_ind_bestlambda]
    coefs_ipf <- coefs_ipf[coefs_ipf != 0]

    # IPF-Lasso 2 (expr + noisy_expr)
    ipflasso.mod2 <- cvr2.ipflasso(
      X = as.matrix(concat2_train), Y = expr.data$y_train, family = "gaussian", type.measure = "mse", standardize = TRUE,
      blocks = blocks_all2, pflist = pflist, nfolds = 5, ncv = 1, plot = FALSE
    )
    ipf_ind_bestlambda2 <- ipflasso.mod2$ind.bestlambda
    coefs_ipf2 <- ipflasso.mod2$coeff[-1,ipf_ind_bestlambda2]
    coefs_ipf2 <- coefs_ipf2[coefs_ipf2 != 0]

    # Separate / Ridge framework
    expr_active <- expr.data$x_train[, colnames(expr.data$x_train) %in% rownames(coefs.expr), drop = FALSE]
    noise_active <- noise.data$x_train[, colnames(noise.data$x_train) %in% rownames(coefs.noise), drop = FALSE]
    sep_active_train <- merge_rows(expr_active, noise_active)
    sep_active.data <- prep_data(train_data = sep_active_train, test_data = concat.data$x_test, drug_data = drug.data, target_drug = drug)

    sep_model <- if (ncol(sep_active.data$x_train) >= 2) cv.glmnet(x = sep_active.data$x_train, y = sep_active.data$y_train, alpha = 0, family = "gaussian", type.measure = "mse", foldid = foldid) else nullmodel
    if (identical(sep_model, nullmodel)) { sep_active.data$x_train <- mock_x_train; sep_active.data$x_test <- mock_x_test }

    # Predictions
    pred.lasso.expr       <- predict(lasso.expr, expr.data$x_test, s = "lambda.min")
    pred.lasso.noise      <- predict(lasso.noise, noise.data$x_test, s = "lambda.min")
    pred.lasso.noisyexpr  <- predict(lasso.noisyexpr, noisy_expr.data$x_test, s = "lambda.min")
    pred.lasso.concat     <- predict(lasso.concat, concat.data$x_test, s = "lambda.min")
    pred.lasso.concat2    <- predict(lasso.concat2, concat2.data$x_test, s = "lambda.min")
    pred.ipflasso         <- ipflasso.predict(object = ipflasso.mod, Xtest = as.matrix(concat.data$x_test))$linpredtest
    pred.ipflasso2        <- ipflasso.predict(object = ipflasso.mod2, Xtest = as.matrix(concat2.data$x_test))$linpredtest
    pred.sep_mods         <- predict(sep_model, sep_active.data$x_test, s = "lambda.min")

    # Metrics
    safe_cor <- function(x, y) suppressWarnings(cor(x, y))

    y_te <- expr.data$y_test
    
    data.frame(
      drug = drug,
      n_samples_train = n_samples_train,
      n_samples_test = n_samples_test,
      prs_expr_lasso = safe_cor(pred.lasso.expr, y_te),
      prs_noise_lasso = safe_cor(pred.lasso.noise, y_te),
      prs_noisyexpr_lasso = safe_cor(pred.lasso.noisyexpr, y_te),
      prs_concat_lasso = safe_cor(pred.lasso.concat, y_te),
      prs_concat2_lasso = safe_cor(pred.lasso.concat2, y_te),
      prs_ipflasso = safe_cor(pred.ipflasso, y_te),
      prs_ipflasso2 = safe_cor(pred.ipflasso2, y_te),
      prs_sep_mods = safe_cor(pred.sep_mods, sep_active.data$y_test),
      check.names = FALSE
    )
  }

  output_df <- as.data.frame(res_list, stringsAsFactors = FALSE)
  output_df <- output_df[, output_cols]
  
  out_file <- file.path(out_results_dir, paste0("Set_", set_idx, "_random_noise_results.csv"))
  fwrite(output_df, file = out_file, sep = ",", row.names = TRUE, col.names = TRUE)
}

parallel::stopCluster(cl)

#-----------------------------
# Averaging CV Results
#-----------------------------
csv_files <- list.files(out_results_dir, pattern = "^Set_[0-9]+_.*\\.csv$", full.names = TRUE)
combined_data <- rbindlist(lapply(csv_files, fread), fill = TRUE)
combined_data[is.na(combined_data)] <- 0

group_cols <- names(combined_data)[!sapply(combined_data, is.numeric)]
numeric_cols <- names(combined_data)[sapply(combined_data, is.numeric)]

averaged_data <- combined_data[, lapply(.SD, mean, na.rm = TRUE), by = group_cols, .SDcols = numeric_cols]
fwrite(averaged_data, file.path(out_results_dir, "Average_CV_results.csv"), sep = ",")

#-----------------------------
# Visualization
#-----------------------------

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

get_metric_label <- function(metric) {
  switch(metric, prs = "Pearson corr", mse = "Mean Squared Error", rmse = "Root Mean Squared Error")
}

is_better_target <- function(target_val, baseline_val, metric) {
  if (metric == "prs") target_val > baseline_val else target_val < baseline_val
}

should_label_drug <- function(target_val, baseline_val, metric) {
  if (metric == "prs") {
    baseline_val <= 0 & target_val > baseline_val
  } else {
    target_val < baseline_val
  }
}

plot_model_comparison <- function(data, target_model, baseline_model, metric, out_dir) {
  target_col <- paste0(metric, "_", target_model)
  baseline_col <- paste0(metric, "_", baseline_model)
  
  if (!(target_col %in% colnames(data)) || !(baseline_col %in% colnames(data))) {
    warning(paste("Columns not found in data. Skipping plot for:", target_model, "vs", baseline_model))
    return(NULL)
  }
  
  plot_data <- data.frame(
    drug = data$drug,
    target = as.numeric(data[[target_col]]),
    baseline = as.numeric(data[[baseline_col]])
  )
  
  plot_data <- plot_data %>%
    mutate(
      position = ifelse(is_better_target(target, baseline, metric), paste("Better", target_model), paste("Better", baseline_model)),
      drug_label = ifelse(should_label_drug(target, baseline, metric), drug, NA)
    )
  
  p_comp <- ggplot(plot_data, aes(x = target, y = baseline, color = position)) +
    geom_point(alpha = 0.6) +
    geom_text_repel(aes(label = drug_label), size = 2.5, min.segment.length = 0, max.overlaps = Inf, seed = 42, color = "black", show.legend = FALSE) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = setNames(c("orange", "darkgrey"), c(paste("Better", target_model), paste("Better", baseline_model)))) +
    labs(title = paste(toupper(target_model), "vs", toupper(baseline_model), "(", toupper(metric), ")"),
         x = paste(target_model, "(", get_metric_label(metric), ")"),
         y = paste(baseline_model, "(", get_metric_label(metric), ")"), color = "Performance") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  
  if (metric == "prs") {
    p_comp <- p_comp + xlim(-0.5, 1) + ylim(-0.5, 1)
  }
  p_comp
  ggsave(file.path(out_dir, paste0("comp_", target_model, "_vs_", baseline_model, "_", metric, ".pdf")), plot = p_comp, width = 8, height = 8)
}

# Example Plot Generation using Averaged Data
plot_model_comparison(averaged_data %>% as.data.frame() %>% clean_names(), target_model = "ipflasso", baseline_model = "expr_lasso", metric = "prs", out_dir = "figures")
plot_model_comparison(averaged_data %>% as.data.frame() %>% clean_names(), target_model = "ipflasso2", baseline_model = "expr_lasso", metric = "prs", out_dir = "figures")
plot_model_comparison(averaged_data %>% as.data.frame() %>% clean_names(), target_model = "noise_lasso", baseline_model = "expr_lasso", metric = "prs", out_dir = "figures")




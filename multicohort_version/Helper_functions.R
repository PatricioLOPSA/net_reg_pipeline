
#Function to merge by row names without adding an extra rowname column
merge_rows <- function(x, y) {
  merged_df <- merge(x, y, by = "row.names")
  rownames(merged_df) <- merged_df[, 1]
  merged_df[, 1] <- NULL
  return(merged_df)
}

# Keep top n features (columns) with highest variance.
# df: data.frame (rows = samples, cols = features)
# n: number of top variable features to keep
filter_top_variable_features <- function(df, n = 500) {
  if (!is.data.frame(df) && !is.matrix(df)) stop("df must be a data.frame or matrix")
  if (n <= 0) return(df[, 0, drop = FALSE])
  # compute column variances (handle NA)
  col_vars <- apply(df, 2, var, na.rm = TRUE)
  # if n >= number of features, return original
  if (n >= length(col_vars)) return(as.data.frame(df))
  top_names <- names(sort(col_vars, decreasing = TRUE))[seq_len(n)]
  return(as.data.frame(df)[, top_names, drop = FALSE])
}

#Function to filter by zero prevalence
filter.zeros.means <- function(m1,prctg_threshold,mean_reads) {
  print("Input (dim): ")
  print(dim(m1))
  #if gene is present in more than x percent of samples, keep
  threshold <- round(dim(m1)[2]*prctg_threshold) 
  print(paste0("Threshold: ",threshold))
  m1 <- m1[rowSums(m1 != 0) >= threshold, ] 

#if gene has on average, more or equal 10 reads, keep
  print(paste0("Rows After Zeros (dim): ",dim(m1)[1]))
  m1 <- m1[rowMeans(m1) >= mean_reads, ]   
  print(paste0("Rows After means (dim): ",dim(m1)[1]))
  return(m1)
}

#Get sd-mean relationship
compute_mean_sd <- function(matrix_data) {
  df <- data.frame(mean = rowMeans(matrix_data), sd = apply(matrix_data, 1, sd))
  df$ranked_mean <- rank(df$mean)  # Rank mean in ascending order
  return(df)
}

# Function to create a mean-SD plot with density contours
plot_mean_sd <- function(df, title, sample_size = 5000) {
  # Sample a subset of the data to reduce computational load
  df_sample <- df[sample(nrow(df), min(sample_size, nrow(df))), ]
  
  ggplot(df_sample, aes(x = ranked_mean, y = sd)) +
    geom_point(alpha = 0.06, size = 1.5, color= "blue") +  # Lightweight scatter plot
    geom_density_2d(color = "lightgrey", alpha = 0.8) +  # Contour density
    geom_smooth(method = NULL, color = "red", se = FALSE, linewidth = 1) +  # Linear trend line
    scale_y_log10() +
    labs(title = title, x = "Ranked Mean Expression", y = "Standard Deviation") +
    theme_minimal()
}


#Function to load and format data with gene expression or indegree shape
#TODO generalize so that there isnt a requirement for column name
read_expr_data <- function(file_path) {
  fread(file_path) %>%
    as.data.frame() %>%
    rename(Target = 1) %>%
    column_to_rownames("Target") %>%
    t() %>%
    as.data.frame()
}

#Function to read and format TF outdegree data
read_TF_data <- function(file_path) {
  fread(file_path) %>%
    as.data.frame() %>%
    rename(TF = 1) %>%
    column_to_rownames("TF") %>%
    t() %>%
    as.data.frame()
}

#Function to read drug response matrices. 
read_drug_response <- function(file_path) {
 print(paste("DEBUG: Attempting to read drug data from path:", file_path)) # Add this line
  print(paste("DEBUG: Does R think this file exists?", file.exists(file_path))) # Add this line

  fread(file_path, sep = "\t") %>%
    as.data.frame() %>%
    rename(compound = 1) %>%
    column_to_rownames("compound") %>%
    t() %>%
    as.data.frame()
}

#Function to extract coefficients from a glmnet object. If alpha > 0, only active predictos will be selected
get_coefs <- function(model, lambda){
coefs_active <- coef.glmnet(model, s=lambda)
coefs_sign_active <- coefs_active[coefs_active[,1] != 0,,drop=F] #get active features != 0
return(coefs_sign_active)
}

#Function to prepare data for regression for a specific drug
prep_data <- function(train_data, test_data, drug_data, target_drug){

# From drug data, filter samples containing drug of interest and filter x_train, x_test, y_train and y_test for selected drug 
    drug_name <- target_drug
    single_drug <- drug_data %>% select(all_of(drug_name)) %>% na.omit() %>% rename(drug_response = 1)

    filt_train <- merge_rows(single_drug, train_data)
    filt_test <- merge_rows(single_drug, test_data) %>%
    select(any_of(colnames(filt_train)))

    missing_genes <- setdiff(colnames(filt_train), colnames(filt_test))

 # if p_i in x_train is missing in x_test, impute p_i in x_test using mean from x_train.
    if (length(missing_genes) > 0) {
    missing_means <- filt_train %>%
    select(all_of(missing_genes)) %>%
    map_dfc(~ mean(.x, na.rm = TRUE))

    n_test <- nrow(filt_test)
    missing_matrix <- map_dfc(missing_means, ~ rep(.x, n_test))
    colnames(missing_matrix) <- names(missing_means)

    filt_test <- filt_test %>%
    bind_cols(missing_matrix)

    x_test <- filt_test %>%
    select(-all_of("drug_response")) %>%
    as.matrix()
    } else {
    x_test <- filt_test %>%
    select(-all_of("drug_response")) %>%
    as.matrix()
    }

    # Define training sets
    x_train <- filt_train %>%
     select(-all_of("drug_response")) %>%
      as.matrix()

    y_train <- filt_train %>%
     pull(all_of("drug_response"))

    y_test <- filt_test %>%
     pull(all_of("drug_response"))

    x_test <- x_test[,colnames(x_train), drop = FALSE]

    #Remove features with zero variance in training and test data
    nzv_cols <- apply(x_train, 2, var) > 0
    x_train <- x_train[, nzv_cols, drop = FALSE]
    x_test <- x_test[, nzv_cols, drop = FALSE]
    #output named list
    list(
        x_train = x_train,
        y_train = y_train,
        x_test = x_test,
        y_test = y_test
    )

}


# generates list of lists of all possible combinations of penalty factors (PFs) for IPFLasso,
# The first modality is always taken as the reference modality that is penalized with a PF of 1
# n_modalities correspond to number of different omics or modalitites
# pen_range corresponds to the range penalty values that grow exponentially (powers of 2) from 1 up to 2^pen_range.
generate_PFlist <- function(n_modalities,pen_range){

penalty_range <- 2^(seq(0, pen_range, by = 1))

# Generate all combinations
pf_combinations <- expand.grid(rep(list(penalty_range), n_modalities))
raw_pflist <- split(pf_combinations, seq(nrow(pf_combinations)))
raw_pflist <- lapply(raw_pflist, function(row) as.numeric(unlist(row)))

# Normalize to first modality
normalize_to_first <- function(v) {
  v / v[1]
}

# Apply normalization
normalized_pflist <- lapply(raw_pflist, normalize_to_first)

# Keep unique normalized vectors
unique_pflist <- unique(normalized_pflist)

return(unique_pflist)
}


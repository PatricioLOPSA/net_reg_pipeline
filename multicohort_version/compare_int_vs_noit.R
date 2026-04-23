#This script will compare the average performance across folds for
# integrated vs non-integrated cases

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# Set exploratory parameters
model_name <- "indeg_lasso" # Model to compare (e.g., ipflasso, prioritylasso, indeg_lasso)
metric <- "prs"          # Performance metric: prs, mse, or rmse
out_dir <- "figures"     # Output directory for plots
show_labels <- FALSE     # Option to show/hide drug labels in the plot

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# paths to read dataframe of average results
noint_path <- '/storage/kuijjerarea/plos/net_reg_pipeline/results_5cv/Average_CV_results.csv'
int_path <- '/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/5cv_results/Average_CV_results.csv'

# read in dataframes
noint_df <- fread(noint_path, sep = ',') %>% as.data.frame()
int_df <- fread(int_path, sep = ',') %>% as.data.frame()

# Clean names to match conventions
noint_df <- janitor::clean_names(noint_df)
int_df <- janitor::clean_names(int_df)

common_cols <- intersect(colnames(noint_df), colnames(int_df))
noint_df <- noint_df %>% select(all_of(common_cols))
int_df <- int_df %>% select(all_of(common_cols))

#--- Helper Functions ---

get_metric_label <- function(metric) {
  switch(metric, prs = "Pearson corr", mse = "Mean Squared Error", rmse = "Root Mean Squared Error")
}

is_better_target <- function(int_val, noint_val, metric) {
  if (metric == "prs") int_val > noint_val else int_val < noint_val
}

should_label_drug <- function(int_val, noint_val, metric) {
  if (metric == "prs") {
    noint_val <= 0 & int_val > noint_val
  } else {
    int_val < noint_val
  }
}

plot_int_vs_noint <- function(noint_data, int_data, model_name, metric, out_dir, show_labels = TRUE) {
  col_name <- paste0(metric, "_", model_name)
  
  if (!(col_name %in% colnames(noint_data)) || !(col_name %in% colnames(int_data))) {
    stop(paste("Model column", col_name, "not found in dataframes."))
  }
  
  df_comb <- data.frame(
    drug = noint_data$drug,
    Non_Integrated = noint_data[[col_name]],
    Integrated = int_data[[col_name]]
  )
  
  df_comb <- df_comb %>%
    mutate(
      position = ifelse(is_better_target(Integrated, Non_Integrated, metric), "Better Integrated", "Better Non-Integrated"),
      drug_label = ifelse(should_label_drug(Integrated, Non_Integrated, metric), drug, NA)
    )
  
  p_comp <- ggplot(df_comb, aes(x = Integrated, y = Non_Integrated, color = position)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = c("Better Integrated" = "orange", "Better Non-Integrated" = "darkgrey")) +
    labs(title = paste("Integrated vs Non-Integrated:", toupper(model_name), "(", toupper(metric), ")"),
         x = paste("Integrated (", get_metric_label(metric), ")"),
         y = paste("Non-Integrated (", get_metric_label(metric), ")"), color = "Performance") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  
  if (show_labels) {
    p_comp <- p_comp + geom_text_repel(aes(label = drug_label), size = 2.5, min.segment.length = 0, max.overlaps = Inf, seed = 42, color = "black", show.legend = FALSE)
  }
  
  if (metric == "prs") {
    p_comp <- p_comp + xlim(-0.5, 1) + ylim(-0.5, 1)
  }
  p_comp
  #ggsave(file.path(out_dir, paste0("int_vs_noint_", model_name, "_", metric, ".pdf")), plot = p_comp, width = 8, height = 8)
}

# Generate Plot
plot_int_vs_noint(noint_df, int_df, model_name, metric, out_dir, show_labels)




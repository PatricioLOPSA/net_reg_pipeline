#========================================================================================
# Visualize_cv.R
#
# Purpose:
#   Visualize the main results behind the regression benchmarks script. Generates
#   performance boxplots and comparison scatter plots based on a chosen performance metric.
#
# Usage:
#   Rscript Visualize_cv.R --input_dir <dir> --avg_file <file> --metric <prs|mse|rmse>
#
# Options:
#   -i, --input_dir   Character. Directory containing the Set_*_CV_results.csv files [required]
#   -a, --avg_file    Character. Path to the Average_CV_results.csv file [required]
#   -o, --out_dir     Character. Output directory for the plots [default: ./figures]
#   -m, --metric      Character. Metric to plot: prs, mse, or rmse [default: prs]
#   -h, --help        Show help message and exit
#========================================================================================

library(data.table)
library(janitor)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL, help = "Directory containing the Set CV results [required]"),
  make_option(c("-a", "--avg_file"), type = "character", default = NULL, help = "Path to Average_CV_results.csv [required]"),
  make_option(c("-o", "--out_dir"), type = "character", default = "figures", help = "Output directory for plots [default: %default]"),
  make_option(c("-m", "--metric"), type = "character", default = "prs", help = "Performance metric: prs, mse, or rmse [default: %default]")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
opts <- parse_args(parser)

if (is.null(opts$input_dir) || is.null(opts$avg_file)) {
  print_help(parser)
  stop("Missing required arguments: --input_dir and/or --avg_file", call. = FALSE)
}

if (!dir.exists(opts$out_dir)) dir.create(opts$out_dir, recursive = TRUE)

#--- Helper Functions ---

get_metric_label <- function(metric) {
  switch(metric, prs = "Pearson corr", mse = "Mean Squared Error", rmse = "Root Mean Squared Error")
}

is_better_target <- function(target_val, expr_val, metric) {
  if (metric == "prs") target_val > expr_val else target_val < expr_val
}

should_label_drug <- function(target_val, expr_val, metric) {
  if (metric == "prs") {
    expr_val <= 0 & target_val > expr_val
  } else {
    target_val < expr_val
  }
}

plot_boxplots <- function(data, metric, out_dir) {
  metric_prefix <- paste0(metric, "_")
  perf_long <- data %>%
    select(set, observation_id, starts_with(metric_prefix)) %>%
    pivot_longer(cols = starts_with(metric_prefix), names_to = "full_metric", values_to = "value") %>%
    mutate(model = sub(paste0("^", metric_prefix), "", full_metric),
           model_set = paste(model, set, sep = "__"),
           model_set = fct_reorder(model_set, value, .fun = median))
  
  p_perf <- ggplot(perf_long, aes(x = value, y = model_set)) +
    geom_point(aes(group = observation_id), alpha = 0.5, color = "steelblue", size = 1.2,
               position = position_jitter(height = 0.05, seed = 42)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "white", color = "black", linewidth = 0.5) +
    scale_y_discrete(labels = function(x) sub("__.*$", "", x)) +
    facet_wrap(~ set, ncol = 3, scales = "free_y") +
    labs(title = paste("Model Performance Across Benchmark Sets (", toupper(metric), ")"),
         x = get_metric_label(metric), y = "Model") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"), axis.text = element_text(color = "gray20"))
  
  ggsave(file.path(out_dir, paste0("model_performance_boxplot_", metric, ".pdf")), plot = p_perf, width = 12, height = 10)
}

plot_scatter_sets <- function(data, metric, out_dir) {
  metric_prefix <- paste0(metric, "_")
  comp_long <- data %>%
    select(set, observation_id, drug, 
           expression = !!sym(paste0(metric_prefix, "expr_lasso")),
           Indegree = !!sym(paste0(metric_prefix, "indeg_lasso")), 
           Outdegree = !!sym(paste0(metric_prefix, "outdeg_lasso")), 
           IPFLasso = !!sym(paste0(metric_prefix, "ipflasso"))) %>%
    pivot_longer(cols = c(Indegree, Outdegree, IPFLasso), names_to = "model", values_to = "value") %>%
    mutate(position = ifelse(is_better_target(value, expression, metric), "Better Target Model", "Better Expression"),
           drug_label = ifelse(should_label_drug(value, expression, metric), drug, NA))
  
  p_comp <- ggplot(comp_long, aes(x = value, y = expression, color = position)) +
    geom_point(alpha = 0.6) + geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = c("Better Target Model" = "orange", "Better Expression" = "darkgrey")) +
    facet_grid(model ~ set) +
    labs(title = paste("Target Models vs Expression Baseline across Benchmark Sets (", toupper(metric), ")"),
         x = paste("Target Model (", get_metric_label(metric), ")"),
         y = paste("Expression Baseline (", get_metric_label(metric), ")"), color = "Performance") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  
  if (metric == "prs") {
    p_comp <- p_comp + xlim(-0.5, 1) + ylim(-0.5, 1)
  }
  
  ggsave(file.path(out_dir, paste0("model_vs_expression_scatter_grid_", metric, ".pdf")), plot = p_comp, width = 14, height = 8)
  
  p_comp_labelled <- p_comp +
    geom_text_repel(aes(label = drug_label), size = 2.5, min.segment.length = 0, max.overlaps = Inf, seed = 42, color = "black", show.legend = FALSE) +
    labs(title = paste("Target Models vs Expression Baseline Across Sets (Highlights -", toupper(metric), ")"))
  ggsave(file.path(out_dir, paste0("model_vs_expression_grid_labelled_", metric, ".pdf")), plot = p_comp_labelled, width = 14, height = 8)
}

plot_scatter_avg <- function(avg_data, metric, out_dir) {
  metric_prefix <- paste0(metric, "_")
  comp_avg <- avg_data %>%
    select(drug, 
           expression = !!sym(paste0(metric_prefix, "expr_lasso")),
           Indegree = !!sym(paste0(metric_prefix, "indeg_lasso")), 
           Outdegree = !!sym(paste0(metric_prefix, "outdeg_lasso")),
           Indeg_Outdeg_lasso = !!sym(paste0(metric_prefix, "sep_inout")), 
           IPFLasso = !!sym(paste0(metric_prefix, "ipflasso")),
           PriorityLasso = !!sym(paste0(metric_prefix, "prioritylasso")),
           Separate_models = !!sym(paste0(metric_prefix, "sep_allmods"))) %>%
    pivot_longer(cols = c(Indegree, Outdegree, Indeg_Outdeg_lasso, IPFLasso, PriorityLasso, Separate_models), names_to = "model", values_to = "value") %>%
    mutate(position = ifelse(is_better_target(value, expression, metric), "Better Target Model", "Better Expression"),
           drug_label = ifelse(should_label_drug(value, expression, metric), drug, NA))
  
  p_comp_avg <- ggplot(comp_avg, aes(x = value, y = expression, color = position)) +
    geom_point(alpha = 0.6) +
    geom_text_repel(aes(label = drug_label), size = 2.5, min.segment.length = 0, max.overlaps = Inf, seed = 42, color = "black", show.legend = FALSE) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = c("Better Target Model" = "orange", "Better Expression" = "darkgrey")) +
    facet_wrap(~ model) +
    labs(title = paste("Average Target Models vs Expression Baseline across Benchmark (", toupper(metric), ")"),
         x = paste("Target Model (", get_metric_label(metric), ")"),
         y = paste("Expression Baseline (", get_metric_label(metric), ")"), color = "Performance") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  
  if (metric == "prs") {
    p_comp_avg <- p_comp_avg + xlim(-0.5, 1) + ylim(-0.5, 1)
  }
  
  ggsave(file.path(out_dir, paste0("model_vs_expression_avg_labelled_", metric, ".pdf")), plot = p_comp_avg, width = 12, height = 8)
}


#--- Main Execution ---

# Combine Sets 1 through 5 into a single dataframe
files <- list.files(opts$input_dir, pattern = "^Set_[0-9]+_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No Set_* files found in input_dir.")

data_list <- lapply(seq_along(files), function(i) {
  dt <- fread(files[i], sep = ",") %>% as.data.frame() %>% clean_names()
  dt[is.na(dt)] <- 0
  dt %>% mutate(observation_id = row_number(), set = paste("Set", i))
})
data <- bind_rows(data_list)

averaged_data <- fread(opts$avg_file, sep = ",") %>% as.data.frame() %>% clean_names()
averaged_data[is.na(averaged_data)] <- 0

# Generate Outputs
plot_boxplots(data, opts$metric, opts$out_dir)
plot_scatter_sets(data, opts$metric, opts$out_dir)
plot_scatter_avg(averaged_data, opts$metric, opts$out_dir)







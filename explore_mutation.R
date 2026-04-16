#Exploration of mutation data in beatAML dataset. Not used in final model
library(dplyr)
library(data.table)
library(tidyverse)


setwd('net_reg_pipeline/raw_data')

mutation_data <- fread('beataml_wes_wv1to4_mutations_dbgap.txt', sep = '\t') %>%
as.data.frame()

mutation_data %>% str()
mutation_data$variant_classification %>% table()

mutation_data %>% select(symbol) %>% table() %>% sort(decreasing = TRUE) %>% head(20)
mutation_data %>% select(symbol) %>% table() %>% sort(decreasing = TRUE) %>% length()


# --- Clean  ---
# Replace empty strings and '""' with NA
mutation_data <- mutation_data %>%
  mutate(across(everything(), ~ ifelse(. == '""' | . == "", NA_character_, .)))

# Quick look
mutation_data %>%
  select(dbgap_sample_id, symbol) %>%
  head()


mut_binary <- mutation_data %>%
    distinct(dbgap_sample_id, symbol) %>%
    mutate(mutated = 1L) %>%
    pivot_wider(
      names_from  = symbol,
      values_from = mutated,
      values_fill = 0L,
      names_prefix = ""
    ) %>%
    column_to_rownames("dbgap_sample_id")
  
  colnames(mut_binary) <- paste0(colnames(mut_binary), "_mutated")

#Mutation burden per sample
mutation_burden <- mutation_data %>%
  group_by(dbgap_sample_id) %>%
  summarise(mutation_burden = n()) 
mutation_burden

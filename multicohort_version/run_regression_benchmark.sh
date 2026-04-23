#!/bin/bash

for i in {1..5}; do
  Rscript Regression_benchmark_par.R \
    --train_indeg "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort_nets/Set_${i}/train/lioness_indegree.csv" \
    --test_indeg "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort_nets/Set_${i}/test/lioness_indegree.csv" \
    --train_outdeg "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort_nets/Set_${i}/train/lioness_outdegree.csv" \
    --test_outdeg "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort_nets/Set_${i}/test/lioness_outdegree.csv" \
    --train_expr "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort/Set_${i}/preprocessing_output/Train_Set_Log2norm_Counts.tsv" \
    --test_expr "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/kFold_splits_multicohort/Set_${i}/preprocessing_output/Test_Set_Log2norm_Counts.tsv" \
    --drug_data "/storage/kuijjerarea/plos/net_reg_pipeline_multicohort/cleaned_data/AUC_matrix_complete.tsv" \
    --out_file "Set_${i}_CV_results.csv" \
    --cores 60

    sleep 10
done

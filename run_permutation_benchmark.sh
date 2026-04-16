#!/bin/bash

for i in {1..5}; do
  Rscript permutation_test.R \
    --train_indeg "/storage/kuijjerarea/plos/net_reg_pipeline/kFold_splits_ProtCod_nets/Set_${i}/train/lioness_indegree.csv" \
    --test_indeg "/storage/kuijjerarea/plos/net_reg_pipeline/kFold_splits_ProtCod_nets/Set_${i}/test/lioness_indegree.csv" \
    --train_outdeg "/storage/kuijjerarea/plos/net_reg_pipeline/kFold_splits_ProtCod_nets/Set_${i}/train/lioness_outdegree.csv" \
    --test_outdeg "/storage/kuijjerarea/plos/net_reg_pipeline/kFold_splits_ProtCod_nets/Set_${i}/test/lioness_outdegree.csv" \
    --train_expr "/storage/kuijjerarea/plos/net_reg_pipeline/kFold_splits_ProtCod/Set_${i}/preprocessing_output/Train_Set_VSTnorm_Counts.tsv" \
    --test_expr "/storage/kuijjerarea/plos/net_reg_pipeline/kFold_splits_ProtCod/Set_${i}/preprocessing_output/Test_Set_VSTnorm_Counts.tsv" \
    --drug_data "/storage/kuijjerarea/plos/net_reg_pipeline/cleaned_data/AUC_matrix_complete.tsv" \
    --out_file "Set_${i}_CV_PERM_testonly_results.csv" \
    --cores 32

    sleep 10
done

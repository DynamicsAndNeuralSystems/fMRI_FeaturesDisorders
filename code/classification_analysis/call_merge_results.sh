#!/usr/bin/env bash

data_path=/headnode1/abry4213/data/TS_feature_manuscript/classification_results
classifier_type=Linear_SVM_sklearn

# dataset_ID=UCLA_CNP
# for disorder in SCZ BP ADHD; do
#     python3 merge_main_and_null_results.py --dataset_ID $dataset_ID --disorder $disorder --classifier_type $classifier_type --data_path $data_path
# done

dataset_ID=ABIDE
for disorder in ASD; do
    python3 merge_main_and_null_results.py --dataset_ID $dataset_ID --disorder $disorder --classifier_type $classifier_type --data_path $data_path
done
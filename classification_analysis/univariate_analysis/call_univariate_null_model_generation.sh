#!/bin/bash
# Define directories
project_dir=/headnode1/abry4213/
github_dir=${project_dir}/github/fMRI_FeaturesDisorders/
mkdir -p $github_dir/cluster_output/

export num_null_permutations=1000

# Run each null iteration script
null_perm_scripts=$(find ${github_dir}/classification_analysis/univariate_analysis/*ABIDE* -name "null_iter_*.pbs")
for script in $null_perm_scripts; do
  echo "Now submitting $script"
  qsub $script
done

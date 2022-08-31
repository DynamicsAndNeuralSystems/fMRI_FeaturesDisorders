##########################################################################################
export github_dir=/headnode1/abry4213/github/fMRI_FeaturesDisorders/


cd $github_dir/data_prep_and_QC/

# Prep univariate data
for run_number in 1 2 3 4 5; do
  qsub -v run_number=$run_number -N prepare_univariate_data${run_number} \
  -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data${run_number}_out.txt \
  call_prepare_univariate_data.pbs
done

# Prep pairwise data


# # Prep pairwise data
# # Get data into .npy files
# qsub data_prep_and_QC/call_prepare_pairwise_data.pbs 
# # Run pyspi-distribute
# bash data_prep_and_QC/call_run_pyspi_distribute.sh 
# # Integrate results from pyspi-distribute
# qsub data_prep_and_QC/call_clean_pairwise_data.pbs 
# 
# # Merge subjects with univariate + pairwise data
# qsub data_prep_and_QC/call_merge_samples_univariate_pairwise.pbs

##########################################################################################
# # Univariate linear SVM
# qsub univariate_analysis/call_univariate_classification.pbs 
# 
# # Generate null model fits
# bash univariate_analysis/call_univariate_null_model_generation.sh
# 
# # Integrate null model fits and calculate p-values
# qsub univariate_analysis/call_univariate_null_model_analysis.pbs

##########################################################################################
# # Pairwise linear SVM
# qsub pairwise_analysis/call_pairwise_classification.pbs 

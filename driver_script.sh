##########################################################################################
export github_dir=/headnode1/abry4213/github/fMRI_FeaturesDisorders/
export univariate_feature_set="catch22"
export pairwise_feature_set="pyspi14"

# UCLA Schizophrenia
export dataset_ID="UCLA_Schizophrenia"
export data_path=${project_dir}/data/${dataset_ID}/
export sample_metadata_file=${dataset_ID}_sample_metadata.Rds
export brain_region_lookup="Brain_Region_info.csv"
export noise_procs="AROMA+2P\ AROMA+2P+GMR\ AROMA+2P+DiCER"

# ABIDE ASD
# export dataset_ID="ABIDE_ASD"
# export data_path=${project_dir}/data/${dataset_ID}/
# export sample_metadata_file=${dataset_ID}_sample_metadata.Rds
# export brain_region_lookup="Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"
# export noise_procs="FC1000"

cd $github_dir/data_prep_and_QC/

# Prep univariate data
for run_number in 1; do #2 3 4 5; do
  cmd="qsub -v run_number=$run_number,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,sample_metadata_file=$sample_metadata_file,brain_region_lookup=$brain_region_lookup,noise_procs=$noise_procs \
  -N prepare_univariate_data_${dataset_ID}${run_number} \
  -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data_${dataset_ID}${run_number}_out.txt \
  -m a \
  call_prepare_univariate_data.pbs"
  echo $cmd
  $cmd
done

# # Prep pairwise data

# # Get data into .npy files
# qsub data_prep_and_QC/call_prepare_pairwise_data.pbs 

# # Run pyspi-distribute
# bash data_prep_and_QC/call_run_pyspi_distribute.sh 

# # Integrate results from pyspi-distribute
# for run_number in 1 2 3 4 5; do
#   qsub -v run_number=$run_number -N clean_pairwise_data${run_number} \
#   -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/clean_pairwise_data${run_number}_out.txt \
#   -m a \
#   call_clean_pairwise_data.pbs
# done

# Merge subjects with univariate + pairwise data
# qsub data_prep_and_QC/call_merge_samples_univariate_pairwise.pbs
# for run_number in 1 2 3 4 5; do
#   qsub -v run_number=$run_number -N merge_samples_univariate_pairwise${run_number} \
#   -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/merge_samples_univariate_pairwise${run_number}_out.txt \
#   -m a \
#   call_merge_samples_univariate_pairwise.pbs
# done


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

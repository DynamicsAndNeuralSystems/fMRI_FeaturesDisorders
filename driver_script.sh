##########################################################################################
export github_dir=/headnode1/abry4213/github/
export univariate_feature_set="catch22"
export pairwise_feature_set="pyspi14"
export email="abry4213@uni.sydney.edu.au"
export python_to_use=/headnode1/abry4213/.conda/envs/pyspi/bin/python3

cd $github_dir/fMRI_FeaturesDisorders/data_prep_and_QC/

# UCLA Schizophrenia
export dataset_ID="UCLA_Schizophrenia"
export data_path=/headnode1/abry4213/data/${dataset_ID}/
export sample_metadata_file=${dataset_ID}_sample_metadata.Rds
export brain_region_lookup="Brain_Region_info.csv"
export noise_procs="AROMA+2P;AROMA+2P+GMR;AROMA+2P+DiCER"
export main_noise_proc="AROMA+2P+GMR"

# ABIDE ASD
# export dataset_ID="ABIDE_ASD"
# export data_path=${project_dir}/data/${dataset_ID}/
# export sample_metadata_file=${dataset_ID}_sample_metadata.Rds
# export brain_region_lookup="Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"
# export noise_procs="FC1000"

# # Prep univariate data
# for run_number in 1 2 3 4 5; do
#   qsub -v run_number=$run_number,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,sample_metadata_file=$sample_metadata_file,brain_region_lookup=$brain_region_lookup,noise_procs=$noise_procs \
#   -N prepare_univariate_data_${dataset_ID}${run_number} \
#   -o $github_dir/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data_${dataset_ID}${run_number}_out.txt \
#   -m a \
#   call_prepare_univariate_data.pbs
# done

# # Prep pairwise data
# # Get data into .npy files
# cmd="qsub -N prepare_pairwise_data_${dataset_ID} -q yossarian -j oe \
# -v github_dir=$github_dir,data_path=$data_path,noise_procs=$noise_procs,dataset_ID=$dataset_ID,sample_metadata_file=$sample_metadata_file \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/prepare_pairwise_data_${dataset_ID}.txt \
# -l select=1:ncpus=1:mem=20GB -l walltime=4:00:00 -M $email -m a -V \
# prepare_pairwise_data.sh"
# $cmd

# Run pyspi-distribute
# bash call_run_pyspi_distribute.sh \
# $github_dir \
# ${github_dir}/fMRI_FeaturesDisorders/data_prep_and_QC/pyspi14_config.yaml \
# $email \
# $dataset_ID \
# $data_path \
# $noise_procs \
# 8

# Integrate results from pyspi-distribute
# for run_number in 1 2 3 4 5; do
#   qsub -v run_number=$run_number,github_dir=$github_dir,data_path=$data_path,python_to_use=$python_to_use,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,brain_region_lookup=$brain_region_lookup,noise_procs=$noise_procs,main_noise_proc=$main_noise_proc,dataset_ID=$dataset_ID \
#   -N clean_pairwise_data_${dataset_ID}${run_number} \
#   -o ${github_dir}/fMRI_FeaturesDisorders/cluster_output/clean_pairwise_data_${dataset_ID}${run_number}_out.txt \
#   -m a -M $email \
#   call_clean_pairwise_data.pbs
# done

# Merge subjects with univariate + pairwise data
# for run_number in 1 2 3 4 5; do
#   qsub -v run_number=$run_number,github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set \
#   -N merge_samples_univariate_pairwise_${dataset_ID}${run_number} \
#   -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/merge_samples_univariate_pairwise_${dataset_ID}${run_number}_out.txt \
#   -m a -M $email \
#   call_merge_samples_univariate_pairwise.pbs
# done


##########################################################################################
cd $github_dir/fMRI_FeaturesDisorders/univariate_analysis/

# Univariate linear SVM
# for run_number in 1 2 3 4 5; do
#   qsub -v run_number=$run_number,github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,noise_procs=$noise_procs,main_noise_proc=$main_noise_proc \
#   -N run_univariate_classification_${dataset_ID}${run_number} \
#   -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_classification_${dataset_ID}${run_number}_out.txt \
#   -m a -M $email \
#   call_univariate_classification.pbs 
# done

# # Generate null model fits
# null_perm_scripts=$(find ${github_dir}/fMRI_FeaturesDisorders/univariate_analysis/null_pbs_scripts/*ROI* -name "null_iter_*.pbs")
# for script in $null_perm_scripts; do
#   echo "Now submitting $script"
#   qsub $script
# done

# Integrate null model fits and calculate p-values
for run_number in 1 2 3 4 5; do
  qsub -v run_number=$run_number,github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,sample_metadata_file=$sample_metadata_file,main_noise_proc=$main_noise_proc \
  -N run_univariate_null_model_analysis${dataset_ID}${run_number} \
  -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_null_model_analysis_${dataset_ID}${run_number}_out.txt \
  -m a -M $email \
  call_univariate_null_model_analysis.pbs 
done

##########################################################################################
# # Pairwise linear SVM
# qsub pairwise_analysis/call_pairwise_classification.pbs 

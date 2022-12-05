##########################################################################################
export github_dir=/headnode1/abry4213/github/
export univariate_feature_set="catch22"
export pairwise_feature_set="pyspi14"
export email="abry4213@uni.sydney.edu.au"
export conda_env="pyspi_annie"
export pyspi_ncpus=1
export pyspi_mem=40
export python_to_use=/headnode1/abry4213/.conda/envs/${conda_env}/bin/python3

cd $github_dir/fMRI_FeaturesDisorders/data_prep_and_QC/

# Define data paths
export dataset_ID="UCLA_CNP_ABIDE_ASD"
export data_path=/headnode1/abry4213/data/${dataset_ID}/
export UCLA_sample_metadata_file=UCLA_sample_metadata.Rds
export ABIDE_sample_metadata_file=ABIDE_sample_metadata.Rds
export noise_proc_UCLA="AROMA+2P+GMR"
export noise_proc_ABIDE="FC1000"
export label_vars="Diagnosis"
export pkl_file="calc_pyspi14.pkl"

##########################################################################################
cd $github_dir/fMRI_FeaturesDisorders/classification_analysis/univariate_analysis/

# Univariate linear SVM
# qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,noise_procs=$noise_procs,main_noise_proc=$main_noise_proc \
# -N run_univariate_classification_${dataset_ID} \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_classification_${dataset_ID}_out.txt \
# -m a -M $email \
# call_univariate_classification.pbs 

# # Generate null model fits
# null_perm_scripts=$(find ${github_dir}/fMRI_FeaturesDisorders/classification_analysis/univariate_analysis/null_pbs_scripts/*${dataset_ID}*${univariate_feature_set}_inv_prob_null_model_fits* -name "null_iter_*.pbs")
# for script in $null_perm_scripts; do
#   echo "Now submitting $script"
#   qsub $script
# done

# # Integrate null model fits and calculate p-values
# qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,sample_metadata_file=$sample_metadata_file,main_noise_proc=$main_noise_proc \
# -N run_univariate_null_model_analysis${dataset_ID} \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_null_model_analysis_${dataset_ID}_out.txt \
# -m a -M $email \
# call_univariate_null_model_analysis.pbs 

##########################################################################################
# Pairwise analysis
cd $github_dir/fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/

# # Pairwise linear SVM
# qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,main_noise_proc=$main_noise_proc  \
# -N pairwise_classification${dataset_ID} \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/pairwise_classification_${dataset_ID}_out.txt \
# -m a -M $email \
# call_pairwise_classification.pbs 

# # Generate null model fits
# null_perm_scripts=$(find ${github_dir}/fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/null_pbs_scripts/*${dataset_ID}*${pairwise_feature_set}_inv_prob_null_model_fits* -name "null_iter*.pbs")
# for script in $null_perm_scripts; do
#   echo "Now submitting $script"
#   qsub $script
# done

# # Integrate null model fits and calculate p-values
# qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,main_noise_proc=$main_noise_proc \
# -N pairwise_null_model_analysis${dataset_ID} \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/pairwise_null_model_analysis_${dataset_ID}_out.txt \
# -m a -M $email \
# call_pairwise_null_model_analysis.pbs 


##########################################################################################
# Combined univariate + pairwise analysis
cd $github_dir/fMRI_FeaturesDisorders/classification_analysis/combined_univariate_pairwise/

# qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,main_noise_proc=$main_noise_proc,email=$email  \
# -N combined_univariate_pairwise_classification_${dataset_ID} \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/combined_univariate_pairwise_classification_${dataset_ID}_out.txt \
# -m a -M $email \
# call_combined_univariate_pairwise_classification.pbs 

# # Generate null model fits
# null_perm_scripts=$(find ${github_dir}/fMRI_FeaturesDisorders/classification_analysis/combined_univariate_pairwise/null_pbs_scripts/*${dataset_ID}*${univariate_feature_set}_pairwise_${pairwise_feature_set}_inv_prob_null_model_fits* -name "null_iter_*.pbs")
# for script in $null_perm_scripts; do
#   echo "Now submitting $script"
#   qsub $script
# done

# # Integrate null model fits and calculate p-values
# qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,main_noise_proc=$main_noise_proc \
# -N uni_pairwise_null_model_analysis_${dataset_ID} \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/combined_univariate_pairwise_null_model_analysis_${dataset_ID}_out.txt \
# -m a -M $email \
# call_combined_univariate_pairwise_null_model_analysis.pbs 
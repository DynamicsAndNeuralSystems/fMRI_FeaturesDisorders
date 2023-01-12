##########################################################################################
export github_dir=/headnode1/abry4213/github/
export pairwise_feature_set="pyspi14"
export email="abry4213@uni.sydney.edu.au"
export conda_env="pyspi"
export pyspi_ncpus=1
export pyspi_mem=40
export python_to_use=/headnode1/abry4213/.conda/envs/${conda_env}/bin/python3

#module load Anaconda3-5.1.0
/usr/physics/Modules/3.2.8/bin/modulecmd bash load Anaconda3-5.1.0 --silent
source /usr/physics/python/anaconda3/etc/profile.d/conda.sh 

# Activate the given conda environment
conda activate $conda_env

# Define data paths
export dataset_ID="UCLA_CNP_ABIDE_ASD"
export data_path=/headnode1/abry4213/data/${dataset_ID}/
export UCLA_sample_metadata_file=UCLA_sample_metadata.Rds
export ABIDE_sample_metadata_file=ABIDE_sample_metadata.Rds
export sample_metadata_file=UCLA_CNP_ABIDE_ASD_sample_metadata.Rds
export noise_proc_UCLA_CNP="AROMA+2P+GMR"
export noise_label_UCLA_CNP="AROMA_2P_GMR"
export noise_proc_ABIDE_ASD="FC1000"
export noise_label_ABIDE_ASD="FC1000"
export label_vars="Diagnosis"
export pkl_file="calc_pyspi14.pkl"

# ##########################################################################################
# # Prepare the UCLA CNP data
# bash $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/dataset_specific_files/prepare_UCLA_data.sh

# # Prepare the ABIDE ASD data
# bash $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/dataset_specific_files/prepare_ABIDE_data.sh

# ##########################################################################################
# # Merge the UCLA CNP and ABIDE ASD datasets that were pre-processed separately
# Rscript $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/dataset_specific_files/merge_UCLA_CNP_ABIDE_ASD_data.R

##########################################################################################
# Univariate analysis
# cd $github_dir/fMRI_FeaturesDisorders/classification_analysis/univariate_analysis/

# # Univariate linear SVM
# # for univariate_feature_set in catch22 catch24 catch2; do
# for univariate_feature_set in catch22; do
#   qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,input_feature_data_file=UCLA_CNP_${noise_label_UCLA_CNP}_ABIDE_ASD_${noise_label_ABIDE_ASD}_${univariate_feature_set}_filtered_zscored.Rds,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file \
#   -N run_univariate_classification_${dataset_ID} \
#   -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_classification_${dataset_ID}_out.txt \
#   -m a -M $email \
#   call_univariate_classification.pbs 
# done

# Generate null model fits
# for univariate_feature_set in catch22 catch24 catch2; do
#   null_perm_scripts=$(find ${github_dir}/fMRI_FeaturesDisorders/classification_analysis/univariate_analysis/null_pbs_scripts/*${dataset_ID}*${univariate_feature_set}_inv_prob_null_model_fits* -name "null_iter_*.pbs")
#   for script in $null_perm_scripts; do
#     echo "Now submitting $script"
#     qsub $script
#   done
# done

# Integrate null model fits and calculate p-values
# for univariate_feature_set in catch22 catch24 catch2; do
#   qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,feature_set=$univariate_feature_set,sample_metadata_file=$sample_metadata_file \
#   -N run_univariate_null_model_analysis_${dataset_ID}_${univariate_feature_set} \
#   -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_null_model_analysis_${dataset_ID}_${univariate_feature_set}_out.txt \
#   -m a -M $email \
#   call_univariate_null_model_analysis.pbs
# done

##########################################################################################
# Pairwise analysis
cd $github_dir/fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/
export univariate_feature_set="catch22"

# Pairwise linear SVM
qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,input_feature_data_file=UCLA_CNP_${noise_label_UCLA_CNP}_ABIDE_ASD_${noise_label_ABIDE_ASD}_${pairwise_feature_set}_filtered_zscored.Rds,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_set=$univariate_feature_set,sample_metadata_file=$sample_metadata_file  \
-N pairwise_classification${dataset_ID} \
-o $github_dir/fMRI_FeaturesDisorders/cluster_output/pairwise_classification_${dataset_ID}_out.txt \
-m a -M $email \
call_pairwise_classification.pbs 

# # Generate null model fits
# null_perm_scripts=$(find ${github_dir}/fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/null_pbs_scripts/*${dataset_ID}*${pairwise_feature_set}_inv_prob_null_model_fits* -name "null_iter*.pbs")
# for script in $null_perm_scripts; do
#   echo "Now submitting $script"
#   qsub $script
# done

# # Integrate null model fits and calculate p-values
# qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file \
# -N pairwise_null_model_analysis${dataset_ID} \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/pairwise_null_model_analysis_${dataset_ID}_out.txt \
# -m a -M $email \
# call_pairwise_null_model_analysis.pbs 


##########################################################################################
# Combined univariate + pairwise analysis
cd $github_dir/fMRI_FeaturesDisorders/classification_analysis/combined_univariate_pairwise/

# qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,email=$email  \
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
# qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file \
# -N uni_pairwise_null_model_analysis_${dataset_ID} \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/combined_univariate_pairwise_null_model_analysis_${dataset_ID}_out.txt \
# -m a -M $email \
# call_combined_univariate_pairwise_null_model_analysis.pbs 

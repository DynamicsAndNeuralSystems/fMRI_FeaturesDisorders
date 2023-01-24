##########################################################################################
export github_dir=/headnode1/abry4213/github/
export pairwise_feature_set="pyspi14"
export email="abry4213@uni.sydney.edu.au"
export conda_env="pyspi"
export pyspi_ncpus=1
export pyspi_mem=40
export pyspi_config=${github_dir}/fMRI_FeaturesDisorders/prep_data_and_QC/pyspi14_config.yaml
export sample_yaml="sample.yaml"
export python_to_use=/headnode1/abry4213/.conda/envs/${conda_env}/bin/python3
export pkl_file="calc_pyspi14.pkl"
export label_vars="Diagnosis"
export num_null_iters=1000

#module load Anaconda3-5.1.0
/usr/physics/Modules/3.2.8/bin/modulecmd bash load Anaconda3-5.1.0 --silent
source /usr/physics/python/anaconda3/etc/profile.d/conda.sh 

# Activate the given conda environment
conda activate $conda_env

##########################################################################################
# UCLA CNP

# Define data paths
export dataset_ID="UCLA_CNP"
export data_path=/headnode1/abry4213/data/UCLA_CNP/
export sample_metadata_file="UCLA_CNP_sample_metadata.feather"
export noise_proc="AROMA+2P+GMR"
export brain_region_lookup="UCLA_CNP_Brain_Region_Lookup.feather"
export pyspi_walltime_hrs=2

# # Run feature extraction
# bash $github_dir/fMRI_FeaturesDisorders/run_feature_extraction.sh 

# # Run feature cleaning
# bash $github_dir/fMRI_FeaturesDisorders/run_feature_cleaning.sh

# Run univariate linear SVM
for comparison_group in Schizophrenia Bipolar ADHD; do
    for univariate_feature_set in catch2 catch22 catch24; do
        qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,num_null_iters=$num_null_iters \
        -N ${dataset_ID}_${comparison_group}_${feature_set} \
        -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_classification_${dataset_ID}_${comparison_group}_${feature_set}_out.txt \
        -m a \
        $github_dir/fMRI_FeaturesDisorders/classification_analysis/univariate_analysis/call_univariate_classification.pbs
    done
done
##########################################################################################
# # ABIDE ASD

# # Define data paths
# export dataset_ID="ABIDE_ASD"
# export data_path=/headnode1/abry4213/data/ABIDE_ASD/
# export sample_metadata_file=ABIDE_ASD_sample_metadata.feather
# export noise_proc="FC1000"
# export brain_region_lookup="ABIDE_ASD_Brain_Region_Lookup.feather"
# export pyspi_walltime_hrs=3

# # Run feature extraction
# bash $github_dir/fMRI_FeaturesDisorders/run_feature_extraction.sh 

# # Run feature cleaning
# bash $github_dir/fMRI_FeaturesDisorders/run_feature_cleaning.sh
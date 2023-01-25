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
export noise_label="AROMA_2P_GMR"
export brain_region_lookup="UCLA_CNP_Brain_Region_Lookup.feather"
export pyspi_walltime_hrs=2

# Run feature extraction
bash $github_dir/fMRI_FeaturesDisorders/data_prep.sh 

##########################################################################################
# ABIDE ASD

# Define data paths
export dataset_ID="ABIDE_ASD"
export data_path=/headnode1/abry4213/data/ABIDE_ASD/
export sample_metadata_file=ABIDE_ASD_sample_metadata.feather
export noise_proc="FC1000"
export noise_label="FC1000"
export brain_region_lookup="ABIDE_ASD_Brain_Region_Lookup.feather"
export pyspi_walltime_hrs=3

# Run feature extraction
bash $github_dir/fMRI_FeaturesDisorders/data_prep.sh 

##########################################################################################
export github_dir=/headnode1/abry4213/github/
export univariate_feature_set="catch22"
export pairwise_feature_set="pyspi14"
export email="abry4213@uni.sydney.edu.au"
export conda_env="pyspi"
export pyspi_ncpus=1
export pyspi_mem=40
export python_to_use=/headnode1/abry4213/.conda/envs/${conda_env}/bin/python3

cd $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/

# ABIDE ASD
export dataset_ID="ABIDE_ASD"
export data_path=/headnode1/abry4213/data/UCLA_CNP_ABIDE_ASD/
export sample_metadata_file=${dataset_ID}_sample_metadata.Rds
export brain_region_lookup="ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"
export sample_yaml="sample.yaml"
export pyspi_config=${github_dir}/fMRI_FeaturesDisorders/prep_data_and_QC/pyspi14_config.yaml
export noise_proc="FC1000"
export label_vars="Diagnosis"
export pyspi_walltime_hrs=1
export pkl_file="calc_pyspi14.pkl"

# # Prep univariate data
# qsub -v dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,sample_metadata_file=$sample_metadata_file,brain_region_lookup=$brain_region_lookup,noise_proc=$noise_proc \
# -N prepare_univariate_data_${dataset_ID} \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data_${dataset_ID}_out.txt \
# -m a \
# call_prepare_univariate_data.pbs

# # Prep pairwise data
# # Get data into .npy files
# cmd="qsub -N prepare_pairwise_data_${dataset_ID} -q yossarian -j oe \
# -v github_dir=$github_dir,data_path=$data_path,noise_proc=$noise_proc,dataset_ID=$dataset_ID,sample_metadata_file=$sample_metadata_file,label_vars=$label_vars \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/prepare_pairwise_data_${dataset_ID}.txt \
# -l select=1:ncpus=1:mem=20GB -l walltime=4:00:00 -M $email -m a -V \
# prepare_pairwise_data.sh"
# $cmd

# # Run pyspi-distribute
# # First run for subjects
# bash $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/call_run_pyspi_distribute.sh \
# $github_dir \
# $pyspi_config \
# $email \
# $dataset_ID \
# $data_path/raw_data/${dataset_ID}/numpy_files \
# $noise_proc \
# $pyspi_walltime_hrs \
# $pyspi_mem \
# $pyspi_ncpus \
# $pkl_file \
# $sample_yaml \
# $conda_env

# Integrate results from pyspi-distribute
qsub -v github_dir=$github_dir,data_path=$data_path,pkl_file=$pkl_file,python_to_use=$python_to_use,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,brain_region_lookup=$brain_region_lookup,noise_proc=$noise_proc,dataset_ID=$dataset_ID \
-N clean_pairwise_data_${dataset_ID} \
-o ${github_dir}/fMRI_FeaturesDisorders/cluster_output/clean_pairwise_data_${dataset_ID}_out.txt \
-m a -M $email \
call_clean_pairwise_data.pbs
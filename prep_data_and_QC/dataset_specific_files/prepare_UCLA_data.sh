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

# UCLA CNP
export dataset_ID="UCLA_CNP"
export data_path=/headnode1/abry4213/data/${dataset_ID}/
export sample_metadata_file=${dataset_ID}_sample_metadata.Rds
export brain_region_lookup="Brain_Region_info.csv"
export sample_yaml="sample.yaml"
export pyspi_config=${github_dir}/fMRI_FeaturesDisorders/prep_data_and_QC/pyspi14_config.yaml
export noise_procs="AROMA+2P+GMR"
export main_noise_proc="AROMA+2P+GMR"
export label_vars="Diagnosis"
export pyspi_walltime_hrs=1
export pkl_file="calc_pyspi14.pkl"

# # Prep univariate data
# qsub -v dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,sample_metadata_file=$sample_metadata_file,brain_region_lookup=$brain_region_lookup,noise_procs=$noise_procs \
# -N prepare_univariate_data_${dataset_ID} \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data_${dataset_ID}_out.txt \
# -m a \
# call_prepare_univariate_data.pbs

# # Prep pairwise data
# # Get data into .npy files
# cmd="qsub -N prepare_pairwise_data_${dataset_ID} -q yossarian -j oe \
# -v github_dir=$github_dir,data_path=$data_path,noise_procs=$noise_procs,dataset_ID=$dataset_ID,sample_metadata_file=$sample_metadata_file,label_vars=$label_vars \
# -o $github_dir/fMRI_FeaturesDisorders/cluster_output/prepare_pairwise_data_${dataset_ID}.txt \
# -l select=1:ncpus=1:mem=20GB -l walltime=2:00:00 -M $email -m a -V \
# prepare_pairwise_data.sh"
# $cmd

# Run pyspi-distribute
bash $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/call_run_pyspi_distribute.sh \
$github_dir \
$pyspi_config \
$email \
$dataset_ID \
$data_path \
$noise_procs \
$pyspi_walltime_hrs \
$pyspi_mem \
$pyspi_ncpus \
$pkl_file \
$sample_yaml \
$conda_env

# # Integrate results from pyspi-distribute
# qsub -v github_dir=$github_dir,data_path=$data_path,pkl_file=$pkl_file,python_to_use=$python_to_use,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,brain_region_lookup=$brain_region_lookup,noise_procs=$noise_procs,main_noise_proc=$main_noise_proc,dataset_ID=$dataset_ID \
# -N clean_pairwise_data_${dataset_ID} \
# -o ${github_dir}/fMRI_FeaturesDisorders/cluster_output/clean_pairwise_data_${dataset_ID}_out.txt \
# -m a -M $email \
# call_clean_pairwise_data.pbs

# # Merge subjects with univariate + pairwise data
# qsub -v github_dir=$github_dir,data_path=$data_path,dataset_ID=$dataset_ID,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set \
# -N merge_samples_univariate_pairwise_${dataset_ID} \
# -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/merge_samples_univariate_pairwise_${dataset_ID}_out.txt \
# -m a -M $email \
# call_merge_samples_univariate_pairwise.pbs

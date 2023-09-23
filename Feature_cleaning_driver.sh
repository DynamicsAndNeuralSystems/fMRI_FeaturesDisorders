##########################################################################################
export github_dir=/headnode1/abry4213/github/
export pairwise_feature_set="pyspi14"
export univariate_feature_set="catch25"
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
cd $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/

##########################################################################################
# UCLA CNP

# Define data paths
export dataset_ID="UCLA_CNP"
export data_path=/headnode1/abry4213/data/UCLA_CNP/
export sample_metadata_file="UCLA_CNP_sample_metadata.feather"
export noise_proc="AROMA+2P+GMR"
export noise_label="AROMA_2P_GMR"
export brain_region_lookup="UCLA_CNP_Brain_Region_Lookup.feather"

# Merge pyspi data from pkl files
qsub -N merge_${dataset_ID}_${pairwise_feature_set} -q yossarian -j oe \
-v github_dir=$github_dir,conda_env=$conda_env,data_path=$data_path,noise_proc=$noise_proc,dataset_ID=$dataset_ID,pairwise_feature_set=$pairwise_feature_set,univariate_feature_set=$univariate_feature_set,brain_region_lookup=$brain_region_lookup \
-o $github_dir/fMRI_FeaturesDisorders/cluster_output/merge_pairwise_data_${dataset_ID}_${pairwise_feature_set}.txt \
-l select=1:ncpus=4:mem=120GB:mpiprocs=4 -l walltime=1:00:00 -M $email -m a -V \
prep_data_and_QC/call_merge_pairwise_data.pbs

# Run feature cleaning for pairwise data
qsub -N ${dataset_ID}_${pairwise_feature_set}_${univariate_feature_set} -q yossarian -j oe \
-v github_dir=$github_dir,conda_env=$conda_env,data_path=$data_path,noise_proc=$noise_proc,dataset_ID=$dataset_ID,pairwise_feature_set=$pairwise_feature_set,univariate_feature_set=$univariate_feature_set,brain_region_lookup=$brain_region_lookup \
-o $github_dir/fMRI_FeaturesDisorders/cluster_output/all_feature_merge_${dataset_ID}_${univariate_feature_set}_${pairwise_feature_set}.txt \
-l select=1:ncpus=1:mem=20GB:mpiprocs=1 -l walltime=1:00:00 -M $email -m a -V call_final_feature_merges.pbs

##########################################################################################
# # ABIDE ASD

# Define data paths
export dataset_ID="ABIDE_ASD"
export data_path=/headnode1/abry4213/data/ABIDE_ASD/
export sample_metadata_file=ABIDE_ASD_sample_metadata.feather
export noise_proc="FC1000"
export noise_label="FC1000"
export brain_region_lookup="ABIDE_ASD_Brain_Region_Lookup.feather"

# Merge pyspi data from pkl files
qsub -N merge_${dataset_ID}_${pairwise_feature_set} -q yossarian -j oe \
-v github_dir=$github_dir,conda_env=$conda_env,data_path=$data_path,noise_proc=$noise_proc,dataset_ID=$dataset_ID,pairwise_feature_set=$pairwise_feature_set,univariate_feature_set=$univariate_feature_set,brain_region_lookup=$brain_region_lookup \
-o $github_dir/fMRI_FeaturesDisorders/cluster_output/merge_pairwise_data_${dataset_ID}_${pairwise_feature_set}.txt \
-l select=1:ncpus=4:mem=120GB:mpiprocs=4 -l walltime=1:00:00 -M $email -m a -V \
prep_data_and_QC/call_merge_pairwise_data.pbs

# Run feature cleaning for pairwise data
qsub -N ${dataset_ID}_${pairwise_feature_set}_${univariate_feature_set} -q yossarian -j oe \
-v github_dir=$github_dir,conda_env=$conda_env,data_path=$data_path,noise_proc=$noise_proc,dataset_ID=$dataset_ID,pairwise_feature_set=$pairwise_feature_set,univariate_feature_set=$univariate_feature_set,brain_region_lookup=$brain_region_lookup \
-o $github_dir/fMRI_FeaturesDisorders/cluster_output/all_feature_merge_${dataset_ID}_${univariate_feature_set}_${pairwise_feature_set}.txt \
-l select=1:ncpus=1:mem=20GB:mpiprocs=1 -l walltime=1:00:00 -M $email -m a -V call_final_feature_merges.pbs
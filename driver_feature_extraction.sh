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
export dataset_ID="UCLA_CNP"
export data_path=/headnode1/abry4213/data/${dataset_ID}/
export sample_metadata_file=UCLA_sample_metadata.Rds
export noise_proc="AROMA+2P+GMR"
export noise_label="AROMA_2P_GMR"
export label_vars="Diagnosis"
export pkl_file="calc_pyspi14.pkl"

# ##########################################################################################
# Read metadata into R and python
echo "Now preparing metadata"
Rscript prep_data_and_QC/dataset_specific_files/parse_metadata_${dataset_ID}.R

# ##########################################################################################
# Prep univariate data

# Prep univariate data in R
echo "Now preparing univariate data"
qsub -v dataset_ID=$dataset_ID,data_path=$data_path,univariate_feature_set=$univariate_feature_set,sample_metadata_file=$sample_metadata_file,brain_region_lookup=$brain_region_lookup,noise_proc=$noise_proc \
-N prepare_univariate_data_${dataset_ID} \
-o $github_dir/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data_${dataset_ID}_out.txt \
-m a \
call_prepare_univariate_data.pbs

# ########################################################################################
# Prep pairwise data
# Get data into .npy files
cmd="qsub -N prepare_pairwise_data_${dataset_ID} -q yossarian -j oe \
-v github_dir=$github_dir,conda_env=$conda_env,data_path=$data_path/raw_data/${dataset_ID},noise_proc=$noise_proc,dataset_ID=$dataset_ID,sample_metadata_file=$sample_metadata_file,label_vars=$label_vars \
-o $github_dir/fMRI_FeaturesDisorders/cluster_output/prepare_pairwise_data_${dataset_ID}.txt \
-l select=1:ncpus=1:mem=20GB -l walltime=2:00:00 -M $email -m a -V \
prepare_pairwise_data.sh"
echo "Now preparing pairwise data"
$cmd

# Run pyspi-distribute
echo "Now submitting pyspi-distribute jobs"
bash $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/call_run_pyspi_distribute.sh \
$github_dir \
$pyspi_config \
$email \
$dataset_ID \
$data_path/raw_data/numpy_files \
$noise_proc \
$pyspi_walltime_hrs \
$pyspi_mem \
$pyspi_ncpus \
$pkl_file \
$sample_yaml \
$conda_env


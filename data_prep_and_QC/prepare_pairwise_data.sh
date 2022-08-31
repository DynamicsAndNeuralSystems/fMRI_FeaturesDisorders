#!/usr/bin/env bash

# Load python
#module load Anaconda3-5.1.0
/usr/physics/Modules/3.2.8/bin/modulecmd bash load Anaconda3-5.1.0 --silent
source /usr/physics/python/anaconda3/etc/profile.d/conda.sh 
conda activate pyspi

# Define paths
export fmri_github_dir=${github_dir}/fMRI_FeaturesDisorders/

noise_procs_list=$(echo $noise_procs | sed "s/;/ /g")

###############################################################################

# Split data into numpy files
cmd="python3 ${fmri_github_dir}/helper_functions/data_prep_and_QC/split_MTS_into_npy.py \
--github_dir ${fmri_github_dir} \
--data_path ${data_path} \
--noise_procs "$noise_procs_list" \
--dataset_ID ${dataset_ID}"
$cmd

# Create sample.yaml file for this dataset
for noise_proc in $noise_procs_list
do
    noise_label=$(echo $noise_proc | sed "s/\+/_/g")
    cmd="Rscript $github_dir/pyspi-distribute/create_yaml_for_samples.R \
    --data_dir ${data_path}/raw_data/numpy_files/${noise_label}/ \
    --sample_metadata_file ${data_path}/${sample_metadata_file} \
    --ID_var Sample_ID \
    --label_vars Diagnosis \
    --dim_order ps \
    --overwrite"
    $cmd
done
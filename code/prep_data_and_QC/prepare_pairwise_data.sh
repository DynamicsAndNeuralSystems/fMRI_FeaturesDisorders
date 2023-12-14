#!/usr/bin/env bash

# Load python
#module load Anaconda3-5.1.0
/usr/physics/Modules/3.2.8/bin/modulecmd bash load Anaconda3-5.1.0 --silent
source /usr/physics/python/anaconda3/etc/profile.d/conda.sh 
conda activate ${conda_env}

# Define paths
export fmri_github_dir=${github_dir}/fMRI_FeaturesDisorders/

###############################################################################

# Split data into numpy files
cmd="python3 ${fmri_github_dir}/helper_functions/data_prep_and_QC/split_MTS_into_npy.py \
--github_dir ${fmri_github_dir} \
--data_path ${data_path} \
--noise_proc "$noise_proc" \
--dataset_ID ${dataset_ID}"
echo $cmd
$cmd

# Create sample.yaml file for this dataset
noise_label=$(echo $noise_proc | sed "s/\+/_/g")

cmd="Rscript $github_dir/pyspi-distribute/create_yaml_for_samples.R \
--data_dir ${data_path}/raw_data/numpy_files/${noise_label}/ \
--sample_metadata_file ${data_path}/study_metadata/${sample_metadata_file} \
--ID_var Sample_ID \
--label_vars $label_vars \
--dim_order ps"
echo $cmd
$cmd
#!/usr/bin/env bash

# Activate given conda env
conda activate ${conda_env}
export pyspi_script_dir=${github_dir}/pyspi-distribute/

noise_label=$(echo $noise_proc | sed "s/\+/_/g")
echo $noise_label
cmd="python $pyspi_script_dir/distribute_jobs.py \
--data_dir ${data_path}/raw_data/numpy_files/${noise_label}/ \
--calc_file_name $pkl_file \
--compute_file $pyspi_script_dir/pyspi_compute.py \
--template_pbs_file $pyspi_script_dir/template.pbs \
--pyspi_config $pyspi_config \
--sample_yaml ${data_path}/raw_data/numpy_files/${noise_label}/${sample_yaml} \
--pbs_notify a \
--email $email \
--conda_env $conda_env \
--queue $queue \
--walltime_hrs $pyspi_walltime_hrs \
--cpu $pyspi_ncpus \
--mem $pyspi_mem \
--table_only"
$cmd
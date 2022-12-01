#!/usr/bin/env bash

#module load Anaconda3-5.1.0
/usr/physics/Modules/3.2.8/bin/modulecmd bash load Anaconda3-5.1.0 --silent
source /usr/physics/python/anaconda3/etc/profile.d/conda.sh 
conda activate pyspi

export github_dir="$1"
export config_file="$2"
export email="$3"
export dataset_ID="$4"
export data_path="$5"
export noise_procs="$6"
export walltime_hrs="$7"
export mem="$8"
export ncpus="$9"
export calc_file_name="${10}"
export sample_yaml="${11}"

noise_procs_list=$(echo $noise_procs | sed "s/;/ /g")
export pyspi_script_dir=${github_dir}/pyspi-distribute/

for noise_proc in $noise_procs_list
do
    noise_label=$(echo $noise_proc | sed "s/\+/_/g")
    echo $noise_label
    python $pyspi_script_dir/distribute_jobs.py \
    --data_dir ${data_path}/raw_data/numpy_files/${noise_label}/ \
    --calc_file_name $calc_file_name \
    --compute_file $pyspi_script_dir/pyspi_compute.py \
    --template_pbs_file $pyspi_script_dir/template.pbs \
    --pyspi_config $config_file \
    --sample_yaml ${data_path}/raw_data/numpy_files/${noise_label}/${sample_yaml} \
    --pbs_notify a \
    --email $email \
    --walltime_hrs $walltime_hrs \
    --cpu $ncpus \
    --mem $mem \
    --table_only \
    --overwrite_pkl
done
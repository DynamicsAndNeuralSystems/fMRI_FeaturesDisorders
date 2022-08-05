export study=/headnode1/abry4213
export github_dir=${study}/github/
export config_file=${github_dir}/fMRI_FeaturesDisorders/data_prep_and_QC/pyspi19_config.yaml
export pyspi_script_dir=${study}/github/pyspi-distribute/

# UCLA Schizophrenia
dataset_ID="UCLA_Schizophrenia"
for noise_label in AROMA_2P AROMA_2P_GMR AROMA_2P_DiCER
do
    python $pyspi_script_dir/distribute_jobs.py --data_dir ${study}/data/${dataset_ID}/pydata/${noise_label}/ \
    --compute_file $pyspi_script_dir/pyspi_compute.py \
    --template_pbs_file $pyspi_script_dir/template.pbs \
    --pyspi_config $config_file \
    --sample_yaml ${study}/data/${dataset_ID}/pydata/${noise_label}/sample.yaml \
    --pbs_notify a \
    --email abry4213@uni.sydney.edu.au --walltime_hrs 2 --cpu 2 --mem 8 --table_only
done


# ABIDE ASD
dataset_ID="ABIDE_ASD"
for noise_label in FC1000
do
    python $pyspi_script_dir/distribute_jobs.py --data_dir ${study}/data/${dataset_ID}/pydata/${noise_label}/ \
    --compute_file $pyspi_script_dir/pyspi_compute.py \
    --template_pbs_file $pyspi_script_dir/template.pbs \
    --pyspi_config $config_file \
    --sample_yaml ${study}/data/${dataset_ID}/pydata/${noise_label}/sample_onesubj.yaml \
    --pbs_notify abe \
    --email abry4213@uni.sydney.edu.au --walltime_hrs 6 --cpu 2 --mem 8 --table_only
done
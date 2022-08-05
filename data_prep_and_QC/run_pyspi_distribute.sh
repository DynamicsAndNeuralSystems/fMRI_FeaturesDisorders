export study=/headnode1/abry4213
export github_dir=${study}/github/
export config_file=${github_dir}/fMRI_FeaturesDisorders/data_prep_and_QC/pyspi19_config.yaml
export pyspi_script_dir=${study}/github/pyspi-distribute/

# UCLA Schizophrenia
for noise_label in AROMA_2P AROMA_2P_GMR AROMA_2P_DiCER
do
    python $pyspi_script_dir/distribute_jobs.py --data_dir ${study}/data/UCLA_Schizophrenia/pydata/${noise_label}/ \
    --compute_file $pyspi_script_dir/pyspi_compute.py \
    --template_pbs_file $pyspi_script_dir/template.pbs \
    --pyspi_config $config_file \
    --sample_yaml ${study}/data/UCLA_Schizophrenia/pydata/${noise_label}/sample.yaml \
    --pbs_notify a --overwrite_pkl \
    --email abry4213@uni.sydney.edu.au --walltime_hrs 12 --cpu 2 --mem 8 --table_only
done
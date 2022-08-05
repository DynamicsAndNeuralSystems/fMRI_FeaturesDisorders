export study=/headnode1/abry4213
export config_file=$study/scripts/scz/pyspi/myconfig.yaml
export pyspi_script_dir=${study}/github/pyspi-distribute

# UCLA SCZ AROMA+2P dataset 
python $pyspi_script_dir/distribute_jobs.py --data_dir ${study}/data/scz/UCLA/pydata/AROMA_2P/ \
--compute_file $pyspi_script_dir/pyspi_compute.py \
--template_pbs_file $pyspi_script_dir/template.pbs \
--pyspi_config $config_file \
--sample_yaml ${study}/data/scz/UCLA/pydata/AROMA_2P/sample_AROMA_2P.yaml \
--pbs_notify a \
--email abry4213@uni.sydney.edu.au --walltime_hrs 12 --cpu 2 --mem 8 --table_only
# Write calc.table pkl files to CSVs
pkl_to_csv_cmd=f"python3 {github_dir}/fMRI_FeaturesDisorders/data_prep_and_QC/pyspi_pickle_to_csv --github_dir {github_dir} --data_path {data_path} --noise_procs {noise_procs_cl} --brain_region_lookup {brain_region_lookup} --parcellation_name {parcellation_name} --dataset_ID {dataset_ID}"
os.system(pkl_to_csv_cmd)  